#!/usr/bin/env python3
"""
seams_check.py: the whole-document checks that a section-by-section pass cannot make.
Not a prose check. Looks only for the defects that survive editing: unresolved references,
unbalanced delimiters, duplicate or missing numbering, orphan fragments, and numbers that
disagree with themselves.
USAGE  python3 seams_check.py
"""
import re, collections
import sections as S, si as SI, refs as R, spine as SP, dedup as DD, drafted as D

def _strings(obj, depth=0):
    """Every string anywhere inside a module-level container, whatever its shape."""
    if depth > 7:
        return
    if isinstance(obj, str):
        yield obj
    elif isinstance(obj, (list, tuple)):
        for x in obj:
            yield from _strings(x, depth + 1)
    elif isinstance(obj, dict):
        for x in obj.values():
            yield from _strings(x, depth + 1)


def paras():
    """Collect the prose the BUILD uses.

    This used to walk drafted.py's P_* and CAP_* lists. The Cell Chem Biol restructure of
    19 October moved the prose into spine.py, build.py was repointed and this checker was
    not, so from that day it reported on a draft nothing builds from. It passed clean while
    the real manuscript went unchecked. Walking every container of the modules build.py
    imports, rather than naming the containers, is what stops that recurring: a container
    that moves or is renamed is still found.

    drafted.py IS walked, after checking what build.py actually does with it. build.py loads
    it through importlib and uses D.P_METHODS for the whole Methods section and D.P_2_<i> /
    D.CAP_2_<i> for the per-section Results, which spine.py then overrides where SPINE
    carries that section. So it is half-live, not dead, and the original checker's mistake
    was never that it read drafted.py; it was that it read ONLY drafted.py and so never saw
    spine.py at all. Over-inclusion is the safe direction here: a retired drafted paragraph
    that spine replaced can make a reference look cited when it is not, which section 7
    catches against the built file.
    """
    seen, out = set(), []
    for mod in (SP, S, SI, DD, D):
        for n in dir(mod):
            if n.startswith("__"):
                continue
            v = getattr(mod, n)
            if not isinstance(v, (str, list, tuple, dict)):
                continue
            for t in _strings(v):
                if len(t) > 40 and t not in seen:
                    seen.add(t)
                    out.append(("%s.%s" % (mod.__name__, n), t))
    return out


import os
OUTDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "results"))

P=paras(); issues=0

_k = set()
for _n, _t in P:
    for _m in re.findall(r"\[\[([a-z0-9,]+)\]\]", _t):
        _k |= set(_m.split(","))
if len(_k) < 20:
    raise SystemExit(
        "seams_check: only %d citation keys found in the prose sources, which is far below "
        "what this manuscript cites. The checker is almost certainly pointed at the wrong "
        "module. Refusing to report on a document it cannot see." % len(_k))

def bad(tag,msg):
    global issues; issues+=1; print("  [%s] %s" % (tag,msg))

print("="*74); print("SEAMS CHECK"); print("="*74)

# 1 delimiters
print("\n1. Unbalanced delimiters")
for n,t in P:
    for o,c,name in (("(",")","paren"),("[","]","bracket")):
        if t.count(o)!=t.count(c):
            bad(n,"%s %d open / %d close :: %s" % (name,t.count(o),t.count(c),t[:80]))
if not issues: print("   none")

# 2 figure references resolve
print("\n2. Figure references")
main = {str(i) for i in range(1,9)}
sifigs = {str(f["no"]).upper().replace("S","") for f in SI.SI_FIGS}
sifigs = {str(f["no"]) for f in SI.SI_FIGS}
cited_main=collections.Counter(); cited_si=collections.Counter(); before=issues
for n,t in P:
    for m in re.findall(r"Figure\s+(\d+)", t): cited_main[m]+=1
    for m in re.findall(r"figure\s+(S\d+)", t): cited_si[m.upper()]+=1
for k in sorted(cited_main):
    if k not in main: bad("main","Figure %s cited but no such main figure" % k)
known_si={s.upper() for s in sifigs}
for k in sorted(cited_si):
    if k not in known_si: bad("SI","figure %s cited but not in SI_FIGS (have %s)" % (k, ", ".join(sorted(known_si))))
for k in sorted(main):
    if k not in cited_main: bad("main","Figure %s never cited in prose" % k)
for k in sorted(known_si):
    if k not in cited_si: bad("SI","figure %s never cited in prose" % k)
if issues==before: print("   all resolve, all cited")

# 3 reference keys
print("\n3. Reference keys")
before=issues; used=set()
for n,t in P:
    for m in re.findall(r"\[\[([a-z0-9,]+)\]\]", t): used|=set(m.split(","))
for k in sorted(used-set(R.REFS)): bad("refs","cited key not in refs.py: %s" % k)
for k in sorted(set(R.REFS)-used): bad("refs","defined but never cited: %s" % k)
if issues==before: print("   all cited keys defined, all defined keys cited")

# 4 SI note numbering
print("\n4. Supplementary note numbering")
before=issues
nums=[int(e[0][1:]) for e in SI.SI_NOTES]
dupes=[k for k,v in collections.Counter(nums).items() if v>1]
for d in dupes: bad("SI","note N%d defined more than once" % d)
if nums!=list(range(1,len(nums)+1)): bad("SI","note numbering not contiguous 1..%d: %s" % (len(nums),nums))
refd=set()
for n,t in P:
    for m in re.findall(r"note\s+(N\d+)", t): refd.add(int(m[1:]))
for r_ in sorted(refd-set(nums)): bad("SI","prose points at note N%d which does not exist" % r_)
if issues==before: print("   contiguous, no duplicates, all pointers resolve")

# 5 orphan fragments
print("\n5. Orphan fragments and stranded headings")
before=issues
for n,t in P:
    # Table cells are fragments by design: a site name, an evidence tier, a citation list. They are
    # in P so that checks 2 and 3 see their figure and reference citations, not so that prose rules
    # are applied to them.
    # The generic walker names this source "si.SI_TABLES" where the old one named it
    # "SI_tab_<no>", so match both rather than let a rename silently re-enable the rule.
    if n.startswith("SI_tab_") or "SI_TABLES" in n: continue
    s=t.strip()
    if len(s.split())<=4 and not s.endswith((".","!","?",'."')) and not s.startswith("**"):
        bad(n,"fragment with no terminal punctuation :: %r" % s[:70])
    if s.endswith(":") and len(s.split())<8: bad(n,"stranded lead-in :: %r" % s[:70])
    if re.search(r"\b(TBD|TODO|XXX|FIXME|\[ *\])", s): bad(n,"editing marker left in :: %r" % s[:70])
if issues==before: print("   none")

# 6 numbers that must agree with themselves
print("\n6. Recurring quantities")
before=issues
KEY=[("295,411","distinct acquired-lysine sites"),("8,470","site-in-cell-line observations"),
     ("0.58%","pooled identification FDP"),("0.83%","distinct-peptide FDP"),
     ("4.70%","FDP at margin 0"),("7.04%","KGAIN FDP at margin 0"),
     ("375","CCLE cell lines"),("387","Nusinow FIRST PAGE, must never be a cell-line count"),
     ("42","CCLE plexes"),
     # Figure 2's detection census has two thresholds and they must never be swapped: 53 / 0.63%
     # is the search's reported 1% FDR, 30 / 0.35% is the margin >= 10 subset section 2.3
     # measures the error rate over. Deposited as c_total_found / c_total_found_margin in
     # results/panel_data/fig2_head.json.
     ("53","detections at the reported 1% filter"),
     ("0.63%","detection rate at the reported 1% filter"),
     ("30","detections at margin >= 10"),
     ("0.35%","detection rate at margin >= 10")]
allt=" ".join(t for _,t in P)
for tok,what in KEY:
    c=allt.count(tok)
    print("   %-10s %-38s %d occurrence(s)" % (tok, what, c))
for a,b_,ctx in [("375","372","cell-line count"),("8,470","8,055","site unit"),
                 ("0.63%","0.35%","detection rate: reported 1% filter vs margin >= 10")]:
    if a in allt and b_ in allt:
        print("   NOTE both %s and %s present (%s), confirm each is in its correct unit" % (a,b_,ctx))

# 7 the built documents
# Every check above reads source modules, which is how this file came to spend two weeks
# reporting on a draft nothing builds from. This section reads the built files instead, so it
# cannot go stale: it is the artifact that ships. One clean document serves both venues, so a
# tag left in it is submitted to the journal and posted to bioRxiv at the same time.
#
# [CUT, [WITHDRAWN], [REWRITTEN] and [NOT YET WRITTEN] are deliberate status labels on
# supplementary notes, printed so a dropped note is on the record rather than silently gone.
# READING COPY belongs only in the reading copy, never in the submission document.
import glob as _glob
import docx as _docx

_OK_TAGS = ("CUT", "WITHDRAWN", "REWRITTEN", "NOT YET WRITTEN")


def _tags(path):
    out = []
    for i, p_ in enumerate(_docx.Document(path).paragraphs):
        for m in re.finditer(r"\[([A-Z][A-Z0-9 _/-]{3,})\b", p_.text):
            if m.group(1).strip() not in _OK_TAGS:
                out.append((i, m.group(1).strip(), p_.text[:70]))
    return out


def _refs(path):
    """printed reference entries against the superscript numbers actually cited"""
    ps = _docx.Document(path).paragraphs
    entries, cited = {}, set()
    for p_ in ps:
        # "[1] " since 10 November; "1. " before that. Accept either, or every citation
        # reads as pointing at nothing.
        m = re.match(r"\s*\[?(\d{1,3})[\].]\s+(\S.{0,70})", p_.text)
        if m and len(p_.text) > 60:
            entries[int(m.group(1))] = m.group(2)
    for p_ in ps:
        if re.match(r"\s*\[?\d{1,3}[\].]\s", p_.text) and len(p_.text) > 60:
            continue
        for r_ in p_.runs:
            if r_.font.superscript and re.search(r"\d", r_.text):
                for tok in re.findall(r"\d+(?:\s*[\u2013-]\s*\d+)?", r_.text):
                    tok = tok.replace(" ", "")
                    if "\u2013" in tok or "-" in tok:
                        a_, b_ = re.split(r"[\u2013-]", tok)
                        cited.update(range(int(a_), int(b_) + 1))
                    else:
                        cited.add(int(tok))
    return entries, cited


print("\n7. The built documents")
before = issues

# The submission document is the one that matters: it goes to Cell Chemical Biology AND to
# bioRxiv, so a tag left in it is both submitted and posted permanently. The reading copy is
# expected to carry apparatus and is reported without being counted against the build.
_sub = [c for c in sorted(_glob.glob(os.path.join(OUTDIR, "Manuscript_20??-??-??.docx")))
        if "_tracked" not in c and "_MIKAIL" not in c]
_rc = sorted(_glob.glob(os.path.join(OUTDIR, "Reading_copy_*.docx")))
if not _sub:
    bad("built", "no Manuscript_<date>.docx in %s. Run build.py; the default build is the "
                 "submission document and --draft gives the reading copy" % OUTDIR)
else:
    _name = os.path.basename(_sub[-1])
    _t = _tags(_sub[-1])
    if _t:
        for _i, _k, _ctx in _t:
            bad("built", "SUBMISSION %s paragraph %d carries [%s], which would be submitted "
                         "and posted :: %s" % (_name, _i, _k, _ctx))
    else:
        print("   %s: no bracketed tags" % _name)
if _rc:
    print("   %s: %d bracketed tag(s), expected in a reading copy"
          % (os.path.basename(_rc[-1]), len(_tags(_rc[-1]))))

# The bibliography check reads the SUBMISSION document, always, not whichever file happened
# to exist. Folding it under the reading-copy branch would have skipped it entirely on a
# machine that had never run --draft.
if _sub:
    _e, _c = _refs(_sub[-1])
    for _n in sorted(set(_e) - _c):
        bad("built", "reference %d printed but never cited: %s" % (_n, _e[_n][:60]))
    for _n in sorted(_c - set(_e)):
        bad("built", "citation %d points at no printed reference" % _n)
    if issues == before:
        print("   %d references, all printed and all cited" % len(_e))

print("\n" + "=" * 74)
print("%d issue(s)" % issues)
print("=" * 74)
