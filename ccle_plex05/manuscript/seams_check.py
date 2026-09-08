#!/usr/bin/env python3
"""
seams_check.py — the whole-document checks that a section-by-section pass cannot make.
Not a prose check. Looks only for the defects that survive editing: unresolved references,
unbalanced delimiters, duplicate or missing numbering, orphan fragments, and numbers that
disagree with themselves.
USAGE  python3 seams_check.py
"""
import re, collections
import drafted as D, sections as S, si as SI, refs as R

def paras():
    out=[]
    for n in dir(D):
        v=getattr(D,n)
        if n.startswith(("P_","CAP_")) and isinstance(v,list):
            for p in v:
                t = p[0] if isinstance(p,tuple) else p
                if isinstance(t,str): out.append((n,t))
    for e in SI.SI_NOTES:
        if len(e)>3:
            for b in e[3]: out.append(("SI_"+e[0], b))
    for f in SI.SI_FIGS:
        for c in f["cap"]: out.append(("SI_cap_"+str(f["no"]), c))
    # Supplementary tables carry citations too; without this the checker reports every reference
    # used only in a table as defined-but-never-cited.
    for t_ in getattr(SI, "SI_TABLES", []):
        for r in t_["rows"]:
            for c in r: out.append(("SI_tab_"+str(t_["no"]), c))
        for b in t_.get("foot", []): out.append(("SI_tab_"+str(t_["no"]), b))
    return out

P=paras(); issues=0
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
    # Table cells are fragments by design — a site name, an evidence tier, a citation list. They are
    # in P so that checks 2 and 3 see their figure and reference citations, not so that prose rules
    # are applied to them.
    if n.startswith("SI_tab_"): continue
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
     ("375","CCLE cell lines"),("387","Nusinow FIRST PAGE — must never be a cell-line count"),
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
        print("   NOTE both %s and %s present (%s) — confirm each is in its correct unit" % (a,b_,ctx))

print("\n" + "="*74)
print("%d issue(s)" % issues)
print("="*74)
