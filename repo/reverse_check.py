#!/usr/bin/env python3
"""
reverse_check.py -- the dangerous direction. Find numbers in the prose that NO figure script
produces.

verify.py checks ledger -> prose: every number a script computes, does it still appear on the
page? A miss there is usually harmless, because scripts compute more than the paper quotes.

This checks prose -> ledger, which is the direction that hides errors. A quantity in the text
that no script emits is one of four things:

  1. legitimately not from a figure  -- a count from Methods, a reference year, a PDB code,
                                        a residue number inside a mutation label (E17K)
  2. from an analysis script that has no ledger  -- step reports, tables
  3. transcribed by hand from an earlier run, and now stale
  4. never computed at all

Only 3 and 4 are defects, but 1 and 2 cannot be told apart from them automatically, so this
script's job is to shrink the haystack: it filters out the shapes that are structurally not
results, and prints what is left for a human to walk.

WHAT COUNTS AS A CANDIDATE RESULT
  - any decimal (0.617, 4.38, 1.68)
  - any thousands-separated integer (8,470, 245,684)
  - any integer >= 100 that is not inside a mutation label, a PDB id, a year, or a figure/
    section number
Integers below 100 are dropped: they match too much to be evidence either way.
"""
import os, sys, json, glob, re, importlib.util

# Resolved from this file's own location since 3 November 2026. It was hardcoded to
# "~/mnt/DepMap", which is one sandbox's mount path: it does not exist on the author's own
# machine, and certainly not in a stranger's clone, so every path below failed for anyone
# who actually took the deposit at its word. The repository root is the parent of repo/.
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)


def _firstdir(*cands):
    for c in cands:
        if c and os.path.isdir(c):
            return c
    return cands[0]


TOP = _firstdir(_ROOT,
                os.path.expanduser("~/mnt/DepMap"),
                os.path.expanduser("~/Desktop/DepMap"))
MAN = _firstdir(os.path.join(TOP, "ccle", "manuscript"),
                os.path.join(TOP, "ccle_plex05", "manuscript"))
LED = _firstdir(os.path.join(TOP, "panel_data"),
                os.path.join(TOP, "results", "panel_data"))

sys.path.insert(0, MAN)
import sections as S, si as SI
spec = importlib.util.spec_from_file_location("drafted", os.path.join(MAN, "drafted.py"))
D = importlib.util.module_from_spec(spec); spec.loader.exec_module(D)

MINUS = "−"


def blocks():
    """(where, text) so a finding can be pointed at its section rather than the whole paper."""
    # Walked only drafted.py's P_* and CAP_* lists until 3 November 2026. The Cell Chem Biol
    # restructure of 19 October moved the prose into spine.py, so from that day this check
    # read a draft nothing builds from and reported on numbers the paper no longer contains.
    # seams_check.py had the identical bug. Walking every container of every live prose module,
    # rather than naming containers in one module, is what stops it recurring.
    out = []
    def _str(o, d=0):
        if d > 7:
            return
        if isinstance(o, str):
            yield o
        elif isinstance(o, (list, tuple)):
            for x in o:
                yield from _str(x, d + 1)
        elif isinstance(o, dict):
            for x in o.values():
                yield from _str(x, d + 1)
    import importlib as _il
    _seen = set()
    for _mn in ("spine", "sections", "si", "dedup", "drafted"):
        try:
            _m = _il.import_module(_mn)
        except Exception:
            continue
        for _n in sorted(dir(_m)):
            if _n.startswith("__"):
                continue
            _v = getattr(_m, _n)
            if not isinstance(_v, (str, list, tuple, dict)):
                continue
            for _t in _str(_v):
                if len(_t) > 40 and _t not in _seen:
                    _seen.add(_t)
                    out.append(("%s.%s" % (_mn, _n), _t))
    if len(out) < 40:
        raise SystemExit("reverse_check: only %d prose blocks found. The checker is probably "
                         "pointed at the wrong modules. Refusing to report a clean paper."
                         % len(out))
    for e in SI.SI_NOTES:
        if len(e) > 3:
            for i, t in enumerate(e[3]):
                out.append(("note %s[%d]" % (e[0], i), t))
    for f in SI.SI_FIGS:
        for i, t in enumerate(f["cap"]):
            out.append(("%s cap[%d]" % (f["no"], i), t))
    for t_ in getattr(SI, "SI_TABLES", []):
        for r in t_["rows"]:
            out.append(("Table %s row" % t_["no"], " | ".join(str(c) for c in r)))
    for i, t in enumerate(S.HIGHLIGHTS["body"]):
        out.append(("Highlight[%d]" % i, t))
    return out


def ledger_forms():
    """Every spelling of every ledger value, so a prose numeral can be matched against it."""
    forms = set()

    def add(x):
        if isinstance(x, bool):
            return
        if isinstance(x, (int, float)):
            f = float(x)
            for s_ in ("%d" % round(f), "%.1f" % f, "%.2f" % f, "%.3f" % f, "%.4f" % f,
                       "%g" % f):
                forms.add(s_.lstrip("+"))
                forms.add(s_.lstrip("+").replace("-", MINUS))
            forms.add("{:,}".format(int(round(abs(f)))))
            # The prose sets a negative as "-0.307" but NUM captures the digits alone, so the
            # unsigned spellings have to be in the set or every negative effect size in the
            # paper reads as a miss. This was 20 false positives across figures 4 and S3.
            for s_ in ("%.0f" % abs(f), "%.1f" % abs(f), "%.2f" % abs(f), "%.3f" % abs(f),
                       "%.4f" % abs(f), "%g" % abs(f)):
                forms.add(s_)
            # scientific notation as the manuscript sets it: 1.6 x 10^-4
            if f != 0 and (abs(f) < 1e-3 or abs(f) >= 1e5):
                m, e = ("%e" % abs(f)).split("e")
                forms.add("%.1f" % float(m))
                forms.add("%.2f" % float(m))
        elif isinstance(x, str):
            forms.add(x.strip())
            m = re.match(r"^([\d.]+)e([+-]\d+)$", x.strip())
            if m:
                forms.add(m.group(1))
        elif isinstance(x, dict):
            for v in x.values():
                add(v)
        elif isinstance(x, (list, tuple)):
            for v in x:
                add(v)

    n = 0
    for p in sorted(glob.glob(os.path.join(LED, "*_head.json"))):
        for v in json.load(open(p))["numbers"].values():
            add(v); n += 1
    return forms, n


# numerals that are structurally not results
MUT = re.compile(r"[A-Z]\d+[A-Z*]")                     # E17K, V600K, K677
# Cell-line and reagent names carry digits that are not quantities. Without this every mention
# of MDA-MB-231 reads as an unsourced "231" -- 18 of the first 143 findings were exactly that.
NAMES = re.compile(r"MDA-MB-\d+|NCI-?H\d+|KU-?19-?19|CVCL_\w+|SW\d+|HCT-?\d+|CAL-?\d+|"
                   r"KYSE\d+|HEC\d+|A\d{3,4}|MCF\d|REH|MDST\d+|IGR\d+|8-?UVY|dTB-STP")
PDB = re.compile(r"\b\d[A-Za-z0-9]{3}\b")               # 8UVY, 6CMS
CITE = re.compile(r"\b(19|20)\d{2}\b")                  # years
FIGREF = re.compile(r"(?:[Ff]igure|[Ff]ig\.?|[Ss]ection|[Tt]able|note|panel|PDB|RRID|CVCL|"
                    r"[Vv]ersion|v)\s*S?\d")
NUM = re.compile(r"(?<![\w.])(\d{1,3}(?:,\d{3})+|\d+\.\d+|\d+)(?![\w])")


def candidates(text):
    """Numerals in one block that could be a result."""
    masked = MUT.sub(" ", NAMES.sub(" ", text))
    out = []
    for m in NUM.finditer(masked):
        s_ = m.group(1)
        lo = max(0, m.start() - 24)
        ctx = masked[lo:m.end() + 6]
        if FIGREF.search(ctx):
            continue
        if CITE.fullmatch(s_):
            continue
        if PDB.fullmatch(s_) and not s_.isdigit():
            continue
        if "," in s_ or "." in s_:
            out.append(s_)
        elif int(s_) >= 100:
            out.append(s_)
    return out


def main():
    forms, nled = ledger_forms()
    print("=" * 78)
    print("REVERSE CHECK -- prose numbers that no figure ledger produces")
    print("=" * 78)
    print("ledger values loaded: %d, in %d spellings" % (nled, len(forms)))
    print()
    miss, seen = [], set()
    for where, text in blocks():
        for s_ in candidates(text):
            v = s_.replace(",", "")
            hit = (s_ in forms or v in forms or
                   ("{:,}".format(int(v)) in forms if v.isdigit() else False))
            if not hit and (where, s_) not in seen:
                seen.add((where, s_))
                i = text.find(s_)
                miss.append((where, s_, text[max(0, i - 70):i + 40].replace("\n", " ")))
    by = {}
    for where, s_, ctx in miss:
        by.setdefault(where, []).append((s_, ctx))
    for where in sorted(by):
        print("%s" % where)
        for s_, ctx in by[where]:
            print("    %-12s ...%s..." % (s_, ctx.strip()))
        print()
    print("=" * 78)
    print("%d prose numerals in %d blocks match no ledger value." % (len(miss), len(by)))
    print("Most will be Methods counts, step-report figures or reference numbers. Walk them.")


if __name__ == "__main__":
    main()
