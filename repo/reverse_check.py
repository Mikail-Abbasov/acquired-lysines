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

TOP = os.path.expanduser("~/mnt/DepMap")
MAN = os.path.join(TOP, "ccle_plex05", "manuscript")
LED = os.path.join(TOP, "results", "panel_data")

sys.path.insert(0, MAN)
import sections as S, si as SI
spec = importlib.util.spec_from_file_location("drafted", os.path.join(MAN, "drafted.py"))
D = importlib.util.module_from_spec(spec); spec.loader.exec_module(D)

MINUS = "−"


def blocks():
    """(where, text) so a finding can be pointed at its section rather than the whole paper."""
    out = []
    for n in sorted(dir(D)):
        v = getattr(D, n)
        if n.startswith(("P_", "CAP_")) and isinstance(v, list):
            for i, x in enumerate(v):
                t = x[0] if isinstance(x, tuple) else x
                if isinstance(t, str):
                    out.append(("%s[%d]" % (n, i), t))
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
