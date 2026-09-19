#!/usr/bin/env python3
"""
verify.py -- check that every number a figure script computes still appears in the manuscript.

THE CLAIM THIS TESTS
Methods: "every number in this paper is printed by the script that produces the corresponding panel
rather than transcribed from the artwork." Until now that was unverifiable: the scripts printed free
text and nothing compared it to the prose. Each figure script now writes results/panel_data/
<stem>_head.json, and this reads those ledgers and searches the assembled manuscript and SI for
each value.

WHAT A MISS MEANS, IN EACH DIRECTION -- they are not symmetric
  ledger value not in the prose : usually fine. Scripts compute more than the paper quotes.
                                  Worth eyeballing once, then ignoring.
  prose number not in any ledger: THE DANGEROUS ONE. A quantity in the text that no script
                                  produces is either transcribed by hand or stale. This script
                                  cannot find those on its own -- it would need every number in
                                  the prose enumerated first -- so it reports coverage instead:
                                  what share of the ledger reaches the page, per figure.

FORMATTING
A value is searched in the spellings the manuscript actually uses: bare (617), thousands-separated
(8,470), one and two decimals (7.28, 4.38), and percent-suffixed. A hit on any spelling counts.
Small integers below 10 are skipped -- "5" matches everything and proves nothing.
"""
import os, sys, json, glob, re, importlib.util

TOP = os.path.expanduser("~/mnt/DepMap")
MAN = os.path.join(TOP, "ccle", "manuscript")
LED = os.path.join(TOP, "results", "panel_data")

sys.path.insert(0, MAN)
import sections as S, si as SI
spec = importlib.util.spec_from_file_location("drafted", os.path.join(MAN, "drafted.py"))
D = importlib.util.module_from_spec(spec); spec.loader.exec_module(D)

def prose():
    out = []
    for n in dir(D):
        v = getattr(D, n)
        if n.startswith(("P_", "CAP_")) and isinstance(v, list):
            out += [x[0] if isinstance(x, tuple) else x for x in v if isinstance(
                    x[0] if isinstance(x, tuple) else x, str)]
    for e in SI.SI_NOTES:
        if len(e) > 3: out += e[3]
    for f in SI.SI_FIGS: out += f["cap"]
    for t in getattr(SI, "SI_TABLES", []):
        out += [c for r in t["rows"] for c in r] + list(t.get("foot", []))
    out += list(S.HIGHLIGHTS["body"])
    return " ".join(out)

TEXT = prose()

def flatten(key, v):
    """A HEAD value may be a dict of sub-values or a list. The manuscript quotes the members, not
    the container, so a dict compared whole can never match and would be reported as a false miss.
    Flatten one level into (key.subkey, scalar) pairs; strings that are prose rather than numbers
    are dropped, because a sentence is not a checkable quantity."""
    out = []
    if isinstance(v, dict):
        for k2, v2 in v.items():
            out += flatten("%s.%s" % (key, k2), v2)
    elif isinstance(v, (list, tuple)):
        for i, v2 in enumerate(v):
            out += flatten("%s[%d]" % (key, i), v2)
    elif isinstance(v, str):
        # A bare label ("L", "SKCM") is not a number. Neither is a formatted sentence kept in HEAD
        # for a human reading the log -- those hold several numbers and match nothing as a whole,
        # so counting them produces false misses. The scalars beside them are what gets checked.
        if len(v) > 60 or v.count(" ") > 6:
            return out
        if re.search(r"\d", v):
            out.append((key, v))
    else:
        out.append((key, v))
    return out


MINUS = "\u2212"   # the manuscript sets negatives with a true minus sign, not a hyphen

def _both_signs(forms):
    out = []
    for f in forms:
        out.append(f)
        if f.startswith("-"):
            out.append(MINUS + f[1:])
    return out


def spellings(v):
    if isinstance(v, bool): return []
    if isinstance(v, int):
        if abs(v) < 10: return []
        return _both_signs([str(v), "{:,}".format(v)])
    if isinstance(v, float):
        out = []
        for d in (0, 1, 2, 3):
            s = ("%%.%df" % d) % v
            out += [s, s + "%"]
        if v != 0 and abs(v) < 1e-3:
            out += ["%.1e" % v, ("%.1e" % v).replace("e-0", " × 10⁻")]
        return _both_signs(out)
    return [str(v)]

print("=" * 78)
print("VERIFY -- figure ledgers against the assembled manuscript")
print("=" * 78)
led = sorted(glob.glob(os.path.join(LED, "*_head.json")))
if not led:
    sys.exit("no ledgers in %s -- run the figure scripts first" % LED)

grand_hit = grand_tot = 0
for path in led:
    d = json.load(open(path))
    nums = d["numbers"]
    hits, missing = 0, []
    checked = 0
    flat = []
    for k, v in nums.items():
        flat += flatten(k, v)
    for k, v in flat:
        sp = spellings(v)
        if not sp:
            continue
        checked += 1
        if any(s in TEXT for s in sp):
            hits += 1
        else:
            missing.append("%s=%s" % (k, v))
    grand_hit += hits; grand_tot += checked
    print("\n%-14s %d of %d checkable values appear in the prose"
          % (d["script"], hits, checked))
    if missing:
        print("   not quoted: %s" % ", ".join(missing))

print("\n" + "=" * 78)
print("coverage: %d of %d ledger values reach the page (%.0f%%)"
      % (grand_hit, grand_tot, 100.0 * grand_hit / max(grand_tot, 1)))
print("ledgers present: %d of 18 figures" % len(led))
print("A value not quoted is usually fine. A number in the prose that NO ledger produces is the")
print("real risk, and needs every figure wired before it can be detected.")
