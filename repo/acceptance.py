#!/usr/bin/env python3
"""
acceptance.py -- does the deposit actually reproduce the paper's numbers?

The claim a deposit makes is that someone else can get the same answers. This tests it the only
way that counts: run each figure script inside the assembled release/ tree and compare the ledger
it writes against the ledger deposited with the paper.

THREE OUTCOMES, and the third is the one worth having
  MATCH    every deposited value reproduced
  PARTIAL  reproduced values agree, but the run produced FEWER -- a panel was skipped because a
           regenerable input has not been built yet in a fresh tree. Figure scripts skip such
           panels silently; this makes the skip visible and names the missing values.
  DIFFER   a value changed. That is a real failure and the run prints which key.

Usage:  python3 acceptance.py [figure_stem ...]      (default: all with a deposited ledger)
Run it from repo/, after build_release.py.
"""
import os, sys, json, subprocess, glob

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
REL = os.path.join(TOP, "release")
DEPOSITED = os.path.join(REL, "panel_data")
PRODUCED = os.path.join(REL, "results", "panel_data")

SCRIPT = {
    "fig1": "ccle/fig1_v3.py",   "fig2": "ccle/fig2_v3.py",
    "fig3": "ccle/fig3_v3.py",   "fig4": "ccle/fig4_v3.py",
    # Keyed on the LEDGER stem, which is the source figure number. fig5 and fig6 are the two
    # sources that merge into published Figure 5; fig7 and fig8 are published 6 and 7.
    "fig5": "ccle/fig5_v2.py",   "fig6": "ccle/fig6_v2.py",
    "fig7": "ccle/fig7_v2.py",   "fig8": "ccle/fig8_v2.py",
    "figS1": "ccle/figS1_v3.py", "figS2": "ccle/figS2_v3.py",
    "figS3": "ccle/figS3_v3.py", "figS4": "ccle/figS4_v3.py",
    "figS5": "ccle/figS5_v3.py", "figS6": "ccle/figS6_v1.py",
    "figS7": "ccle/figS7_v1.py", "figS8": "ccle/figS8_v1.py",
    "figS9": "ccle/figS9_v1.py", "figS10": "ccle/figS10_v1.py",
}


def compare(stem):
    dep = os.path.join(DEPOSITED, "%s_head.json" % stem)
    got = os.path.join(PRODUCED, "%s_head.json" % stem)
    if not os.path.exists(dep):
        return "NO LEDGER", []
    if not os.path.exists(got):
        return "NOT RUN", []
    a = json.load(open(dep))["numbers"]
    b = json.load(open(got))["numbers"]
    differ = [k for k in a if k in b and json.dumps(a[k], sort_keys=True, default=str)
              != json.dumps(b[k], sort_keys=True, default=str)]
    missing = [k for k in a if k not in b]
    if differ:
        return "DIFFER", differ
    if missing:
        return "PARTIAL", missing
    return "MATCH", []


def main():
    stems = sys.argv[1:] or sorted(SCRIPT)
    env = dict(os.environ)
    env["PYTHONPATH"] = os.path.join(TOP, "figtools", "py")
    rows = []
    for stem in stems:
        rel = SCRIPT.get(stem)
        if rel is None:
            continue
        d = os.path.join(REL, os.path.dirname(rel))
        try:
            subprocess.run([sys.executable, os.path.basename(rel)], cwd=d, env=env,
                           stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, timeout=900)
        except subprocess.TimeoutExpired:
            rows.append((stem, "TIMEOUT", [])); continue
        rows.append((stem,) + compare(stem))

    print("=" * 78)
    print("ACCEPTANCE -- release tree reproduces the deposited ledgers")
    print("=" * 78)
    for stem, verdict, keys in rows:
        print("  %-8s %-8s %s" % (stem, verdict,
              ("" if not keys else "%d: %s%s" % (len(keys), ", ".join(keys[:4]),
                                                " ..." if len(keys) > 4 else ""))))
    n = {v: sum(1 for _, x, _ in rows if x == v) for v in
         {x for _, x, _ in rows}}
    print("\n  " + " | ".join("%s %d" % (k, v) for k, v in sorted(n.items())))
    return 1 if any(v == "DIFFER" for _, v, _ in rows) else 0


if __name__ == "__main__":
    sys.exit(main())
