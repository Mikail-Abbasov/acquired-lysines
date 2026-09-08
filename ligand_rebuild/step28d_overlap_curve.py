#!/usr/bin/env python3
"""
step28d_overlap_curve.py — the overlap result with NO chosen threshold.

WHY
Figure 8 panel c originally reported the fraction of each population inside the region holding
90% of engaged lysines. 90% was a convention, not a derived quantity, and a reader is entitled
to ask why that number. This emits the whole curve instead: for every radius, the fraction of
each population within it. The 90% figure is then one readable point on a curve rather than a
choice the figure has to defend.

Distance is Euclidean in the nine-feature space, standardised on the ENGAGED population and
measured from the engaged centroid -- so the axis is "how far from a typical engaged lysine",
in engaged standard deviations, and it does not move when the union gets bigger.

INPUT   union_distribution.csv, measured_cb.csv, step28_engaged_keys.csv, features_expanded.csv
OUTPUT  step28d_curve.csv, step28d_report.txt
"""
import numpy as np, pandas as pd
from sklearn.preprocessing import StandardScaler

FEATS = ["crowding_8A", "atoms_10A", "hydrophobic_8A", "aromatic_8A", "polar_8A",
         "anionic_8A", "cationic_8A", "net_charge_8A", "n_lys_8A"]
log = []
def say(s=""):
    print(s, flush=True); log.append(s)

ek = set(map(tuple, pd.read_csv("step28_engaged_keys.csv").values))
feat = pd.read_csv("features_expanded.csv", low_memory=False).drop_duplicates(["uniprot", "position"])
key_eng = {(u, int(p)) for u, p in zip(feat.uniprot, feat.position) if (u, int(p)) in ek}

u = pd.read_csv("union_distribution.csv")
m = pd.read_csv("measured_cb.csv")
m["engaged"] = [(x, int(p)) in key_eng for x, p in zip(m.uniprot, m.position)]

def grab(df, pre=""):
    return df[[pre + f for f in FEATS]].apply(pd.to_numeric, errors="coerce").dropna().values

U = grab(u); E = grab(m[m.engaged], "cb_"); N = grab(m[~m.engaged], "cb_")
sc = StandardScaler().fit(E)
c = sc.transform(E).mean(axis=0)
dU = np.linalg.norm(sc.transform(U) - c, axis=1)
dE = np.linalg.norm(sc.transform(E) - c, axis=1)
dN = np.linalg.norm(sc.transform(N) - c, axis=1)

say("=" * 78)
say("STEP 28d - OVERLAP WITHOUT A CHOSEN THRESHOLD")
say("=" * 78)
say("populations: union %d, engaged %d, labelled-not-engaged %d" % (len(U), len(E), len(N)))
say("distance from the engaged centroid, in engaged SD units")
say("  %-26s %8s %8s %8s %8s" % ("population", "median", "q25", "q75", "q90"))
for lab, d in (("acquired union", dU), ("engaged", dE), ("labelled, not engaged", dN)):
    say("  %-26s %8.2f %8.2f %8.2f %8.2f"
        % (lab, np.median(d), np.quantile(d, .25), np.quantile(d, .75), np.quantile(d, .90)))

grid = np.linspace(0, np.quantile(dU, 0.995), 400)
out = pd.DataFrame(dict(radius=grid,
                        union=[(dU <= r).mean() for r in grid],
                        engaged=[(dE <= r).mean() for r in grid],
                        not_engaged=[(dN <= r).mean() for r in grid]))
out.to_csv("step28d_curve.csv", index=False)

say()
say("THE SINGLE NUMBER THE OLD PANEL REPORTED, AND ITS NEIGHBOURS")
say("  %-8s %-10s %10s %10s" % ("engaged", "radius", "union in", "not-eng in"))
for q in (0.50, 0.75, 0.90, 0.95):
    r = np.quantile(dE, q)
    say("  %-8.0f%% %-10.2f %9.1f%% %9.1f%%" % (100 * q, r, 100 * (dU <= r).mean(),
                                                100 * (dN <= r).mean()))
say()
say("THRESHOLD-FREE SUMMARIES, for the legend")
med = np.median(dE)
say("  median acquired lysine sits at the %.0fth percentile of the engaged distribution"
    % (100 * (dE < np.median(dU)).mean()))
say("  median not-engaged lysine sits at the %.0fth percentile"
    % (100 * (dE < np.median(dN)).mean()))
say("  P(a random acquired lysine is closer to the engaged centre than a random engaged one)")
rng = np.random.default_rng(0)
a = rng.choice(dU, 200000); b = rng.choice(dE, 200000)
say("    = %.3f   (0.5 would mean the two are interchangeable)" % (a < b).mean())
a = rng.choice(dN, 200000)
say("    same for labelled-not-engaged = %.3f" % (a < b).mean())

open("step28d_report.txt", "w").write("\n".join(log) + "\n")
print("\nwrote step28d_curve.csv, step28d_report.txt")
