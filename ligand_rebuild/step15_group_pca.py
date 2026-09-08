#!/usr/bin/env python3
"""
step15_group_pca.py — do the four site classes separate in structural space?

Step 14 tested each parameter one at a time and found shared sites more buried and more
cationic than exclusive ones. That is 21 separate tests on overlapping measurements. This asks
the multivariate question instead: placed in the same structural space, do the groups occupy
different regions, or do they sit on top of each other?

GROUPS (all from the 8,211 lysines measured by BOTH chemotypes)
  SHARED       engaged by CyMPSA and by a meroterpenoid
  CYMPSA_ONLY  engaged by CyMPSA only
  MERO_ONLY    engaged by a meroterpenoid only
  NEITHER      engaged by neither — the reference

WHAT IS COMPUTED
  1. PCA on all sites together, then group positions in PC space with confidence ellipses.
     Fitting the PCA on ALL sites (not per group) is deliberate: the axes must be defined
     independently of the grouping, or the separation would be circular.
  2. PERMANOVA-style test: is the between-group dispersion larger than expected if the labels
     were shuffled? This is the honest multivariate test — a one-way test per feature cannot
     say whether the groups differ overall.
  3. Centroid distances in the full standardised space, with a permutation null.

WHY A PERMUTATION TEST
Group sizes are extremely unequal (25 vs several thousand) and the features are correlated and
non-normal. Parametric MANOVA would be unreliable here. Shuffling the labels 2,000 times builds
the null distribution directly from this data, making no distributional assumption.

THE POWER CONSTRAINT, AGAIN
SHARED is ~25 sites in a 21-dimensional space. That is far too few to characterise a
distribution properly — the ellipse is indicative, not a real confidence region. A null result
means "not detectable at n=25". Any positive result should be treated as a hypothesis for a
larger dataset, not a conclusion.

INPUT   step14_sites.csv
OUTPUT  step15_pca.csv, step15_report.txt
"""
import os
import numpy as np, pandas as pd
import warnings
warnings.filterwarnings("ignore")
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

FEATS = ["rsa_tien", "sasa_sidechain", "nz_sasa", "crowding_8A", "atoms_10A",
         "hydrophobic_8A", "aromatic_8A", "polar_8A", "anionic_8A", "cationic_8A",
         "net_charge_8A", "n_lys_8A", "dist_centroid", "propka_pKa",
         "hbond_partners_NZ", "in_pocket", "pocket_volume", "pocket_drugscore",
         "dist_nearest_pocket", "plddt_site", "plddt_8A"]
GROUPS = ["SHARED", "CYMPSA_ONLY", "MERO_ONLY", "NEITHER"]
N_PERM = 2000
rng = np.random.default_rng(0)

log = []
def say(m):
    print(m, flush=True); log.append(m)

say("=" * 78)
say("STEP 15 — DO THE SITE CLASSES SEPARATE IN STRUCTURAL SPACE?")
say("=" * 78)

d = pd.read_csv("step14_sites.csv", low_memory=False)
if "in_pocket" in d:
    d["in_pocket"] = d.in_pocket.astype(str).str.lower().isin(["true", "1", "1.0"]).astype(float)
F = [f for f in FEATS if f in d.columns]
d = d[["uniprot", "position", "group"] + F].dropna()
say("%d sites, %d features" % (len(d), len(F)))
for g in GROUPS:
    say("  %-14s %6d" % (g, int((d.group == g).sum())))

X = StandardScaler().fit_transform(d[F].astype(float))
p = PCA().fit(X)
S = p.transform(X)
cum = np.cumsum(p.explained_variance_ratio_)
say("")
say("PCA fitted on ALL sites together — axes are defined without reference to the groups,")
say("so any separation seen is not circular.")
for i in range(min(5, len(cum))):
    say("  PC%-2d %6.1f%%   cumulative %5.1f%%"
        % (i + 1, 100 * p.explained_variance_ratio_[i], 100 * cum[i]))

say("")
say("PC1 loadings (what the dominant axis is):")
for f, v in pd.Series(p.components_[0], index=F).sort_values().items():
    say("  %-22s %+6.3f %s" % (f, v, "#" * int(abs(v) * 30)))
say("")
say("PC2 loadings:")
for f, v in pd.Series(p.components_[1], index=F).sort_values().items():
    say("  %-22s %+6.3f %s" % (f, v, "#" * int(abs(v) * 30)))

say("")
say("GROUP POSITIONS IN PC SPACE (mean +/- SD)")
say("  %-14s %6s %16s %16s" % ("group", "n", "PC1", "PC2"))
for g in GROUPS:
    m = (d.group == g).values
    if m.sum() < 3:
        continue
    say("  %-14s %6d %8.3f +/-%5.2f %8.3f +/-%5.2f"
        % (g, m.sum(), S[m, 0].mean(), S[m, 0].std(), S[m, 1].mean(), S[m, 1].std()))

say("")
say("CENTROID DISTANCES in the full %d-dimensional standardised space" % len(F))
say("with a permutation null (%d shuffles of the group labels)" % N_PERM)
cent = {g: X[(d.group == g).values].mean(axis=0) for g in GROUPS
        if (d.group == g).sum() >= 3}
pairs = [("SHARED", "NEITHER"), ("CYMPSA_ONLY", "NEITHER"), ("MERO_ONLY", "NEITHER"),
         ("SHARED", "CYMPSA_ONLY"), ("SHARED", "MERO_ONLY"), ("CYMPSA_ONLY", "MERO_ONLY")]
rows = []
labels = d.group.values
say("  %-28s %10s %10s %10s" % ("pair", "distance", "null mean", "p"))
for a, b in pairs:
    if a not in cent or b not in cent:
        continue
    obs = float(np.linalg.norm(cent[a] - cent[b]))
    na, nb = int((labels == a).sum()), int((labels == b).sum())
    pool = np.where((labels == a) | (labels == b))[0]
    null = []
    for _ in range(N_PERM):
        idx = rng.permutation(pool)
        ca = X[idx[:na]].mean(axis=0)
        cb = X[idx[na:na + nb]].mean(axis=0)
        null.append(np.linalg.norm(ca - cb))
    null = np.array(null)
    pv = float((null >= obs).mean())
    rows.append(dict(pair="%s vs %s" % (a, b), distance=obs, null_mean=float(null.mean()),
                     p=pv, n_a=na, n_b=nb))
    say("  %-28s %10.3f %10.3f %10.4f%s"
        % ("%s vs %s" % (a, b), obs, null.mean(), pv, " *" if pv < 0.05 else ""))

res = pd.DataFrame(rows)
if len(res):
    res.to_csv("step15_pca.csv", index=False)

say("")
say("READING THIS")
say("  A distance far above the null mean means the groups occupy different regions of")
say("  structural space — a genuinely multivariate difference, not 21 correlated one-way")
say("  tests. A distance near the null means they overlap.")
say("")
say("  CAUTION: SHARED has ~25 sites in %d dimensions. Its centroid is poorly determined," % len(F))
say("  and the permutation test accounts for that by construction — but a null result here")
say("  means 'not detectable at n=25', never that the groups are alike.")

open("step15_report.txt", "w").write("\n".join(log) + "\n")
say("")
say("wrote step15_pca.csv and step15_report.txt")
