#!/usr/bin/env python3
"""
step10_composite.py — do the features act together as one physical situation?

THE CLAIM BEING TESTED
Step 9 reported seven or eight features that each distinguish engaged lysines, and it is
tempting to rank them. The better reading is that they are facets of a single situation: a
lysine sitting in a pocket is simultaneously buried, crowded, ringed by hydrophobic and
aromatic neighbours, nearer the protein core, and — because burial shifts local electrostatics
— lower in predicted pKa. Not seven findings. One situation measured seven ways.

That claim is testable, and this script tests it three ways.

  1. HOW MANY INDEPENDENT DIMENSIONS ARE THERE?
     Correlation matrix and PCA. If eight features collapse onto two or three components, the
     honest statement is "two or three properties", not "eight". The Part 1 plan raised this
     ("items 1/2/5 and 3/4 are near-duplicates") and it has never been quantified here.

  2. IS THE JOINT EFFECT LARGER THAN ANY SINGLE FEATURE?
     A composite score is built and its separation compared against the best single feature.
     If "taken together" means anything, the composite must beat its own best component by
     more than fitting noise would explain.

  3. IS THE GAIN REAL OR JUST EXTRA PARAMETERS?
     **Cross-validated.** Any model with more parameters fits its own training data better;
     that is not evidence. The composite is fitted on one set of PROTEINS and evaluated on
     proteins it has never seen. Splitting by protein rather than by site matters: lysines on
     the same protein are correlated, so a random site-level split would leak information
     between train and test and inflate the result.

REPORTED HONESTLY
  - AUC with a confidence interval from repeated splits, not a single number.
  - The null: a well-powered AUC near 0.5-0.65 is what the Part 1 ceiling analysis predicted
    is achievable given label reproducibility. Beating 0.65 would be surprising.
  - Per-feature loadings, so the composite is interpretable rather than a black box.

WHAT THIS CANNOT SHOW
A composite that separates groups is still not a predictor of individual sites. Separation of
distributions and prediction of cases are different things, and with effect sizes near 0.5 SD
the distributions overlap heavily. This measures how much the features share, not whether
ligandability can be predicted — step 7 already answered that (pseudo-R^2 < 0.03).

INPUT   features_expanded.csv, occ_long.csv
OUTPUT  step10_composite.csv, step10_report.txt
"""
import os
import numpy as np, pandas as pd
import warnings
warnings.filterwarnings("ignore")
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score
from scipy import stats as sstats

TIER = int(os.environ.get("TIER", 60))
N_SPLITS = 20
FEATS = ["rsa_tien", "sasa_sidechain", "nz_sasa", "crowding_8A", "atoms_10A",
         "hydrophobic_8A", "aromatic_8A", "polar_8A", "dist_centroid",
         "net_charge_8A", "anionic_8A", "cationic_8A", "propka_pKa", "plddt_site"]

log = []
def say(m):
    print(m, flush=True); log.append(m)

say("=" * 78)
say("STEP 10 — DO THE FEATURES ACT TOGETHER?")
say("engaged = blocked >=%d%% by at least one compound" % TIER)
say("=" * 78)

feat = pd.read_csv("features_expanded.csv", low_memory=False)
occ = pd.read_csv("occ_long.csv", low_memory=False)
eng = set(map(tuple, occ[occ.occupancy_pct >= TIER][["uniprot", "position"]]
              .drop_duplicates().values))
d = feat.drop_duplicates(["uniprot", "position"]).copy()
d["engaged"] = [(u, int(p)) in eng for u, p in zip(d.uniprot, d.position)]
F = [f for f in FEATS if f in d.columns]
d = d[["uniprot", "position", "engaged"] + F].dropna()
say("%d distinct lysines | %d engaged (%.2f%%) | %d proteins | %d features"
    % (len(d), d.engaged.sum(), 100 * d.engaged.mean(), d.uniprot.nunique(), len(F)))

X = StandardScaler().fit_transform(d[F].astype(float))
y = d.engaged.values

# ---------------------------------------------------------------- 1. dimensionality
say("")
say("1. HOW MANY INDEPENDENT DIMENSIONS?")
p = PCA().fit(X)
cum = np.cumsum(p.explained_variance_ratio_)
say("  variance explained by successive components:")
for i in range(min(6, len(cum))):
    say("    PC%-2d %6.1f%%   cumulative %5.1f%%"
        % (i + 1, 100 * p.explained_variance_ratio_[i], 100 * cum[i]))
n80 = int(np.searchsorted(cum, 0.80) + 1)
n90 = int(np.searchsorted(cum, 0.90) + 1)
say("  %d components carry 80%% of the variance; %d carry 90%%" % (n80, n90))
say("  -> %d features, but effectively ~%d independent properties" % (len(F), n80))

say("")
say("  PC1 loadings (what the dominant axis actually is):")
l1 = pd.Series(p.components_[0], index=F).sort_values()
for k, v in l1.items():
    bar = "#" * int(abs(v) * 40)
    say("    %-18s %+6.3f %s" % (k, v, bar))

say("")
say("  strongest feature pairs (|r| > 0.5):")
C = pd.DataFrame(X, columns=F).corr()
seen = set()
for a in F:
    for b in F:
        if a >= b or (b, a) in seen:
            continue
        seen.add((a, b))
        if abs(C.loc[a, b]) > 0.5:
            say("    %-18s %-18s r = %+.3f" % (a, b, C.loc[a, b]))

# ---------------------------------------------------------------- 2/3. composite, CV
say("")
say("2-3. IS THE JOINT EFFECT LARGER, AND IS IT REAL?")
say("  Cross-validated by PROTEIN: the model never sees a test protein during fitting.")
say("  Splitting by site instead would leak information, since lysines on one protein")
say("  are correlated.")

rng = np.random.default_rng(0)
prots = d.uniprot.unique()
single_auc = {f: [] for f in F}
comp_auc, comp_train = [], []
coefs = []
for it in range(N_SPLITS):
    test_p = set(rng.choice(prots, size=max(1, len(prots) // 5), replace=False))
    te = d.uniprot.isin(test_p).values
    tr = ~te
    if y[tr].sum() < 20 or y[te].sum() < 10:
        continue
    sc = StandardScaler().fit(d.loc[tr, F].astype(float))
    Xtr, Xte = sc.transform(d.loc[tr, F].astype(float)), sc.transform(d.loc[te, F].astype(float))
    m = LogisticRegression(max_iter=2000, C=1.0).fit(Xtr, y[tr])
    comp_auc.append(roc_auc_score(y[te], m.decision_function(Xte)))
    comp_train.append(roc_auc_score(y[tr], m.decision_function(Xtr)))
    coefs.append(m.coef_[0])
    for j, f in enumerate(F):
        single_auc[f].append(roc_auc_score(y[te], Xte[:, j]))

def ci(v):
    v = np.array(v)
    return v.mean(), np.percentile(v, 2.5), np.percentile(v, 97.5)

say("")
say("  %-20s %8s %18s" % ("model", "AUC", "95% CI over splits"))
best_f, best_v = None, 0
for f in F:
    if not single_auc[f]:
        continue
    m_, lo, hi = ci(single_auc[f])
    v = max(m_, 1 - m_)                      # direction-agnostic
    if v > best_v:
        best_f, best_v = f, v
for f in sorted(F, key=lambda z: -max(np.mean(single_auc[z]), 1 - np.mean(single_auc[z]))
                if single_auc[z] else 0)[:5]:
    m_, lo, hi = ci(single_auc[f])
    say("  %-20s %8.3f     [%.3f, %.3f]" % (f[:20], max(m_, 1 - m_), min(lo, 1 - hi), max(hi, 1 - lo)))
if comp_auc:
    m_, lo, hi = ci(comp_auc)
    tm, tlo, thi = ci(comp_train)
    say("  %-20s %8.3f     [%.3f, %.3f]   <- all %d together, held-out proteins"
        % ("COMPOSITE (test)", m_, lo, hi, len(F)))
    say("  %-20s %8.3f     [%.3f, %.3f]   (training fit, for reference only)"
        % ("composite (train)", tm, tlo, thi))
    say("")
    say("  gain over best single feature (%s): %+.3f AUC" % (best_f, m_ - best_v))
    say("  optimism (train minus test): %+.3f" % (tm - m_))
    if tm - m_ > 0.05:
        say("  !! substantial optimism — the training fit overstates real performance")

    say("")
    say("  mean coefficients across splits (direction and weight in the composite):")
    cf = pd.Series(np.mean(coefs, axis=0), index=F).sort_values(key=abs, ascending=False)
    for k, v in cf.items():
        say("    %-18s %+7.3f" % (k, v))

say("")
say("INTERPRETATION GUIDE")
say("  The Part 1 ceiling analysis put the maximum achievable AUC for any site-level")
say("  classifier at 0.55-0.65, given how reproducibly these labels replicate. An AUC in")
say("  that band means the model is at the ceiling, not that it is weak.")
say("  A composite that beats its best component shows the features carry complementary")
say("  information — which is the 'taken together' claim. It does NOT show ligandability")
say("  is predictable: step 7 measured that directly (pseudo-R^2 < 0.03).")

pd.DataFrame(dict(feature=F, mean_coef=np.mean(coefs, axis=0) if coefs else np.nan,
                  single_auc=[np.mean(single_auc[f]) if single_auc[f] else np.nan for f in F],
                  pc1_loading=p.components_[0])).to_csv("step10_composite.csv", index=False)
open("step10_report.txt", "w").write("\n".join(log) + "\n")
say("")
say("wrote step10_composite.csv and step10_report.txt")
