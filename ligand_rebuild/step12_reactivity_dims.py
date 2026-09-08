#!/usr/bin/env python3
"""
step12_reactivity_dims.py — are the four surviving reactivity features four things or one?

WHY THIS IS NEEDED
Step 11 gated on structural confidence and found that four features survive and strengthen,
while the shape features collapse:

    anionic_8A      d -0.307 -> -0.454
    salt_bridge     d -0.221 -> -0.327
    propka_pKa      d -0.263 -> -0.286
    net_charge_8A   d +0.347 -> +0.352

But these four are chained by construction:
  - net_charge_8A = cationic_8A - anionic_8A, so it contains anionic_8A
  - a salt bridge IS an anionic residue within 4 A, so it is a stricter anionic_8A
  - PROPKA computes pKa partly FROM those same nearby charges

So "four independent effects" may be one effect counted four times. Step 10 showed exactly this
problem for the shape features: 14 measurements collapsed to about 5 real properties, with
pairs correlating up to r = 0.96. The same check must be applied here before any claim about
what governs reactivity.

DIFFERENCES FROM STEP 10, both deliberate
  1. Run WITHIN THE GATE (pLDDT >= 70). These effects only exist there; including disordered
     sites would mix in the geometry we just showed is unreliable.
  2. Run on the SURVIVING FEATURES ONLY. Including the collapsed shape features would let them
     dominate PC1 and hide whatever structure exists among the electrostatic ones.

A CAUTION ON METHOD
With four variables PCA can only say so much — you cannot extract more structure than the
variables contain. The correlation matrix is likely more informative than the components, so
both are reported, and partial correlations are given as well: does each feature still relate
to labelling once the others are held constant? That is the direct test of whether they are
separate effects.

INPUT   step9_all_lysines.csv
OUTPUT  step12_reactivity_dims.csv, step12_report.txt
"""
import os
import numpy as np, pandas as pd
import warnings
warnings.filterwarnings("ignore")
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from scipy import stats as sstats

GATE = int(os.environ.get("PLDDT_GATE", 70))
SURVIVORS = ["anionic_8A", "salt_bridge", "propka_pKa", "net_charge_8A", "cationic_8A"]

log = []
def say(m):
    print(m, flush=True); log.append(m)

say("=" * 78)
say("STEP 12 — ARE THE SURVIVING REACTIVITY FEATURES INDEPENDENT?")
say("gated at pLDDT >= %d, surviving features only" % GATE)
say("=" * 78)

allK = pd.read_csv("step9_all_lysines.csv", low_memory=False)
allK["salt_bridge"] = allK.salt_bridge.astype(str).str.lower().isin(["true", "1", "1.0"]).astype(float)
d = allK[allK.plddt_site >= GATE].copy()
F = [f for f in SURVIVORS if f in d.columns]
d = d[["uniprot", "position", "labelled"] + F].dropna()
say("%d lysines above the gate | labelled %d (%.1f%%) | %d proteins"
    % (len(d), d.labelled.sum(), 100 * d.labelled.mean(), d.uniprot.nunique()))

X = StandardScaler().fit_transform(d[F].astype(float))
Xd = pd.DataFrame(X, columns=F)

say("")
say("1. CORRELATION MATRIX — how much do they overlap?")
C = Xd.corr()
say("  %-16s %s" % ("", "  ".join("%14s" % f[:14] for f in F)))
for a in F:
    say("  %-16s %s" % (a[:16], "  ".join("%14.3f" % C.loc[a, b] for b in F)))
say("")
say("  pairs with |r| > 0.5:")
found = False
for i, a in enumerate(F):
    for b in F[i + 1:]:
        if abs(C.loc[a, b]) > 0.5:
            say("    %-16s %-16s r = %+.3f" % (a, b, C.loc[a, b])); found = True
if not found:
    say("    none — the features are less redundant than their construction suggests")

say("")
say("2. PCA — how many dimensions are really there?")
p = PCA().fit(X)
cum = np.cumsum(p.explained_variance_ratio_)
for i in range(len(F)):
    say("  PC%-2d %6.1f%%   cumulative %5.1f%%"
        % (i + 1, 100 * p.explained_variance_ratio_[i], 100 * cum[i]))
n80 = int(np.searchsorted(cum, 0.80) + 1)
say("  %d of %d components carry 80%% of the variance" % (n80, len(F)))
say("")
for k in range(min(2, len(F))):
    say("  PC%d loadings:" % (k + 1))
    l = pd.Series(p.components_[k], index=F).sort_values()
    for f, v in l.items():
        say("    %-16s %+6.3f %s" % (f, v, "#" * int(abs(v) * 30)))

say("")
say("3. PARTIAL EFFECTS — does each feature survive holding the others constant?")
say("   This is the direct test. Univariate effects can all be the same effect seen")
say("   through different variables; a multivariate fit separates them.")
say("")
uni, multi = {}, {}
y = d.labelled.values.astype(int)
for f in F:
    m = LogisticRegression(max_iter=2000).fit(Xd[[f]], y)
    uni[f] = float(m.coef_[0][0])
mm = LogisticRegression(max_iter=2000).fit(Xd[F], y)
for f, c in zip(F, mm.coef_[0]):
    multi[f] = float(c)
say("  %-16s %12s %12s %10s" % ("feature", "alone", "adjusted", "retained"))
for f in F:
    keep = "yes" if abs(multi[f]) > 0.3 * abs(uni[f]) and np.sign(multi[f]) == np.sign(uni[f]) else "NO"
    say("  %-16s %12.3f %12.3f %10s" % (f[:16], uni[f], multi[f], keep))
say("")
say("  'retained' = the adjusted coefficient keeps its sign and at least a third of its")
say("  magnitude. A feature that collapses when the others are present was not an")
say("  independent effect.")

say("")
say("4. WITHIN-PROTEIN CHECK on the adjusted model")
say("   The comparison that matters is within protein, so the multivariate fit is repeated")
say("   with protein-level differencing on the strongest surviving feature set.")
keep = [f for f in F if abs(multi[f]) > 0.3 * abs(uni[f]) and np.sign(multi[f]) == np.sign(uni[f])]
say("   features retained: %s" % (", ".join(keep) if keep else "none"))
if keep:
    diffs = {f: [] for f in keep}
    for _, g in d.groupby("uniprot"):
        a = g[g.labelled]
        b = g[~g.labelled]
        if len(a) >= 2 and len(b) >= 2:
            for f in keep:
                diffs[f].append(a[f].astype(float).mean() - b[f].astype(float).mean())
    say("   %-16s %10s %8s %10s" % ("feature", "proteins", "d", "p"))
    for f in keep:
        v = np.array(diffs[f])
        if len(v) < 30:
            continue
        t, pv = sstats.ttest_1samp(v, 0)
        say("   %-16s %10d %8.3f %10.3g" % (f[:16], len(v), v.mean() / v.std(), pv))

pd.DataFrame(dict(feature=F, uni_coef=[uni[f] for f in F],
                  multi_coef=[multi[f] for f in F],
                  pc1_loading=p.components_[0])).to_csv("step12_reactivity_dims.csv", index=False)
open("step12_report.txt", "w").write("\n".join(log) + "\n")
say("")
say("wrote step12_reactivity_dims.csv and step12_report.txt")
say("")
say("READ: if one feature retains its effect and the others collapse, reactivity is governed")
say("by ONE property measured several ways — report it as one, not four.")
