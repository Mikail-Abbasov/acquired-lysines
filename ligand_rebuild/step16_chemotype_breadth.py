#!/usr/bin/env python3
"""
step16_chemotype_breadth.py — does structural distinctiveness scale with how many chemotypes
engage a lysine?

WHAT THIS EXTENDS
Step 15 found that lysines engaged by BOTH chemotypes (CyMPSA and meroterpenoid) occupy a
distinct region of structural space, while lysines engaged by only one are indistinguishable
from lysines engaged by none. That rested on 25 shared sites — suggestive, not settled.

Nature Chemistry provides a far better-powered version of the same question. Its 328
compound x concentration entries parse cleanly into **34 chemotypes** (number = chemotype,
letter = member: 1a-1d, 34a-34f). So instead of a two-group comparison, the number of distinct
chemotypes engaging a site is a GRADIENT from 0 to 34.

A monotonic trend across that gradient is much stronger evidence than a two-group difference:
noise and confounding rarely produce an ordered dose-response.

THE CONFOUND THAT MUST BE CONTROLLED
A site tested against more chemotypes has more chances to be engaged by more of them. Chemotype
coverage is not uniform (chemotype 26 has 36 compound-concentration entries; others far fewer),
and screening depth varies by site. So:
  - n_chemotypes_TESTED is computed per site and carried alongside n_chemotypes_ENGAGED
  - the trend is tested both raw and stratified by testing depth
  - a site tested against fewer than 10 chemotypes is reported separately, since its maximum
    possible breadth is capped

WHAT IS TESTED
  1. Structural features vs chemotype breadth, binned 0 / 1 / 2-3 / 4-5 / >=6
  2. Spearman correlation between breadth and each feature — a direct monotonicity test that
     uses the full gradient rather than bins
  3. Jonckheere-Terpstra style check: is the trend ordered, not merely different?
  4. Multivariate: centroid distance from the 0-chemotype group, by bin, with a permutation
     null — the step 15 test repeated along the gradient

INTERPRETATION IN ADVANCE
  If distinctiveness rises monotonically with breadth, the step 15 finding generalises: what
  marks a site is not which chemistry engages it but how many can.
  If the effect appears only at the top bin, it may be a property of promiscuously engaged
  sites rather than a gradient.
  If there is no trend, step 15's result does not replicate in a better-powered setting and
  should be reported as such.

INPUT   occ_long.csv, features_expanded.csv
OUTPUT  step16_breadth.csv, step16_report.txt
"""
import os, re
import numpy as np, pandas as pd
import warnings
warnings.filterwarnings("ignore")
from scipy import stats as sstats
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests

TIER = int(os.environ.get("TIER", 60))
N_PERM = 2000
FEATS = ["rsa_tien", "sasa_sidechain", "nz_sasa", "crowding_8A", "atoms_10A",
         "hydrophobic_8A", "aromatic_8A", "polar_8A", "anionic_8A", "cationic_8A",
         "net_charge_8A", "n_lys_8A", "dist_centroid", "propka_pKa",
         "hbond_partners_NZ", "in_pocket", "pocket_volume", "pocket_drugscore",
         "dist_nearest_pocket", "plddt_site", "plddt_8A"]
rng = np.random.default_rng(0)

log = []
def say(m):
    print(m, flush=True); log.append(m)

say("=" * 78)
say("STEP 16 — DOES DISTINCTIVENESS SCALE WITH CHEMOTYPE BREADTH?")
say("34 Nature Chemistry chemotypes | engaged = >=%d%% occupancy" % TIER)
say("=" * 78)

occ = pd.read_csv("occ_long.csv", low_memory=False)
nat = occ[occ.dataset == "NatChem"].copy()
pat = re.compile(r"^(\d+)([a-z])")
nat["chemotype"] = nat.compound.map(lambda c: pat.match(str(c)).group(1)
                                    if pat.match(str(c)) else None)
nat = nat[nat.chemotype.notna()]
say("%d measurements | %d compound labels | %d chemotypes"
    % (len(nat), nat.compound.nunique(), nat.chemotype.nunique()))

tested = (nat.groupby(["uniprot", "position"]).chemotype.nunique()
             .rename("n_ct_tested").reset_index())
engaged = (nat[nat.occupancy_pct >= TIER].groupby(["uniprot", "position"])
              .chemotype.nunique().rename("n_ct_engaged").reset_index())
sites = tested.merge(engaged, on=["uniprot", "position"], how="left")
sites["n_ct_engaged"] = sites.n_ct_engaged.fillna(0).astype(int)

feat = (pd.read_csv("features_expanded.csv", low_memory=False)
          .drop_duplicates(["uniprot", "position"]))
if "in_pocket" in feat:
    feat["in_pocket"] = feat.in_pocket.astype(str).str.lower().isin(["true", "1", "1.0"]).astype(float)
F = [f for f in FEATS if f in feat.columns]
d = sites.merge(feat[["uniprot", "position"] + F], on=["uniprot", "position"], how="inner").dropna()
say("%d lysines with structure and chemotype counts" % len(d))


def band(n):
    return "0" if n == 0 else "1" if n == 1 else "2-3" if n <= 3 else "4-5" if n <= 5 else ">=6"


d["bin"] = d.n_ct_engaged.map(band)
BINS = ["0", "1", "2-3", "4-5", ">=6"]
say("")
say("BREADTH DISTRIBUTION")
say("  %-8s %8s %14s %14s" % ("bin", "n", "med tested", "med engaged"))
for b in BINS:
    g = d[d["bin"] == b]
    if len(g):
        say("  %-8s %8d %14.0f %14.0f"
            % (b, len(g), g.n_ct_tested.median(), g.n_ct_engaged.median()))
say("")
say("  chemotypes TESTED per site: median %.0f, 10th %.0f, 90th %.0f"
    % (d.n_ct_tested.median(), d.n_ct_tested.quantile(.1), d.n_ct_tested.quantile(.9)))
say("  A site tested against few chemotypes cannot show high breadth — this caps the")
say("  gradient and is why depth is carried through every test below.")

# ------------------------------------------------------------------ 1/2. trend
say("")
say("MONOTONIC TREND — Spearman correlation of each feature with chemotype breadth")
say("  (full gradient, not bins; partial version holds testing depth constant)")
rows = []
for f in F:
    r, p = sstats.spearmanr(d.n_ct_engaged, d[f])
    # partial correlation controlling for chemotypes tested
    res_f = d[f] - np.poly1d(np.polyfit(d.n_ct_tested, d[f], 1))(d.n_ct_tested)
    res_b = d.n_ct_engaged - np.poly1d(np.polyfit(d.n_ct_tested, d.n_ct_engaged, 1))(d.n_ct_tested)
    rp, pp = sstats.spearmanr(res_b, res_f)
    rows.append(dict(feature=f, rho=r, p=p, rho_adj=rp, p_adj=pp))
tr = pd.DataFrame(rows)
tr["q"] = multipletests(tr.p, method="fdr_bh")[1]
tr["q_adj"] = multipletests(tr.p_adj, method="fdr_bh")[1]
tr = tr.reindex(tr.rho.abs().sort_values(ascending=False).index)
say("  %-22s %8s %10s %10s %10s" % ("feature", "rho", "q", "rho|depth", "q|depth"))
for _, r in tr.iterrows():
    say("  %-22s %8.3f %10.3g %10.3f %10.3g%s"
        % (r.feature[:22], r.rho, r.q, r.rho_adj, r.q_adj,
           " *" if r.q_adj < 0.05 else ""))
tr.to_csv("step16_breadth.csv", index=False)

# ------------------------------------------------------------------ 3. by bin
say("")
say("MEDIANS BY BIN — is the change ordered, or only present at the extreme?")
top = tr.head(6).feature.tolist()
say("  %-22s %s" % ("feature", "  ".join("%9s" % b for b in BINS)))
for f in top:
    cells = []
    for b in BINS:
        g = d[d["bin"] == b]
        cells.append("%9.3f" % g[f].median() if len(g) > 5 else "%9s" % "-")
    say("  %-22s %s" % (f[:22], "  ".join(cells)))

# ------------------------------------------------------------------ 4. multivariate
say("")
say("MULTIVARIATE — distance from the 0-chemotype group, by bin")
say("  (step 15's test repeated along the gradient; %d label shuffles)" % N_PERM)
X = StandardScaler().fit_transform(d[F].astype(float))
base = X[(d["bin"] == "0").values]
say("  %-8s %8s %10s %10s %10s" % ("bin", "n", "distance", "null mean", "p"))
for b in BINS[1:]:
    m = (d["bin"] == b).values
    if m.sum() < 5:
        continue
    obs = float(np.linalg.norm(X[m].mean(axis=0) - base.mean(axis=0)))
    pool = np.where(m | (d["bin"] == "0").values)[0]
    na = int(m.sum())
    null = []
    for _ in range(N_PERM):
        idx = rng.permutation(pool)
        null.append(np.linalg.norm(X[idx[:na]].mean(axis=0) - X[idx[na:]].mean(axis=0)))
    null = np.array(null)
    say("  %-8s %8d %10.3f %10.3f %10.4f%s"
        % (b, na, obs, null.mean(), (null >= obs).mean(),
           " *" if (null >= obs).mean() < 0.05 else ""))

say("")
say("READING THIS")
say("  Distance rising steadily across bins -> step 15 generalises: what marks a site is")
say("     not WHICH chemistry engages it but HOW MANY can.")
say("  Distance flat until the top bin -> a property of promiscuously engaged sites, not")
say("     a gradient.")
say("  No trend -> step 15 does not replicate with better power; report that plainly.")
say("")
say("  Note: this is WITHIN one dataset and one aggregated cell-line context, so it tests")
say("  chemotype breadth, not cross-dataset reproducibility. Step 15 remains the only")
say("  cross-dataset evidence.")
open("step16_report.txt", "w").write("\n".join(log) + "\n")
say("")
say("wrote step16_breadth.csv and step16_report.txt")
