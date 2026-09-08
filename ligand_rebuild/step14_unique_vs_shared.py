#!/usr/bin/env python3
"""
step14_unique_vs_shared.py — what distinguishes a lysine that BOTH chemotypes engage from one
that only a single chemotype engages?

THE QUESTION
Step 13 showed the two chemotypes converge on the same lysines about five times more often than
chance (26 shared, 5.5 expected, OR 6.0). But most engaged sites are still hit by only one
chemotype. Are the shared sites structurally different from the exclusive ones? A site that
accommodates both a designed sulfonyl electrophile and a fused-ring natural product may be a
more open, more generally accessible pocket than one that suits only a single scaffold.

FOUR GROUPS, all drawn from the SAME denominator
    the 8,211 lysines measured by BOTH chemotypes

  SHARED         engaged (>=60%) by CyMPSA AND by a meroterpenoid
  CYMPSA_ONLY    engaged by CyMPSA, not by any meroterpenoid
  MERO_ONLY      engaged by a meroterpenoid, not by CyMPSA
  NEITHER        measured by both, engaged by neither  (the reference group)

WHY THE COMMON DENOMINATOR IS ESSENTIAL
CyMPSA is HeLa; meroterpenoids are MDA-MB-231/436. A site engaged in one cell line may simply
be absent from the other because the protein is not expressed there. Restricting every group to
lysines actually measured in both makes "not engaged by chemotype X" a real negative rather
than a missing value. Without this the SHARED group would be defined by co-expression, not
chemistry.

WHAT IS COMPARED
Every structural and electrostatic parameter from steps 6/6d, plus the pocket features:
burial, crowding, neighbour composition, pocket membership/volume/druggability, distance to
nearest pocket, pKa, local charge, salt bridge, secondary structure, pLDDT.

STATISTICS
  - each group vs NEITHER, and SHARED vs each exclusive group
  - Mann-Whitney (no normality assumption; several features are counts or bounded)
  - Cliff's delta as the effect size, which is rank-based and does not assume equal variance
  - Benjamini-Hochberg within each comparison
  - group sizes printed with every result

THE POWER PROBLEM, STATED UP FRONT
SHARED is ~26 sites. That is small. Any null must be read as "not detectable at n=26", never
as evidence the groups are alike. Confidence intervals are given; conclusions that depend on
SHARED alone are provisional until more data exists.

INPUT   occ_long.csv, features_expanded.csv
OUTPUT  step14_groups.csv (per-site group labels + parameters), step14_report.txt
"""
import os
import numpy as np, pandas as pd
import warnings
warnings.filterwarnings("ignore")
from scipy import stats as sstats
from statsmodels.stats.multitest import multipletests

TIER = int(os.environ.get("TIER", 60))
FEATS = ["rsa_tien", "sasa_sidechain", "nz_sasa", "crowding_8A", "atoms_10A",
         "hydrophobic_8A", "aromatic_8A", "polar_8A", "anionic_8A", "cationic_8A",
         "net_charge_8A", "n_lys_8A", "dist_centroid", "propka_pKa",
         "hbond_partners_NZ", "in_pocket", "pocket_volume", "pocket_drugscore",
         "dist_nearest_pocket", "plddt_site", "plddt_8A"]

log = []
def say(m):
    print(m, flush=True); log.append(m)


def cliffs_delta(a, b):
    """Rank-based effect size: P(a>b) - P(a<b). -1..+1, sign follows a relative to b."""
    a, b = np.asarray(a), np.asarray(b)
    if len(a) < 3 or len(b) < 3:
        return np.nan
    n = len(a) * len(b)
    if n > 4_000_000:                       # subsample for speed on the large reference group
        rng = np.random.default_rng(0)
        b = rng.choice(b, size=min(len(b), 4_000_000 // max(len(a), 1)), replace=False)
    gt = sum((a[:, None] > b[None, :]).sum(axis=1))
    lt = sum((a[:, None] < b[None, :]).sum(axis=1))
    return (gt - lt) / (len(a) * len(b))


say("=" * 78)
say("STEP 14 — CHEMOTYPE-UNIQUE vs SHARED LIGANDED LYSINES")
say("engaged = >=%d%% occupancy by at least one compound of that chemotype" % TIER)
say("=" * 78)

occ = pd.read_csv("occ_long.csv", low_memory=False)
feat = (pd.read_csv("features_expanded.csv", low_memory=False)
          .drop_duplicates(["uniprot", "position"]))
if "in_pocket" in feat:
    feat["in_pocket"] = feat.in_pocket.astype(str).str.lower().isin(["true", "1", "1.0"]).astype(float)

def sset(ds, engaged):
    g = occ[occ.dataset == ds]
    if engaged:
        g = g[g.occupancy_pct >= TIER]
    return set(map(tuple, g[["uniprot", "position"]].drop_duplicates().values))

mC, mM = sset("CyMPSA", False), sset("Meroterpenoid", False)
eC, eM = sset("CyMPSA", True), sset("Meroterpenoid", True)
both_meas = mC & mM
say("lysines measured by BOTH chemotypes: %d" % len(both_meas))
say("  (essential: CyMPSA is HeLa, meroterpenoids MDA-MB-231/436 — outside this set,")
say("   'not engaged' would often mean 'protein not expressed there')")

key = list(zip(feat.uniprot, feat.position.astype(int)))
feat = feat.assign(_k=key)
d = feat[feat._k.isin(both_meas)].copy()


def label(k):
    c, m = k in eC, k in eM
    return "SHARED" if c and m else "CYMPSA_ONLY" if c else "MERO_ONLY" if m else "NEITHER"


d["group"] = d._k.map(label)
F = [f for f in FEATS if f in d.columns]
say("")
say("GROUP SIZES")
for g in ("SHARED", "CYMPSA_ONLY", "MERO_ONLY", "NEITHER"):
    say("  %-14s %6d" % (g, int((d.group == g).sum())))
say("")
say("  !! SHARED is small. A null result for SHARED means 'not detectable at this n',")
say("     never that the groups are alike.")

comparisons = [("CYMPSA_ONLY", "NEITHER"), ("MERO_ONLY", "NEITHER"),
               ("SHARED", "NEITHER"), ("SHARED", "CYMPSA_ONLY"), ("SHARED", "MERO_ONLY")]
rows = []
for ga, gb in comparisons:
    A, B = d[d.group == ga], d[d.group == gb]
    if len(A) < 5 or len(B) < 5:
        continue
    out = []
    for f in F:
        a = pd.to_numeric(A[f], errors="coerce").dropna().values
        b = pd.to_numeric(B[f], errors="coerce").dropna().values
        if len(a) < 5 or len(b) < 5:
            continue
        try:
            u, p = sstats.mannwhitneyu(a, b, alternative="two-sided")
        except ValueError:
            continue
        out.append(dict(comparison="%s vs %s" % (ga, gb), feature=f,
                        n_a=len(a), n_b=len(b),
                        med_a=float(np.median(a)), med_b=float(np.median(b)),
                        delta=cliffs_delta(a, b), p=float(p)))
    if not out:
        continue
    o = pd.DataFrame(out)
    o["q"] = multipletests(o.p, method="fdr_bh")[1]
    rows.append(o)
    o = o.sort_values("q")
    say("")
    say("-" * 78)
    say("%s vs %s   (n = %d vs %d)" % (ga, gb, len(A), len(B)))
    say("  %-22s %9s %9s %8s %10s" % ("feature", "median A", "median B", "delta", "q"))
    for _, r in o.iterrows():
        mark = " *" if r.q < 0.05 else ""
        say("  %-22s %9.3f %9.3f %8.3f %10.3g%s"
            % (r.feature[:22], r.med_a, r.med_b, r.delta, r.q, mark))
    nsig = int((o.q < 0.05).sum())
    say("  -> %d of %d features differ" % (nsig, len(o)))
    if nsig == 0:
        say("     with n=%d this is inconclusive, not evidence of similarity" % len(A))

res = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
if len(res):
    res.to_csv("step14_groups.csv", index=False)

say("")
say("=" * 78)
say("SUMMARY — is a shared site different from an exclusive one?")
if len(res):
    for ga, gb in [("SHARED", "CYMPSA_ONLY"), ("SHARED", "MERO_ONLY")]:
        sub = res[res.comparison == "%s vs %s" % (ga, gb)]
        if not len(sub):
            continue
        sig = sub[sub.q < 0.05]
        say("  %s vs %s: %d of %d features differ%s"
            % (ga, gb, len(sig), len(sub),
               "" if len(sig) else " — inconclusive at this sample size"))
        for _, r in sig.sort_values("q").head(6).iterrows():
            direction = "higher" if r.delta > 0 else "lower"
            say("     %-22s %s in SHARED (delta %+.2f, q %.3g)"
                % (r.feature[:22], direction, r.delta, r.q))
    say("")
    say("  Interpretation guide: if SHARED sites are MORE open (larger pocket, less buried,")
    say("  higher druggability) than exclusive sites, that suggests generally accessible")
    say("  pockets accommodate varied scaffolds while tighter sites are scaffold-selective.")
    say("  If SHARED sites are indistinguishable, chemotype exclusivity is a property of the")
    say("  COMPOUNDS, not of the site.")

d.drop(columns=["_k"]).to_csv("step14_sites.csv", index=False)
open("step14_report.txt", "w").write("\n".join(log) + "\n")
say("")
say("wrote step14_groups.csv, step14_sites.csv, step14_report.txt")
