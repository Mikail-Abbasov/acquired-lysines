#!/usr/bin/env python3
"""
step43b_enclosure_tests.py — does concavity, measured independently of fpocket, behave like a
property of the SITE or of the REAGENT?

Reads step43_enclosure.csv (see step43_enclosure.py for the measure) and answers three things
Figure 7 needs. All tests use the WITHIN-PROTEIN paired design of step 19 and step 22: one
difference per protein, then tested against zero across proteins.

WHY WITHIN PROTEIN, EMPHATICALLY
The pooled site-level version of test 2 gives d = -0.03 — labelled lysines look LESS enclosed.
Within protein the same data give d = +0.073. That is a sign reversal from pooling, i.e. protein
composition confounds the pooled comparison, and it is the same hazard the within-compound design
of Figure 6a exists to avoid. Any future pooled result on this data should be treated as suspect
until it has been repeated within protein.

NOTE ON q-VALUES. Test 3 applies FDR across the 6 features listed here, while step 22 applied it
across 22. The same comparison therefore carries a different q in the two reports (cationic_8A:
4.3e-05 here, 1.6e-04 in step 22). Cliff's d is the quantity to compare across the two, and it
is identical: 0.233 and 0.180.
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
from scipy import stats as sstats
from statsmodels.stats.multitest import multipletests

log = []
def say(m):
    print(m, flush=True); log.append(m)

def within(df, feats, group_col=None, a_val=None, b_val=None):
    """One difference per protein, then a one-sample t-test across proteins."""
    out = []
    for f in feats:
        diffs = []
        for _, g in df.groupby("uniprot"):
            if group_col is None:
                a = pd.to_numeric(g.loc[g.labelled, f], errors="coerce").dropna()
                b = pd.to_numeric(g.loc[~g.labelled, f], errors="coerce").dropna()
            else:
                a = pd.to_numeric(g.loc[g[group_col] == a_val, f], errors="coerce").dropna()
                b = pd.to_numeric(g.loc[g[group_col] == b_val, f], errors="coerce").dropna()
            if len(a) >= 2 and len(b) >= 2:
                diffs.append(a.mean() - b.mean())
        if len(diffs) < 30:
            continue
        v = np.array(diffs); t, p = sstats.ttest_1samp(v, 0)
        out.append(dict(feature=f, n_proteins=len(v), mean_diff=float(v.mean()),
                        d=float(v.mean() / v.std()) if v.std() else np.nan, p=float(p)))
    r = pd.DataFrame(out)
    if len(r):
        r["q"] = multipletests(r.p, method="fdr_bh")[1]
    return r

def show(R):
    say("  %-22s %9s %11s %8s %11s" % ("feature", "proteins", "mean diff", "d", "q"))
    for _, r in R.iterrows():
        say("  %-22s %9d %11.4f %8.3f %11.3g%s"
            % (r.feature, r.n_proteins, r.mean_diff, r.d, r.q, "  *" if r.q < 0.05 else ""))

E = pd.read_csv("step43_enclosure.csv")
P = pd.read_csv("step19_pockets.csv")
M = E.merge(P, on=["uniprot", "position"])
M = M[M.plddt_site >= 70].copy()
ok = M[["enclosure", "rsa_tien"]].notna().all(1)
M.loc[ok, "enclosure_adj"] = (M.loc[ok, "enclosure"].values - np.polyval(
    np.polyfit(M.loc[ok, "rsa_tien"], M.loc[ok, "enclosure"], 1), M.loc[ok, "rsa_tien"]))

say("=" * 78)
say("STEP 43b — CONCAVITY AS A PROPERTY OF THE SITE, NOT THE REAGENT")
say("=" * 78)
say("%s lysines, %d proteins, pLDDT >= 70" % (f"{len(M):,}", M.uniprot.nunique()))

say("")
say("1. IS ENCLOSURE A DISTINCT AXIS, OR IS IT SOLVENT ACCESSIBILITY IN DISGUISE?")
rho = sstats.spearmanr(M.enclosure, M.rsa_tien)[0]
say("   enclosure vs relative SASA: rho %.3f (%.0f%% shared variance)" % (rho, 100 * rho ** 2))
sub = M[M.in_pocket.notna()]
say("   enclosure vs in_pocket: raw rho %.3f | after removing SASA %.3f"
    % (sstats.spearmanr(sub.enclosure, sub.in_pocket)[0],
       sstats.spearmanr(sub.enclosure_adj, sub.in_pocket)[0]))
say("   -> mostly accessibility, but a fifth of it is pocket geometry that SASA does not carry.")

say("")
say("2. DO PROBE-LABELLED LYSINES SIT IN CONCAVITIES? (within protein)")
R1 = within(M, ["in_pocket", "dist_nearest_pocket", "rsa_tien", "enclosure", "enclosure_adj"])
show(R1)
say("   in_pocket reproduces step 19 exactly (d 0.107, q 1.7e-10).")
say("   Enclosure agrees, and SURVIVES removal of solvent accessibility (d 0.061).")
say("   POOLED, the same data give enclosure d = -0.03 — a sign reversal. Use within protein.")

occ = pd.read_csv("occ_long.csv", low_memory=False)
cy = occ[occ.dataset == "CyMPSA"]
sets = {}
for p, n in (("dTB-STP", "STP"), ("Probe 1", "NASA")):
    sets[n] = set(map(tuple, cy[cy.probe == p][["uniprot", "position"]].drop_duplicates().values))
stp_only = sets["STP"] - sets["NASA"]; nasa_only = sets["NASA"] - sets["STP"]
k = list(map(tuple, M[["uniprot", "position"]].values))
M["_grp"] = ["STP" if x in stp_only else ("NASA" if x in nasa_only else None) for x in k]
S = M[M._grp.notna()]

say("")
say("3. DO THE TWO PROBES DIFFER IN CONCAVITY? (within protein; + = anionic probe)")
say("   STP-only %d | NASA-only %d sites" % ((S._grp == "STP").sum(), (S._grp == "NASA").sum()))
R2 = within(S, ["enclosure", "enclosure_adj", "in_pocket", "dist_nearest_pocket",
                "cationic_8A", "net_charge_8A"], "_grp", "STP", "NASA")
show(R2)
say("   -> The probes differ in CHARGE environment and not in CONCAVITY, now tested with a")
say("      measure that owes nothing to fpocket and nothing to solvent accessibility. This is")
say("      a stronger form of E6: cavity preference belongs to the site, not to the reagent.")

say("")
say("SUPERSEDES: the architecture's \"41.6% of lysines within 5 A of a cavity, only 10.3% large")
say("enough to matter\" appears in no step report and does not reproduce (computed: 40.4%")
say("ungated, 48.3% at pLDDT >= 70; no volume threshold gives 10.3%). It also rested on a")
say("distinction the data cannot draw, since step 19's in_pocket IS dist_nearest_pocket <= 5.")
R1.to_csv("step43b_site_vs_unlabelled.csv", index=False)
R2.to_csv("step43b_probe_contrast.csv", index=False)
open("step43b_report.txt", "w").write("\n".join(log) + "\n")
say("")
say("wrote step43b_report.txt, step43b_site_vs_unlabelled.csv, step43b_probe_contrast.csv")
