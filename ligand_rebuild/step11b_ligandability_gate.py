#!/usr/bin/env python3
"""
step11b_ligandability_gate.py — the missing half of step 11.

WHY THIS EXISTS
Step 11 gates Comparison A (reactivity) on pLDDT and shows that its shape effects collapse
while its charge effects hold. Figure 4 then compares those GATED reactivity effects against
UNGATED ligandability effects, because a gated ligandability contrast had never been run.
That asymmetry is a real weakness: if AlphaFold confidence distorts one contrast it can
distort the other, and an ungated comparator cannot rule that out.

This script runs Comparison B (engaged vs all probe-labelled, at >=60% occupancy) through the
identical gate ladder, with the identical paired within-protein test, so that both sides of
Figure 4 can be read at the same confidence.

THE SAME TWO CAUTIONS AS STEP 11 APPLY, AND ARE HANDLED THE SAME WAY

  1. THE GATE IS APPLIED TO BOTH GROUPS. Engaged and unengaged labelled sites are filtered at
     the same pLDDT. Gating one side would manufacture a difference.

  2. RANGE RESTRICTION IS REPORTED, NOT IGNORED. The SD of every feature is printed at every
     gate, so a fall in d can be checked against a fall in spread before it is called
     artefact removal.

ONE CAUTION THAT IS SHARPER HERE THAN IN STEP 11
Comparison B starts from 3,047 engaged sites, not 30,805 labelled ones, and the paired test
needs at least two sites per group per protein. Tightening the gate therefore costs proteins
much faster than it does in step 11. n_proteins is reported at every gate and every feature;
where it falls below 30 the feature is dropped rather than reported on thin evidence.

INPUT   features_expanded.csv, occ_long.csv
OUTPUT  step11b_ligandability_gates.csv, step11b_report.txt
"""
import numpy as np, pandas as pd
import warnings
warnings.filterwarnings("ignore")
from scipy import stats as sstats
from statsmodels.stats.multitest import multipletests

GATES = [(0, "none"), (50, ">=50 low-conf cut"), (70, ">=70 confident"), (90, ">=90 very high")]
# salt_bridge added 6 September 2026: step 11 measures it for reactivity and this list had
# simply omitted it, so the two ladders were not the identical feature set this script exists
# to provide. It is a charge feature, so the paired set becomes 15 (6 charge / 9 shape).
FEATS = ["propka_pKa", "net_charge_8A", "anionic_8A", "cationic_8A", "n_lys_8A", "salt_bridge",
         "crowding_8A", "atoms_10A", "hydrophobic_8A", "aromatic_8A", "polar_8A",
         "rsa_tien", "sasa_sidechain", "nz_sasa", "dist_centroid"]
TIER = 60
MIN_PER_PROT = 2
MIN_PROTEINS = 30

log = []
def say(m):
    print(m, flush=True); log.append(m)

say("=" * 78)
say("STEP 11b — pLDDT GATE ON COMPARISON B (LIGANDABILITY)")
say("the symmetric partner of step 11, so Figure 4 can gate both contrasts")
say("=" * 78)

feat = pd.read_csv("features_expanded.csv", low_memory=False)
occ  = pd.read_csv("occ_long.csv", low_memory=False)
sites = feat.drop_duplicates(["uniprot", "position"]).copy()
# salt_bridge lives in step9_all_lysines.csv, the same table step 11 measures reactivity on,
# and not in features_expanded.csv. Merged on (uniprot, position) so both ladders carry the
# identical feature set rather than differing by one charge term.
_sb = pd.read_csv("step9_all_lysines.csv", low_memory=False,
                  usecols=["uniprot", "position", "salt_bridge"]).drop_duplicates(
                      ["uniprot", "position"])
sites = sites.merge(_sb, on=["uniprot", "position"], how="left")
eng = set(map(tuple, occ[occ.occupancy_pct >= TIER][["uniprot", "position"]]
              .drop_duplicates().values))
sites["engaged"] = [(u, int(p)) in eng for u, p in zip(sites.uniprot, sites.position)]
sites["plddt_site"] = pd.to_numeric(sites.plddt_site, errors="coerce")

say("%d distinct probe-labelled lysines on %d proteins | engaged at >=%d%%: %d (%.2f%%)"
    % (len(sites), sites.uniprot.nunique(), TIER, sites.engaged.sum(),
       100 * sites.engaged.mean()))
say("pLDDT: engaged median %.1f | unengaged median %.1f"
    % (sites[sites.engaged].plddt_site.median(), sites[~sites.engaged].plddt_site.median()))
say("salt_bridge merged from step9_all_lysines.csv: %d of %d sites have a value (%.1f%%)"
    % (sites.salt_bridge.notna().sum(), len(sites),
       100 * sites.salt_bridge.notna().mean()))
say("")


def within_protein(df, feats):
    """One difference per protein (engaged minus unengaged), then test against zero —
    the same estimator step 9 and step 11 use, so the numbers are comparable."""
    out = []
    for f in feats:
        diffs, sds = [], []
        for _, g in df.groupby("uniprot"):
            a = pd.to_numeric(g.loc[g.engaged, f], errors="coerce").dropna()
            b = pd.to_numeric(g.loc[~g.engaged, f], errors="coerce").dropna()
            if len(a) >= MIN_PER_PROT and len(b) >= MIN_PER_PROT:
                diffs.append(a.mean() - b.mean())
        if len(diffs) < MIN_PROTEINS:
            out.append(dict(feature=f, n_proteins=len(diffs), mean_diff=np.nan,
                            d=np.nan, p=np.nan, sd_all=np.nan, dropped=True))
            continue
        v = np.array(diffs)
        t, p = sstats.ttest_1samp(v, 0)
        out.append(dict(feature=f, n_proteins=len(v), mean_diff=float(v.mean()),
                        d=float(v.mean() / v.std()) if v.std() else np.nan,
                        p=float(p), sd_all=float(pd.to_numeric(df[f], errors="coerce").std()),
                        dropped=False))
    return pd.DataFrame(out)


rows = []
for g, label in GATES:
    sub = sites[sites.plddt_site >= g] if g else sites
    n_e, n_u = int(sub.engaged.sum()), int((~sub.engaged).sum())
    say("GATE %s" % label)
    say("  engaged %d (%d%% of all engaged) | unengaged %d | proteins %d"
        % (n_e, round(100 * n_e / max(sites.engaged.sum(), 1)), n_u, sub.uniprot.nunique()))
    r = within_protein(sub, FEATS)
    ok = r[~r.dropped]
    if len(ok):
        ok = ok.copy(); ok["q"] = multipletests(ok.p, method="fdr_bh")[1]
        r = r.merge(ok[["feature", "q"]], on="feature", how="left")
    else:
        r["q"] = np.nan
    say("  %-18s %8s %10s %8s %10s %8s" % ("feature", "proteins", "mean diff", "d", "q", "SD"))
    for _, x in r.sort_values("d", key=abs, ascending=False).iterrows():
        if x.dropped:
            say("  %-18s %8d   dropped — fewer than %d proteins with both classes"
                % (x.feature[:18], x.n_proteins, MIN_PROTEINS))
        else:
            say("  %-18s %8d %10.4f %8.3f %10.3g %8.3f%s"
                % (x.feature[:18], x.n_proteins, x.mean_diff, x.d, x.q, x.sd_all,
                   " *" if x.q < 0.05 else ""))
    r["gate"] = g
    rows.append(r)
    say("")

out = pd.concat(rows, ignore_index=True)
out.to_csv("step11b_ligandability_gates.csv", index=False)

say("=" * 78)
say("WHAT MOVED")
g0 = out[out.gate == 0].set_index("feature")
gL = out[out.gate == GATES[-1][0]].set_index("feature")
common = [f for f in FEATS if f in g0.index and f in gL.index
          and not g0.loc[f, "dropped"] and not gL.loc[f, "dropped"]]
say("  %-18s %10s %10s %10s %10s" % ("feature", "no gate", ">=90", "change", "SD ratio"))
for f in common:
    say("  %-18s %10.3f %10.3f %10.3f %10.2f"
        % (f[:18], g0.loc[f, "d"], gL.loc[f, "d"], gL.loc[f, "d"] - g0.loc[f, "d"],
           gL.loc[f, "sd_all"] / g0.loc[f, "sd_all"] if g0.loc[f, "sd_all"] else np.nan))
say("")
say("  %d of %d features measurable at both ends of the ladder" % (len(common), len(FEATS)))
if common:
    weaker = sum(abs(gL.loc[f, "d"]) < abs(g0.loc[f, "d"]) for f in common)
    say("  %d weaker under the gate, %d stronger" % (weaker, len(common) - weaker))

open("step11b_report.txt", "w").write("\n".join(log) + "\n")
say("")
say("wrote step11b_ligandability_gates.csv, step11b_report.txt")
