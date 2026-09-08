#!/usr/bin/env python3
"""
step11_plddt_gate.py — is Comparison A real chemistry, or unreliable geometry?

THE PROBLEM THIS TESTS
Comparison A (step 9) compares probe-labelled lysines against every other lysine on the same
protein. The workflow justifies that comparison: probe labels, protein is digested,
streptavidin pulls down labelled peptides, then MS. An unlabelled lysine did not react — it is
not a detection artefact.

But the largest single difference between the groups was pLDDT (d = 0.47), roughly double any
other feature. The chemical reading is not obvious. The alternative is that our FEATURES are
unreliable in disordered regions: a lysine in a floppy loop has no fixed environment, so its
computed burial, crowding, neighbour counts and pKa are meaningless. If labelled lysines are
concentrated in ordered regions, the other effects could partly reflect that.

THE TEST
Restrict both groups to sites above a pLDDT threshold and see what happens to the effects.

  effects WEAKEN as the gate tightens  -> they were driven by unreliable geometry
  effects HOLD or STRENGTHEN           -> they are real, and disorder was adding noise

Gates: none, >=50 (AlphaFold's "low confidence" boundary), >=70 ("confident"), >=90
("very high"). Published thresholds, so the choice is defensible rather than arbitrary.

TWO THINGS THIS SCRIPT IS CAREFUL ABOUT

  1. THE GATE IS APPLIED TO BOTH GROUPS. Labelled and unlabelled sites are filtered at the
     same pLDDT. Gating only one side would compare well-modelled labelled lysines against
     everything else and manufacture a difference.

  2. RANGE RESTRICTION IS REPORTED, NOT IGNORED. pLDDT correlates with the features themselves
     — ordered regions are more buried and more crowded almost by construction. Tightening the
     gate therefore narrows the range of every feature, which shrinks effect sizes for a purely
     statistical reason. The SD of each feature is printed at every gate so that a decline in d
     can be checked against a decline in spread before it is called artefact removal.

ALSO REPORTED
  what each gate discards, on both sides. If >=90 removes most labelled sites, the surviving
  comparison is on a narrow and unusual subset, and that must be stated.

INPUT   step9_all_lysines.csv (every lysine on every relevant protein, labelled or not)
OUTPUT  step11_gates.csv, step11_report.txt
"""
import os
import numpy as np, pandas as pd
import warnings
warnings.filterwarnings("ignore")
from scipy import stats as sstats
from statsmodels.stats.multitest import multipletests

GATES = [(0, "none"), (50, ">=50 low-conf cut"), (70, ">=70 confident"), (90, ">=90 very high")]
FEATS = ["propka_pKa", "salt_bridge", "net_charge_8A", "anionic_8A", "cationic_8A",
         "crowding_8A", "atoms_10A", "hydrophobic_8A", "aromatic_8A", "polar_8A",
         "rsa_tien", "sasa_sidechain", "nz_sasa", "dist_centroid", "n_lys_8A"]
MIN_PER_PROT = 2

log = []
def say(m):
    print(m, flush=True); log.append(m)

say("=" * 78)
say("STEP 11 — pLDDT GATE: is Comparison A chemistry or unreliable geometry?")
say("=" * 78)

allK = pd.read_csv("step9_all_lysines.csv", low_memory=False)
allK["salt_bridge"] = allK.salt_bridge.astype(str).str.lower().isin(["true", "1", "1.0"]).astype(float)
say("%d lysines on %d proteins | labelled %d (%.1f%%)"
    % (len(allK), allK.uniprot.nunique(), allK.labelled.sum(), 100 * allK.labelled.mean()))
say("")
say("pLDDT distribution: labelled median %.1f | unlabelled median %.1f"
    % (allK[allK.labelled].plddt_site.median(), allK[~allK.labelled].plddt_site.median()))


def within_protein(df, feats):
    """One difference per protein (labelled minus unlabelled), then test against zero."""
    out = []
    for f in feats:
        diffs = []
        for _, g in df.groupby("uniprot"):
            a = pd.to_numeric(g.loc[g.labelled, f], errors="coerce").dropna()
            b = pd.to_numeric(g.loc[~g.labelled, f], errors="coerce").dropna()
            if len(a) >= MIN_PER_PROT and len(b) >= MIN_PER_PROT:
                diffs.append(a.mean() - b.mean())
        if len(diffs) < 30:
            continue
        v = np.array(diffs)
        t, p = sstats.ttest_1samp(v, 0)
        out.append(dict(feature=f, n_proteins=len(v), mean_diff=float(v.mean()),
                        d=float(v.mean() / v.std()) if v.std() else np.nan, p=float(p)))
    r = pd.DataFrame(out)
    if len(r):
        r["q"] = multipletests(r.p, method="fdr_bh")[1]
    return r


rows = []
for cut, label in GATES:
    sub = allK if cut == 0 else allK[allK.plddt_site >= cut]
    nl, nu = int(sub.labelled.sum()), int((~sub.labelled).sum())
    say("")
    say("-" * 78)
    say("GATE %s" % label)
    say("  retained: labelled %d (%.0f%% of all labelled) | unlabelled %d (%.0f%%)"
        % (nl, 100 * nl / max(allK.labelled.sum(), 1),
           nu, 100 * nu / max((~allK.labelled).sum(), 1)))
    if nl < 200:
        say("  too few labelled sites; skipped")
        continue
    r = within_protein(sub, [f for f in FEATS if f in sub.columns])
    if not len(r):
        say("  too few proteins with both groups")
        continue
    r["gate"] = cut
    # feature spread at this gate — the range-restriction check
    for f in r.feature:
        r.loc[r.feature == f, "sd_all"] = pd.to_numeric(sub[f], errors="coerce").std()
    rows.append(r)
    r = r.sort_values("q")
    say("  %-18s %8s %10s %8s %9s %8s" % ("feature", "proteins", "mean diff", "d", "q", "SD"))
    for _, x in r.iterrows():
        say("  %-18s %8d %10.4f %8.3f %9.3g %8.3f%s"
            % (x.feature[:18], x.n_proteins, x.mean_diff, x.d, x.q, x.sd_all,
               " *" if x.q < 0.05 else ""))

res = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
if len(res):
    res.to_csv("step11_gates.csv", index=False)
    say("")
    say("=" * 78)
    say("DOSE-RESPONSE ACROSS GATES — effect size (d) and feature spread (SD)")
    say("=" * 78)
    say("%-18s %s" % ("feature", "  ".join("%14s" % g[1].split()[0] for g in GATES)))
    for f in FEATS:
        g = res[res.feature == f]
        if len(g) < 2:
            continue
        cells = []
        for cut, _ in GATES:
            row = g[g.gate == cut]
            cells.append("%6.3f/%-6.2f" % (row.d.iloc[0], row.sd_all.iloc[0])
                         if len(row) else "%13s" % "-")
        say("%-18s %s" % (f[:18], "  ".join(cells)))
    say("")
    say("  each cell is  effect-size / feature-SD  at that gate")
    say("")
    say("HOW TO READ THIS")
    say("  d falls AND SD falls proportionally -> range restriction, not artefact removal.")
    say("     Tightening the gate narrows every feature's spread because ordered regions are")
    say("     more uniformly buried; smaller spread mechanically shrinks d.")
    say("  d falls while SD holds -> the effect really was driven by disordered sites where")
    say("     the geometry is unreliable.")
    say("  d holds or rises -> real chemistry; disorder was adding noise, not signal.")
    say("")
    base = res[res.gate == 0].set_index("feature")
    tight = res[res.gate == 90].set_index("feature")
    if len(tight):
        say("SUMMARY, no gate vs >=90:")
        for f in FEATS:
            if f in base.index and f in tight.index:
                d0, d9 = base.loc[f, "d"], tight.loc[f, "d"]
                s0, s9 = base.loc[f, "sd_all"], tight.loc[f, "sd_all"]
                dd = abs(d9) - abs(d0)
                ds = 100 * (s9 - s0) / s0 if s0 else np.nan
                verdict = ("holds/strengthens" if dd > -0.02 else
                           "weakens (SD also fell %.0f%%)" % ds if ds < -10 else
                           "weakens")
                say("  %-18s d %+.3f -> %+.3f   %s" % (f[:18], d0, d9, verdict))

open("step11_report.txt", "w").write("\n".join(log) + "\n")
say("")
say("wrote step11_gates.csv and step11_report.txt")
