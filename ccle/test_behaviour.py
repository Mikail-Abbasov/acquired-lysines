#!/usr/bin/env python3
"""
test_behaviour.py — do reactivity and ligandability really behave differently, or is the
apparent difference a property of how each was measured?

THE PROBLEM
Across all 16 shared features, the reactivity effect and the ligandability effect correlate at
r = 0.941. The same features push the same way for both. Figure 4's title claims "distinct
determinants", and that correlation says otherwise.

What remains is a claim about BEHAVIOUR: ligandability shows a monotonic gradient across
occupancy tiers, and its shape features collapse under a structural-confidence gate while its
electrostatic features strengthen. Reactivity does neither.

But there are two explanations for that, and they are not the same claim:

  (1) MEASUREMENT. Occupancy tiers order sites by how extremely they were engaged, so effects
      grow with tier for the trivial reason that higher tiers are more extreme selections of the
      same population. The pLDDT gate orders nothing, so no gradient can appear. On this reading
      the two properties are one property measured two ways.

  (2) REAL. Under gating, electrostatic and shape features move in OPPOSITE directions. Simple
      selection on extremity cannot produce that.

TEST 1 — tier reactivity the same way.
If labelling can be tiered by intensity and shows the same gradient, the gradient is a property
of the design, not of ligandability.

TEST 2 — does the correlation survive gating?
If the same features matter but their relative weights shift under a pLDDT gate, that is a
difference in kind. If the ordering is unchanged, the two properties differ only in scale.

OUTPUT  test_behaviour_report.txt
"""
import os, re, sys
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
from scipy import stats as sstats

log = []
def say(m):
    print(m, flush=True); log.append(m)

def find(n):
    for d in (".", "..", os.path.expanduser("~/Desktop/DepMap/ligand_rebuild")):
        p = os.path.join(d, n)
        if os.path.exists(p): return p
    return None

say("=" * 78)
say("DO THE TWO PROPERTIES BEHAVE DIFFERENTLY, OR ARE THEY MEASURED DIFFERENTLY?")
say("=" * 78)

S9 = pd.read_csv(find("step9_results.csv"))
G = pd.read_csv(find("step11_gates.csv"))
A = S9[S9.comparison == "A_reactivity"][["feature", "cohens_d"]].set_index("feature").cohens_d
B = S9[S9.comparison == "B_ligandability"]
tiers = sorted(B.tier.dropna().unique()) if "tier" in B.columns else []

# ---------------------------------------------------------------- baseline
say("")
say("BASELINE — the correlation that prompted this")
for t in tiers:
    b = B[B.tier == t].set_index("feature").cohens_d
    common = A.index.intersection(b.index)
    if len(common) > 3:
        r = np.corrcoef(A[common], b[common])[0, 1]
        rho = sstats.spearmanr(A[common], b[common]).correlation
        say("  ligandability tier %-6s vs reactivity:  r = %.3f, rho = %.3f, n = %d"
            % (t, r, rho, len(common)))

# ---------------------------------------------------------------- test 1
say("")
say("=" * 78)
say("TEST 1 — CAN REACTIVITY BE TIERED THE SAME WAY?")
say("=" * 78)
say("If the ligandability gradient comes from ordering sites by extremity, then ordering")
say("labelled sites by labelling intensity should produce the same gradient.")
say("")
OCC = find("occ_long.csv")
FEA = find("features_final_full.csv") or find("features_expanded.csv") or find("features.csv")
if not OCC or not FEA:
    say("  occ_long.csv or the feature table not found — cannot run")
else:
    O = pd.read_csv(OCC, low_memory=False, usecols=["uniprot", "position", "occupancy_pct"])
    F = pd.read_csv(FEA, low_memory=False)
    fcols = [c for c in ("anionic_8A", "net_charge_8A", "propka_pKa", "crowding_8A",
                         "hydrophobic_8A", "aromatic_8A", "rsa_tien") if c in F.columns]
    if not fcols:
        say("  no shared feature columns in %s" % os.path.basename(FEA))
    else:
        # labelling intensity per site: a labelled lysine's evidence of being labelled is that it
        # was quantified at all, so the closest available analogue is how many compounds and
        # experiments detected it. Tier on that.
        depth = O.groupby(["uniprot", "position"]).size().rename("n_obs").reset_index()
        M = depth.merge(F.drop_duplicates(["uniprot", "position"]), on=["uniprot", "position"])
        M = M.dropna(subset=fcols)
        M["tier"] = pd.qcut(M.n_obs, 10, labels=False, duplicates="drop")
        base = M[M.tier == M.tier.min()]
        say("  %s labelled sites tiered by detection depth into %d bins"
            % (f"{len(M):,}", M.tier.nunique()))
        say("")
        say("  %-16s %s" % ("feature", "  ".join("t%d" % i for i in sorted(M.tier.unique()))))
        slopes = {}
        for f in fcols:
            sd = M[f].std()
            vals = [(M.loc[M.tier == t, f].mean() - base[f].mean()) / sd
                    for t in sorted(M.tier.unique())]
            slopes[f] = np.polyfit(range(len(vals)), vals, 1)[0]
            say("  %-16s %s" % (f, "  ".join("%+.2f" % v for v in vals)))
        say("")
        say("  slope per tier:")
        for f, s_ in sorted(slopes.items(), key=lambda kv: -abs(kv[1])):
            say("    %-18s %+.4f" % (f, s_))
        say("")
        mx = max(abs(v) for v in slopes.values())
        if mx < 0.03:
            say("  -> reactivity does NOT show a gradient when tiered. The ligandability")
            say("     gradient is therefore a property of ligandability, not of tiering.")
        else:
            say("  -> reactivity DOES show a gradient when tiered (max slope %.3f). The" % mx)
            say("     ligandability gradient may be a property of the design, not the property.")

# ---------------------------------------------------------------- test 2
say("")
say("=" * 78)
say("TEST 2 — DOES THE CORRELATION SURVIVE STRUCTURAL GATING?")
say("=" * 78)
say("If the same features matter in the same order at every gate, the two properties differ")
say("only in scale. If the ordering changes, they differ in kind.")
say("")
gates = sorted(G.gate.unique())
say("  reactivity effects by gate, and their correlation with the no-gate ordering:")
g0 = G[G.gate == gates[0]].set_index("feature").d
for gt in gates:
    gg = G[G.gate == gt].set_index("feature").d
    common = g0.index.intersection(gg.index)
    r = np.corrcoef(g0[common], gg[common])[0, 1]
    rho = sstats.spearmanr(g0[common], gg[common]).correlation
    say("    gate %-5s r = %.3f, rho = %.3f, n = %d" % (gt, r, rho, len(common)))
say("")
say("  and per feature, effect at no gate vs the tightest gate:")
gL = G[G.gate == gates[-1]].set_index("feature").d
common = g0.index.intersection(gL.index)
tab = pd.DataFrame({"no_gate": g0[common], "gated": gL[common]})
tab["change"] = tab.gated - tab.no_gate
for r_ in tab.reindex(tab.change.abs().sort_values(ascending=False).index).head(8).itertuples():
    say("    %-18s %+.3f -> %+.3f   (%+.3f)" % (r_.Index, r_.no_gate, r_.gated, r_.change))
say("")
opp = ((tab.no_gate > 0) & (tab.change < -abs(tab.no_gate))) | \
      ((tab.no_gate < 0) & (tab.change > abs(tab.no_gate)))
say("  features that REVERSE or vanish under gating: %d of %d" % (int(opp.sum()), len(tab)))
say("  features that strengthen: %d" % int(((tab.gated.abs() > tab.no_gate.abs())).sum()))
say("")
say("  If some strengthen while others collapse, selection on extremity cannot be the")
say("  explanation, because that would move everything the same way.")

open("test_behaviour_report.txt", "w").write("\n".join(log) + "\n")
say("")
say("wrote test_behaviour_report.txt")
