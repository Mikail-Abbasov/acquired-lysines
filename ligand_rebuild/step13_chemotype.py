#!/usr/bin/env python3
"""
step13_chemotype.py — is ligandability a property of the SITE, or of the site-and-scaffold pair?

WHY IT MATTERS
Step 9 found that engaged lysines are more buried, more crowded and more aromatic, and step 10
combined those into a composite. Both treated "ligandable" as a property of the lysine. If the
two chemotypes in this data prefer different sites, that framing is wrong: ligandability would
be a property of the pair, each scaffold has its own preferences, and a general composite score
overstates what has been shown.

The two chemotypes are about as different as this data allows:
  CyMPSA        designed acylating electrophiles on a sulfonyl scaffold, 15 compounds, HeLa
  Meroterpenoid natural products, fused ring systems, aldehyde reactivity, 10 compounds,
                MDA-MB-231 and MDA-MB-436

THREE TESTS

  1. DO THE SAME FEATURES PREDICT ENGAGEMENT FOR BOTH?
     Comparison B run within each chemotype, then the difference in effect sizes tested
     directly rather than eyeballed. Two effects that "look similar" may still differ
     significantly, and two that look different may not.

  2. DO THEY ENGAGE THE SAME SITES?
     The sharper test. Among lysines measured by BOTH chemotypes, is a site engaged by one
     more likely than chance to be engaged by the other? If not, ligandability is
     scaffold-specific whatever the feature analysis says.

  3. WITHIN meroterpenoids, do ENANTIOMERS agree?
     (+)- and (-)-guadial B are chemically identical and differ only in three-dimensional
     shape. If even they engage different sites, site preference is shape-specific at the
     finest possible grain — and that is the cleanest evidence that binding, not reactivity,
     determines engagement.

TWO CONSTRAINTS THAT MUST BE REPORTED WITH ANY RESULT

  POWER. Meroterpenoids have ~311 engaged sites against CyMPSA's ~1,354. A failure to detect a
  difference may reflect sample size, not agreement. Confidence intervals are given for every
  comparison and a null is never reported as evidence of equivalence.

  DIFFERENT CELL LINES. CyMPSA is HeLa; meroterpenoids are MDA-MB-231/436. A site can be
  engaged in one and simply ABSENT from the other because the protein is not expressed there.
  Test 2 is therefore restricted to lysines actually measured in both, and the number lost to
  that restriction is reported.

INPUT   occ_long.csv, features_expanded.csv
OUTPUT  step13_chemotype.csv, step13_report.txt
"""
import os
import numpy as np, pandas as pd
import warnings
warnings.filterwarnings("ignore")
from scipy import stats as sstats
from statsmodels.stats.multitest import multipletests

TIER = int(os.environ.get("TIER", 60))
FEATS = ["rsa_tien", "sasa_sidechain", "nz_sasa", "crowding_8A", "atoms_10A",
         "hydrophobic_8A", "aromatic_8A", "polar_8A", "dist_centroid",
         "anionic_8A", "net_charge_8A", "propka_pKa", "plddt_site"]

log = []
def say(m):
    print(m, flush=True); log.append(m)

say("=" * 78)
say("STEP 13 — IS LIGANDABILITY CHEMOTYPE-SPECIFIC?")
say("engaged = >=%d%% occupancy by at least one compound of that chemotype" % TIER)
say("=" * 78)

occ = pd.read_csv("occ_long.csv", low_memory=False)
feat = (pd.read_csv("features_expanded.csv", low_memory=False)
          .drop_duplicates(["uniprot", "position"]))
F = [f for f in FEATS if f in feat.columns]

meas, eng = {}, {}
for ds in ("CyMPSA", "Meroterpenoid", "NatChem"):
    g = occ[occ.dataset == ds]
    meas[ds] = set(map(tuple, g[["uniprot", "position"]].drop_duplicates().values))
    eng[ds] = set(map(tuple, g[g.occupancy_pct >= TIER][["uniprot", "position"]]
                      .drop_duplicates().values))
    say("  %-14s measured %6d lysines | engaged %5d (%.2f%%)"
        % (ds, len(meas[ds]), len(eng[ds]), 100 * len(eng[ds]) / max(len(meas[ds]), 1)))

# ------------------------------------------------------------------ test 1
say("")
say("TEST 1 — do the same features predict engagement in each chemotype?")
say("  within-protein differences, then the two chemotypes compared directly")


def within(df, feats):
    out = []
    for f in feats:
        diffs = []
        for _, g in df.groupby("uniprot"):
            a = pd.to_numeric(g.loc[g.engaged, f], errors="coerce").dropna()
            b = pd.to_numeric(g.loc[~g.engaged, f], errors="coerce").dropna()
            if len(a) >= 1 and len(b) >= 2:
                diffs.append(a.mean() - b.mean())
        if len(diffs) < 25:
            continue
        v = np.array(diffs)
        out.append(dict(feature=f, n=len(v), d=float(v.mean() / v.std()) if v.std() else np.nan,
                        mean=float(v.mean()), sem=float(v.std() / np.sqrt(len(v)))))
    return pd.DataFrame(out)


res = {}
for ds in ("CyMPSA", "Meroterpenoid"):
    sub = feat[feat.apply(lambda r: (r.uniprot, int(r.position)) in meas[ds], axis=1)].copy()
    sub["engaged"] = [(u, int(p)) in eng[ds] for u, p in zip(sub.uniprot, sub.position)]
    res[ds] = within(sub, F)
    say("    %-14s %d lysines, %d engaged, %d proteins contribute"
        % (ds, len(sub), int(sub.engaged.sum()), res[ds].n.max() if len(res[ds]) else 0))

if all(len(res[k]) for k in res):
    a, b = res["CyMPSA"].set_index("feature"), res["Meroterpenoid"].set_index("feature")
    common = [f for f in F if f in a.index and f in b.index]
    rows = []
    for f in common:
        z = (a.loc[f, "mean"] - b.loc[f, "mean"]) / np.sqrt(a.loc[f, "sem"] ** 2 + b.loc[f, "sem"] ** 2)
        rows.append(dict(feature=f, d_cympsa=a.loc[f, "d"], d_mero=b.loc[f, "d"],
                         z=z, p=2 * sstats.norm.sf(abs(z))))
    cmp = pd.DataFrame(rows)
    cmp["q"] = multipletests(cmp.p, method="fdr_bh")[1]
    say("")
    say("  %-18s %10s %10s %9s %10s" % ("feature", "d CyMPSA", "d Mero", "same sign", "q (differ)"))
    for _, r in cmp.sort_values("q").iterrows():
        same = "yes" if np.sign(r.d_cympsa) == np.sign(r.d_mero) else "NO"
        say("  %-18s %10.3f %10.3f %9s %10.3g%s"
            % (r.feature[:18], r.d_cympsa, r.d_mero, same, r.q, " *" if r.q < 0.05 else ""))
    say("")
    say("  features with the same direction: %d of %d"
        % (sum(np.sign(cmp.d_cympsa) == np.sign(cmp.d_mero)), len(cmp)))
    say("  features differing significantly : %d of %d" % ((cmp.q < 0.05).sum(), len(cmp)))
    cmp.to_csv("step13_chemotype.csv", index=False)

# ------------------------------------------------------------------ test 2
say("")
say("TEST 2 — do the two chemotypes engage the SAME sites?")
both_meas = meas["CyMPSA"] & meas["Meroterpenoid"]
say("  lysines measured by BOTH chemotypes: %d" % len(both_meas))
say("  (CyMPSA is HeLa; meroterpenoids are MDA-MB-231/436, so a site can be engaged in one")
say("   and simply absent from the other — this restriction is essential)")
if both_meas:
    ec = eng["CyMPSA"] & both_meas
    em = eng["Meroterpenoid"] & both_meas
    both_eng = ec & em
    n = len(both_meas)
    exp = len(ec) * len(em) / n
    say("  engaged by CyMPSA        : %d (%.2f%%)" % (len(ec), 100 * len(ec) / n))
    say("  engaged by Meroterpenoid : %d (%.2f%%)" % (len(em), 100 * len(em) / n))
    say("  engaged by BOTH          : %d   (expected by chance %.1f)" % (len(both_eng), exp))
    a_ = len(both_eng); b_ = len(ec) - a_; c_ = len(em) - a_; d_ = n - a_ - b_ - c_
    orr, p = sstats.fisher_exact([[a_, b_], [c_, d_]])
    say("  odds ratio %.2f, p = %.3g" % (orr, p))
    if orr > 1 and p < 0.05:
        say("  -> the same sites tend to be engaged by both: ligandability is largely a")
        say("     SITE property, and a general composite is defensible")
    elif p >= 0.05:
        say("  -> no detectable overlap beyond chance. With %d and %d engaged sites this may")
        say("     be POWER, not evidence of independence — do not report as equivalence"
            % (len(ec), len(em)))
    else:
        say("  -> sites engaged by one are LESS likely to be engaged by the other:")
        say("     ligandability is scaffold-specific")

# ------------------------------------------------------------------ test 3
say("")
say("TEST 3 — do enantiomers engage the same sites?")
say("  same chemistry, opposite three-dimensional shape: the finest-grain test available")
me = occ[occ.dataset == "Meroterpenoid"]
pairs = []
for c in sorted(me.compound.unique()):
    if c.startswith("(+)-"):
        for minus in ("(\u2013)-" + c[4:], "(-)-" + c[4:]):
            if minus in set(me.compound):
                pairs.append((c, minus))
for plus, minus in pairs:
    w = me[me.compound.isin([plus, minus])].pivot_table(
        index=["uniprot", "position"], columns="compound", values="occupancy_pct").dropna()
    if plus not in w or minus not in w or len(w) < 50:
        continue
    p_e = w[plus] >= TIER; m_e = w[minus] >= TIER
    both_ = int((p_e & m_e).sum()); only_p = int((p_e & ~m_e).sum()); only_m = int((~p_e & m_e).sum())
    n = len(w)
    exp = p_e.sum() * m_e.sum() / n
    say("  %-24s vs %-24s  n=%5d" % (plus[:24], minus[:24], n))
    say("      both %3d (expected %.1f) | only(+) %3d | only(-) %3d | r = %.3f"
        % (both_, exp, only_p, only_m, w[plus].corr(w[minus])))
    if only_p + only_m > 0:
        say("      %d sites engaged by ONE enantiomer only — cannot be an intrinsic"
            " reactivity effect" % (only_p + only_m))

open("step13_report.txt", "w").write("\n".join(log) + "\n")
say("")
say("wrote step13_chemotype.csv and step13_report.txt")
