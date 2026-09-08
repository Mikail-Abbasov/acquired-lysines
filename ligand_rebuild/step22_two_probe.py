#!/usr/bin/env python3
"""
step22_two_probe.py — is the charge effect lysine nucleophilicity, or electrostatic steering of
an anionic reagent?

WHAT IS AT STAKE
A49 concluded that probe labelling is governed by the local charge environment: labelled lysines
have FEWER acidic neighbours (d = -0.46 at pLDDT >=70), MORE cationic ones (+0.40), fewer salt
bridges (-0.32) and lower predicted pKa (-0.30). I read that as nucleophilicity — fewer
carboxylates stabilising the ammonium, so more neutral amine available to attack.

There is a second reading. **An anionic reagent is electrostatically steered toward cationic
surfaces and repelled from acidic ones.** On that account, "labelled lysines sit in cationic
environments" describes where the PROBE goes, not which lysines are reactive. This is the
bound-anionic-reagent picture the lab has argued supersedes simple pKa depression.

The two readings make opposite predictions about probe identity, and CyMPSA can separate them.

THE TWO PROBES — same desthiobiotin-hexanoyl core, same proteome, same cell line, differing
only in the acyl leaving group:

  dTB-STP    4-sulfo-2,3,5,6-tetrafluorophenyl ester.
             **ANIONIC** — aromatic sulfonate, negatively charged at physiological pH.
  Probe 1    N-(cyanomethyl)-N-(4-nitrobenzenesulfonyl) amide ("NASA").
             **NEUTRAL** overall. Carries enolizable alpha-cyano hydrogens alpha to an
             electron-withdrawing acylsulfonamide; enolisation may occur at physiological pH
             but the pKa is not known, so this probe may acquire partial anionic character.
             That uncertainty weakens the contrast and must be stated with any result.

PREDICTIONS, fixed before computing
  REAGENT STEERING: the anionic STP probe shows a STRONGER preference for cationic, acid-poor
      environments than the neutral NASA probe. The charge effect in A49 is then partly about
      reagent electrostatics, not lysine nucleophilicity, and A49's wording must change.
  LYSINE PROPERTY: both probes show the SAME charge preference. The effect belongs to the site
      and A49 stands as written.
  If NASA enolises substantially, both probes are effectively anionic and the test is
      underpowered — a null result would then be uninformative rather than supportive.

ALSO TESTED — the leaving-group bulk hypothesis (A56)
  STP is a medium aromatic fluorinated ring; NASA is a bulkier nitrophenyl sulfonamide. If
  pocket preference (A55: in_pocket d = 0.107 at pLDDT >=70) reflects accommodation of the
  leaving group at the transition state, the two probes should differ there too.

DESIGN
  Three site sets, all from the same HeLa proteome:
      STP_ONLY     labelled by dTB-STP, not by Probe 1      (6,294 sites)
      NASA_ONLY    labelled by Probe 1, not by dTB-STP      (5,851)
      BOTH         labelled by both                          (12,298)
  Each compared against unlabelled lysines on the same proteins, within-protein, at pLDDT >=70.
  Then STP_ONLY vs NASA_ONLY directly — the discriminating contrast, since these differ only in
  which probe reached them.

INPUT   occ_long.csv, step9_all_lysines.csv, step19_pockets.csv
OUTPUT  step22_two_probe.csv, step22_report.txt
"""
import os
import numpy as np, pandas as pd
import warnings
warnings.filterwarnings("ignore")
from scipy import stats as sstats
from statsmodels.stats.multitest import multipletests

GATE = float(os.environ.get("PLDDT_GATE", 70))
# All available features are tested. Restricting to a subset would let the choice of which
# features to examine shape the conclusion, and n_lys_8A in particular is directly relevant:
# if an anionic reagent is steered toward cationic surfaces, nearby LYSINES are that positive
# charge.
CHARGE = ["anionic_8A", "cationic_8A", "net_charge_8A", "propka_pKa", "salt_bridge",
          "hbond_partners_NZ", "n_lys_8A"]
POCKET = ["in_pocket", "dist_nearest_pocket", "pocket_volume", "pocket_score",
          "pocket_drugscore"]
SHAPE = ["crowding_8A", "atoms_10A", "hydrophobic_8A", "aromatic_8A", "polar_8A",
         "rsa_tien", "sasa_sidechain", "nz_sasa", "dist_centroid"]
CONTEXT = ["plddt_site", "plddt_8A"]

log = []
def say(m):
    print(m, flush=True); log.append(m)

say("=" * 78)
say("STEP 22 — REAGENT STEERING OR LYSINE REACTIVITY?")
say("dTB-STP (anionic leaving group) vs Probe 1 / NASA (neutral, possibly enolisable)")
say("=" * 78)

occ = pd.read_csv("occ_long.csv", low_memory=False)
cy = occ[occ.dataset == "CyMPSA"]
sets = {}
for p, name in (("dTB-STP", "STP"), ("Probe 1", "NASA")):
    sets[name] = set(map(tuple, cy[cy.probe == p][["uniprot", "position"]]
                         .drop_duplicates().values))
stp_only = sets["STP"] - sets["NASA"]
nasa_only = sets["NASA"] - sets["STP"]
both = sets["STP"] & sets["NASA"]
say("  STP-only  %6d | NASA-only %6d | both %6d" % (len(stp_only), len(nasa_only), len(both)))

allK = pd.read_csv("step9_all_lysines.csv", low_memory=False)
if os.path.exists("step19_pockets.csv"):
    pk = pd.read_csv("step19_pockets.csv", low_memory=False)
    cols = [c for c in ("uniprot", "position") + tuple(POCKET) if c in pk.columns]
    allK = allK.merge(pk[cols].drop_duplicates(["uniprot", "position"]),
                      on=["uniprot", "position"], how="left", suffixes=("", "_pk"))
for c in ("salt_bridge", "in_pocket"):
    if c in allK:
        allK[c] = pd.to_numeric(allK[c].replace({True: 1.0, False: 0.0}), errors="coerce")
allK = allK[allK.plddt_site >= GATE].copy()
allK["_k"] = list(zip(allK.uniprot, allK.position.astype(int)))
FEATS = [f for f in CHARGE + POCKET + SHAPE + CONTEXT if f in allK.columns]
say("  %d lysines above pLDDT %.0f | %d features" % (len(allK), GATE, len(FEATS)))


def within(mask_a, mask_b, feats):
    """One difference per protein, then tested against zero."""
    out = []
    d = allK[mask_a | mask_b].copy()
    d["_g"] = mask_a[mask_a | mask_b]
    for f in feats:
        diffs = []
        for _, g in d.groupby("uniprot"):
            a = pd.to_numeric(g.loc[g._g, f], errors="coerce").dropna()
            b = pd.to_numeric(g.loc[~g._g, f], errors="coerce").dropna()
            if len(a) >= 2 and len(b) >= 2:
                diffs.append(a.mean() - b.mean())
        if len(diffs) < 25:
            continue
        v = np.array(diffs)
        t, p = sstats.ttest_1samp(v, 0)
        out.append(dict(feature=f, n_proteins=len(v), mean_diff=float(v.mean()),
                        d=float(v.mean() / v.std()) if v.std() else np.nan, p=float(p)))
    r = pd.DataFrame(out)
    if len(r):
        r["q"] = multipletests(r.p, method="fdr_bh")[1]
    return r


labelled_any = allK._k.isin(sets["STP"] | sets["NASA"]).values
unlab = ~labelled_any
rows = []
for name, S in (("STP_ONLY", stp_only), ("NASA_ONLY", nasa_only), ("BOTH", both)):
    m = allK._k.isin(S).values
    if m.sum() < 200:
        say("  %s: too few sites above the gate" % name); continue
    r = within(m, unlab, FEATS)
    if not len(r):
        continue
    r["group"] = name
    rows.append(r)
    say("")
    say("-" * 78)
    say("%s vs unlabelled  (n = %d vs %d)" % (name, int(m.sum()), int(unlab.sum())))
    say("  %-22s %9s %10s %8s %10s" % ("feature", "proteins", "mean diff", "d", "q"))
    for _, x in r.sort_values("q").iterrows():
        blk = "chg" if x.feature in CHARGE else ("pkt" if x.feature in POCKET else "shp")
        say("  %-22s %9d %10.4f %8.3f %10.3g %s%s"
            % (x.feature[:22], x.n_proteins, x.mean_diff, x.d, x.q, blk,
               " *" if x.q < 0.05 else ""))

say("")
say("=" * 78)
say("THE DISCRIMINATING CONTRAST — STP-only vs NASA-only")
say("  These sites differ only in which probe reached them, in the same proteome and cells.")
mS = allK._k.isin(stp_only).values
mN = allK._k.isin(nasa_only).values
direct = within(mS, mN, FEATS)
if len(direct):
    say("  positive d = higher in STP-only (the ANIONIC probe)")
    say("  %-22s %9s %10s %8s %10s" % ("feature", "proteins", "mean diff", "d", "q"))
    for _, x in direct.sort_values("q").iterrows():
        blk = ("chg" if x.feature in CHARGE else "pkt" if x.feature in POCKET
               else "shp" if x.feature in SHAPE else "ctx")
        say("  %-22s %9d %10.4f %8.3f %10.3g %s%s"
            % (x.feature[:22], x.n_proteins, x.mean_diff, x.d, x.q, blk,
               " *" if x.q < 0.05 else ""))
    direct["group"] = "STP_ONLY vs NASA_ONLY"
    rows.append(direct)
    say("")
    chg = direct[direct.feature.isin(CHARGE)]
    pkt = direct[direct.feature.isin(POCKET)]
    say("  charge features differing between probes: %d of %d"
        % (int((chg.q < 0.05).sum()), len(chg)))
    say("  pocket features differing between probes: %d of %d"
        % (int((pkt.q < 0.05).sum()), len(pkt)))
    shp = direct[direct.feature.isin(SHAPE)]
    say("  shape  features differing between probes: %d of %d"
        % (int((shp.q < 0.05).sum()), len(shp)))
    say("  ALL    features differing between probes: %d of %d"
        % (int((direct.q < 0.05).sum()), len(direct)))
    say("")
    if len(chg) and (chg.q < 0.05).any():
        say("  -> REAGENT STEERING is supported: the two probes prefer different charge")
        say("     environments. A49's charge effect is then partly about reagent")
        say("     electrostatics, not lysine nucleophilicity. A49 must be reworded.")
    else:
        say("  -> No detectable difference in charge preference between an anionic and a")
        say("     neutral leaving group. Consistent with the effect belonging to the SITE.")
        say("     CAVEAT: if the NASA probe enolises at physiological pH it is also anionic,")
        say("     and this test is then underpowered rather than negative. The alpha-cyano")
        say("     pKa is unknown, so a null here cannot fully exclude reagent steering.")

if rows:
    pd.concat(rows, ignore_index=True).to_csv("step22_two_probe.csv", index=False)
open("step22_report.txt", "w").write("\n".join(log) + "\n")
say("")
say("wrote step22_two_probe.csv and step22_report.txt")
