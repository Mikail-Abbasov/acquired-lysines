#!/usr/bin/env python3
"""
step44_charge_decomposition.py — how much of A49's charge effect is the reagent?

WHY
A57 concluded that A49's cationic-environment finding is "at least partly" electrostatic
steering of the anionic reagent. "At least partly" is a hedge, and the obvious way to sharpen it
is wrong: noting that STP-vs-unlabelled is significant for cationic_8A while NASA-vs-unlabelled
is not, and concluding the effect is all reagent. A difference between a significant and a
non-significant result is not itself a significant difference. That inference was drafted in
conversation on 2026-09-03 and is not used here.

WHAT IS ACTUALLY TESTED
Three quantities per charge feature, all on step 22's within-protein design (one difference per
protein, then across proteins), with bootstrap confidence intervals over PROTEINS — the unit of
independence, since sites on one protein share an environment:

    STP-only  vs unlabelled     what the anionic probe prefers
    NASA-only vs unlabelled     what the neutral probe prefers
    STP-only  vs NASA-only      the difference, which is the formal reagent test

The claim "the cationic preference is the reagent" is only supportable if the NASA interval
brackets zero AND excludes the STP estimate. If the NASA interval is wide, "at least partly"
stands and this script says so rather than reaching.

The corresponding claim for anionic_8A runs the other way: if both probes avoid acidic
neighbourhoods by a similar amount and the direct contrast is null, that half of A49 is a
property of the site and survives A57 untouched.
"""
import os, sys
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
from scipy import stats as sstats

GATE = float(os.environ.get("PLDDT_GATE", 70))
NBOOT = int(os.environ.get("NBOOT", 4000))
RNG = np.random.default_rng(0)
CHARGE = ["cationic_8A", "net_charge_8A", "anionic_8A", "propka_pKa", "salt_bridge", "n_lys_8A"]

log = []
def say(m):
    print(m, flush=True); log.append(m)

say("=" * 78)
say("STEP 44 — HOW MUCH OF A49'S CHARGE EFFECT BELONGS TO THE REAGENT?")
say("=" * 78)

occ = pd.read_csv("occ_long.csv", low_memory=False)
cy = occ[occ.dataset == "CyMPSA"]
sets = {}
for p, n in (("dTB-STP", "STP"), ("Probe 1", "NASA")):
    sets[n] = set(map(tuple, cy[cy.probe == p][["uniprot", "position"]].drop_duplicates().values))
stp_only = sets["STP"] - sets["NASA"]; nasa_only = sets["NASA"] - sets["STP"]

allK = pd.read_csv("step9_all_lysines.csv", low_memory=False)
pk = pd.read_csv("step19_pockets.csv", low_memory=False)
for c in ("salt_bridge",):
    if c in allK: allK[c] = pd.to_numeric(allK[c].replace({True: 1.0, False: 0.0}), errors="coerce")
allK = allK[allK.plddt_site >= GATE].copy()
allK["_k"] = list(zip(allK.uniprot, allK.position.astype(int)))
mS = allK._k.isin(stp_only).values
mN = allK._k.isin(nasa_only).values
mU = ~allK._k.isin(sets["STP"] | sets["NASA"]).values
say("lysines above pLDDT %.0f: %s | STP-only %d | NASA-only %d | unlabelled %s"
    % (GATE, f"{len(allK):,}", mS.sum(), mN.sum(), f"{mU.sum():,}"))

def per_protein(mask_a, mask_b, feat):
    """One difference per protein: mean(a) - mean(b), both groups >= 2 sites."""
    d = allK[mask_a | mask_b].copy()
    d["_g"] = mask_a[mask_a | mask_b]
    out = []
    for _, g in d.groupby("uniprot"):
        a = pd.to_numeric(g.loc[g._g, feat], errors="coerce").dropna()
        b = pd.to_numeric(g.loc[~g._g, feat], errors="coerce").dropna()
        if len(a) >= 2 and len(b) >= 2:
            out.append(a.mean() - b.mean())
    return np.array(out)

def boot_d(v):
    """Standardised effect (mean/sd across proteins) with a percentile CI, resampling proteins."""
    if len(v) < 25 or v.std() == 0: return np.nan, np.nan, np.nan, len(v)
    idx = RNG.integers(0, len(v), size=(NBOOT, len(v)))
    bs = v[idx]
    ds = bs.mean(1) / bs.std(1)
    return float(v.mean()/v.std()), float(np.percentile(ds, 2.5)), float(np.percentile(ds, 97.5)), len(v)

CMP = [("STP-only vs unlabelled", mS, mU), ("NASA-only vs unlabelled", mN, mU),
       ("STP-only vs NASA-only", mS, mN)]
rows = []
for feat in CHARGE:
    for name, a, b in CMP:
        v = per_protein(a, b, feat)
        d, lo, hi, n = boot_d(v)
        rows.append(dict(feature=feat, comparison=name, n_proteins=n, d=d, lo=lo, hi=hi))
R = pd.DataFrame(rows)

for feat in CHARGE:
    say("")
    say("  %s" % feat)
    say("    %-26s %8s %8s %20s" % ("comparison", "proteins", "d", "95% CI"))
    for _, r in R[R.feature == feat].iterrows():
        say("    %-26s %8d %8.3f   [%+.3f, %+.3f]"
            % (r.comparison, r.n_proteins, r.d, r.lo, r.hi))

say("")
say("-" * 78)
say("THE DECOMPOSITION")
def get(f, c):
    x = R[(R.feature == f) & (R.comparison == c)]
    return x.iloc[0] if len(x) else None
for feat, story in (("cationic_8A", "cationic attraction"), ("anionic_8A", "acidic avoidance"),
                    ("net_charge_8A", "net local charge")):
    s_, n_, dcon = get(feat, CMP[0][0]), get(feat, CMP[1][0]), get(feat, CMP[2][0])
    if s_ is None or n_ is None or dcon is None: continue
    nasa_null = n_.lo <= 0 <= n_.hi
    excl = not (n_.lo <= s_.d <= n_.hi)
    contrast_sig = not (dcon.lo <= 0 <= dcon.hi)
    say("")
    say("  %s (%s)" % (feat, story))
    say("    neutral probe's CI brackets zero:            %s" % ("YES" if nasa_null else "no"))
    say("    neutral probe's CI excludes the anionic d:   %s" % ("YES" if excl else "no"))
    say("    direct contrast excludes zero:               %s" % ("YES" if contrast_sig else "no"))
    if contrast_sig and nasa_null and excl:
        say("    -> attributable to the REAGENT. A49 may state this, not 'at least partly'.")
    elif contrast_sig:
        say("    -> the reagent CONTRIBUTES, but the neutral probe's own effect is not")
        say("       bounded away from the anionic one. 'At least partly' STANDS.")
    else:
        say("    -> no difference between probes: consistent with a property of the SITE.")

R.to_csv("step44_charge_decomposition.csv", index=False)
open("step44_report.txt", "w").write("\n".join(log) + "\n")
say("")
say("wrote step44_charge_decomposition.csv, step44_report.txt  (%d bootstrap resamples)" % NBOOT)
