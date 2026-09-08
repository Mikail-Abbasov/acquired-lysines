#!/usr/bin/env python3
"""
step31_compound_site_matching.py — do compound properties correspond to the microenvironments of
the lysines they engage?

THE QUESTION, and why it is different from everything before it
Every analysis so far has treated the compound as a means of measurement and asked what
distinguishes an engaged SITE. That framing found little: protein structural features explain
0.6-2.5% of variance in engagement rate, and a composite classifier reaches AUC 0.617 against a
label-reproducibility ceiling of at least 0.72.

The PI's question is a matching question rather than a counting one: do compounds with
particular properties engage lysines with particular microenvironments? A bulky compound finding
open sites, an anionic compound finding cationic environments, a lipophilic compound finding
hydrophobic pockets. That is a correspondence between two sets of features, and it cannot be
detected by any analysis that collapses across compounds — which is what all the previous work
did.

THE PRINCIPLE IS ALREADY DEMONSTRATED AT n = 2
Step 22 compared two probes differing only in leaving group. The anionic dTB-STP engaged lysines
in more cationic environments than the neutral NASA probe: 2 of 22 site features differed, both
charge, none of shape. That is exactly this analysis with two compounds. This scales it to 212.

DESIGN
  For each compound: compute molecular descriptors from SMILES.
  For each compound: compute the MEAN microenvironment of the lysines it engages (>=60% occupancy).
  Then correlate the two across compounds.

  A significant correlation between, say, compound formal charge and mean engaged-site
  net_charge_8A means compounds are finding sites that suit them — a compound-site matching
  effect invisible to the site-only analyses.

WHAT MUST BE CONTROLLED, and why the naive version fails
  1. **Compounds engaging few sites give noisy means.** A compound with three engaged sites has
     an unreliable mean environment. Minimum site count enforced, and results reported at
     several thresholds.
  2. **Different compounds were screened against different site sets.** A compound tested only
     on abundant proteins will show their environments regardless of its own properties. Mean
     engaged-site features are therefore expressed RELATIVE to the mean of the sites that
     compound was actually tested against — a within-compound baseline.
  3. **Dataset differences.** Compounds from different papers were run in different cell lines
     and sample states. Correlations are computed within dataset as well as pooled.
  4. **Multiple testing.** Descriptors x site features is a large grid; Benjamini-Hochberg
     across all pairs, and the number expected by chance is stated.

INPUT   /path/to/SMILES.xlsx (three sheets), occ_long.csv, features_expanded.csv
OUTPUT  step31_matching.csv, step31_matching_report.txt
        (the report was called step31_report.txt until 2026-09-07, which
        step31_crystal_validation.py also writes -- it silently overwrote this
        one, and the within-dataset arm of this analysis was lost with it.)
"""
import os, re
import numpy as np, pandas as pd
import warnings
warnings.filterwarnings("ignore")
from scipy import stats as sstats
from statsmodels.stats.multitest import multipletests

SMILES_XLSX = os.environ.get("SMILES_XLSX", "SMILES.xlsx")
TIER = int(os.environ.get("TIER", 60))
# Meroterpenoids are only 11 compounds and several engage fewer than 10 sites; a threshold of
# 10 excluded the dataset entirely. That is the wrong trade here — the meroterpenoids carry a
# measured formal charge of -1, span distinct scaffolds, and include an enantiomer pair known to
# engage different sites, so they are the dataset where compound-site matching is most likely to
# be visible. The threshold is lowered and its effect reported.
MIN_SITES = int(os.environ.get("MIN_SITES", 5))

SITE_FEATS = ["rsa_tien", "nz_sasa", "crowding_8A", "atoms_10A", "hydrophobic_8A",
              "aromatic_8A", "polar_8A", "anionic_8A", "cationic_8A", "net_charge_8A",
              "n_lys_8A", "dist_centroid", "propka_pKa", "in_pocket", "pocket_volume"]

log = []
def say(m):
    print(m, flush=True); log.append(m)

say("=" * 78)
say("STEP 31 — DO COMPOUND PROPERTIES MATCH ENGAGED-SITE MICROENVIRONMENTS?")
say("=" * 78)

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors, Crippen
except ImportError:
    say("RDKit not installed:  pip install rdkit --break-system-packages")
    raise SystemExit(1)

# ---- compound descriptors
x = pd.ExcelFile(SMILES_XLSX)
sheet_ds = {"NatChem": "NatChem", "Meroterpenoids": "Meroterpenoid", "CyMPSA": "CyMPSA"}
rows = []
for sh, ds in sheet_ds.items():
    if sh not in x.sheet_names:
        continue
    d = x.parse(sh).dropna(subset=["SMILES"])
    for r in d.itertuples():
        m = Chem.MolFromSmiles(str(r.SMILES))
        if m is None:
            continue
        rows.append(dict(dataset=ds, compound=str(r.Compounds).strip(),
                         MW=Descriptors.MolWt(m),
                         logP=Crippen.MolLogP(m),
                         TPSA=Descriptors.TPSA(m),
                         HBD=rdMolDescriptors.CalcNumHBD(m),
                         HBA=rdMolDescriptors.CalcNumHBA(m),
                         rings=rdMolDescriptors.CalcNumRings(m),
                         aromatic_rings=rdMolDescriptors.CalcNumAromaticRings(m),
                         rotatable=rdMolDescriptors.CalcNumRotatableBonds(m),
                         heavy=m.GetNumHeavyAtoms(),
                         formal_charge=Chem.GetFormalCharge(m),
                         fsp3=rdMolDescriptors.CalcFractionCSP3(m)))
cmp = pd.DataFrame(rows)
DESC = ["MW", "logP", "TPSA", "HBD", "HBA", "rings", "aromatic_rings", "rotatable",
        "heavy", "formal_charge", "fsp3"]
say("compounds with descriptors: %d" % len(cmp))
for ds, g in cmp.groupby("dataset"):
    say("  %-14s %3d" % (ds, len(g)))

# ---- site features
feat = pd.read_csv("features_expanded.csv", low_memory=False).drop_duplicates(["uniprot", "position"])
if "in_pocket" in feat:
    feat["in_pocket"] = feat.in_pocket.astype(str).str.lower().isin(["true", "1", "1.0"]).astype(float)
SF = [f for f in SITE_FEATS if f in feat.columns]
fidx = feat.set_index(["uniprot", "position"])[SF]

occ = pd.read_csv("occ_long.csv", low_memory=False)
occ["base"] = [re.match(r"^(\d+[a-z])", str(c)).group(1)
               if (d == "NatChem" and re.match(r"^(\d+[a-z])", str(c))) else str(c)
               for c, d in zip(occ.compound, occ.dataset)]

# ---- per compound: mean engaged-site features RELATIVE to sites it was tested against
say("")
say("computing per-compound engaged-site environments (relative to each compound's own")
say("tested-site baseline, so differences in which sites were screened do not leak in)...")
recs = []
for (ds, base), g in occ.groupby(["dataset", "base"]):
    eng = g[g.occupancy_pct >= TIER][["uniprot", "position"]].drop_duplicates()
    tst = g[["uniprot", "position"]].drop_duplicates()
    if len(eng) < MIN_SITES:
        continue
    ei = fidx.reindex(pd.MultiIndex.from_frame(eng)).dropna(how="all")
    ti = fidx.reindex(pd.MultiIndex.from_frame(tst)).dropna(how="all")
    if len(ei) < MIN_SITES or len(ti) < 50:
        continue
    rec = dict(dataset=ds, compound=base, n_engaged=len(ei), n_tested=len(ti))
    for f in SF:
        a = pd.to_numeric(ei[f], errors="coerce").dropna()
        b = pd.to_numeric(ti[f], errors="coerce").dropna()
        if len(a) >= 5 and len(b) >= 20 and b.std() > 0:
            rec["d_" + f] = float((a.mean() - b.mean()) / b.std())   # standardised shift
    recs.append(rec)
env = pd.DataFrame(recs)
say("  %d compounds with >=%d engaged sites" % (len(env), MIN_SITES))
if len(env) < 15:
    say("  too few compounds to correlate — lower MIN_SITES or check the compound-name join")
    raise SystemExit(0)

m = env.merge(cmp, on=["dataset", "compound"], how="inner")
say("  %d matched to a structure" % len(m))
unmatched = set(zip(env.dataset, env.compound)) - set(zip(m.dataset, m.compound))
if unmatched:
    say("  UNMATCHED (name mismatch between SMILES file and data): %d" % len(unmatched))
    for u in sorted(unmatched)[:8]:
        say("     %s / %s" % u)

# ---- correlate
DS = [c for c in DESC if c in m.columns and m[c].std() > 0]
SFD = [c for c in m.columns if c.startswith("d_")]
out = []
for dcol in DS:
    for scol in SFD:
        sub = m[[dcol, scol]].dropna()
        if len(sub) < 15:
            continue
        r, p = sstats.spearmanr(sub[dcol], sub[scol])
        out.append(dict(descriptor=dcol, site_feature=scol[2:], n=len(sub), rho=r, p=p))
res = pd.DataFrame(out)
res["q"] = multipletests(res.p, method="fdr_bh")[1]
res = res.reindex(res.rho.abs().sort_values(ascending=False).index)
res.to_csv("step31_matching.csv", index=False)

say("")
say("COMPOUND DESCRIPTOR vs ENGAGED-SITE ENVIRONMENT (pooled across datasets)")
say("  %d descriptor x site-feature pairs tested; %.1f expected significant by chance at q<0.05"
    % (len(res), 0.05 * len(res)))
say("")
say("  %-16s %-18s %5s %8s %10s" % ("compound property", "site feature", "n", "rho", "q"))
for _, r in res.head(20).iterrows():
    say("  %-16s %-18s %5d %8.3f %10.3g%s"
        % (r.descriptor, r.site_feature[:18], r.n, r.rho, r.q, " *" if r.q < 0.05 else ""))
sig = int((res.q < 0.05).sum())
say("")
say("  significant pairs: %d of %d" % (sig, len(res)))

say("")
say("WITHIN DATASET (cell line and sample state held constant)")
for ds, g in m.groupby("dataset"):
    if len(g) < 6:
        say("  %-14s %d compounds — too few to correlate" % (ds, len(g))); continue
    if len(g) < 15:
        say("  %-14s %d compounds — UNDERPOWERED, correlations shown without correction"
            % (ds, len(g)))
    o = []
    for dcol in DS:
        for scol in SFD:
            sub = g[[dcol, scol]].dropna()
            if len(sub) < 12 or sub[dcol].std() == 0:
                continue
            r, p = sstats.spearmanr(sub[dcol], sub[scol])
            o.append(dict(descriptor=dcol, site_feature=scol[2:], rho=r, p=p, n=len(sub)))
    if not o:
        continue
    oo = pd.DataFrame(o)
    oo["q"] = multipletests(oo.p, method="fdr_bh")[1]
    top = oo.reindex(oo.rho.abs().sort_values(ascending=False).index).head(5)
    say("  %-14s %d compounds | %d of %d pairs significant" % (ds, len(g), int((oo.q < 0.05).sum()), len(oo)))
    for _, r in top.iterrows():
        say("     %-16s %-18s rho %+.3f  q %.3g%s"
            % (r.descriptor, r.site_feature[:18], r.rho, r.q, " *" if r.q < 0.05 else ""))

# ---- sensitivity to the site-count threshold
say("")
say("SENSITIVITY TO THE MINIMUM ENGAGED-SITE THRESHOLD")
say("  %-6s %10s %12s" % ("min", "compounds", "sig pairs"))
for thr in (3, 5, 10, 20):
    keep = m[m.n_engaged >= thr]
    if len(keep) < 15:
        say("  %-6d %10d %12s" % (thr, len(keep), "too few")); continue
    o = []
    for dcol in DS:
        for scol in SFD:
            sub = keep[[dcol, scol]].dropna()
            if len(sub) < 15 or sub[dcol].std() == 0:
                continue
            r, p = sstats.spearmanr(sub[dcol], sub[scol])
            o.append(p)
    if o:
        q = multipletests(o, method="fdr_bh")[1]
        say("  %-6d %10d %12d" % (thr, len(keep), int((q < 0.05).sum())))

# ---- enantiomers: identical 2D descriptors, so a real difference is stereochemical
say("")
say("ENANTIOMER CHECK — same 2D descriptors, different engaged sites?")
mm = m.copy()
mm["stem"] = mm.compound.str.replace(r"^\(\+\)-|^\(–\)-|^\(-\)-", "", regex=True)
for stem, g in mm.groupby("stem"):
    if len(g) < 2:
        continue
    say("  %s: %d forms" % (stem, len(g)))
    for f in SFD[:6]:
        vals = ["%s %+.2f" % (c[:12], v) for c, v in zip(g.compound, g[f]) if v == v]
        if len(vals) >= 2:
            say("     %-18s %s" % (f[2:], " | ".join(vals)))
    say("     (2D descriptors are identical for these; any difference above is stereochemical")
    say("      and cannot be captured by the descriptor correlations)")

say("")
say("READING THIS")
say("  A correlation here means compounds are engaging sites that SUIT them — a matching effect")
say("  invisible to every site-only analysis in this work, all of which collapsed across")
say("  compounds. Step 22 showed the effect exists at n = 2 (anionic probe -> cationic sites).")
say("")
say("  Interpret with care:")
say("    - descriptors are correlated with each other (MW with heavy atoms, rings with logP),")
say("      so significant pairs are not independent findings")
say("    - site-feature shifts are relative to each compound's own tested-site baseline, which")
say("      controls for screening coverage but not for potency or concentration")
say("    - a null result does not exclude matching: it may require pairing at the level of")
say("      individual site-compound observations rather than compound means")
open("step31_matching_report.txt", "w").write("\n".join(log) + "\n")
say("")
say("wrote step31_matching.csv and step31_matching_report.txt")
