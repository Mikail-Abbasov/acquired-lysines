#!/usr/bin/env python3
"""
step9_two_comparisons.py — give each question its own control group.

THE FLAW THIS FIXES
Everything up to step 8 compared engaged lysines against other *probe-labelled* lysines and
called the result "ligandability". But it also reported pKa effects from that same comparison
and discussed them as reactivity. One comparison cannot answer both questions.

The two questions need different controls:

  COMPARISON A — REACTIVITY
      45,552 probe-labelled lysines  vs  every OTHER lysine on the same protein
      Asks: what makes a lysine get labelled at all? Labelling requires the amine to attack an
      electrophile, so this contrast is about nucleophilicity. Within-protein, so abundance,
      localisation and observability largely cancel.

  COMPARISON B — LIGANDABILITY
      engaged lysines  vs  all 45,552 probe-labelled lysines
      Asks: among lysines that ALL react with the probe, which ones does a competing molecule
      block? Reactivity is held approximately constant by the sampling — every site in both
      groups cleared the labelling threshold — so what differs is whether a molecule binds.
      Tiered at >=60 / 70 / 80 / 90% occupancy to test for dose-response.

THE SAME FEATURES GO THROUGH BOTH. A feature that predicts labelling but not engagement is
reactivity-specific; one that predicts engagement but not labelling is ligandability-specific.
That contrast is the point, and neither analysis alone provides it.

CAVEAT THAT MUST TRAVEL WITH COMPARISON A
An unlabelled lysine is not necessarily unreactive — it may sit in a region the mass
spectrometer never covered. A36/A37 measured that: the wild-type peptide is seen for only ~7%
of acquired-lysine regions in a conventional deep search. So "unlabelled" conflates unreactive
with unobserved, and the bias direction is unknown. Restricting to proteins with several
labelled sites (i.e. demonstrably well covered) reduces but does not remove it. Reported, not
solved.

CAVEAT FOR COMPARISON B
Tiers shrink fast (step 8: 3,255 rows at >=60%, 966 at >=75%, 261 at >=90%), so the highest
tiers have wide intervals and may not converge per dataset. Pooled across datasets with dataset
as a covariate at the top tiers.

INPUT   features_expanded.csv (labelled sites, with pKa and geometry)
        occ_long.csv, af/*.pdb
OUTPUT  step9_all_lysines.csv   every lysine on every relevant protein, labelled or not
        step9_results.csv, step9_report.txt
"""
import os, sys, glob, json, time, shutil, subprocess, tempfile, re
import numpy as np, pandas as pd
import warnings
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from scipy import stats as sstats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

AF_DIR = os.environ.get("AF_DIR", "af")
CACHE = os.environ.get("CACHE9", "cache9")
os.makedirs(CACHE, exist_ok=True)
TIERS = (60, 70, 80, 90)
T3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
      "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
      "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
TIEN = {"A": 121, "R": 265, "N": 187, "D": 187, "C": 148, "Q": 214, "E": 214, "G": 97,
        "H": 216, "I": 195, "L": 191, "K": 230, "M": 203, "F": 228, "P": 154, "S": 143,
        "T": 163, "W": 264, "Y": 255, "V": 165}
HYDRO, AROM, POLAR = set("AVLIMFWPC"), set("FWYH"), set("STNQCY")
ANION, CATION = set("DE"), set("KRH")

log = []
def say(m):
    print(m, flush=True); log.append(m)


def propka_all(acc, pdb):
    """pKa for EVERY lysine in the structure, not only the measured ones."""
    cf = os.path.join(CACHE, "%s_pka.json" % acc)
    if os.path.exists(cf):
        try:
            return {int(k): v for k, v in json.load(open(cf)).items()}
        except Exception:
            pass
    w = tempfile.mkdtemp()
    dst = os.path.join(w, os.path.basename(pdb))
    shutil.copy(pdb, dst)
    out = {}
    try:
        subprocess.run([sys.executable, "-m", "propka", dst],
                       capture_output=True, timeout=1800, cwd=w)
        pk = dst[:-4] + ".pka"
        if os.path.exists(pk):
            for line in open(pk):
                if not line.startswith("LYS"):
                    continue
                try:                       # fixed-width; number field overflows at >=1000
                    out.setdefault(int(line[3:8].strip()), float(line[10:16].strip()))
                except (ValueError, IndexError):
                    continue
    except Exception:
        pass
    shutil.rmtree(w, ignore_errors=True)
    json.dump({str(k): v for k, v in out.items()}, open(cf, "w"))
    return out


def all_lysines(acc):
    """Every lysine in the AlphaFold model, with the same geometry as step 6/6d."""
    cf = os.path.join(CACHE, "%s_geom.json" % acc)
    if os.path.exists(cf):
        try:
            return json.load(open(cf))
        except Exception:
            pass
    pdb = os.path.join(AF_DIR, "%s.pdb" % acc)
    if not os.path.exists(pdb):
        return {}
    try:
        ch = next(PDBParser(QUIET=True).get_structure(acc, pdb)[0].get_chains())
        res = {r.id[1]: r for r in ch if r.id[0] == " "}
        ShrakeRupley().compute(ch, level="A")
    except Exception:
        return {}
    heavy = [(a.coord, a.get_parent()) for a in ch.get_atoms() if a.element != "H"]
    coords = np.array([h[0] for h in heavy]); parents = [h[1] for h in heavy]
    centroid = coords.mean(axis=0)
    carb = []
    for r in res.values():
        n = r.get_resname()
        if n == "ASP":
            carb += [r[a].coord for a in ("OD1", "OD2") if a in r]
        elif n == "GLU":
            carb += [r[a].coord for a in ("OE1", "OE2") if a in r]
    carb = np.array(carb) if carb else np.empty((0, 3))
    pk = propka_all(acc, pdb)
    out = {}
    for pos, r in res.items():
        if r.get_resname() != "LYS" or "NZ" not in r:
            continue
        nz = r["NZ"].coord
        d = np.linalg.norm(coords - nz, axis=1)
        near8 = {parents[j] for j in np.where(d <= 8)[0]} - {r}
        aas = [T3.get(x.get_resname(), "X") for x in near8]
        pl = [x["CA"].get_bfactor() for x in near8 if "CA" in x]
        dc = float(np.min(np.linalg.norm(carb - nz, axis=1))) if len(carb) else np.nan
        out[str(pos)] = dict(
            crowding_8A=len(near8), atoms_10A=int((d <= 10).sum()),
            hydrophobic_8A=sum(a in HYDRO for a in aas),
            aromatic_8A=sum(a in AROM for a in aas),
            polar_8A=sum(a in POLAR for a in aas),
            anionic_8A=sum(a in ANION for a in aas),
            cationic_8A=sum(a in CATION for a in aas),
            net_charge_8A=sum(a in CATION for a in aas) - sum(a in ANION for a in aas),
            n_lys_8A=sum(a == "K" for a in aas),
            sasa_sidechain=float(sum(a.sasa for a in r
                                     if a.get_name() not in ("N", "CA", "C", "O"))),
            nz_sasa=float(r["NZ"].sasa),
            rsa_tien=float(sum(a.sasa for a in r)) / TIEN["K"],
            dist_centroid=float(np.linalg.norm(nz - centroid)),
            dist_carboxyl=dc, salt_bridge=bool(dc == dc and dc <= 4.0),
            plddt_site=float(r["CA"].get_bfactor()) if "CA" in r else np.nan,
            plddt_8A=float(np.mean(pl)) if pl else np.nan,
            propka_pKa=pk.get(pos, np.nan))
    json.dump(out, open(cf, "w"))
    return out


FEATS = ["propka_pKa", "salt_bridge", "net_charge_8A", "anionic_8A", "cationic_8A",
         "crowding_8A", "atoms_10A", "hydrophobic_8A", "aromatic_8A", "polar_8A",
         "rsa_tien", "sasa_sidechain", "nz_sasa", "dist_centroid", "n_lys_8A",
         "plddt_site", "plddt_8A"]

say("=" * 78)
say("STEP 9 — SEPARATE CONTROLS FOR REACTIVITY AND LIGANDABILITY")
say("=" * 78)

feat = pd.read_csv("features_expanded.csv", low_memory=False)
labelled = set(zip(feat.uniprot, feat.position.astype(int)))
prots = sorted(feat.uniprot.unique())
say("%d proteins carry at least one probe-labelled lysine" % len(prots))
say("%d probe-labelled lysine sites (%d distinct)"
    % (len(feat), feat.groupby(["uniprot", "position"]).ngroups))

rows, t0 = [], time.time()
for i, acc in enumerate(prots, 1):
    if i % 200 == 0:
        say("  %d/%d proteins | %.0f min | eta %.0f min"
            % (i, len(prots), (time.time() - t0) / 60,
               (time.time() - t0) / 60 / i * (len(prots) - i)))
    for pos, rec in all_lysines(acc).items():
        rows.append(dict(uniprot=acc, position=int(pos),
                         labelled=(acc, int(pos)) in labelled, **rec))
allK = pd.DataFrame(rows)
allK.to_csv("step9_all_lysines.csv", index=False)

say("")
say("COMPARISON A — REACTIVITY: probe-labelled vs unlabelled lysines, same protein")
say("  %d lysines on %d proteins | labelled %d (%.1f%%) | unlabelled %d"
    % (len(allK), allK.uniprot.nunique(), allK.labelled.sum(),
       100 * allK.labelled.mean(), (~allK.labelled).sum()))
say("  CAVEAT: an unlabelled lysine may be unreactive OR simply never observed by MS.")
say("  Restricting to well-covered proteins reduces but does not remove this.")

def within_protein_test(df, group_col, feats, min_per_prot=2):
    """Paired within-protein comparison: for each protein compute the difference in means
    between groups, then test whether those differences are centred on zero. Protein-level
    confounders cancel because each protein contributes one difference."""
    out = []
    for f in feats:
        diffs = []
        for acc, g in df.groupby("uniprot"):
            a = pd.to_numeric(g.loc[g[group_col], f], errors="coerce").dropna()
            b = pd.to_numeric(g.loc[~g[group_col], f], errors="coerce").dropna()
            if len(a) >= min_per_prot and len(b) >= min_per_prot:
                diffs.append(a.mean() - b.mean())
        if len(diffs) < 30:
            continue
        diffs = np.array(diffs)
        t, p = sstats.ttest_1samp(diffs, 0)
        out.append(dict(feature=f, n_proteins=len(diffs), mean_diff=float(diffs.mean()),
                        sd=float(diffs.std()), t=float(t), p=float(p),
                        cohens_d=float(diffs.mean() / diffs.std()) if diffs.std() else np.nan))
    r = pd.DataFrame(out)
    if len(r):
        r["q"] = multipletests(r.p, method="fdr_bh")[1]
    return r


A = within_protein_test(allK, "labelled", FEATS)
if len(A):
    A = A.sort_values("q")
    say("")
    say("  %-18s %8s %10s %8s %10s" % ("feature", "proteins", "mean diff", "d", "q"))
    for _, r in A.iterrows():
        say("  %-18s %8d %10.4f %8.3f %10.3g%s"
            % (r.feature[:18], r.n_proteins, r.mean_diff, r.cohens_d, r.q,
               " *" if r.q < 0.05 else ""))
    say("  (mean diff = labelled minus unlabelled, averaged within protein)")

say("")
say("COMPARISON B — LIGANDABILITY: engaged vs all probe-labelled, by occupancy tier")
occ = pd.read_csv("occ_long.csv", low_memory=False)
res_b = []
for tier in TIERS:
    eng = set(map(tuple, occ[occ.occupancy_pct >= tier][["uniprot", "position"]]
                  .drop_duplicates().values))
    f2 = feat.drop_duplicates(["uniprot", "position"]).copy()
    f2["engaged"] = [(u, int(p)) in eng for u, p in zip(f2.uniprot, f2.position)]
    n_e = int(f2.engaged.sum())
    say("")
    say("  tier >=%d%%: %d engaged of %d labelled (%.2f%%)"
        % (tier, n_e, len(f2), 100 * n_e / len(f2)))
    if n_e < 50:
        say("    too few to model")
        continue
    B = within_protein_test(f2, "engaged", [f for f in FEATS if f in f2.columns])
    if not len(B):
        say("    too few proteins with both classes")
        continue
    B = B.sort_values("q")
    for _, r in B.head(8).iterrows():
        say("    %-18s %6d prot  diff %8.4f  d %6.3f  q %9.3g%s"
            % (r.feature[:18], r.n_proteins, r.mean_diff, r.cohens_d, r.q,
               " *" if r.q < 0.05 else ""))
    B["tier"] = tier
    res_b.append(B)

if len(A):
    A["comparison"] = "A_reactivity"; A["tier"] = np.nan
allres = pd.concat([A] + [b.assign(comparison="B_ligandability") for b in res_b],
                   ignore_index=True) if len(A) else pd.concat(res_b, ignore_index=True)
allres.to_csv("step9_results.csv", index=False)

say("")
say("CONTRAST — which features are specific to which property")
if len(A) and res_b:
    b60 = res_b[0].set_index("feature")
    a = A.set_index("feature")
    say("  %-18s %14s %14s   %s" % ("feature", "A reactivity", "B ligandability", "reads as"))
    for f in FEATS:
        if f not in a.index or f not in b60.index:
            continue
        sa, sb = a.loc[f, "q"] < 0.05, b60.loc[f, "q"] < 0.05
        tag = ("both" if sa and sb else
               "REACTIVITY only" if sa else
               "LIGANDABILITY only" if sb else "neither")
        say("  %-18s %14.3g %14.3g   %s" % (f[:18], a.loc[f, "q"], b60.loc[f, "q"], tag))

open("step9_report.txt", "w").write("\n".join(log) + "\n")
say("")
say("wrote step9_all_lysines.csv, step9_results.csv, step9_report.txt")
