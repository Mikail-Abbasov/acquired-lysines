#!/usr/bin/env python3
"""
step21_union_characterise.py — place the acquired-lysine union in the structural space learned
from measured chemoproteomics.

WHAT THIS IS, AND WHAT IT IS NOT
This CHARACTERISES the union — it does not prioritise it. The composite reaches AUC 0.617
cross-validated by protein, which sits at the ceiling that label reproducibility permits
(Part 1 estimated 0.55-0.65) but still means the engaged and unengaged distributions overlap
heavily. A score at this level describes a population; it cannot identify individual targetable
sites. Any output framed as a ranked target list would contradict A47 (structure plus predicted
reactivity explain under 1% of engagement-rate variance) and A48 (the decision not to prioritise
the union computationally). Say "these sites resemble / do not resemble measured engaged
lysines", never "these sites are ligandable".

WHY NO MUTATION MODELLING — the justification that must appear in the methods
An acquired lysine is not a lysine in the reference structure: position 545 of PIK3CA is a
glutamate, with no NZ atom, no lysine side chain, and no residue for PROPKA. Modelling each
mutation would restore those features but does not scale — PROPKA runs per STRUCTURE, so
295,411 mutants at 10-60 s each is 34-200 days, and any repacking tool adds an error term that
cannot be quantified at that scale.

Step 20 measured what those features are actually worth. Cross-validated by protein:

    FULL         23 features   AUC 0.6223  [0.5945, 0.6428]
    NO_PKA       22 features   AUC 0.6203  [0.6001, 0.6449]
    ENVIRONMENT  17 features   AUC 0.6163  [0.5945, 0.6434]

Dropping predicted pKa costs 0.002 AUC; dropping every residue-dependent feature costs 0.006,
with the paired 95% CI spanning zero. **The signal is carried by the surroundings, and the
surroundings do not change when one residue mutates.** So the analysis uses environment features
computed at the reference position, no mutation is modelled, and there is no rotamer choice or
repacking error to defend. What is lost — predicted pKa and NZ-specific geometry — is stated
rather than approximated.

Residual assumption, to state: the mutation is taken not to perturb local folding. Nothing here
tests that, and it is the honest limit of the approach.

DESIGN
  1. AlphaFold structures for the union's genes (UniProt accessions resolved from gene symbols)
  2. pLDDT >= 70 gate at the site, per A48 — below that the geometry is unreliable and step 11
     showed shape features collapse there
  3. fpocket once per protein
  4. 17 environment features at the reference position
  5. the step 20 ENVIRONMENT model, trained on the measured chemoproteomic data, applied to
     score how much each union site RESEMBLES a measured engaged lysine
  6. distribution comparison against the measured engaged and unengaged sets — the output is a
     distribution, not a ranked list

INPUT   union_sites.tsv, features_expanded.csv, occ_long.csv
OUTPUT  union_scored.csv, step21_report.txt, step21_methods.txt
"""
import os, sys, re, glob, json, time, shutil, subprocess, tempfile, urllib.request
import numpy as np, pandas as pd
import warnings
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

UNION = os.environ.get("UNION", "union_sites.tsv")
AF_DIR = os.environ.get("AF_DIR2", "af_union")
CACHE = os.environ.get("CACHE21", "cache21")
PLDDT_GATE = float(os.environ.get("PLDDT_GATE", 70))
TIER = int(os.environ.get("TIER", 60))
AF_VERSIONS = (6, 5, 4)
AF_URL = "https://alphafold.ebi.ac.uk/files/AF-%s-F1-model_v%d.pdb"
ENV_FEATS = ["crowding_8A", "atoms_10A", "hydrophobic_8A", "aromatic_8A", "polar_8A",
             "anionic_8A", "cationic_8A", "net_charge_8A", "n_lys_8A", "dist_centroid",
             "in_pocket", "pocket_volume", "pocket_score", "pocket_drugscore",
             "dist_nearest_pocket", "plddt_site", "plddt_8A"]
T3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
      "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
      "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
HYDRO, AROM, POLAR = set("AVLIMFWPC"), set("FWYH"), set("STNQCY")
ANION, CATION = set("DE"), set("KRH")
os.makedirs(AF_DIR, exist_ok=True); os.makedirs(CACHE, exist_ok=True)

log = []
def say(m):
    print(m, flush=True); log.append(m)

say("=" * 78)
say("STEP 21 — CHARACTERISING THE ACQUIRED-LYSINE UNION")
say("characterisation, NOT prioritisation (see docstring and A48)")
say("=" * 78)

u = pd.read_csv(UNION, sep="\t")
say("union: %d sites, %d genes" % (len(u), u.gene.nunique()))

# ---- gene -> UniProt, from the measured data where possible, else UniProt REST
map_file = os.path.join(CACHE, "gene2acc.json")
g2a = json.load(open(map_file)) if os.path.exists(map_file) else {}
if not g2a:
    say("resolving gene symbols to UniProt accessions...")
    src = pd.read_csv("features_expanded.csv", low_memory=False)
    if "gene" in src.columns:
        g2a.update(src.dropna(subset=["gene"]).groupby("gene").uniprot.first().to_dict())
    missing = sorted(set(u.gene) - set(g2a))
    say("  %d from existing data, %d to fetch from UniProt" % (len(g2a), len(missing)))
    for i in range(0, len(missing), 100):
        chunk = missing[i:i + 100]
        q = "+OR+".join("gene_exact:%s" % g for g in chunk)
        url = ("https://rest.uniprot.org/uniprotkb/search?query=(%s)+AND+"
               "organism_id:9606+AND+reviewed:true&fields=accession,gene_primary"
               "&format=tsv&size=500" % q)
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                for line in r.read().decode().splitlines()[1:]:
                    p = line.split("\t")
                    if len(p) >= 2 and p[1] and p[1] not in g2a:
                        g2a[p[1]] = p[0]
        except Exception:
            pass
        if i % 1000 == 0:
            say("    %d/%d" % (i, len(missing)))
    json.dump(g2a, open(map_file, "w"))
u["uniprot"] = u.gene.map(g2a)
say("  resolved %d of %d sites (%.1f%%)" % (u.uniprot.notna().sum(), len(u),
                                            100 * u.uniprot.notna().mean()))
u = u[u.uniprot.notna()]

def fetch(acc):
    p = os.path.join(AF_DIR, "%s.pdb" % acc)
    if os.path.exists(p) and os.path.getsize(p) > 1000:
        return p
    for v in AF_VERSIONS:
        try:
            urllib.request.urlretrieve(AF_URL % (acc, v), p)
            if os.path.getsize(p) > 1000:
                return p
        except Exception:
            pass
        if os.path.exists(p):
            os.remove(p)
    return None


def pockets(acc, pdb):
    cf = os.path.join(CACHE, "%s_pk.json" % acc)
    if os.path.exists(cf):
        try:
            return {int(k): v for k, v in json.load(open(cf)).items()}
        except Exception:
            pass
    out, work = {}, tempfile.mkdtemp()
    dst = os.path.join(work, "%s.pdb" % acc); shutil.copy(pdb, dst)
    try:
        subprocess.run(["fpocket", "-f", dst], capture_output=True, timeout=600)
        odir = dst[:-4] + "_out"; info = os.path.join(odir, "%s_info.txt" % acc)
        pk, cur = [], None
        if os.path.exists(info):
            for line in open(info):
                m = re.match(r"Pocket\s+(\d+)", line)
                if m:
                    cur = {"id": int(m.group(1))}; pk.append(cur)
                elif cur is not None and ":" in line:
                    k, v = [x.strip() for x in line.split(":", 1)]
                    try: cur[k] = float(v)
                    except ValueError: pass
        cen = {}
        for f in glob.glob(os.path.join(odir, "pockets", "pocket*_vert.pqr")):
            pid = int(re.search(r"pocket(\d+)_", os.path.basename(f)).group(1)) + 1
            pts = [[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                   for l in open(f) if l.startswith(("ATOM", "HETATM"))]
            if pts: cen[pid] = np.array(pts)
        out = {"_pockets": [{k: v for k, v in p.items() if k != "centres"} for p in pk],
               "_centres": {str(k): v.tolist() for k, v in cen.items()}}
    except Exception:
        pass
    shutil.rmtree(work, ignore_errors=True)
    json.dump(out, open(cf, "w"))
    return out


parser, sr = PDBParser(QUIET=True), ShrakeRupley()
rows, t0, miss = [], time.time(), 0
by_prot = u.groupby("uniprot").position.apply(list).to_dict()
accs = sorted(by_prot)
say("")
say("computing environment features at %d positions across %d proteins" % (len(u), len(accs)))
for i, acc in enumerate(accs, 1):
    if i % 250 == 0:
        say("  %d/%d proteins | %d sites | %.0f min" % (i, len(accs), len(rows),
                                                        (time.time() - t0) / 60))
    pdb = fetch(acc)
    if not pdb:
        miss += 1; continue
    try:
        ch = next(parser.get_structure(acc, pdb)[0].get_chains())
        res = {r.id[1]: r for r in ch if r.id[0] == " "}
    except Exception:
        miss += 1; continue
    pk = pockets(acc, pdb)
    plist = pk.get("_pockets", []) if isinstance(pk, dict) else []
    cen = {int(k): np.array(v) for k, v in (pk.get("_centres", {}) or {}).items()} if isinstance(pk, dict) else {}
    heavy = [(a.coord, a.get_parent()) for a in ch.get_atoms() if a.element != "H"]
    if not heavy:
        continue
    co = np.array([h[0] for h in heavy]); par = [h[1] for h in heavy]
    centroid = co.mean(axis=0)
    for pos in set(int(p) for p in by_prot[acc]):
        r = res.get(pos)
        if r is None:
            continue
        pl = r["CA"].get_bfactor() if "CA" in r else np.nan
        if not (pl == pl) or pl < PLDDT_GATE:
            continue
        c = r["CB"].coord if "CB" in r else (r["CA"].coord if "CA" in r else None)
        if c is None:
            continue
        d = np.linalg.norm(co - c, axis=1)
        near = {par[j] for j in np.where(d <= 8)[0]} - {r}
        aas = [T3.get(x.get_resname(), "X") for x in near]
        pls = [x["CA"].get_bfactor() for x in near if "CA" in x]
        best, dmin = None, np.inf
        for p in plist:
            cc = cen.get(p.get("id"))
            if cc is None: continue
            dd = float(np.min(np.linalg.norm(cc - c, axis=1)))
            if dd < dmin: dmin, best = dd, p
        rows.append(dict(uniprot=acc, position=pos,
            crowding_8A=len(near), atoms_10A=int((d <= 10).sum()),
            hydrophobic_8A=sum(a in HYDRO for a in aas),
            aromatic_8A=sum(a in AROM for a in aas),
            polar_8A=sum(a in POLAR for a in aas),
            anionic_8A=sum(a in ANION for a in aas),
            cationic_8A=sum(a in CATION for a in aas),
            net_charge_8A=sum(a in CATION for a in aas) - sum(a in ANION for a in aas),
            n_lys_8A=sum(a == "K" for a in aas),
            dist_centroid=float(np.linalg.norm(c - centroid)),
            in_pocket=float(np.isfinite(dmin) and dmin <= 5.0),
            dist_nearest_pocket=float(dmin) if np.isfinite(dmin) else np.nan,
            pocket_volume=best.get("Volume") if best else np.nan,
            pocket_score=best.get("Score") if best else np.nan,
            pocket_drugscore=best.get("Druggability Score") if best else np.nan,
            plddt_site=float(pl), plddt_8A=float(np.mean(pls)) if pls else np.nan))

U = pd.DataFrame(rows)
say("")
say("union sites with usable geometry: %d (%.1f%% of resolved)" % (len(U), 100 * len(U) / max(len(u), 1)))
say("  structures unavailable for %d proteins" % miss)
say("  sites below pLDDT %.0f were excluded at source (A48: geometry unreliable there)" % PLDDT_GATE)

# ---- train the ENVIRONMENT model on the measured data
feat = pd.read_csv("features_expanded.csv", low_memory=False).drop_duplicates(["uniprot", "position"])
if "in_pocket" in feat:
    feat["in_pocket"] = feat.in_pocket.astype(str).str.lower().isin(["true", "1", "1.0"]).astype(float)
occ = pd.read_csv("occ_long.csv", low_memory=False)
eng = set(map(tuple, occ[occ.occupancy_pct >= TIER][["uniprot", "position"]].drop_duplicates().values))
feat["engaged"] = [(u_, int(p)) in eng for u_, p in zip(feat.uniprot, feat.position)]
F = [f for f in ENV_FEATS if f in feat.columns and f in U.columns]
tr = feat[["engaged"] + F].dropna()
sc = StandardScaler().fit(tr[F].astype(float))
mdl = LogisticRegression(max_iter=3000).fit(sc.transform(tr[F].astype(float)), tr.engaged.astype(int))
say("")
say("model trained on %d measured lysines (%d engaged), %d environment features"
    % (len(tr), int(tr.engaged.sum()), len(F)))

Uc = U.dropna(subset=F).copy()
Uc["score"] = mdl.decision_function(sc.transform(Uc[F].astype(float)))
tr2 = tr.copy(); tr2["score"] = mdl.decision_function(sc.transform(tr2[F].astype(float)))
Uc.to_csv("union_scored.csv", index=False)

say("")
say("DISTRIBUTION COMPARISON — the output is a distribution, not a ranked list")
say("  %-28s %8s %9s %9s %9s" % ("set", "n", "median", "25th", "75th"))
for lab, v in (("union acquired lysines", Uc.score),
               ("measured, engaged", tr2.loc[tr2.engaged, "score"]),
               ("measured, not engaged", tr2.loc[~tr2.engaged, "score"])):
    say("  %-28s %8d %9.3f %9.3f %9.3f" % (lab, len(v), v.median(), v.quantile(.25), v.quantile(.75)))

thr = tr2.loc[tr2.engaged, "score"].median()
say("")
say("  union sites scoring above the median of MEASURED ENGAGED lysines: %d (%.1f%%)"
    % ((Uc.score > thr).sum(), 100 * (Uc.score > thr).mean()))
say("  for reference, %.1f%% of measured UNENGAGED lysines also exceed that threshold"
    % (100 * (tr2.loc[~tr2.engaged, "score"] > thr).mean()))
say("  -> the overlap is the point: this separates distributions, it does not identify sites")

m = []
m.append("Characterisation of the acquired-lysine union. Structural features were computed at "
         "each acquired-lysine position in the AlphaFold model of the reference protein, "
         "restricted to positions with pLDDT >= %.0f." % PLDDT_GATE)
m.append("An acquired lysine is not a lysine in the reference structure, so features requiring "
         "the lysine side chain (predicted pKa, NZ solvent accessibility, side-chain SASA, "
         "relative accessibility, hydrogen-bond partners at NZ, salt-bridge geometry) cannot be "
         "computed without modelling the mutation. Modelling each of the 295,411 mutations was "
         "not attempted: PROPKA operates per structure, making per-site pKa prediction "
         "infeasible at this scale, and side-chain repacking would introduce an error term that "
         "could not be quantified across so many sites.")
m.append("Instead we measured what those features contribute. In the chemoproteomic training "
         "data, a composite trained on all 23 features reached a cross-validated AUC of 0.622 "
         "(95%% CI 0.595-0.643, cross-validated by protein); removing predicted pKa cost 0.002 "
         "AUC, and removing all residue-dependent features cost 0.006 (paired 95%% CI spanning "
         "zero). The discriminating signal therefore resides in the residue's surroundings, "
         "which are unchanged by substitution at the site itself, and the analysis proceeds "
         "with the %d environment features alone. Predicted pKa and NZ-specific geometry are "
         "unavailable for the union and are not estimated." % len(F))
m.append("The analysis assumes that the substitution does not perturb local folding; this is "
         "not tested and is a limitation of the approach.")
m.append("This is a characterisation, not a prioritisation. The composite's cross-validated AUC "
         "of ~0.62 lies at the ceiling implied by the reproducibility of the underlying "
         "engagement measurements, and the engaged and unengaged distributions overlap "
         "substantially. Scores describe how closely a population of sites resembles measured "
         "engaged lysines; they do not identify individual ligandable sites, and are not used "
         "to rank targets.")
open("step21_methods.txt", "w").write("\n\n".join(m) + "\n")
open("step21_report.txt", "w").write("\n".join(log) + "\n")
say("")
say("wrote union_scored.csv, step21_methods.txt, step21_report.txt")
