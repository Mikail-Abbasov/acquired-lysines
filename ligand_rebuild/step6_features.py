#!/usr/bin/env python3
"""
step6_features.py — download AlphaFold models and compute per-lysine structural features.

BUILT FROM SCRATCH. Reuses no previously generated structure or feature file.

TWO JOBS, and the first matters more than the second
----------------------------------------------------
1. VERIFY RESIDUE NUMBERING. Steps 2-3 showed the three datasets agree with each other, but
   agreement is not correctness — they could share an upstream numbering source and be
   consistently wrong. The AlphaFold model carries the UniProt canonical sequence, so the
   test is direct: is the residue at the reported position actually a lysine? A site that
   fails is either mis-numbered or on a non-canonical isoform, and must not enter the model.
   This check is reported before any feature is used.

2. COMPUTE FEATURES around the lysine NZ atom — the nitrogen that does the chemistry. All
   geometry is measured from NZ rather than from CA or the residue centroid, because
   reachability of NZ is the physical question.

FEATURES
  crowding_8A          residues with any heavy atom within 8 A of NZ
  atoms_10A            heavy atoms within 10 A of NZ
  hydrophobic_8A       AVLIMFWPC residues within 8 A
  aromatic_8A          FWYH residues within 8 A
  polar_8A             STNQCY residues within 8 A
  anionic_8A / cationic_8A / net_charge_8A     DE / KRH within 8 A
  sasa_total           solvent accessible surface area of the whole residue (Shrake-Rupley)
  sasa_sidechain       SASA of side-chain atoms only
  rsa                  sasa_total / 200.0, the Tien et al. theoretical maximum for lysine
  nz_sasa              SASA of the NZ atom alone — the most direct measure of reachability
  dist_centroid        NZ to protein centre of mass
  dist_nearest_cys     NZ to nearest cysteine SG
  plddt_site           AlphaFold confidence at this residue (B-factor column)
  plddt_8A             mean confidence of residues within 8 A — local order
  n_lys_8A             other lysines within 8 A (competition for the same probe)

NO FEATURE IS INTERPRETED HERE. This step produces measurements; step 5's model consumes them.

CAVEATS THAT MUST TRAVEL WITH THE OUTPUT
  - AlphaFold omits proteins above ~2,700 residues. Those failures are non-random, biased
    against giant repeat and coiled-coil proteins (AHNAK, plectin, desmoplakin, dystonin,
    DNA-PK). Reported, not silently dropped.
  - A predicted structure is not an experimental one. Low-pLDDT regions are disordered or
    uncertain; plddt_site and plddt_8A are provided so this can be conditioned on rather
    than assumed away.
  - SASA is computed on the isolated chain. Sites buried by a binding partner in vivo will
    appear exposed here.

USAGE
    python3 step6_features.py            # download + compute
    AF_DIR=/path/to/af python3 step6_features.py    # reuse an existing structure directory
"""
import os, sys, gzip, time, math, urllib.request, urllib.error
import numpy as np, pandas as pd
import warnings
warnings.filterwarnings("ignore")

from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

AF_DIR = os.environ.get("AF_DIR", "af")
# AlphaFold DB retired v4 (all v4 URLs now 404) and serves v6. Versions are tried newest
# first so this does not silently return zero structures the next time EBI increments.
AF_VERSIONS = (6, 5, 4)
AF_URL = "https://alphafold.ebi.ac.uk/files/AF-%s-F1-model_v%d.pdb"
HYDRO = set("AVLIMFWPC")
AROM = set("FWYH")
POLAR = set("STNQCY")
ANION = set("DE")
CATION = set("KRH")
THREE2ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
             "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
             "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
             "TYR": "Y", "VAL": "V"}

os.makedirs(AF_DIR, exist_ok=True)
log = []
def say(m):
    print(m, flush=True); log.append(m)


def fetch(acc):
    """Return path to the AlphaFold model, downloading if absent. None if unavailable."""
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


sites = pd.read_csv("sites.csv", low_memory=False)
# Numbering convention differs BETWEEN DATASETS for the same protein. Verified on beta-actin
# (P60709): NatChem is UniProt-canonical at 13/13 sites, while both JACS datasets sit one
# residue higher (their 51/62/192/316/327/329 = canonical 50/61/191/315/326/328). Pooling
# datasets before choosing an offset therefore guarantees that neither convention resolves.
# Offsets are resolved per (dataset, protein).
want = (sites.groupby(["dataset", "uniprot", "position"]).size().reset_index()
             [["dataset", "uniprot", "position"]])
by_ds_prot = want.groupby(["dataset", "uniprot"]).position.apply(list).to_dict()
by_prot = want.groupby("uniprot").position.apply(list).to_dict()
accs = sorted(by_prot)
say("=" * 74)
say("STEP 6 — ALPHAFOLD STRUCTURES AND PER-SITE GEOMETRY")
say("=" * 74)
say("%d proteins, %d sites" % (len(accs), len(want)))

parser = PDBParser(QUIET=True)
sr = ShrakeRupley()
rows, no_struct, t0 = [], [], time.time()
offsets = {}          # uniprot -> (offset, n_ok_at_offset, n_ok_at_zero, n_sites, verdict)
OFFSET_RANGE = range(-3, 4)


def read_chain(acc, path):
    st = parser.get_structure(acc, path)
    chain = next(st[0].get_chains())
    return chain, {r.id[1]: r for r in chain if r.id[0] == " "}


def resolve_offset(res, positions):
    """Choose the position offset that makes the most reported sites lysines.

    Lysine is 5.7% of residues, so a single site can match a wrong offset by luck; testing
    7 offsets on one site would give a ~30% false-correction rate. The offset is therefore
    scored against ALL of the protein's sites at once and accepted only if it strictly
    beats offset 0. A protein whose sites are resolved by different offsets is not
    systematically shifted - those sites are genuinely mis-numbered and are dropped.
    """
    def n_lys(off):
        return sum(1 for p in positions
                   if THREE2ONE.get(res[p + off].get_resname(), "X") == "K"
                   if (p + off) in res)
    score = {o: sum(1 for p in positions
                    if (p + o) in res and
                    THREE2ONE.get(res[p + o].get_resname(), "X") == "K")
             for o in OFFSET_RANGE}
    z = score[0]
    best = max(score, key=lambda o: (score[o], -abs(o)))
    n = len(positions)
    if score[best] <= z:
        return 0, z, z, n, ("clean" if z == n else "partial_no_offset_helps")
    if score[best] == n:
        return best, score[best], z, n, "systematic_offset_all_resolved"
    return best, score[best], z, n, "offset_improves_but_incomplete"


for i, acc in enumerate(accs, 1):
    if i % 250 == 0:
        say("  %d/%d proteins | %d sites done | %.0f min elapsed"
            % (i, len(accs), len(rows), (time.time() - t0) / 60))
    path = fetch(acc)
    if not path:
        no_struct.append(acc)
        continue
    try:
        chain, res = read_chain(acc, path)
        sr.compute(chain, level="A")
    except Exception:
        no_struct.append(acc)
        continue

    heavy = [(a.coord, a.get_parent()) for a in chain.get_atoms() if a.element != "H"]
    coords = np.array([h[0] for h in heavy])
    parents = [h[1] for h in heavy]
    centroid = coords.mean(axis=0)
    cys_sg = np.array([r["SG"].coord for r in res.values()
                       if r.get_resname() == "CYS" and "SG" in r]) if res else np.empty((0, 3))

    for ds in [d for (d, a) in by_ds_prot if a == acc]:
        positions = [int(p) for p in by_ds_prot[(ds, acc)]]
        off, n_ok, n_zero, n_tot, verdict = resolve_offset(res, positions)
        offsets[(ds, acc)] = (off, n_ok, n_zero, n_tot, verdict)

        for pos in positions:
            rec = {"dataset": ds, "uniprot": acc, "position": int(pos),
                   "offset_applied": off, "offset_verdict": verdict}
            r = res.get(int(pos) + off)
            if r is None:
                rec["qc"] = "position_absent"; rows.append(rec); continue
            aa = THREE2ONE.get(r.get_resname(), "X")
            rec["af_residue"] = aa
            if aa != "K":
                rec["qc"] = "not_lysine"; rows.append(rec); continue
            if "NZ" not in r:
                rec["qc"] = "no_NZ"; rows.append(rec); continue

            nz = r["NZ"].coord
            d = np.linalg.norm(coords - nz, axis=1)
            near8 = {parents[j] for j in np.where(d <= 8)[0]} - {r}
            aas = [THREE2ONE.get(x.get_resname(), "X") for x in near8]
            plddts = [x["CA"].get_bfactor() for x in near8 if "CA" in x]

            rec.update(
                qc="ok",
                crowding_8A=len(near8),
                atoms_10A=int((d <= 10).sum()),
                hydrophobic_8A=sum(a in HYDRO for a in aas),
                aromatic_8A=sum(a in AROM for a in aas),
                polar_8A=sum(a in POLAR for a in aas),
                anionic_8A=sum(a in ANION for a in aas),
                cationic_8A=sum(a in CATION for a in aas),
                net_charge_8A=sum(a in CATION for a in aas) - sum(a in ANION for a in aas),
                n_lys_8A=sum(a == "K" for a in aas),
                sasa_total=float(sum(a.sasa for a in r)),
                sasa_sidechain=float(sum(a.sasa for a in r
                                         if a.get_name() not in ("N", "CA", "C", "O"))),
                nz_sasa=float(r["NZ"].sasa),
                rsa=float(sum(a.sasa for a in r)) / 200.0,
                dist_centroid=float(np.linalg.norm(nz - centroid)),
                dist_nearest_cys=float(np.min(np.linalg.norm(cys_sg - nz, axis=1)))
                                 if len(cys_sg) else np.nan,
                plddt_site=float(r["CA"].get_bfactor()) if "CA" in r else np.nan,
                plddt_8A=float(np.mean(plddts)) if plddts else np.nan,
            )
            rows.append(rec)

feat = pd.DataFrame(rows)
feat.to_csv("features.csv", index=False)

od = pd.DataFrame([(k[0], k[1]) + v for k, v in offsets.items()],
                  columns=["dataset", "uniprot", "offset", "n_lys_at_offset",
                           "n_lys_at_zero", "n_sites", "verdict"])
od.to_csv("numbering_offsets.csv", index=False)
say("")
say("NUMBERING CONVENTION PER PROTEIN")
say("  %d (dataset, protein) combinations with a structure" % len(od))
say("  offset by dataset:")
for ds, g in od.groupby("dataset"):
    say("    %-14s %s" % (ds, g.offset.value_counts().sort_index().to_dict()))
for v, g in od.groupby("verdict"):
    say("    %-34s %5d proteins, %6d sites" % (v, len(g), g.n_sites.sum()))
shifted = od[od.offset != 0]
if len(shifted):
    say("  proteins needing an offset: %d (%.1f%%)" % (len(shifted), 100*len(shifted)/len(od)))
    say("    offsets used: %s" % shifted.offset.value_counts().to_dict())
    say("    A -1 offset is mature-protein numbering: the initiator methionine is cleaved,")
    say("    so literature positions run one lower than UniProt canonical (beta-actin is the")
    say("    classic case). Applied only where it resolves MORE of the protein's sites than")
    say("    offset 0, evaluated across all sites at once.")
    say("    examples: %s" % ", ".join("%s%+d (%d/%d)" % (r.uniprot, r.offset,
        r.n_lys_at_offset, r.n_sites) for r in shifted.head(8).itertuples()))
say("")
say("NUMBERING VERIFICATION — the test steps 2-3 could not perform")
n = len(feat)
for k, lbl in [("ok", "position IS a lysine in the AlphaFold sequence"),
               ("not_lysine", "position is NOT a lysine  <-- MIS-NUMBERED or non-canonical isoform"),
               ("position_absent", "position beyond the modelled sequence"),
               ("no_NZ", "lysine present but side chain incomplete")]:
    c = int((feat.get("qc") == k).sum())
    say("  %-52s %6d  %5.2f%%" % (lbl, c, 100 * c / max(n, 1)))
bad = feat[feat.qc == "not_lysine"]
if len(bad):
    say("  observed residues at mis-numbered positions: %s"
        % bad.af_residue.value_counts().head(6).to_dict())
    say("  examples: %s" % ", ".join("%s:%d=%s" % (r.uniprot, r.position, r.af_residue)
                                     for r in bad.head(6).itertuples()))

say("")
say("STRUCTURES UNAVAILABLE: %d of %d proteins (%.1f%%)"
    % (len(no_struct), len(accs), 100 * len(no_struct) / max(len(accs), 1)))
if no_struct:
    say("  AlphaFold DB omits proteins above ~2,700 residues; this loss is NON-RANDOM,")
    say("  biased against giant repeat and coiled-coil proteins. Report as a limitation.")
    say("  first 12: %s" % ", ".join(no_struct[:12]))
    pd.Series(no_struct).to_csv("no_structure.csv", index=False, header=["uniprot"])

ok = feat[feat.qc == "ok"]
say("")
say("USABLE SITES: %d" % len(ok))
if len(ok):
    say("%-20s %8s %8s %8s" % ("feature", "median", "10th", "90th"))
    for c in ["crowding_8A", "aromatic_8A", "hydrophobic_8A", "rsa", "nz_sasa",
              "dist_nearest_cys", "plddt_site"]:
        if c in ok:
            q = ok[c].dropna()
            if len(q):
                say("%-20s %8.2f %8.2f %8.2f"
                    % (c, q.median(), q.quantile(.1), q.quantile(.9)))

open("step6_report.txt", "w").write("\n".join(log) + "\n")
say("")
say("wrote features.csv (%d rows) and step6_report.txt" % len(feat))
say("Sites with qc != 'ok' MUST be excluded from modelling — see the numbering table above.")
