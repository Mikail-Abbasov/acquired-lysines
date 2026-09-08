#!/usr/bin/env python3
"""
step24_cb_reference.py — recompute the MEASURED lysine features from CB, so the union
comparison uses the same reference atom on both sides.

THE ERROR THIS FIXES
Step 21 compared acquired lysines against measured lysines and found enormous differences:
crowding 12 vs 7, atom density 117 vs 75, centroid distances of 2.3-2.5 against 0.58 between
the two measured populations. That was an artefact, and it was mine.

  measured features (steps 6, 6d, 9) were computed from **NZ** — the terminal amine of the
      lysine side chain, roughly 6.4 A out from the backbone, projecting into solvent
  union features (step 21) were computed from **CB** — because at an acquired-lysine position
      the reference residue is a glutamate, valine, etc., and there is no NZ to measure from

CB sits against the protein body; NZ reaches away from it. Counting neighbours within 8 A of
one and comparing against the other measures the difference between two reference points, not
between two populations of lysines. Every feature that showed a difference — crowding, atom
density, hydrophobic and polar contacts — is exactly what that displacement would produce.

THE FIX
Recompute every measured lysine's environment from **CB**, using identical code and cutoffs to
step 21. Both sides then share a backbone-anchored reference and no side chain is assumed on
either side, which is also the honest position given that the union's residues are not lysines.

WHAT IS LOST, and it should be stated
NZ is the atom that does the chemistry, so NZ-anchored geometry is the more meaningful
description of a reactive site. The CB comparison answers a narrower question — is the
*backbone environment* of acquired-lysine positions like that of measured lysine positions —
and that is the only question the union can support. All NZ-based results (A45, A49, A50, A55,
A57) are unaffected: they compare measured lysines with measured lysines, all from NZ.

CHECK BUILT IN
The script reports the CB-vs-NZ difference within the measured set itself. That quantifies the
artefact directly and gives the magnitude to cite when explaining why step 21's first run was
discarded.

INPUT   features_expanded.csv (which lysines), af/*.pdb
OUTPUT  measured_cb.csv, step24_report.txt
"""
import os, glob, json, time
import numpy as np, pandas as pd
import warnings
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser

AF_DIR = os.environ.get("AF_DIR", "af")
CACHE = os.environ.get("CACHE24", "cache24")
os.makedirs(CACHE, exist_ok=True)
T3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
      "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
      "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
HYDRO, AROM, POLAR = set("AVLIMFWPC"), set("FWYH"), set("STNQCY")
ANION, CATION = set("DE"), set("KRH")

log = []
def say(m):
    print(m, flush=True); log.append(m)

say("=" * 78)
say("STEP 24 — MEASURED LYSINES RECOMPUTED FROM CB")
say("so the union comparison uses the same reference atom on both sides")
say("=" * 78)

feat = pd.read_csv("features_expanded.csv", low_memory=False).drop_duplicates(["uniprot", "position"])
by_prot = feat.groupby("uniprot").position.apply(lambda s: sorted(set(int(x) for x in s))).to_dict()
accs = sorted(by_prot)
say("%d measured lysines across %d proteins" % (len(feat), len(accs)))

parser = PDBParser(QUIET=True)
rows, t0 = [], time.time()
for i, acc in enumerate(accs, 1):
    if i % 500 == 0:
        say("  %d/%d proteins | %.0f min" % (i, len(accs), (time.time() - t0) / 60))
    cf = os.path.join(CACHE, "%s.json" % acc)
    if os.path.exists(cf):
        try:
            for r in json.load(open(cf)):
                rows.append(r)
            continue
        except Exception:
            pass
    pdb = os.path.join(AF_DIR, "%s.pdb" % acc)
    if not os.path.exists(pdb):
        continue
    try:
        ch = next(parser.get_structure(acc, pdb)[0].get_chains())
        res = {r.id[1]: r for r in ch if r.id[0] == " "}
    except Exception:
        continue
    heavy = [(a.coord, a.get_parent()) for a in ch.get_atoms() if a.element != "H"]
    if not heavy:
        continue
    co = np.array([h[0] for h in heavy]); par = [h[1] for h in heavy]
    centroid = co.mean(axis=0)
    got = []
    for pos in by_prot[acc]:
        r = res.get(pos)
        if r is None:
            continue
        cb = r["CB"] if "CB" in r else (r["CA"] if "CA" in r else None)
        if cb is None:
            continue
        c = cb.coord
        d = np.linalg.norm(co - c, axis=1)
        near = {par[j] for j in np.where(d <= 8)[0]} - {r}
        aas = [T3.get(x.get_resname(), "X") for x in near]
        pl = [x["CA"].get_bfactor() for x in near if "CA" in x]
        rec = dict(uniprot=acc, position=int(pos),
                   cb_crowding_8A=len(near), cb_atoms_10A=int((d <= 10).sum()),
                   cb_hydrophobic_8A=sum(a in HYDRO for a in aas),
                   cb_aromatic_8A=sum(a in AROM for a in aas),
                   cb_polar_8A=sum(a in POLAR for a in aas),
                   cb_anionic_8A=sum(a in ANION for a in aas),
                   cb_cationic_8A=sum(a in CATION for a in aas),
                   cb_net_charge_8A=sum(a in CATION for a in aas) - sum(a in ANION for a in aas),
                   cb_n_lys_8A=sum(a == "K" for a in aas),
                   cb_dist_centroid=float(np.linalg.norm(c - centroid)),
                   cb_plddt_site=float(r["CA"].get_bfactor()) if "CA" in r else np.nan,
                   cb_plddt_8A=float(np.mean(pl)) if pl else np.nan)
        got.append(rec); rows.append(rec)
    json.dump(got, open(cf, "w"))

cb = pd.DataFrame(rows).drop_duplicates(["uniprot", "position"])
cb.to_csv("measured_cb.csv", index=False)
say("")
say("recomputed %d measured lysines from CB" % len(cb))

m = feat.merge(cb, on=["uniprot", "position"], how="inner")
say("")
say("HOW BIG WAS THE ARTEFACT? — same lysines, NZ-anchored vs CB-anchored")
say("  %-20s %12s %12s %10s" % ("feature", "from NZ", "from CB", "CB - NZ"))
pairs = [("crowding_8A", "cb_crowding_8A"), ("atoms_10A", "cb_atoms_10A"),
         ("hydrophobic_8A", "cb_hydrophobic_8A"), ("aromatic_8A", "cb_aromatic_8A"),
         ("polar_8A", "cb_polar_8A"), ("anionic_8A", "cb_anionic_8A"),
         ("cationic_8A", "cb_cationic_8A"), ("net_charge_8A", "cb_net_charge_8A"),
         ("n_lys_8A", "cb_n_lys_8A"), ("dist_centroid", "cb_dist_centroid")]
for a, b in pairs:
    if a in m and b in m:
        va = pd.to_numeric(m[a], errors="coerce")
        vb = pd.to_numeric(m[b], errors="coerce")
        say("  %-20s %12.2f %12.2f %+10.2f" % (a, va.median(), vb.median(),
                                               vb.median() - va.median()))
say("")
say("  This is the displacement step 21's first run mistook for a biological difference:")
say("  the union was measured from CB, the reference populations from NZ.")
say("  CB sits against the protein body; NZ projects ~6.4 A into solvent.")

open("step24_report.txt", "w").write("\n".join(log) + "\n")
say("")
say("wrote measured_cb.csv and step24_report.txt")
say("NEXT: re-run the union comparison using cb_* features on BOTH sides.")
