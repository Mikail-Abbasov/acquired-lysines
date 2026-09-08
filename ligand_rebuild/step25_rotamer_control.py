#!/usr/bin/env python3
"""
step25_rotamer_control.py — the reviewer's question, answered with data: why not model the
mutation in silico and compute NZ features exactly as for the measured lysines?

THE QUESTION
Step 20 justified using environment features only, on the grounds that per-site PROPKA does not
scale (295,411 mutants x 10-60 s = 34-200 days). A reviewer will reasonably reply: then model
the mutation and skip pKa, since step 20 also showed pKa costs only 0.002 AUC. Placing a lysine
rotamer is cheap. Why not do it?

The scaling argument is insufficient. The real question is whether a modelled side chain
recovers the true environment of the amine. That is testable directly, and it does not require
the union at all — it can be answered on lysines whose real NZ position is known.

THE TEST
For every measured lysine we already have the experimental (AlphaFold) NZ position. So:

  1. take the real backbone (N, CA, CB) and DISCARD the side chain beyond CB
  2. rebuild CG, CD, CE, NZ from ideal internal geometry using the most populated lysine
     rotamer (all-trans, chi1-chi4 = 180 deg), which is what a rotamer-based tool places in the
     absence of any other information
  3. compare the modelled NZ position with the real one
  4. compute the 8 A environment at the modelled NZ and compare with the environment at the
     real NZ

If modelled features track real ones closely, modelling is viable and the union analysis should
use it. If they do not, modelling is refuted **with data** rather than argued away on runtime,
and the environment-at-CB approach is justified positively rather than by default.

WHY ALL-TRANS
Lysine's side chain is the most flexible of the twenty; chi1-chi4 give it a large accessible
volume and no rotamer dominates. All-trans is the single most populated state in the Dunbrack
library (~30-40% for chi1), and it is what a tool must fall back on when the true conformation
is unknown — which is exactly the situation at an acquired-lysine position, where the side chain
has never existed. Using the most common rotamer is therefore the FAIREST possible test of
modelling: a real tool doing energy minimisation in context could do better, but it also has no
experimental side chain to minimise toward.

This is a lower bound on modelling error, not an upper bound. State it that way.

INPUT   features_expanded.csv (which lysines), af/*.pdb
OUTPUT  step25_rotamer.csv, step25_report.txt
"""
import os, time
import numpy as np, pandas as pd
import warnings
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser
from scipy import stats as sstats

AF_DIR = os.environ.get("AF_DIR", "af")
N_PROT = int(os.environ.get("N_PROT", 1500))     # sample of proteins; the test does not need all
T3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
      "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
      "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
HYDRO, AROM, POLAR = set("AVLIMFWPC"), set("FWYH"), set("STNQCY")
ANION, CATION = set("DE"), set("KRH")

log = []
def say(m):
    print(m, flush=True); log.append(m)


def place(a, b, c, bond, angle, dihedral):
    """Natural extension reference frame: position an atom given three predecessors,
    a bond length, a bond angle and a dihedral. Standard internal-coordinate construction."""
    ang, dih = np.radians(angle), np.radians(dihedral)
    bc = c - b
    bc /= np.linalg.norm(bc)
    n = np.cross(b - a, bc)
    n /= np.linalg.norm(n)
    m = np.cross(n, bc)
    d = np.array([-bond * np.cos(ang),
                  bond * np.cos(dih) * np.sin(ang),
                  bond * np.sin(dih) * np.sin(ang)])
    return c + d[0] * bc + d[1] * m + d[2] * n


def build_lys(N, CA, CB, chi=(180., 180., 180., 180.)):
    """Rebuild CG, CD, CE, NZ from the backbone using ideal geometry."""
    CG = place(N, CA, CB, 1.520, 114.1, chi[0])
    CD = place(CA, CB, CG, 1.520, 111.3, chi[1])
    CE = place(CB, CG, CD, 1.520, 111.3, chi[2])
    NZ = place(CG, CD, CE, 1.489, 111.9, chi[3])
    return CG, CD, CE, NZ


def env(coords, parents, centre, exclude):
    d = np.linalg.norm(coords - centre, axis=1)
    near = {parents[j] for j in np.where(d <= 8)[0]} - {exclude}
    aas = [T3.get(x.get_resname(), "X") for x in near]
    return dict(crowding_8A=len(near), atoms_10A=int((d <= 10).sum()),
                hydrophobic_8A=sum(a in HYDRO for a in aas),
                aromatic_8A=sum(a in AROM for a in aas),
                polar_8A=sum(a in POLAR for a in aas),
                anionic_8A=sum(a in ANION for a in aas),
                cationic_8A=sum(a in CATION for a in aas),
                net_charge_8A=sum(a in CATION for a in aas) - sum(a in ANION for a in aas),
                n_lys_8A=sum(a == "K" for a in aas))


say("=" * 78)
say("STEP 25 — CAN A MODELLED SIDE CHAIN RECOVER THE REAL NZ ENVIRONMENT?")
say("=" * 78)

feat = pd.read_csv("features_expanded.csv", low_memory=False).drop_duplicates(["uniprot", "position"])
by_prot = feat.groupby("uniprot").position.apply(lambda s: sorted(set(int(x) for x in s))).to_dict()
accs = sorted(by_prot)
rng = np.random.default_rng(0)
if len(accs) > N_PROT:
    accs = list(rng.choice(accs, N_PROT, replace=False))
say("testing on %d proteins (sampled)" % len(accs))

parser = PDBParser(QUIET=True)
rows, t0 = [], time.time()
for i, acc in enumerate(accs, 1):
    if i % 300 == 0:
        say("  %d/%d | %d lysines | %.0f min" % (i, len(accs), len(rows), (time.time() - t0) / 60))
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
    for pos in by_prot[acc]:
        r = res.get(pos)
        if r is None or r.get_resname() != "LYS":
            continue
        if not all(a in r for a in ("N", "CA", "CB", "NZ")):
            continue
        real_nz = np.array(r["NZ"].coord)
        try:
            _, _, _, mod_nz = build_lys(np.array(r["N"].coord), np.array(r["CA"].coord),
                                        np.array(r["CB"].coord))
        except Exception:
            continue
        e_real = env(co, par, real_nz, r)
        e_mod = env(co, par, mod_nz, r)
        rec = dict(uniprot=acc, position=pos,
                   nz_error=float(np.linalg.norm(mod_nz - real_nz)),
                   plddt=float(r["CA"].get_bfactor()) if "CA" in r else np.nan)
        for k in e_real:
            rec["real_" + k] = e_real[k]
            rec["mod_" + k] = e_mod[k]
        rows.append(rec)

d = pd.DataFrame(rows)
d.to_csv("step25_rotamer.csv", index=False)
say("")
say("%d lysines tested" % len(d))
say("")
say("HOW FAR IS THE MODELLED AMINE FROM THE REAL ONE?")
say("  median %.2f A | 25th %.2f | 75th %.2f | 90th %.2f | max %.2f"
    % (d.nz_error.median(), d.nz_error.quantile(.25), d.nz_error.quantile(.75),
       d.nz_error.quantile(.9), d.nz_error.max()))
say("  within 1 A: %.1f%% | within 2 A: %.1f%% | beyond 4 A: %.1f%%"
    % (100 * (d.nz_error <= 1).mean(), 100 * (d.nz_error <= 2).mean(),
       100 * (d.nz_error > 4).mean()))
say("  for scale, the CB-to-NZ distance in a lysine is ~6.4 A")

say("")
say("DO MODELLED FEATURES MATCH REAL ONES?")
say("  %-20s %10s %10s %10s %10s" % ("feature", "real med", "model med", "Pearson r", "Spearman"))
for k in ("crowding_8A", "atoms_10A", "hydrophobic_8A", "aromatic_8A", "polar_8A",
          "anionic_8A", "cationic_8A", "net_charge_8A", "n_lys_8A"):
    a, b = d["real_" + k], d["mod_" + k]
    r_p = np.corrcoef(a, b)[0, 1]
    r_s = sstats.spearmanr(a, b).statistic
    say("  %-20s %10.2f %10.2f %10.3f %10.3f" % (k, a.median(), b.median(), r_p, r_s))

say("")
say("  agreement on the feature that matters most for reactivity (acidic neighbours):")
a, b = d.real_anionic_8A, d.mod_anionic_8A
say("    exact match: %.1f%% of sites | within 1: %.1f%% | mean |difference|: %.2f residues"
    % (100 * (a == b).mean(), 100 * ((a - b).abs() <= 1).mean(), (a - b).abs().mean()))

say("")
say("  does accuracy improve in well-ordered regions?")
for lo, hi, lab in ((0, 70, "pLDDT < 70"), (70, 90, "70-90"), (90, 101, ">= 90")):
    g = d[(d.plddt >= lo) & (d.plddt < hi)]
    if len(g) > 50:
        say("    %-12s n=%6d | median NZ error %.2f A | anionic_8A r = %.3f"
            % (lab, len(g), g.nz_error.median(),
               np.corrcoef(g.real_anionic_8A, g.mod_anionic_8A)[0, 1]))

say("")
say("INTERPRETATION")
say("  This uses the single most populated lysine rotamer (all-trans), which is what any tool")
say("  must fall back on when the true conformation is unknown — the situation at an acquired-")
say("  lysine position, where the side chain has never existed. A tool doing energy")
say("  minimisation in context could do better, so this is a LOWER BOUND on modelling error.")
say("")
say("  If the modelled amine sits far from the real one and the environments correlate weakly,")
say("  then modelling the mutation would not recover the NZ environment it is meant to supply,")
say("  and the environment-at-CB approach is justified POSITIVELY rather than by runtime.")
open("step25_report.txt", "w").write("\n".join(log) + "\n")
say("")
say("wrote step25_rotamer.csv and step25_report.txt")
