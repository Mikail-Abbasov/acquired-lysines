#!/usr/bin/env python3
"""
step32_akt1_site.py — the drugged acquired lysine, measured from the deposited coordinates.

WHY THIS EXISTS
Figure 8 panel e was mostly words. The structure of the one acquired lysine that has been
drugged is public (PDB 8UVY, "STRUCTURE OF AKT1(E17K) WITH COMPOUND 3"), so the panel can show
the measured geometry instead of describing it. Nothing here is copied from the paper's figures;
every number is computed from the deposited coordinates.

WHAT IT EXTRACTS
  - Lys17 side chain (CB-CG-CD-CE-NZ)
  - the covalent ligand XOO and its bond to NZ
  - the zinc ion and everything coordinating it
  - a 2D projection of that assembly, from a PCA of the atoms themselves, so the figure draws a
    real geometry rather than a cartoon

INPUT   xtal/8UVY.pdb  (downloaded by step31)
OUTPUT  step32_akt1_atoms.csv, step32_akt1_bonds.csv, step32_report.txt
"""
import os
import numpy as np, pandas as pd
import warnings
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser

PID, POS, CH = "8UVY", 17, "A"
log = []
def say(s=""):
    print(s, flush=True); log.append(s)

st = PDBParser(QUIET=True).get_structure(PID, os.path.join("xtal", "%s.pdb" % PID))
ch = st[0][CH]
res = {r.id[1]: r for r in ch if r.id[0] == " "}
K = res[POS]
say("=" * 78)
say("STEP 32 - AKT1 E17K, THE DRUGGED ACQUIRED LYSINE, FROM 8UVY")
say("=" * 78)
say("PDB 8UVY chain %s: residue %d is %s" % (CH, POS, K.get_resname()))
assert K.get_resname() == "LYS"

NZ = np.array(K["NZ"].coord)
lig = zn = None
for r in ch:
    if r.get_resname() == "XOO":
        lig = r
    if r.get_resname() == "ZN":
        zn = r
say("ligand present: %s | zinc present: %s" % (lig is not None, zn is not None))

d_lig = min((float(np.linalg.norm(np.array(a.coord) - NZ)), a.get_id()) for a in lig)
say()
say("THE TWO MEASUREMENTS THAT MAKE THE PANEL")
say("  Lys17 NZ to nearest ligand atom (%s) : %.2f A   <- a C-N bond, not a contact"
    % (d_lig[1], d_lig[0]))
znc = np.array(zn["ZN"].coord) if "ZN" in zn else np.array(list(zn)[0].coord)
say("  Lys17 NZ to zinc                     : %.2f A" % float(np.linalg.norm(znc - NZ)))

say()
say("WHAT ELSE COORDINATES THE ZINC (within 2.6 A)")
coord_atoms = []
for r in ch:
    for a in r:
        if a.element == "H" or r is zn:
            continue
        d = float(np.linalg.norm(np.array(a.coord) - znc))
        if d <= 2.6:
            tag = "%s%s" % (r.get_resname(), r.id[1] if r.id[0] == " " else "")
            coord_atoms.append((tag, a.get_id(), d))
            say("  %-10s %-5s %.2f A" % (tag, a.get_id(), d))
if not coord_atoms:
    say("  none within 2.6 A")

# ---- the atoms the panel draws
rows = []
for a in K:
    if a.get_id() in ("CB", "CG", "CD", "CE", "NZ", "CA", "N", "C", "O"):
        rows.append(dict(part="lys17", name=a.get_id(), element=a.element, *[],
                         x=float(a.coord[0]), y=float(a.coord[1]), z=float(a.coord[2])))
for a in lig:
    if a.element != "H":
        rows.append(dict(part="ligand", name=a.get_id(), element=a.element,
                         x=float(a.coord[0]), y=float(a.coord[1]), z=float(a.coord[2])))
rows.append(dict(part="zinc", name="ZN", element="ZN",
                 x=float(znc[0]), y=float(znc[1]), z=float(znc[2])))
for tag, aid, d in coord_atoms:
    r = [rr for rr in ch if "%s%s" % (rr.get_resname(), rr.id[1] if rr.id[0] == " " else "") == tag]
    if r:
        a = r[0][aid]
        rows.append(dict(part="znligand", name="%s %s" % (tag, aid), element=a.element,
                         x=float(a.coord[0]), y=float(a.coord[1]), z=float(a.coord[2])))
A = pd.DataFrame(rows)

# ---- 2D view: PCA of the assembly itself, so the projection is the data's own best plane
X = A[["x", "y", "z"]].values
Xc = X - X.mean(axis=0)
U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
A["u"], A["v"] = Xc @ Vt[0], Xc @ Vt[1]
say()
say("2D PROJECTION: principal plane of the assembly itself")
say("  variance captured by the two drawn axes: %.1f%%" % (100 * (S[0] ** 2 + S[1] ** 2) / (S ** 2).sum()))

# ---- bonds, by distance, within ligand and within the lysine chain
bonds = []
for i in range(len(A)):
    for j in range(i + 1, len(A)):
        pi, pj = A.part[i], A.part[j]
        if "zinc" in (pi, pj) or "znligand" in (pi, pj):
            continue
        d = float(np.linalg.norm(X[i] - X[j]))
        if d <= 1.75:
            bonds.append(dict(i=i, j=j, d=d, kind="covalent"))
A.to_csv("step32_akt1_atoms.csv", index=False)
pd.DataFrame(bonds).to_csv("step32_akt1_bonds.csv", index=False)
say()
say("drawable: %d atoms, %d covalent bonds" % (len(A), len(bonds)))
say()
say("LIMITS")
say("  One chain of one structure. The zinc is modelled by the depositors; this script reports")
say("  the deposited geometry and does not re-refine it. Distances are single measurements, not")
say("  averages, and carry the coordinate error of a %s structure." % "2.0 A-class crystal")
open("step32_report.txt", "w").write("\n".join(log) + "\n")
print("\nwrote step32_akt1_atoms.csv, step32_akt1_bonds.csv, step32_report.txt")
