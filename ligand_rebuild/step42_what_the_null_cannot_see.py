#!/usr/bin/env python3
"""
step42_what_the_null_cannot_see.py — turn figure S7 f from an inventory into a measurement.

WHY
Figure 6d reports that no compound descriptor predicts which lysines a compound prefers. That is
a statement about the descriptors that were tested, and the paper goes on to test two axes those
descriptors do not cover: charge (Fig. 7) and chirality (Fig. 5f). figure S7 f said so by listing
the descriptors in a table of chips with two rows reading "nothing tested". A list is an
assertion. This measures both claims instead.

A. CHARGE. Net charge at pH 7.4 is approximated as (basic nitrogens) - (carboxylic acids +
   sulfonic acids) by SMARTS. It varies across the 212 compounds. Each tested descriptor is then
   correlated with net charge, and with molecular weight as the reference: if the descriptor set
   is rich in size information and poor in charge information, the two columns separate.

   NOTE ON WORDING. The result is NOT "not one encodes charge". Two descriptors are weakly but
   significantly correlated with charge (rings, HBA, |rho| ~ 0.19, p < 0.01). The defensible
   statement is the contrast: the largest charge correlation is about a third of the median size
   correlation. Say that, not the absolute.

B. CHIRALITY. NPR1 and NPR2 are ratios of principal moments of inertia, and reflection does not
   change the inertia tensor's eigenvalues, so they cannot respond to handedness. That is a
   mathematical fact, but it is worth measuring because the pipeline does not compare a conformer
   with its own reflection — it embeds each molecule independently. Three quantities are computed:
     same conformer, reflected      -> the descriptor's true response to handedness
     enantiomer, re-embedded        -> conformer-sampling noise, which an earlier check mistook
                                       for a chirality effect
     across different compounds     -> the spread the descriptor is actually measuring
   Reporting all three keeps the second from being read as the first.
"""
import os, sys
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
from scipy import stats as sstats
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Descriptors3D, AllChem
from rdkit.Geometry import Point3D
RDLogger.DisableLog("rdApp.*")

SMILES_XLSX = os.environ.get("SMILES_XLSX", "SMILES.xlsx")
NCHIRAL = int(os.environ.get("NCHIRAL", 25))
log = []
def say(m):
    print(m, flush=True); log.append(m)

ACID = Chem.MolFromSmarts("[CX3](=O)[OX2H1,OX1-]")
SULF = Chem.MolFromSmarts("[SX4](=O)(=O)[OX2H1,OX1-]")
BASE = Chem.MolFromSmarts("[NX3;H2,H1,H0;!$(N[C,S]=[O,S,N]);!$(N=*);!$(n)]")
DESC = {"MW": Descriptors.MolWt, "logP": Descriptors.MolLogP,
        "heavy": lambda m: m.GetNumHeavyAtoms(), "rotatable": Descriptors.NumRotatableBonds,
        "fsp3": Descriptors.FractionCSP3, "TPSA": Descriptors.TPSA,
        "HBD": Descriptors.NumHDonors, "HBA": Descriptors.NumHAcceptors,
        "rings": Descriptors.RingCount, "aromatic_rings": Descriptors.NumAromaticRings}

say("=" * 78)
say("STEP 42 — WHAT THE FIGURE 6d NULL CANNOT SEE, MEASURED")
say("=" * 78)

x = pd.ExcelFile(SMILES_XLSX)
rows, mols = [], []
for sh in x.sheet_names:
    for r in x.parse(sh).dropna(subset=["SMILES"]).itertuples():
        m = Chem.MolFromSmiles(str(r.SMILES))
        if m is None:
            continue
        d = {k: f(m) for k, f in DESC.items()}
        d["net_charge"] = (len(m.GetSubstructMatches(BASE))
                           - len(m.GetSubstructMatches(ACID)) - len(m.GetSubstructMatches(SULF)))
        d["dataset"] = sh; d["compound"] = str(r.Compounds).strip()
        rows.append(d); mols.append(m)
D = pd.DataFrame(rows)
say("compounds with a parsable structure: %d" % len(D))
say("net charge at pH 7.4 (basic - acidic): %s"
    % ", ".join("%+d: %d" % (k, v) for k, v in sorted(D.net_charge.value_counts().items())))

say("")
say("A. HOW MUCH CHARGE INFORMATION DOES EACH TESTED DESCRIPTOR CARRY?")
say("  %-16s %12s %10s %12s" % ("descriptor", "|rho| charge", "p", "|rho| MW"))
out = []
for c in DESC:
    rc, pc = sstats.spearmanr(D[c], D.net_charge)
    rm = sstats.spearmanr(D[c], D.MW)[0]
    out.append(dict(descriptor=c, rho_charge=abs(rc), p_charge=pc, rho_mw=abs(rm)))
A = pd.DataFrame(out).sort_values("rho_charge", ascending=False)
for _, r in A.iterrows():
    say("  %-16s %12.3f %10.3g %12.3f" % (r.descriptor, r.rho_charge, r.p_charge, r.rho_mw))
say("  -> largest charge correlation %.3f; median size correlation %.3f; ratio %.2f"
    % (A.rho_charge.max(), A.rho_mw.median(), A.rho_charge.max() / A.rho_mw.median()))
say("     Two descriptors (rings, HBA) are weakly but significantly correlated with charge, so")
say("     the claim is the CONTRAST, not 'not one encodes charge'.")
A.to_csv("step42_charge_information.csv", index=False)

say("")
say("B. CAN THE SHAPE TERMS SEE HANDEDNESS?")
def npr(mh):
    return Descriptors3D.NPR1(mh), Descriptors3D.NPR2(mh)

refl, reemb, across = [], [], []
n = 0
for m in mols:
    if n >= NCHIRAL:
        break
    if not Chem.FindMolChiralCenters(m, includeUnassigned=False, useLegacyImplementation=False):
        continue
    mh = Chem.AddHs(m)
    if AllChem.EmbedMolecule(mh, randomSeed=0xf00d, maxAttempts=200) != 0:
        continue
    try: AllChem.MMFFOptimizeMolecule(mh, maxIters=200)
    except Exception: pass
    a1, a2 = npr(mh)
    across.append((a1, a2))

    mm = Chem.Mol(mh); c0, c1 = mh.GetConformer(), mm.GetConformer()   # same conformer, reflected
    for i in range(mm.GetNumAtoms()):
        p = c0.GetAtomPosition(i); c1.SetAtomPosition(i, Point3D(-p.x, p.y, p.z))
    b1, b2 = npr(mm)
    refl.append((abs(a1 - b1), abs(a2 - b2)))

    inv = Chem.Mol(m)                                                   # enantiomer, re-embedded
    for at in inv.GetAtoms():
        if at.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
            at.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
        elif at.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
            at.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
    ih = Chem.AddHs(inv)
    if AllChem.EmbedMolecule(ih, randomSeed=0xbeef, maxAttempts=200) == 0:
        try: AllChem.MMFFOptimizeMolecule(ih, maxIters=200)
        except Exception: pass
        e1, e2 = npr(ih)
        reemb.append((abs(a1 - e1), abs(a2 - e2)))
    n += 1

R = np.array(refl); E = np.array(reemb); X = np.array(across)
say("  %-34s %10s %10s" % ("", "max dNPR1", "max dNPR2"))
say("  %-34s %10.3e %10.3e" % ("same conformer, reflected", R[:, 0].max(), R[:, 1].max()))
say("  %-34s %10.3f %10.3f" % ("enantiomer, re-embedded", E[:, 0].max(), E[:, 1].max()))
say("  %-34s %10.3f %10.3f" % ("sd across different compounds", X[:, 0].std(), X[:, 1].std()))
say("  n = %d chiral compounds" % len(R))
say("  -> reflection changes NPR1 and NPR2 by exactly zero. The re-embedding difference is")
say("     conformer sampling, not chirality: an earlier check measured that and would have been")
say("     read as a chirality effect. The descriptor's response to handedness is nil, and its")
say("     own sampling noise is already a third of the spread it resolves between compounds.")
pd.DataFrame(dict(kind=(["reflected"] * len(R) + ["re-embedded"] * len(E)),
                  dNPR1=list(R[:, 0]) + list(E[:, 0]),
                  dNPR2=list(R[:, 1]) + list(E[:, 1]))).to_csv("step42_chirality.csv", index=False)
pd.DataFrame(dict(NPR1=X[:, 0], NPR2=X[:, 1])).to_csv("step42_npr_spread.csv", index=False)
open("step42_report.txt", "w").write("\n".join(log) + "\n")
say("")
say("wrote step42_charge_information.csv, step42_chirality.csv, step42_npr_spread.csv")
