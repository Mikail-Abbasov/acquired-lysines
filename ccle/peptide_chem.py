#!/usr/bin/env python3
"""
peptide_chem.py — every peptide species built as a molecule, labelled as the workflow labels it,
and characterised from the resulting structure. Charge included.

CORRECTIONS IN THIS VERSION
  1. TMT is now placed on the peptide N-TERMINUS of every species, including the probe-labelled
     ones. The earlier version left probe peptides entirely bare, which modelled a label-free
     workflow rather than the intended one: probe-label the proteome, enrich, digest, then TMT.
     The probe peptide therefore carries desthiobiotin on the acquired lysine AND TMT on the
     N-terminus and on every other lysine. That is a heavier molecule than modelled before, and
     mass is the property that governs identification.
  2. CHARGE is computed. The first structure-based pass reported mass, logP, TPSA and HBD and
     omitted charge entirely — so the property most sensitive to neutral versus cationic
     labelling was missing from the table that concluded mass was decisive.

HOW CHARGE IS COUNTED, from the structure rather than from a rule
Basic nitrogens are found by substructure at acquisition pH (~2.7): guanidine (Arg), imidazole
(His), primary aliphatic amine (free Lys, peptide N-terminus), and tertiary aliphatic amine
(the TMT reporter). An amide nitrogen is not basic and is not counted — which is exactly the
difference between the two labels:
  TMT           acylates the amine but contributes a tertiary amine of its own -> charge kept
  desthiobiotin acylates the amine and contributes no basic group             -> charge lost
Both fall out of the structure; neither is asserted.

STRUCTURES, supplied by the PI
  desthiobiotin    OC(CCCCC[C@H]([C@H](C)N1)NC1=O)=O            adds 196.1212 Da
  TMT (NHS ester)  CC1N(CC(NCCC(ON2C(CCC2=O)=O)=O)=O)C(C)CCC1   transfers 224.1525 Da
  224.1525 + five heavy-isotope shifts = 229.163, the value in the workflow file. Isotopes do
  not change logP or charge, so the light structure is used for properties.

ORDER OF LABELLING, which matters
  probe workflow: desthiobiotin goes on the acquired lysine FIRST (it labels the intact
  proteome), then TMT goes on whatever amines remain after digestion — the N-terminus and any
  other lysines. Reversing the order would put TMT on the acquired lysine, which the chemistry
  does not do.

USAGE   python3 peptide_chem.py [n_per_species]
"""
import os, sys, time
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, rdMolDescriptors, AllChem
from rdkit.Chem.Descriptors import ExactMolWt
RDLogger.DisableLog("rdApp.*")

DTB = Chem.MolFromSmiles("OC(CCCCC[C@H]([C@H](C)N1)NC1=O)=O")
TMT = Chem.MolFromSmiles("CC1N(CC(NCCC(ON2C(CCC2=O)=O)=O)=O)C(C)CCC1")
RXN_ACID = AllChem.ReactionFromSmarts(
    "[NX3;H2:1][CH2:2].[C:3](=[O:4])[OX2H1]>>[N:1](-[C:3]=[O:4])[CH2:2]")
RXN_NHS = AllChem.ReactionFromSmarts(
    "[NX3;H2:1][CH2:2].[C:3](=[O:4])O[N:5]1[C:6](=O)[C:7][C:8][C:9]1=O"
    ">>[N:1](-[C:3]=[O:4])[CH2:2]")
LYS_NZ = Chem.MolFromSmarts("[NX3;H2][CH2][CH2][CH2][CH2][CH1]")
ANY_PRIMARY = Chem.MolFromSmarts("[NX3;H2][CH2]")

# basic nitrogens at pH ~2.7, found in the structure
BASIC_PATTS = [("guanidine (Arg)", Chem.MolFromSmarts("[NX3][CX3](=[NX2,NX3+])[NX3]")),
               ("imidazole (His)", Chem.MolFromSmarts("c1cnc[nH]1")),
               ("primary amine",   Chem.MolFromSmarts("[NX3;H2;!$(NC=O)]")),
               ("secondary amine", Chem.MolFromSmarts("[NX3;H1;!$(NC=O);!R]")),
               ("tertiary amine",  Chem.MolFromSmarts("[NX3;H0;!$(NC=O);!$(N=*)]([CX4])([CX4])[CX4]"))]

log = []
def say(m):
    print(m, flush=True); log.append(m)


def acylate_once(mol, reagent, rxn, patt):
    """Acylate one amine matching patt. Returns None if there is nothing to acylate."""
    if not mol.GetSubstructMatches(patt):
        return None
    pr = rxn.RunReactants((mol, reagent))
    if not pr:
        return None
    m = pr[0][0]
    try:
        Chem.SanitizeMol(m)
    except Exception:
        return None
    return m


def n_basic(mol):
    """Count protonatable nitrogens at acquisition pH, from the structure."""
    seen, total = set(), 0
    for _, patt in BASIC_PATTS:
        if patt is None:
            continue
        for match in mol.GetSubstructMatches(patt):
            n = [a for a in match if mol.GetAtomWithIdx(a).GetSymbol() == "N"]
            for a in n:
                if a not in seen:
                    seen.add(a); total += 1
    return total


def characterise(seq, workflow):
    """workflow: 'tmt' (digest then TMT) or 'probe' (probe first, then TMT on what is left)."""
    m = Chem.MolFromSequence(seq)
    if m is None:
        return None
    if workflow == "probe":
        # desthiobiotin on the acquired lysine first — the probe labels the intact proteome
        r = acylate_once(m, DTB, RXN_ACID, LYS_NZ)
        if r is None:
            return None
        m = r
    if workflow in ("tmt", "probe"):
        # TMT on every amine that remains: other lysines and the peptide N-terminus
        for _ in range(8):
            r = acylate_once(m, TMT, RXN_NHS, ANY_PRIMARY)
            if r is None:
                break
            m = r
    try:
        return dict(mw=ExactMolWt(m), logp=Crippen.MolLogP(m),
                    tpsa=rdMolDescriptors.CalcTPSA(m), basic=n_basic(m),
                    hbd=rdMolDescriptors.CalcNumHBD(m), heavy=m.GetNumHeavyAtoms())
    except Exception:
        return None


if __name__ == "__main__":
    N_PER = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    P = pd.read_csv("peptide_properties.csv", low_memory=False)
    P = P[P.peptide.str.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", na=False)]

    say("=" * 78)
    say("STRUCTURE-BASED PROPERTIES, WITH N-TERMINAL TMT AND CHARGE")
    say("=" * 78)
    demo = "LAGEGLKR"
    for lab, wf in (("unlabelled", None), ("TMT workflow", "tmt"), ("probe workflow", "probe")):
        r = characterise(demo, wf)
        say("  %-16s %-9s mass %8.2f  logP %+6.2f  basic N %d"
            % (demo, lab, r["mw"], r["logp"], r["basic"]))
    say("")
    say("  the probe peptide keeps TMT on the N-terminus and other lysines;")
    say("  desthiobiotin occupies the acquired lysine and contributes no basic nitrogen")
    say("")
    say("sampling %d peptides per species; all identified ones included" % N_PER)

    JOBS = [("mutant cleaved", "tmt", "mutant, cleaved (TMT)"),
            ("mutant missed cleavage", "tmt", "mutant, missed cleavage (TMT)"),
            ("downstream", "tmt", "downstream (TMT)"),
            ("mutant missed cleavage", "probe", "probe + TMT (enrichment)")]
    rows, t0 = [], time.time()
    for sp, wf, lab in JOBS:
        s = P[P.species == sp]
        keep = pd.concat([s[s.found],
                          s[~s.found].sample(min(N_PER, int((~s.found).sum())), random_state=0)])
        say("")
        say("  %s: %d peptides" % (lab, len(keep)))
        for i, r in enumerate(keep.itertuples(), 1):
            c = characterise(r.peptide, wf)
            if c:
                rows.append(dict(species=lab, workflow=wf, peptide=r.peptide, found=r.found,
                                 L=len(r.peptide), n_lys=r.peptide.count("K"), **c))
            if i % 250 == 0:
                say("    %d/%d  (%.0fs)" % (i, len(keep), time.time() - t0))

    D = pd.DataFrame(rows)
    D.to_csv("peptide_chem.csv", index=False)

    say("")
    say("MEDIAN PROPERTIES BY SPECIES")
    say("  %-32s %6s %8s %8s %8s %7s %7s" % ("", "n", "mass", "logP", "TPSA", "basic", "length"))
    for _, g in D.groupby("species", sort=False):
        say("  %-32s %6d %8.0f %+8.2f %8.0f %7.1f %7.0f"
            % (g.species.iloc[0], len(g), g.mw.median(), g.logp.median(), g.tpsa.median(),
               g.basic.median(), g.L.median()))

    say("")
    say("IDENTIFIED vs NOT, in standard deviations")
    for sp, g in D.groupby("species", sort=False):
        a, b = g[g.found], g[~g.found]
        if len(a) < 5:
            continue
        say("  %s   (%d identified)" % (sp, len(a)))
        for c, nm in (("L", "length"), ("mw", "mass"), ("tpsa", "polar surface"),
                      ("logp", "logP"), ("basic", "basic N")):
            sd = g[c].std()
            say("    %-14s %+6.2f SD   (%.1f vs %.1f)"
                % (nm, (a[c].median() - b[c].median()) / sd if sd else np.nan,
                   a[c].median(), b[c].median()))

    open("peptide_chem_report.txt", "w").write("\n".join(log) + "\n")
    say("")
    say("wrote peptide_chem.csv and peptide_chem_report.txt  (%.0fs)" % (time.time() - t0))
