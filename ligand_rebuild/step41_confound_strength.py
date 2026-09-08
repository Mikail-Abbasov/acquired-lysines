#!/usr/bin/env python3
"""
step41_confound_strength.py — does a descriptor's POWER TO SEPARATE THE THREE SCREENS predict
how many pooled false positives it produced?

WHY THIS EXISTS
Figures 6f and S7a asserted the confound by naming a category: "every one of the 21 pooled hits
involves a SIZE descriptor". Two problems with that.
  1. It requires calling rotatable-bond count and fsp3 "size", which they are not. fsp3 is a
     saturation term.
  2. It is a statement about step 31 (10 descriptors x 15 site features = 150 pairs). The panel
     it was used to justify, Fig. 6d, plots step 37 Q2 (6 shape descriptors x 14 preference
     features = 84 pairs). Those two analyses both happen to return 21 significant pairs, which
     is a coincidence, and step 37's 21 include SIX from logP — not a size descriptor at all.

This replaces the category claim with a measurement. For each compound descriptor, compute how
much of its variance is explained by which screen the compound came from (eta-squared across the
three datasets), and correlate that against the number of significant pairs that descriptor
produced in the pooled step 31 test. If the confound is what generated the positives, descriptors
that separate the datasets should be the descriptors that produced hits.

RESULT  Spearman rho = 0.69, p = 0.027 across the ten descriptors. MW (eta2 0.255) and heavy-atom
count (0.252) are the two strongest separators and produced 7 hits each; rotatable bonds (0.134)
produced 5. TPSA separates moderately (0.124) and produced none, so the relationship is real but
not deterministic — which is the honest way to state it.
"""
import os, sys
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
from scipy import stats as sstats
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
RDLogger.DisableLog('rdApp.*')

SMILES_XLSX = os.environ.get("SMILES_XLSX", "SMILES.xlsx")
log = []
def say(m):
    print(m, flush=True); log.append(m)

DESC = {
    "MW": Descriptors.MolWt, "logP": Descriptors.MolLogP,
    "heavy": lambda m: m.GetNumHeavyAtoms(), "rotatable": Descriptors.NumRotatableBonds,
    "fsp3": Descriptors.FractionCSP3, "TPSA": Descriptors.TPSA,
    "HBD": Descriptors.NumHDonors, "HBA": Descriptors.NumHAcceptors,
    "rings": Descriptors.RingCount, "aromatic_rings": Descriptors.NumAromaticRings,
}

x = pd.ExcelFile(SMILES_XLSX); rows = []
for sh, ds in (('NatChem','NatChem'), ('Meroterpenoids','Meroterpenoid'), ('CyMPSA','CyMPSA')):
    if sh not in x.sheet_names: continue
    for r in x.parse(sh).dropna(subset=['SMILES']).itertuples():
        m = Chem.MolFromSmiles(str(r.SMILES))
        if m is None: continue
        rows.append(dict(dataset=ds, compound=str(r.Compounds).strip(),
                         **{k: f(m) for k, f in DESC.items()}))
D = pd.DataFrame(rows)
say("=" * 78)
say("STEP 41 — DOES DATASET-SEPARATING POWER PREDICT POOLED FALSE POSITIVES?")
say("=" * 78)
say("compounds with a parsable structure: %d  (%s)"
    % (len(D), ", ".join("%s %d" % (k, v) for k, v in D.dataset.value_counts().items())))

M = pd.read_csv("step31_matching.csv")
hits = M[M.q < 0.05].groupby("descriptor").size()
say("step 31 pooled: %d of %d pairs significant" % ((M.q < 0.05).sum(), len(M)))

def eta2(col):
    g = [v.values for _, v in D.groupby("dataset")[col]]
    gm = np.concatenate(g).mean()
    ssb = sum(len(a) * (a.mean() - gm) ** 2 for a in g)
    sst = ((np.concatenate(g) - gm) ** 2).sum()
    return ssb / sst

out = []
for c in DESC:
    kw = sstats.kruskal(*[v.values for _, v in D.groupby("dataset")[c]])
    out.append(dict(descriptor=c, eta2=eta2(c), kw_p=kw[1], hits=int(hits.get(c, 0)),
                    n_pairs=int((M.descriptor == c).sum())))
R = pd.DataFrame(out).sort_values("eta2", ascending=False)
say("")
say("  %-16s %8s %11s %7s %8s" % ("descriptor", "eta2", "KW p", "hits", "of pairs"))
for _, r in R.iterrows():
    say("  %-16s %8.3f %11.2g %7d %8d" % (r.descriptor, r.eta2, r.kw_p, r.hits, r.n_pairs))

rho, p = sstats.spearmanr(R.eta2, R.hits)
say("")
say("  -> Spearman rho %.3f, p %.4f  (n = %d descriptors)" % (rho, p, len(R)))
say("     Descriptors that separate the three screens are the descriptors that produced")
say("     pooled positives. TPSA is the exception: it separates (eta2 %.3f) and produced none."
    % float(R.loc[R.descriptor == "TPSA", "eta2"].iloc[0]))
R.to_csv("step41_confound_strength.csv", index=False)
open("step41_report.txt", "w").write("\n".join(log) + "\n")
say("")
say("wrote step41_confound_strength.csv and step41_report.txt")
