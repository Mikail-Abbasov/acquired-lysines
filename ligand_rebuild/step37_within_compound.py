#!/usr/bin/env python3
"""
step37_within_compound.py — does a compound's SHAPE determine WHICH of the lysines it labels it
also ligands?

WHY THIS DESIGN, AFTER FOUR NULLS
  step 31  compound descriptors vs mean engaged-site environment -> null; the pooled positives
           were a dataset confound (MW encodes which paper a compound came from)
  step 32  compound 3D shape vs site size/enclosure -> null, and the site side had no shape
           descriptor, so shape complementarity was never actually tested
  step 34  cationic environment vs number of liganding compounds -> null both ways
  step 36  compound volume/shape vs fpocket cavity geometry -> invalid. Occupancy fraction came
           out at 0.91 (real complexes are ~0.4-0.6) and a quarter of pairs had the compound
           LARGER than its pocket. fpocket cavities here are small surface depressions, not
           binding sites, so the test could not be performed rather than answered.

This asks the question without any cavity definition, and WITHIN compound rather than across
compounds — which removes the confound that killed step 31.

THE DESIGN
For each compound separately:
    LIGANDED      lysines this compound blocked to >=60% occupancy
    LABELLED-ONLY lysines this compound was TESTED against and did not block
The comparison group is the compound's own tested sites, so screening coverage cannot leak in.

Each compound yields a vector of standardised differences — its "site preference". Then:
    Q1  Do compounds have a consistent site preference at all? (are the vectors non-zero)
    Q2  Does the DIRECTION or MAGNITUDE of that preference depend on compound shape?
        This is the second-order question. Step 9 already showed that pooled across compounds,
        liganded sites are more buried and crowded than labelled ones. The question here is
        whether rod-like and globular compounds differ in HOW they prefer.

EXPECTATION TO SET IN ADVANCE
Q1 will very likely be positive and is largely a re-derivation of step 9's Comparison B per
compound. Q2 is the new question, it is second-order, and it requires the per-compound
difference vectors to be stable — which needs compounds with enough liganded sites.
"""
import os, re
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
from scipy import stats as sstats
from statsmodels.stats.multitest import multipletests
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors3D, Descriptors
RDLogger.DisableLog('rdApp.*')

SMILES_XLSX=os.environ.get("SMILES_XLSX","SMILES.xlsx")
TIER=int(os.environ.get("TIER",60)); GATE=float(os.environ.get("PLDDT_GATE",70))
MIN_LIG=int(os.environ.get("MIN_LIG",15))
SF=["rsa_tien","nz_sasa","sasa_sidechain","crowding_8A","atoms_10A","hydrophobic_8A",
    "aromatic_8A","polar_8A","anionic_8A","cationic_8A","net_charge_8A","n_lys_8A",
    "dist_centroid","propka_pKa"]
log=[]
def say(m):
    print(m,flush=True); log.append(m)

say("="*78); say("STEP 37 — WITHIN-COMPOUND SITE PREFERENCE vs COMPOUND SHAPE"); say("="*78)

def shape(smi):
    m=Chem.MolFromSmiles(str(smi))
    if m is None: return None
    mh=Chem.AddHs(m)
    ids=AllChem.EmbedMultipleConfs(mh,numConfs=3,randomSeed=0xf00d,maxAttempts=200)
    if len(ids)==0: return None
    try: AllChem.MMFFOptimizeMoleculeConfs(mh,maxIters=200)
    except Exception: pass
    n1=[];n2=[];v=[]
    for c in ids:
        try:
            n1.append(Descriptors3D.NPR1(mh,confId=c)); n2.append(Descriptors3D.NPR2(mh,confId=c))
            v.append(AllChem.ComputeMolVolume(mh,confId=c))
        except Exception: pass
    if not n1: return None
    return dict(NPR1=np.mean(n1),NPR2=np.mean(n2),volume=np.mean(v),
                rodness=np.mean(n2)-np.mean(n1),MW=Descriptors.MolWt(m),
                logP=Descriptors.MolLogP(m))

x=pd.ExcelFile(SMILES_XLSX); rows=[]
for sh,ds in (('NatChem','NatChem'),('Meroterpenoids','Meroterpenoid'),('CyMPSA','CyMPSA')):
    if sh not in x.sheet_names: continue
    for r in x.parse(sh).dropna(subset=['SMILES']).itertuples():
        s=shape(r.SMILES)
        if s: rows.append(dict(dataset=ds,compound=str(r.Compounds).strip(),**s))
C=pd.DataFrame(rows); say("compounds with shape: %d"%len(C))

feat=pd.read_csv('features_expanded.csv',low_memory=False).drop_duplicates(['uniprot','position'])
feat=feat[feat.plddt_site>=GATE]
F=[f for f in SF if f in feat.columns]
fidx=feat.set_index(['uniprot','position'])[F]
occ=pd.read_csv('occ_long.csv',low_memory=False)
occ['base']=[re.match(r'^(\d+[a-z])',str(c)).group(1) if (d=='NatChem' and re.match(r'^(\d+[a-z])',str(c))) else str(c) for c,d in zip(occ.compound,occ.dataset)]

recs=[]
for (ds,base),g in occ.groupby(['dataset','base']):
    lig=g[g.occupancy_pct>=TIER][['uniprot','position']].drop_duplicates()
    non=g[g.occupancy_pct<TIER][['uniprot','position']].drop_duplicates()
    if len(lig)<MIN_LIG or len(non)<50: continue
    L=fidx.reindex(pd.MultiIndex.from_frame(lig)).dropna(how='all')
    N=fidx.reindex(pd.MultiIndex.from_frame(non)).dropna(how='all')
    if len(L)<MIN_LIG or len(N)<50: continue
    rec=dict(dataset=ds,compound=base,n_lig=len(L),n_non=len(N))
    for f in F:
        a=pd.to_numeric(L[f],errors='coerce').dropna(); b=pd.to_numeric(N[f],errors='coerce').dropna()
        if len(a)>=10 and len(b)>=30 and b.std()>0:
            rec['p_'+f]=float((a.mean()-b.mean())/b.std())
    recs.append(rec)
E=pd.DataFrame(recs)
say("compounds with >=%d liganded sites above the gate: %d"%(MIN_LIG,len(E)))
if len(E)<15: say("too few"); raise SystemExit(0)

say("")
say("Q1 — DO COMPOUNDS HAVE A CONSISTENT SITE PREFERENCE?")
say("  (mean standardised shift, liganded vs labelled-only, tested against zero across compounds)")
PC=[c for c in E.columns if c.startswith('p_')]
o=[]
for c in PC:
    v=E[c].dropna()
    if len(v)<15: continue
    t,p=sstats.ttest_1samp(v,0)
    o.append(dict(feature=c[2:],mean_shift=v.mean(),sd=v.std(),n=len(v),p=p))
oo=pd.DataFrame(o); oo['q']=multipletests(oo.p,method='fdr_bh')[1]
oo=oo.reindex(oo.mean_shift.abs().sort_values(ascending=False).index)
say("  %-18s %10s %8s %5s %10s"%("feature","mean shift","sd","n","q"))
for _,r in oo.iterrows():
    say("  %-18s %10.3f %8.3f %5d %10.3g%s"%(r.feature[:18],r.mean_shift,r.sd,r.n,r.q," *" if r.q<0.05 else ""))
say("  -> %d of %d features show a consistent preference"%(int((oo.q<0.05).sum()),len(oo)))

say("")
say("Q2 — DOES THAT PREFERENCE DEPEND ON COMPOUND SHAPE?")
m=E.merge(C,on=['dataset','compound'],how='inner')
say("  %d compounds matched to a structure"%len(m))
DS=['rodness','NPR1','NPR2','volume','MW','logP']
o2=[]
for d in DS:
    for c in PC:
        s=m[[d,c]].dropna()
        if len(s)<15 or s[d].std()==0: continue
        r,p=sstats.spearmanr(s[d],s[c]); o2.append(dict(shape=d,pref=c[2:],rho=r,p=p,n=len(s)))
if o2:
    oo2=pd.DataFrame(o2); oo2['q']=multipletests(oo2.p,method='fdr_bh')[1]
    oo2=oo2.reindex(oo2.rho.abs().sort_values(ascending=False).index)
    say("  %-12s %-18s %5s %8s %10s"%("compound","site preference","n","rho","q"))
    for _,r in oo2.head(12).iterrows():
        say("  %-12s %-18s %5d %8.3f %10.3g%s"%(r['shape'],r.pref[:18],r.n,r.rho,r.q," *" if r.q<0.05 else ""))
    say("  -> %d of %d significant (%.1f expected by chance)"%(int((oo2.q<0.05).sum()),len(oo2),0.05*len(oo2)))
    # within dataset
    say("")
    say("  WITHIN DATASET (removes the confound that produced false positives in step 31)")
    for ds,g in m.groupby('dataset'):
        if len(g)<15: say("    %-14s %d compounds — too few"%(ds,len(g))); continue
        o3=[]
        for d in DS:
            for c in PC:
                s=g[[d,c]].dropna()
                if len(s)<12 or s[d].std()==0: continue
                r,p=sstats.spearmanr(s[d],s[c]); o3.append(dict(shape=d,pref=c[2:],rho=r,p=p))
        if not o3: continue
        oo3=pd.DataFrame(o3); oo3['q']=multipletests(oo3.p,method='fdr_bh')[1]
        top=oo3.reindex(oo3.rho.abs().sort_values(ascending=False).index).head(3)
        say("    %-14s %d compounds | %d of %d significant"%(ds,len(g),int((oo3.q<0.05).sum()),len(oo3)))
        for _,r in top.iterrows():
            say("       %-12s %-18s rho %+.3f q %.3g"%(r['shape'],r.pref[:18],r.rho,r.q))
    oo2.to_csv('step37_within_compound.csv',index=False)
open('step37_report.txt','w').write("\n".join(log)+"\n")
