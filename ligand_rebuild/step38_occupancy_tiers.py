#!/usr/bin/env python3
"""
step38_occupancy_tiers.py — does a compound's site preference strengthen with occupancy, and
does that depend on the compound?

WHY TIER
Everything so far split sites at 60% occupancy: liganded or not. That lumps a site a compound
blocks at 62% with one it blocks at 95%, and lumps a site blocked at 55% with one blocked at 2%.
If a compound's shape or size influences WHICH sites it engages, the effect should be strongest
where it engages best — and a binary split averages that away.

TIERS, as specified: 0-10, 10-20, ... 80-90, >=90 percent occupancy.

TWO QUESTIONS
  Q1  Does the site preference strengthen monotonically with occupancy?
      Pooled, step 9 showed burial rising from d = 0.19 at 60% to 0.58 at 90%. Per compound this
      becomes a SLOPE: how much does the preference sharpen as occupancy rises.
  Q2  Does that slope depend on the compound's shape or size?
      This is the tiered version of the question that has now failed five times as a binary
      comparison (steps 31, 32, 34, 36, 37). If shape matters at all it should show here, at the
      sites each compound blocks hardest.

WITHIN NATURE CHEMISTRY ONLY
Every pooled positive so far has turned out to be a dataset confound: compound size encodes
which paper a compound came from (CyMPSA median MW 479, meroterpenoid 372, NatChem 331), and
each dataset was screened against a different site population. NatChem is the only dataset with
enough compounds to test alone, so Q2 is asked there and the pooled version is not reported.

REFERENCE GROUP
Each compound's own tested sites in the 0-10% tier — the sites it demonstrably reached and
demonstrably failed to block. That keeps the comparison within compound and within its own
screening coverage.
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
GATE=float(os.environ.get("PLDDT_GATE",70)); MIN_PER_TIER=int(os.environ.get("MIN_PER_TIER",8))
SF=["rsa_tien","nz_sasa","crowding_8A","atoms_10A","hydrophobic_8A","aromatic_8A",
    "polar_8A","anionic_8A","cationic_8A","net_charge_8A","dist_centroid","propka_pKa"]
BINS=[(0,10),(10,20),(20,30),(30,40),(40,50),(50,60),(60,70),(70,80),(80,90),(90,101)]
LAB=["0-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80","80-90",">=90"]
log=[]
def say(m):
    print(m,flush=True); log.append(m)

say("="*78); say("STEP 38 — SITE PREFERENCE ACROSS OCCUPANCY TIERS"); say("="*78)

def shape(smi):
    m=Chem.MolFromSmiles(str(smi))
    if m is None: return None
    mh=Chem.AddHs(m); ids=AllChem.EmbedMultipleConfs(mh,numConfs=3,randomSeed=0xf00d,maxAttempts=200)
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
    return dict(NPR1=np.mean(n1),NPR2=np.mean(n2),rodness=np.mean(n2)-np.mean(n1),
                volume=np.mean(v),MW=Descriptors.MolWt(m),logP=Descriptors.MolLogP(m))

x=pd.ExcelFile(SMILES_XLSX)
C=[]
for r in x.parse('NatChem').dropna(subset=['SMILES']).itertuples():
    s=shape(r.SMILES)
    if s: C.append(dict(compound=str(r.Compounds).strip(),**s))
C=pd.DataFrame(C); say("NatChem compounds with shape: %d"%len(C))

feat=pd.read_csv('features_expanded.csv',low_memory=False).drop_duplicates(['uniprot','position'])
feat=feat[feat.plddt_site>=GATE]
F=[f for f in SF if f in feat.columns]
fidx=feat.set_index(['uniprot','position'])[F]
occ=pd.read_csv('occ_long.csv',low_memory=False)
occ=occ[occ.dataset=='NatChem'].copy()
occ['base']=[re.match(r'^(\d+[a-z])',str(c)).group(1) if re.match(r'^(\d+[a-z])',str(c)) else str(c) for c in occ.compound]

# ---- pooled tier profile first
say("")
say("Q1a — POOLED across compounds: site features by occupancy tier")
o=occ.merge(feat[['uniprot','position']+F],on=['uniprot','position'],how='inner')
say("  %d site-compound observations above the gate"%len(o))
base=o[(o.occupancy_pct>=0)&(o.occupancy_pct<10)]
say("  %-8s %8s %s"%("tier","n"," ".join("%9s"%f[:9] for f in F[:6])))
for (lo,hi),lab in zip(BINS,LAB):
    s=o[(o.occupancy_pct>=lo)&(o.occupancy_pct<hi)]
    if len(s)<50: continue
    cells=[]
    for f in F[:6]:
        a=pd.to_numeric(s[f],errors='coerce').dropna(); b=pd.to_numeric(base[f],errors='coerce').dropna()
        cells.append("%9.3f"%((a.mean()-b.mean())/b.std()) if len(a)>10 and b.std()>0 else "%9s"%"-")
    say("  %-8s %8d %s"%(lab,len(s)," ".join(cells)))
say("  (values are standardised shifts from the 0-10%% tier)")

# ---- per compound slopes
say("")
say("Q1b — PER COMPOUND: does preference strengthen with occupancy?")
mid=np.array([(lo+hi)/2 for lo,hi in BINS])
recs=[]
for base_c,g in occ.groupby('base'):
    gg=g.merge(feat[['uniprot','position']+F],on=['uniprot','position'],how='inner')
    if len(gg)<200: continue
    ref=gg[(gg.occupancy_pct>=0)&(gg.occupancy_pct<10)]
    if len(ref)<50: continue
    rec=dict(compound=base_c,n_obs=len(gg))
    for f in F:
        xs=[];ys=[]
        b=pd.to_numeric(ref[f],errors='coerce').dropna()
        if len(b)<30 or b.std()==0: continue
        for (lo,hi),mp in zip(BINS,mid):
            s=gg[(gg.occupancy_pct>=lo)&(gg.occupancy_pct<hi)]
            a=pd.to_numeric(s[f],errors='coerce').dropna()
            if len(a)>=MIN_PER_TIER:
                xs.append(mp); ys.append((a.mean()-b.mean())/b.std())
        if len(xs)>=5:
            sl,ic,r,p,se=sstats.linregress(xs,ys)
            rec['slope_'+f]=sl*100          # shift per 100% occupancy
    if len([k for k in rec if k.startswith('slope_')])>=6: recs.append(rec)
S=pd.DataFrame(recs)
say("  compounds with slopes across >=5 tiers: %d"%len(S))
SL=[c for c in S.columns if c.startswith('slope_')]
o2=[]
for c in SL:
    v=S[c].dropna()
    if len(v)<15: continue
    t,p=sstats.ttest_1samp(v,0)
    o2.append(dict(feature=c[6:],mean_slope=v.mean(),sd=v.std(),n=len(v),p=p))
oo=pd.DataFrame(o2); oo['q']=multipletests(oo.p,method='fdr_bh')[1]
oo=oo.reindex(oo.mean_slope.abs().sort_values(ascending=False).index)
say("  %-18s %12s %8s %5s %10s"%("feature","slope/100%","sd","n","q"))
for _,r in oo.iterrows():
    say("  %-18s %12.3f %8.3f %5d %10.3g%s"%(r.feature[:18],r.mean_slope,r.sd,r.n,r.q," *" if r.q<0.05 else ""))

say("")
say("Q2 — DOES THE SLOPE DEPEND ON THE COMPOUND? (NatChem only)")
m=S.merge(C,on='compound',how='inner')
say("  %d compounds with slopes and structures"%len(m))
DS=['rodness','NPR1','NPR2','volume','MW','logP']
o3=[]
for d in DS:
    for c in SL:
        s=m[[d,c]].dropna()
        if len(s)<15 or s[d].std()==0: continue
        r,p=sstats.spearmanr(s[d],s[c]); o3.append(dict(shape=d,slope=c[6:],rho=r,p=p,n=len(s)))
if o3:
    oo3=pd.DataFrame(o3); oo3['q']=multipletests(oo3.p,method='fdr_bh')[1]
    oo3=oo3.reindex(oo3.rho.abs().sort_values(ascending=False).index)
    say("  %-12s %-18s %5s %8s %10s"%("compound","slope in","n","rho","q"))
    for _,r in oo3.head(10).iterrows():
        say("  %-12s %-18s %5d %8.3f %10.3g%s"%(r['shape'],r.slope[:18],r.n,r.rho,r.q," *" if r.q<0.05 else ""))
    say("  -> %d of %d significant (%.1f expected by chance)"%(int((oo3.q<0.05).sum()),len(oo3),0.05*len(oo3)))
    oo3.to_csv('step38_tiers.csv',index=False)
open('step38_report.txt','w').write("\n".join(log)+"\n")
