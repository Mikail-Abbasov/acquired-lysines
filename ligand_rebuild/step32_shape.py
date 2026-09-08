#!/usr/bin/env python3
"""
step32_shape.py — does compound SHAPE match the shape of the sites it engages?

WHY SHAPE RATHER THAN SIZE
Step 31 tested size, lipophilicity and polarity. The pooled size correlations turned out to be a
dataset confound: compounds with >=20 engaged sites split 48 NatChem / 14 CyMPSA / 7
meroterpenoid with mean MW 331 / 479 / 372, so "large compound" largely meant "CyMPSA compound",
and the effect vanished within NatChem (0 of 150 pairs).

Shape is the better-motivated property. Complementarity between a molecule's three-dimensional
form and the geometry of the site it occupies is the actual mechanism of binding, and the
site-side features (burial, crowding, pocket membership and volume) already describe geometry.

DESIGN — DIRECTIONAL PREDICTIONS, NOT A BLIND GRID
Step 31 tested 150 descriptor x feature pairs on 48 compounds and could not survive correction.
Here the predictions are specified in advance, giving three tests rather than 150:

  P1  ROD-LIKE compounds (elongated, high NPR2 - NPR1) should engage sites that are more OPEN
      and less enclosed: higher nz_sasa, lower in_pocket. A groove accommodates a rod.
  P2  SPHERICAL compounds (both NPR ratios near 1) should engage sites that are more ENCLOSED:
      higher in_pocket, higher pocket_volume. A cavity accommodates a globe.
  P3  FLAT/DISC compounds (low NPR1, NPR2 near 1) should engage sites with more AROMATIC
      neighbours, where face-to-face contact is available.

All tests are run WITHIN NatChem only (48 compounds with >=20 engaged sites, one cell line, one
probe, one sample state), because the pooled analysis is confounded by dataset.

WHAT THIS CANNOT DO
Principal-moment descriptors are IDENTICAL for enantiomers: mirror images have the same moments
of inertia. So this cannot explain the eight sites engaged by one hand of guadial B and not the
other. Chirality requires a different treatment and is not addressed here.

Conformer generation is stochastic and these molecules are flexible, so descriptors are averaged
over multiple embedded conformers and the spread is reported.
"""
import os, re
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
from scipy import stats as sstats
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors3D
RDLogger.DisableLog('rdApp.*')

SMILES_XLSX=os.environ.get("SMILES_XLSX","SMILES.xlsx")
TIER=int(os.environ.get("TIER",60)); MIN_SITES=int(os.environ.get("MIN_SITES",20))
N_CONF=int(os.environ.get("N_CONF",5))
log=[]
def say(m):
    print(m,flush=True); log.append(m)

say("="*78); say("STEP 32 — COMPOUND SHAPE vs ENGAGED-SITE GEOMETRY"); say("="*78)

def shape_desc(smi):
    m=Chem.MolFromSmiles(str(smi))
    if m is None: return None
    m=Chem.AddHs(m)
    ok=AllChem.EmbedMultipleConfs(m,numConfs=N_CONF,randomSeed=0xf00d,maxAttempts=200)
    if len(ok)==0: return None
    try: AllChem.MMFFOptimizeMoleculeConfs(m,maxIters=200)
    except Exception: pass
    out=[]
    for cid in ok:
        try:
            out.append(dict(NPR1=Descriptors3D.NPR1(m,confId=cid),
                            NPR2=Descriptors3D.NPR2(m,confId=cid),
                            asphericity=Descriptors3D.Asphericity(m,confId=cid),
                            eccentricity=Descriptors3D.Eccentricity(m,confId=cid),
                            spherocity=Descriptors3D.SpherocityIndex(m,confId=cid),
                            radius_gyr=Descriptors3D.RadiusOfGyration(m,confId=cid),
                            PBF=Descriptors3D.PBF(m,confId=cid)))
        except Exception: pass
    if not out: return None
    d=pd.DataFrame(out)
    r=d.mean().to_dict(); r['NPR1_sd']=d.NPR1.std(); r['NPR2_sd']=d.NPR2.std()
    r['rodness']=r['NPR2']-r['NPR1']          # high = elongated
    r['sphericity']=r['NPR1']+r['NPR2']-1     # high = globular
    return r

x=pd.ExcelFile(SMILES_XLSX)
d=x.parse('NatChem').dropna(subset=['SMILES'])
say("generating %d conformers per compound for %d NatChem compounds..."%(N_CONF,len(d)))
rows=[]
for i,r in enumerate(d.itertuples(),1):
    if i%40==0: say("  %d/%d"%(i,len(d)))
    s=shape_desc(r.SMILES)
    if s: rows.append(dict(compound=str(r.Compounds).strip(),**s))
S=pd.DataFrame(rows)
say("  shape computed for %d compounds"%len(S))
say("  conformer spread: NPR1 sd median %.3f, NPR2 sd median %.3f"%(S.NPR1_sd.median(),S.NPR2_sd.median()))

occ=pd.read_csv('occ_long.csv',low_memory=False)
occ=occ[occ.dataset=='NatChem'].copy()
occ['base']=[re.match(r'^(\d+[a-z])',str(c)).group(1) if re.match(r'^(\d+[a-z])',str(c)) else str(c) for c in occ.compound]
feat=pd.read_csv('features_expanded.csv',low_memory=False).drop_duplicates(['uniprot','position'])
if 'in_pocket' in feat:
    feat['in_pocket']=feat.in_pocket.astype(str).str.lower().isin(['true','1','1.0']).astype(float)
SF=['nz_sasa','rsa_tien','in_pocket','pocket_volume','aromatic_8A','crowding_8A','atoms_10A']
SF=[f for f in SF if f in feat.columns]
fidx=feat.set_index(['uniprot','position'])[SF]

recs=[]
for base,g in occ.groupby('base'):
    e=g[g.occupancy_pct>=TIER][['uniprot','position']].drop_duplicates()
    t=g[['uniprot','position']].drop_duplicates()
    if len(e)<MIN_SITES: continue
    ei=fidx.reindex(pd.MultiIndex.from_frame(e)).dropna(how='all')
    ti=fidx.reindex(pd.MultiIndex.from_frame(t)).dropna(how='all')
    if len(ei)<MIN_SITES or len(ti)<50: continue
    rec=dict(compound=base,n_engaged=len(ei))
    for f in SF:
        a=pd.to_numeric(ei[f],errors='coerce').dropna(); b=pd.to_numeric(ti[f],errors='coerce').dropna()
        if len(a)>=5 and len(b)>=20 and b.std()>0: rec['d_'+f]=float((a.mean()-b.mean())/b.std())
    recs.append(rec)
E=pd.DataFrame(recs)
m=E.merge(S,on='compound',how='inner')
say("  %d compounds with >=%d engaged sites and shape descriptors"%(len(m),MIN_SITES))
if len(m)<20:
    say("  too few — lower MIN_SITES"); raise SystemExit(0)

say("")
say("PRE-SPECIFIED PREDICTIONS (one-sided, 3 tests not 150)")
P=[("P1 rod-like -> more open sites","rodness","d_nz_sasa","greater"),
   ("P1 rod-like -> less enclosed","rodness","d_in_pocket","less"),
   ("P2 spherical -> more enclosed","sphericity","d_in_pocket","greater"),
   ("P2 spherical -> larger pockets","sphericity","d_pocket_volume","greater"),
   ("P3 flat -> aromatic neighbours","PBF","d_aromatic_8A","less")]
for lab,dc,sc,side in P:
    if dc not in m.columns or sc not in m.columns: continue
    sub=m[[dc,sc]].dropna()
    if len(sub)<15: continue
    r,p2=sstats.spearmanr(sub[dc],sub[sc])
    p1=p2/2 if ((r>0 and side=="greater") or (r<0 and side=="less")) else 1-p2/2
    say("  %-34s rho %+.3f  one-sided p %.4f  n=%d%s"%(lab,r,p1,len(sub)," *" if p1<0.05 else ""))

say("")
say("EXPLORATORY (all shape descriptors x site features, corrected)")
from statsmodels.stats.multitest import multipletests
SD=['NPR1','NPR2','rodness','sphericity','asphericity','eccentricity','spherocity','radius_gyr','PBF']
o=[]
for dc in [c for c in SD if c in m.columns]:
    for sc in [c for c in m.columns if c.startswith('d_')]:
        sub=m[[dc,sc]].dropna()
        if len(sub)<15 or sub[dc].std()==0: continue
        r,p=sstats.spearmanr(sub[dc],sub[sc]); o.append(dict(shape=dc,site=sc[2:],rho=r,p=p,n=len(sub)))
oo=pd.DataFrame(o); oo['q']=multipletests(oo.p,method='fdr_bh')[1]
oo=oo.reindex(oo.rho.abs().sort_values(ascending=False).index)
say("  %-14s %-16s %5s %8s %10s"%("shape","site feature","n","rho","q"))
for _,r in oo.head(10).iterrows():
    say("  %-14s %-16s %5d %8.3f %10.3g%s"%(r.shape,r.site[:16],r.n,r.rho,r.q," *" if r.q<0.05 else ""))
say("  significant: %d of %d"%(int((oo.q<0.05).sum()),len(oo)))
oo.to_csv('step32_shape.csv',index=False)
say("")
say("NOTE: principal-moment descriptors are identical for enantiomers, so none of this can")
say("explain the eight sites engaged by one hand of guadial B and not the other.")
open('step32_report.txt','w').write("\n".join(log)+"\n")
