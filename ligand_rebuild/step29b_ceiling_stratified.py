#!/usr/bin/env python3
"""
step29b_ceiling_stratified.py — is measurement reproducibility different in ordered and
disordered regions?

WHY
Step 29 gave a label-reproducibility ceiling of 0.72 (Nature Chemistry, 20 compounds per half),
well above the composite's 0.617. That was computed across all sites regardless of structural
confidence.

The question is whether the ceiling itself depends on structural order. If lysines in disordered
regions are intrinsically noisier to measure — more conformational states, more variable
accessibility — then 0.72 averages two different populations, and the fair comparison for a
structural model is the ceiling in the region where structural features are meaningful at all
(pLDDT >= 70, per the gating established earlier).

NOTE ON A COINCIDENCE
The ceiling value 0.72 and the pLDDT threshold 70 are not related. The ceiling is an AUC, a
unitless probability between 0.5 and 1.0; pLDDT is a 0-100 confidence score. Step 29 used only
occupancy measurements and no structural information whatsoever, so the two cannot influence
each other. This script asks a real question that the coincidence happens to suggest.

METHOD
Identical split-half procedure, stratified by pLDDT band at the site.
"""
import os
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
from sklearn.metrics import roc_auc_score

TIER=int(os.environ.get("TIER",60)); N_REP=int(os.environ.get("N_REP",50))
rng=np.random.default_rng(0)
log=[]
def say(m):
    print(m,flush=True); log.append(m)

say("="*78); say("STEP 29b — CEILING STRATIFIED BY STRUCTURAL CONFIDENCE"); say("="*78)
occ=pd.read_csv("occ_long.csv",low_memory=False)
feat=pd.read_csv("features_expanded.csv",low_memory=False).drop_duplicates(["uniprot","position"])
pl=dict(zip(zip(feat.uniprot,feat.position.astype(int)),feat.plddt_site))

BANDS=[(0,50,"< 50  disordered"),(50,70,"50-70 low"),(70,90,"70-90 confident"),(90,101,">= 90 very high")]

for ds in ("NatChem","CyMPSA","Meroterpenoid"):
    g=occ[occ.dataset==ds]
    m=g.groupby(["uniprot","position","compound"],as_index=False).occupancy_pct.max()
    grouped=list(m.groupby(["uniprot","position"]).occupancy_pct)
    keys=[k for k,_ in grouped]; vals=[v.values for _,v in grouped]
    band=[]
    for u,p in keys:
        x=pl.get((u,int(p)),np.nan)
        band.append(next((lab for lo,hi,lab in BANDS if x==x and lo<=x<hi),"no structure"))
    band=np.array(band)
    say(""); say("-"*78); say("%s"%ds)
    for lo,hi,lab in BANDS+[(None,None,"no structure")]:
        sel=np.where(band==lab)[0]
        if len(sel)<200: continue
        depth=np.array([len(vals[i]) for i in sel])
        best=None
        for half in (2,3,5,10,20):
            need=half*2
            use=[vals[i] for i in sel if len(vals[i])>=need]
            if len(use)<200: continue
            aucs=[]
            for _ in range(N_REP):
                y=np.empty(len(use),dtype=int); s=np.empty(len(use))
                for j,v in enumerate(use):
                    perm=rng.permutation(len(v))
                    y[j]=int(v[perm[:half]].max()>=TIER); s[j]=v[perm[half:need]].max()
                if y.sum()<10 or y.sum()==len(y): continue
                aucs.append(roc_auc_score(y,s))
            if len(aucs)>=5:
                best=(half,len(use),np.mean(aucs),np.percentile(aucs,2.5),np.percentile(aucs,97.5),100*y.mean())
        if best:
            say("  %-18s n=%5d | deepest split %2d/half (%d sites) | ceiling %.3f [%.3f, %.3f] | %.1f%% engaged"
                %(lab,len(sel),best[0],best[1],best[2],best[3],best[4],best[5]))
say("")
say("READING THIS")
say("  If the ceiling is HIGHER at pLDDT >= 70, then 0.72 averaged two populations and the fair")
say("  comparison for a structural model is the confident-region ceiling — which would make the")
say("  gap to 0.617 LARGER, not smaller.")
say("  If the ceiling is FLAT across bands, measurement reproducibility is independent of")
say("  structural order and 0.72 stands as the comparison.")
open("step29b_report.txt","w").write("\n".join(log)+"\n")
