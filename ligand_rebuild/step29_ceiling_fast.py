#!/usr/bin/env python3
"""
step29_ceiling_fast.py — same derivation as step29, restructured for speed.

The original grouped by site inside the replicate loop, repeating the same grouping 50 times per
depth per dataset. This groups once into arrays and permutes indices instead.
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

say("="*78); say("STEP 29 (fast) — LABEL-REPRODUCIBILITY CEILING"); say("="*78)
occ=pd.read_csv("occ_long.csv",low_memory=False)

for ds in sorted(occ.dataset.unique()):
    g=occ[occ.dataset==ds]
    # collapse to one occupancy per (site, compound): max across replicates
    m=g.groupby(["uniprot","position","compound"],as_index=False).occupancy_pct.max()
    # one array of occupancies per site, grouped ONCE
    per_site=[v.values for _,v in m.groupby(["uniprot","position"]).occupancy_pct]
    depth=np.array([len(v) for v in per_site])
    say(""); say("-"*78)
    say("%s — %d sites | compounds per site: median %d (25th %d, 75th %d)"
        %(ds,len(per_site),np.median(depth),np.percentile(depth,25),np.percentile(depth,75)))
    for half in (2,3,5,10,20):
        need=half*2
        idx=np.where(depth>=need)[0]
        if len(idx)<200: continue
        vals=[per_site[i] for i in idx]
        aucs=[]
        for _ in range(N_REP):
            y=np.empty(len(vals),dtype=int); s=np.empty(len(vals))
            for j,v in enumerate(vals):
                p=rng.permutation(len(v))
                y[j]=int(v[p[:half]].max()>=TIER)
                s[j]=v[p[half:need]].max()
            if y.sum()<10 or y.sum()==len(y): continue
            aucs.append(roc_auc_score(y,s))
        if len(aucs)<5: continue
        a=np.array(aucs)
        say("  %2d compounds/half | %5d sites | %4.1f%% engaged | ceiling AUC %.3f [%.3f, %.3f]"
            %(half,len(vals),100*y.mean(),a.mean(),np.percentile(a,2.5),np.percentile(a,97.5)))
say("")
say("COMPOSITE ACHIEVES 0.617 (95% CI 0.594-0.640), cross-validated by protein.")
say("Each half is noisier than the full screen, so these are LOWER bounds on the true ceiling.")
open("step29_report.txt","w").write("\n".join(log)+"\n")
