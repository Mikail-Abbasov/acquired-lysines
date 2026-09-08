#!/usr/bin/env python3
"""
step40_threshold_free.py — does the breadth result depend on the 60% cut-off?

THE PROBLEM
Step 38 showed the site-feature profile shifts smoothly across occupancy: no discontinuity
anywhere, and 60% sits mid-slope. So a lysine blocked at 59% and one blocked at 61% are placed
in different categories on the basis of nothing.

Step 39 compounds this. Counting "how many compounds ligand this site" applies that arbitrary
line once per compound: a site blocked at 58% by ten compounds counts as ZERO, one blocked at
62% by ten counts as TEN. Chemically near-identical, maximally different in the analysis.

TWO INDEPENDENT FIXES, because they answer different objections

  A. THRESHOLD-FREE. Replace the count with continuous measures that never apply a cut:
       mean_occ      mean occupancy across all compounds tested at that site
       max_occ       highest occupancy observed
       sum_occ_norm  total occupancy divided by compounds tested
     Correlate site features against these directly. No threshold anywhere.

  B. THRESHOLD-VARYING. Keep the intuitive "how many compounds" framing but recompute the
     gradient at 40, 50, 60, 70 and 80%. If the shape is the same at every cut, the published
     convention is harmless.

THE CONFOUND IN (A), which must be handled
Total occupancy accumulates with screening depth: a site tested against 200 compounds
accumulates more than one tested against 20, whatever its properties. So sum_occ is normalised
by compounds tested, and every analysis is restricted to sites tested against >=20 compounds.
Normalisation is itself a choice, but a less arbitrary one than a cut-off.
"""
import os
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
from scipy import stats as sstats
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests

GATE=float(os.environ.get("PLDDT_GATE",70)); MIN_TESTED=int(os.environ.get("MIN_TESTED",20))
N_PERM=1000
SF=["rsa_tien","nz_sasa","sasa_sidechain","crowding_8A","atoms_10A","hydrophobic_8A",
    "aromatic_8A","polar_8A","anionic_8A","cationic_8A","net_charge_8A","n_lys_8A",
    "dist_centroid","propka_pKa","in_pocket"]
rng=np.random.default_rng(0)
log=[]
def say(m):
    print(m,flush=True); log.append(m)

say("="*78); say("STEP 40 — IS THE BREADTH RESULT AN ARTEFACT OF THE 60% CUT-OFF?"); say("="*78)

occ=pd.read_csv('occ_long.csv',low_memory=False)
feat=pd.read_csv('features_expanded.csv',low_memory=False).drop_duplicates(['uniprot','position'])
if 'in_pocket' in feat:
    feat['in_pocket']=feat.in_pocket.astype(str).str.lower().isin(['true','1','1.0']).astype(float)
feat=feat[feat.plddt_site>=GATE]
F=[f for f in SF if f in feat.columns]

agg=(occ.groupby(['uniprot','position'])
     .agg(n_tested=('compound','nunique'),mean_occ=('occupancy_pct','mean'),
          max_occ=('occupancy_pct','max'),sum_occ=('occupancy_pct','sum')).reset_index())
agg['sum_occ_norm']=agg.sum_occ/agg.n_tested
d=feat[['uniprot','position']+F].merge(agg,on=['uniprot','position'],how='inner')
d=d[d.n_tested>=MIN_TESTED].dropna()
say("%d lysines above pLDDT %.0f, tested against >=%d compounds"%(len(d),GATE,MIN_TESTED))
say("  mean_occ: median %.1f | max_occ: median %.1f | sum_occ_norm: median %.1f"
    %(d.mean_occ.median(),d.max_occ.median(),d.sum_occ_norm.median()))

say("")
say("="*78)
say("A. THRESHOLD-FREE — site features vs continuous occupancy measures")
o=[]
for meas in ('mean_occ','max_occ','sum_occ_norm'):
    for f in F:
        r,p=sstats.spearmanr(d[meas],d[f])
        o.append(dict(measure=meas,feature=f,rho=r,p=p))
oo=pd.DataFrame(o); oo['q']=multipletests(oo.p,method='fdr_bh')[1]
for meas in ('mean_occ','max_occ','sum_occ_norm'):
    s=oo[oo.measure==meas].reindex(oo[oo.measure==meas].rho.abs().sort_values(ascending=False).index)
    say("")
    say("  %s (n=%d)"%(meas,len(d)))
    say("  %-18s %8s %10s"%("feature","rho","q"))
    for _,r in s.head(8).iterrows():
        say("  %-18s %8.3f %10.3g%s"%(r.feature[:18],r.rho,r.q," *" if r.q<0.05 else ""))
    say("  -> %d of %d features significant"%(int((s.q<0.05).sum()),len(s)))
oo.to_csv('step40_continuous.csv',index=False)

say("")
say("  DECILE PROFILE on max_occ — no threshold, just the observed distribution")
d2=d.copy(); d2['dec']=pd.qcut(d2.max_occ,10,labels=False,duplicates='drop')
say("  %-6s %8s %9s %s"%("decile","n","med max",""))
say("  %-6s %8s %9s %s"%("","","", " ".join("%9s"%f[:9] for f in ['rsa_tien','crowding_8A','atoms_10A','propka_pKa','in_pocket'])))
base=d2[d2.dec==0]
for k in sorted(d2.dec.dropna().unique()):
    s=d2[d2.dec==k]
    cells=[]
    for f in ['rsa_tien','crowding_8A','atoms_10A','propka_pKa','in_pocket']:
        a=pd.to_numeric(s[f],errors='coerce').dropna(); b=pd.to_numeric(base[f],errors='coerce').dropna()
        cells.append("%9.3f"%((a.mean()-b.mean())/b.std()) if b.std()>0 else "%9s"%"-")
    say("  %-6d %8d %9.1f %s"%(k+1,len(s),s.max_occ.median()," ".join(cells)))
say("  (standardised shifts from the lowest decile)")

say("")
say("="*78)
say("B. THRESHOLD-VARYING — the compound-count gradient at five different cuts")
for cut in (40,50,60,70,80):
    nlig=(occ[occ.occupancy_pct>=cut].groupby(['uniprot','position']).compound.nunique()
          .rename('n_lig').reset_index())
    dd=d.merge(nlig,on=['uniprot','position'],how='left')
    dd['n_lig']=dd.n_lig.fillna(0).astype(int)
    X=StandardScaler().fit_transform(dd[F].astype(float))
    base_m=(dd.n_lig==0).values
    if base_m.sum()<100: say("  cut %d%%: too few zero sites"%cut); continue
    say("")
    say("  CUT AT %d%% — %d sites with zero liganding compounds"%(cut,int(base_m.sum())))
    say("  %-8s %7s %10s %10s %8s"%("n compounds","n","distance","null","p"))
    prev=None; mono=True
    for lo,hi,lab in ((1,1,'1'),(2,3,'2-3'),(4,5,'4-5'),(6,9,'6-9'),(10,10000,'>=10')):
        m=((dd.n_lig>=lo)&(dd.n_lig<=hi)).values
        if m.sum()<15: continue
        obs=float(np.linalg.norm(X[m].mean(0)-X[base_m].mean(0)))
        pool=np.where(m|base_m)[0]; na=int(m.sum())
        null=np.array([np.linalg.norm(X[i[:na]].mean(0)-X[i[na:]].mean(0))
                       for i in (rng.permutation(pool) for _ in range(N_PERM))])
        pv=float((null>=obs).mean())
        say("  %-8s %7d %10.3f %10.3f %8.4f%s"%(lab,na,obs,null.mean(),pv," *" if pv<0.05 else ""))
        if prev is not None and obs<prev: mono=False
        prev=obs
    say("  -> %s"%("monotonic" if mono else "NOT monotonic"))

say("")
say("READING THIS")
say("  If A shows the same features trending with continuous occupancy, and B gives the same")
say("  gradient shape at every cut, then the 60% convention is harmless and the breadth result")
say("  does not depend on it.")
open('step40_report.txt','w').write("\n".join(log)+"\n")
