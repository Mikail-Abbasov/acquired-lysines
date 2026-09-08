#!/usr/bin/env python3
"""
step29_ceiling.py — what is the highest AUC any predictor of ligandability could achieve?

WHY THIS IS NEEDED
The composite classifier reaches AUC 0.617 cross-validated by protein (step 10). Whether that
number is disappointing or at the limit depends on how reproducible the engagement labels are.
Part 1 estimated a ceiling of 0.55-0.65, but that analysis was found unreproducible and its
value cannot be cited. This re-derives the ceiling from the measurements directly.

THE LOGIC
"Engaged" is not a fact about a site; it is the outcome of a measurement. If the same site is
assayed twice and the two assays disagree about whether it is engaged, then no predictor —
however good — can agree with both. The reproducibility of the label therefore bounds the
achievable AUC.

So: split the compound measurements at each site into two independent halves. Use one half to
define the label (engaged at >=60% occupancy) and the other half to produce a score (maximum
observed occupancy). The AUC of the second half predicting the first is what a *perfect*
structural model would achieve if structure captured everything the assay captures. That is the
ceiling.

This is an empirical, assumption-light derivation: no model of label noise is imposed, and the
answer comes from the same data the composite was trained on.

TWO THINGS THAT MUST BE REPORTED WITH IT

  SPLIT SIZE MATTERS. Each half has fewer compounds than the full screen, so each half is
  noisier than the real label. The ceiling estimated this way is therefore a LOWER bound on the
  true ceiling. It is computed at several split sizes so the trend is visible: if the ceiling
  rises with more compounds per half, the asymptote is the quantity of interest.

  SCREENING DEPTH IS UNEQUAL. Nature Chemistry sites see a median of 41 compounds; the JACS
  datasets exactly 5. A site with 5 measurements cannot be split into two informative halves, so
  the ceiling is computed per dataset and only where depth allows.

INPUT   occ_long.csv
OUTPUT  step29_ceiling.csv, step29_report.txt
"""
import os
import numpy as np, pandas as pd
import warnings
warnings.filterwarnings("ignore")
from sklearn.metrics import roc_auc_score

TIER = int(os.environ.get("TIER", 60))
N_REP = int(os.environ.get("N_REP", 50))
rng = np.random.default_rng(0)

log = []
def say(m):
    print(m, flush=True); log.append(m)

say("=" * 78)
say("STEP 29 — LABEL-REPRODUCIBILITY CEILING FOR THE COMPOSITE AUC")
say("=" * 78)

occ = pd.read_csv("occ_long.csv", low_memory=False)
say("%d measurements across %d datasets" % (len(occ), occ.dataset.nunique()))

rows = []
for ds in sorted(occ.dataset.unique()):
    g = occ[occ.dataset == ds]
    depth = g.groupby(["uniprot", "position"]).compound.nunique()
    say("")
    say("-" * 78)
    say("%s — %d sites | compounds per site: median %d, 25th %d, 75th %d"
        % (ds, len(depth), depth.median(), depth.quantile(.25), depth.quantile(.75)))

    for half in (2, 3, 5, 10, 20):
        need = half * 2
        ok = set(depth[depth >= need].index)
        if len(ok) < 200:
            continue
        sub = g[list(zip(g.uniprot, g.position)) if False else
                g.set_index(["uniprot", "position"]).index.isin(ok)]
        aucs = []
        for rep in range(N_REP):
            lab, scr = {}, {}
            for key, gg in sub.groupby(["uniprot", "position"]):
                comps = gg.compound.unique()
                if len(comps) < need:
                    continue
                pick = rng.permutation(comps)
                a, b = set(pick[:half]), set(pick[half:half * 2])
                oa = gg[gg.compound.isin(a)].occupancy_pct
                ob = gg[gg.compound.isin(b)].occupancy_pct
                if not len(oa) or not len(ob):
                    continue
                lab[key] = int(oa.max() >= TIER)
                scr[key] = float(ob.max())
            y = np.array([lab[k] for k in lab])
            s = np.array([scr[k] for k in lab])
            if len(y) < 100 or y.sum() < 10 or y.sum() == len(y):
                continue
            aucs.append(roc_auc_score(y, s))
        if len(aucs) < 5:
            continue
        a = np.array(aucs)
        rows.append(dict(dataset=ds, compounds_per_half=half, n_sites=len(lab),
                         positive_rate=float(y.mean()), auc=a.mean(),
                         lo=np.percentile(a, 2.5), hi=np.percentile(a, 97.5)))
        say("  %2d compounds/half | %5d sites | %4.1f%% engaged | ceiling AUC %.3f [%.3f, %.3f]"
            % (half, len(lab), 100 * y.mean(), a.mean(),
               np.percentile(a, 2.5), np.percentile(a, 97.5)))

res = pd.DataFrame(rows)
if len(res):
    res.to_csv("step29_ceiling.csv", index=False)
    say("")
    say("=" * 78)
    say("SUMMARY")
    for ds, g in res.groupby("dataset"):
        g = g.sort_values("compounds_per_half")
        say("  %-14s %s" % (ds, " -> ".join("%d:%.3f" % (r.compounds_per_half, r.auc)
                                            for r in g.itertuples())))
    best = res.groupby("dataset").auc.max()
    say("")
    say("  highest ceiling observed per dataset:")
    for ds, v in best.items():
        say("    %-14s %.3f" % (ds, v))
    say("")
    say("  COMPOSITE ACHIEVES 0.617 (95%% CI 0.594-0.640), cross-validated by protein.")
    say("")
    say("READING THIS")
    say("  Each half is noisier than the full screen, so these are LOWER bounds on the true")
    say("  ceiling. If the value rises with compounds per half, the asymptote is the quantity")
    say("  of interest and should be quoted rather than the smallest split.")
    say("  If the composite's 0.617 sits within the range of these ceilings, the claim that it")
    say("  is at the limit of what the labels permit is supported. If the ceiling is far higher,")
    say("  the claim must be dropped and the AUC reported without qualification.")
open("step29_report.txt", "w").write("\n".join(log) + "\n")
say("")
say("wrote step29_ceiling.csv and step29_report.txt")
