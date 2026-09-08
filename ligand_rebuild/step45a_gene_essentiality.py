#!/usr/bin/env python3
"""
step45a_gene_essentiality.py — classify every gene by dependency, once, and cache it.

CLASSES
  common essential      DepMap's own inferred list (CRISPRInferredCommonEssentials.csv)
  strongly selective    not common essential, but a strong dependency in a minority of lines
  weakly selective      dependent in a minority, but never strongly
  not essential         essentially no line depends on it

WHY NOT normLRT. DepMap's published "strongly selective" flag is a likelihood-ratio statistic
comparing a skewed-t to a normal fit, and that file is not in Datasets/. Rather than approximate
it with a differently-shaped distribution and call it the same thing, this uses a transparent
count-based rule on the Chronos gene effect scores, stated in full so a reader can reproduce or
disagree with it:

  dependent in a line   gene effect < -0.5   (DepMap's conventional dependency threshold)
  strongly dependent    gene effect < -1.0   (the median common-essential effect is about -1)

  strongly selective    strongly dependent in >= 3 lines AND dependent in < 20% of lines
  weakly selective      dependent in >= 3 lines AND dependent in < 20% of lines
  not essential         everything else that is not common essential

The 20% ceiling is what makes it selective rather than broad; the >= 3 line floor keeps single
noisy lines from creating a category. Both are reported alongside the counts so the sensitivity
is visible.
"""
import os, re, sys, time
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")

DS = os.environ.get("DEPMAP_DIR", "../Datasets")
DEP, STRONG = -0.5, -1.0
MIN_LINES, MAX_FRAC = 3, 0.20

log = []
def say(m):
    print(m, flush=True); log.append(m)

say("=" * 78)
say("STEP 45a — GENE DEPENDENCY CLASSES")
say("=" * 78)

t0 = time.time()
CE = pd.read_csv(os.path.join(DS, "CRISPRInferredCommonEssentials.csv"))
ce = set(re.sub(r"\s*\(\d+\)$", "", str(x)).strip() for x in CE.iloc[:, 0].dropna())
say("common essential genes (DepMap inferred): %d" % len(ce))

G = pd.read_csv(os.path.join(DS, "CRISPRGeneEffect.csv"), index_col=0)
G.columns = [re.sub(r"\s*\(\d+\)$", "", c).strip() for c in G.columns]
say("gene effect matrix: %d cell lines x %d genes (%.0f s)" % (G.shape[0], G.shape[1], time.time()-t0))

n_lines = G.notna().sum(0)
n_dep = (G < DEP).sum(0)
n_strong = (G < STRONG).sum(0)
frac_dep = n_dep / n_lines.replace(0, np.nan)

R = pd.DataFrame({"gene": G.columns, "n_lines": n_lines.values, "n_dep": n_dep.values,
                  "n_strong": n_strong.values, "frac_dep": frac_dep.values,
                  "median_effect": G.median(0).values, "min_effect": G.min(0).values})
def klass(r):
    if r.gene in ce: return "common essential"
    if r.n_strong >= MIN_LINES and r.frac_dep < MAX_FRAC: return "strongly selective"
    if r.n_dep >= MIN_LINES and r.frac_dep < MAX_FRAC: return "weakly selective"
    return "not essential"
R["klass"] = R.apply(klass, axis=1)

say("")
say("  %-22s %8s %10s %12s" % ("class", "genes", "% of all", "median effect"))
for k in ["common essential", "strongly selective", "weakly selective", "not essential"]:
    x = R[R.klass == k]
    say("  %-22s %8d %9.1f%% %12.2f" % (k, len(x), 100*len(x)/len(R), x.median_effect.median()))

say("")
say("  sensitivity of the selective classes to the two cut-offs:")
for mf in (0.10, 0.15, 0.20, 0.30):
    n = int(((R.n_strong >= MIN_LINES) & (R.frac_dep < mf) & (~R.gene.isin(ce))).sum())
    say("    dependent in < %2.0f%% of lines -> %d strongly selective" % (100*mf, n))
for ml in (2, 3, 5, 10):
    n = int(((R.n_strong >= ml) & (R.frac_dep < MAX_FRAC) & (~R.gene.isin(ce))).sum())
    say("    strongly dependent in >= %2d lines -> %d strongly selective" % (ml, n))

R.to_csv("step45a_gene_classes.csv", index=False)
open("step45a_report.txt", "w").write("\n".join(log) + "\n")
say("")
say("wrote step45a_gene_classes.csv (%d genes), step45a_report.txt" % len(R))
