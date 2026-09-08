#!/usr/bin/env python3
"""
step4_collapse.py — reduce many compound measurements to one row per site per experiment.

INPUT   occ_long.csv   (step 2: 1,642,330 site x compound occupancy measurements)
OUTPUT  sites.csv      one row per (dataset, experiment, probe, cell_line, uniprot, position)
        step4_report.txt

TWO NUMBERS PER SITE, not one. Agreed before writing, and the reason matters.

  max_occupancy_pct   the single best engagement any compound achieved.
                      This is what defines ligandability: a site is engageable if ANY
                      compound engages it. One compound out of 300 still proves the site
                      can be reached.

  n_compounds_over_X  how many distinct compounds exceeded threshold X.
                      A lysine engaged by 200 of 328 compounds is a different kind of
                      site from one engaged by exactly one. The first is likely
                      hyperreactive — engaged because it is intrinsically nucleophilic,
                      not because anything binds there. The second is the interesting
                      case: selective engagement implies a binding event.
                      Collapsing to the maximum alone makes these indistinguishable and
                      would merge REACTIVITY into LIGANDABILITY, which must not happen.

  Also kept: n_compounds_tested (how many compounds this site was measured against —
  the denominator that makes a negative meaningful) and the median, for reference.

KNOWN ASYMMETRY, to carry forward and state in any writeup
  NatChem tested 328 compound x concentration combinations; CyMPSA 15; Meroterpenoid 10.
  The maximum of many draws is higher than the maximum of few, for noise reasons alone.
  So max_occupancy_pct is NOT comparable across datasets in absolute terms. This is
  handled downstream by stratifying within experiment (a site is only ever compared with
  other sites measured the same way), never by pooling raw maxima across datasets.

WHAT IS STILL NOT DECIDED HERE
  No case/control threshold. That is step 5, kept separate so the choice stays visible.
"""
import os
import numpy as np, pandas as pd
import warnings
warnings.filterwarnings("ignore")

THRESHOLDS = (50, 60, 75, 90)     # counts reported at several cuts; none is privileged here
KEY = ["dataset", "experiment", "probe", "cell_line", "uniprot", "position"]

log = []
def say(m):
    print(m); log.append(m)

say("=" * 72)
say("STEP 4 — COLLAPSE TO ONE ROW PER SITE PER EXPERIMENT")
say("=" * 72)

occ = pd.read_csv("occ_long.csv", low_memory=False)
say("read occ_long.csv: %d measurements" % len(occ))

# a compound may appear more than once per site within an experiment only if the source
# had duplicate columns; average those first so no compound is counted twice
pre = len(occ)
occ = (occ.groupby(KEY + ["compound"], as_index=False)
          .agg(occupancy_pct=("occupancy_pct", "mean"),
               peptide=("peptide", "first"),
               peptide_marked=("peptide_marked", "first"),
               cell_fraction=("cell_fraction", "first")))
if len(occ) != pre:
    say("  collapsed %d duplicate site x compound rows (averaged)" % (pre - len(occ)))

agg = {"max_occupancy_pct": ("occupancy_pct", "max"),
       "median_occupancy_pct": ("occupancy_pct", "median"),
       "n_compounds_tested": ("compound", "nunique"),
       "peptide": ("peptide", "first"),
       "peptide_marked": ("peptide_marked", "first"),
       "cell_fraction": ("cell_fraction", "first")}
sites = occ.groupby(KEY, as_index=False).agg(**agg)

for t in THRESHOLDS:
    c = (occ[occ.occupancy_pct >= t].groupby(KEY, as_index=False)
            .compound.nunique().rename(columns={"compound": f"n_compounds_over_{t}"}))
    sites = sites.merge(c, on=KEY, how="left")
    sites[f"n_compounds_over_{t}"] = sites[f"n_compounds_over_{t}"].fillna(0).astype(int)

say("")
say("collapsed to %d site x experiment rows" % len(sites))
say("  distinct (uniprot, position) sites: %d" % sites.groupby(["uniprot", "position"]).ngroups)
say("")
say("%-14s %8s %8s %10s %10s %10s" % ("dataset", "rows", "sites", "cmpds/site", "med max%", "90th max%"))
for ds, g in sites.groupby("dataset"):
    say("%-14s %8d %8d %10.1f %10.1f %10.1f"
        % (ds, len(g), g.groupby(["uniprot", "position"]).ngroups,
           g.n_compounds_tested.median(), g.max_occupancy_pct.median(),
           g.max_occupancy_pct.quantile(.9)))

say("")
say("PROMISCUITY — how many compounds engage each site (at >=60% occupancy)")
say("%-14s %10s %10s %10s %10s" % ("dataset", "0 cmpds", "exactly 1", "2-5", ">5"))
for ds, g in sites.groupby("dataset"):
    n = g.n_compounds_over_60
    say("%-14s %10d %10d %10d %10d"
        % (ds, (n == 0).sum(), (n == 1).sum(), ((n >= 2) & (n <= 5)).sum(), (n > 5).sum()))
say("")
say("  Sites engaged by MANY compounds are candidates for hyperreactivity rather than")
say("  ligandability. Sites engaged by exactly one, especially where an enantiomer was")
say("  also tested and failed, are the strongest evidence of binding-driven engagement.")

# enantiomer pairs: the cleanest available separation of binding from reactivity
say("")
say("ENANTIOMER PAIRS (Meroterpenoid) — same chemistry, opposite three-dimensional shape")
me = occ[occ.dataset == "Meroterpenoid"]
pairs = {}
for c in me.compound.unique():
    s = str(c)
    for a, b in (("(+)-", "(\u2013)-"), ("(+)-", "(-)-")):
        if s.startswith(a):
            partner = b + s[len(a):]
            if partner in set(me.compound):
                pairs[s] = partner
for plus, minus in pairs.items():
    w = me[me.compound.isin([plus, minus])].pivot_table(
        index=["uniprot", "position"], columns="compound", values="occupancy_pct")
    if plus not in w or minus not in w:
        continue
    w = w.dropna()
    both = ((w[plus] >= 60) & (w[minus] >= 60)).sum()
    only_p = ((w[plus] >= 60) & (w[minus] < 60)).sum()
    only_m = ((w[plus] < 60) & (w[minus] >= 60)).sum()
    say("  %-22s vs %-22s  sites %5d | both %3d | only(+) %3d | only(-) %3d"
        % (plus[:22], minus[:22], len(w), both, only_p, only_m))
say("  Sites engaged by one enantiomer but not the other cannot be explained by")
say("  intrinsic lysine reactivity; the two compounds are chemically identical.")

sites.to_csv("sites.csv", index=False)
open("step4_report.txt", "w").write("\n".join(log) + "\n")
say("")
say("wrote sites.csv (%d rows) and step4_report.txt" % len(sites))
say("STILL NO THRESHOLD APPLIED — cases and controls are defined in step 5.")
