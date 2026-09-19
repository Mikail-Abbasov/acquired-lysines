#!/usr/bin/env python3
"""
figS9_v1.py — Support for Figure 8: the anchoring question, settled.

Figure 8 measures the environment of positions where a lysine does not yet exist, so it has to
choose a reference atom. That choice was made once and never tested. This figure tests it: the
same comparison is run under three legitimate anchorings and under the mismatched one that
produced step 21's original result, with a positive control carried through every arm.

NO SENTENCES INSIDE PANELS. Everything explanatory belongs in the legend.

  a  centroid distances, four ways, with the engaged-vs-unengaged control   step28
  b  the same, expressed against that control                              step28
  c  overlap with the engaged region, four ways                            step28
  d  what one-site-per-protein does to the excess                          step28b
  e  the side-chain build, checked against real lysines                    step27 / step25b
  f  overlap is not an artefact of the 90% cut                             step28b

OUTPUT  figS9.pdf, figS9.png
"""
import os, sys, re
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import abbasov_style as st
import figure_data as fd
st.apply()
ORANGE, BLUE, GREY, LIGHT, INK = st.ORANGE, st.BLUE, st.GREEN, st.LIGHTGREY, st.GREY

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
LR = os.path.join(ROOT, "ligand_rebuild")
R = pd.read_csv(os.path.join(LR, "step28_results.csv"))
RAT = pd.read_csv(os.path.join(LR, "step28b_ratios.csv"))
SUM = pd.read_csv(os.path.join(LR, "step25b_summary.csv"))
MOD = pd.read_csv(os.path.join(LR, "step27_modnz_measured.csv"), usecols=["nz_error"])
TXT = open(os.path.join(LR, "step28b_report.txt")).read()

SHORT = {"A  both sides at CB (the plan)": "CB\nboth sides",
         "B  both sides modelled, at NZ": "modelled NZ\nboth sides",
         "C  MISMATCHED: union CB vs measured NZ": "MISMATCHED\nCB vs NZ",
         "D  both sides at CB, symmetric convention": "CB\nsymmetric"}
SLUG = {"A": "cb", "B": "nz", "C": "mismatched", "D": "cb_symmetric"}
R["slug"] = R.arm.str[0].map(SLUG)
R["short"] = R.arm.map(SHORT)
R["art"] = R.arm.str.startswith("C")
R = R.iloc[[0, 1, 3, 2]].reset_index(drop=True)
x = np.arange(len(R))

fig = plt.figure(figsize=(st.DOUBLE, st.DOUBLE * 0.80))
gs = fig.add_gridspec(2, 3, hspace=0.62, wspace=0.42,
                      left=0.075, right=0.985, top=0.96, bottom=0.085)

# ------------------------------------------------------------------ a
ax = fig.add_subplot(gs[0, 0]); st.panel_label(ax, "a", dx=-0.26, dy=1.09)
ax.set_title("Distances, four ways", fontsize=7.5, pad=6, color=INK)
w = 0.26
ax.bar(x - w, R.d_union_eng, w, color=ORANGE, edgecolor=INK, lw=0.5)
ax.bar(x,     R.d_union_not, w, color=LIGHT,  edgecolor=INK, lw=0.5)
ax.bar(x + w, R.d_control,   w, color=BLUE,   edgecolor=INK, lw=0.5)
fd.panel("figS9", "a",
         [dict(arm=r.arm, slug=r.slug, comparison=nm, distance=float(getattr(r, col)),
               n_union=int(r.n_union), n_engaged=int(r.n_eng), n_not_engaged=int(r.n_not))
          for r in R.itertuples()
          for nm, col in (("acquired vs engaged", "d_union_eng"),
                          ("acquired vs unengaged", "d_union_not"),
                          ("engaged vs unengaged (control)", "d_control"))],
         note="distance between group averages, four anchorings x three comparisons")
ax.set_xticks(x); ax.set_xticklabels(R.short, fontsize=5.9, linespacing=1.3)
ax.set_ylabel("Distance between group averages", fontsize=6.7)
ax.set_ylim(0, 3.05)
for lab, col, yy in (("acquired vs engaged", ORANGE, 2.92),
                     ("acquired vs unengaged", GREY, 2.72),
                     ("engaged vs unengaged (control)", BLUE, 2.52)):
    ax.text(-0.42, yy, lab, color=col, fontsize=5.9, ha="left", va="center")

# ------------------------------------------------------------------ b
ax = fig.add_subplot(gs[0, 1]); st.panel_label(ax, "b", dx=-0.26, dy=1.09)
ax.set_title("Measured against the control", fontsize=7.5, pad=6, color=INK)
ax.bar(x, R.ratio_to_control, 0.55, color=[ORANGE if a else BLUE for a in R.art],
       edgecolor=INK, lw=0.5)
for xi, v in zip(x, R.ratio_to_control):
    ax.text(xi, v + 0.13, "%.2f" % v, ha="center", va="bottom", fontsize=6.6, color=INK)
ax.axhline(1.0, color=INK, lw=0.9)
fd.panel("figS9", "b",
         [dict(arm=r.arm, slug=r.slug, mismatched=bool(r.art),
               d_union_engaged=float(r.d_union_eng), d_control=float(r.d_control),
               ratio_to_control=float(r.ratio_to_control)) for r in R.itertuples()],
         note="each arm's acquired-vs-engaged distance divided by its own control")
ax.set_xticks(x); ax.set_xticklabels(R.short, fontsize=5.9, linespacing=1.3)
ax.set_ylabel("Multiple of the control", fontsize=6.7)
ax.set_ylim(0, 8.0)

# ------------------------------------------------------------------ c
ax = fig.add_subplot(gs[0, 2]); st.panel_label(ax, "c", dx=-0.26, dy=1.09)
ax.set_title("Overlap with the engaged region", fontsize=7.5, pad=6, color=INK)
w = 0.32
ax.bar(x - w / 2, 100 * R.overlap_union, w, color=ORANGE, edgecolor=INK, lw=0.5)
ax.bar(x + w / 2, 100 * R.overlap_not, w, color=GREY, edgecolor=INK, lw=0.5)
ax.axhline(90, color=BLUE, lw=1.0, ls=(0, (3, 2)))
fd.panel("figS9", "c",
         [dict(arm=r.arm, slug=r.slug, population=nm, pct_inside_90_region=100 * float(getattr(r, col)))
          for r in R.itertuples()
          for nm, col in (("acquired union", "overlap_union"),
                          ("labelled, not engaged", "overlap_not"))],
         note="fraction of each population inside the region holding 90% of engaged lysines")
ax.set_xticks(x); ax.set_xticklabels(R.short, fontsize=5.9, linespacing=1.3)
ax.set_ylabel("Inside the 90% engaged region (%)", fontsize=6.7)
ax.set_ylim(0, 118)
ax.text(-0.42, 112, "acquired", color=ORANGE, fontsize=5.9, ha="left", va="center")
ax.text(0.85, 112, "labelled, not engaged", color=GREY, fontsize=5.9, ha="left", va="center")
ax.text(3.45, 90, "engaged", color=BLUE, fontsize=5.9, ha="right", va="bottom")

# ------------------------------------------------------------------ d
ax = fig.add_subplot(gs[1, 0]); st.panel_label(ax, "d", dx=-0.26, dy=1.09)
ax.set_title("One site per protein", fontsize=7.5, pad=6, color=INK)
xx = np.arange(2); w = 0.32
ax.bar(xx - w / 2, RAT.ratio_site, w, color=LIGHT, edgecolor=INK, lw=0.5)
ax.bar(xx + w / 2, RAT.ratio_prot, w, color=BLUE, edgecolor=INK, lw=0.5)
for xi, a_, b_ in zip(xx, RAT.ratio_site, RAT.ratio_prot):
    ax.text(xi - w / 2, a_ + 0.035, "%.2f" % a_, ha="center", va="bottom", fontsize=6.4, color=INK)
    ax.text(xi + w / 2, b_ + 0.035, "%.2f" % b_, ha="center", va="bottom", fontsize=6.4, color=BLUE)
ax.axhline(1.0, color=INK, lw=0.9)
fd.panel("figS9", "d",
         [dict(arm=r.arm, design=nm, ratio_to_control=float(getattr(r, col)))
          for r in RAT.itertuples()
          for nm, col in (("every site", "ratio_site"), ("one site per protein", "ratio_prot"))],
         note="what one-site-per-protein does to the excess; 20 draws, step 28b")
ax.set_xticks(xx); ax.set_xticklabels(["CB\nboth sides", "modelled NZ\nboth sides"],
                                      fontsize=5.9, linespacing=1.3)
ax.set_ylabel("Multiple of the control", fontsize=6.7)
ax.set_ylim(0, 2.1); ax.set_xlim(-0.62, 1.62)
ax.text(-0.55, 1.95, "every site", color=GREY, fontsize=5.9, ha="left", va="center")
ax.text(0.42, 1.95, "one per protein", color=BLUE, fontsize=5.9, ha="left", va="center")

# ------------------------------------------------------------------ e
ax = fig.add_subplot(gs[1, 1]); st.panel_label(ax, "e", dx=-0.26, dy=1.09)
ax.set_title("The built side chain, checked", fontsize=7.5, pad=6, color=INK)
v = MOD.nz_error.dropna().values
_cnt, _edges, _ = ax.hist(v, bins=60, color=LIGHT, edgecolor=INK, linewidth=0.25)
med = float(np.median(v))
ax.axvline(med, color=ORANGE, lw=1.2)
ax.text(med + 0.15, ax.get_ylim()[1] * 0.92, "median %.2f Å" % med, color=ORANGE,
        fontsize=6.2, ha="left", va="top")
CB_NZ = 6.4   # extended all-trans lysine, Cbeta to Nzeta; declared, not measured
fd.panel("figS9", "e",
         [dict(bin_left=float(_edges[i]), bin_right=float(_edges[i + 1]), measured_lysines=int(_cnt[i]))
          for i in range(len(_cnt))],
         note=("modelled all-trans amine vs real amine, %s lysines, %d bins; median %.4f A, "
               "CB-to-NZ reference %.1f A" % (format(len(v), ","), len(_cnt), med, 6.4)))
ax.axvline(CB_NZ, color=BLUE, lw=1.0, ls=(0, (3, 2)))
ax.text(6.25, ax.get_ylim()[1] * 0.62, "CB\u2192NZ\n%.1f \u00c5" % CB_NZ, color=BLUE, fontsize=6.0,
        ha="right", va="top", linespacing=1.3)
ax.set_xlabel("Modelled amine vs real amine (Å)", fontsize=6.7)
ax.set_ylabel("Measured lysines (n = %s)" % format(len(v), ","), fontsize=6.7)
ax.set_xlim(0, 9)

# ------------------------------------------------------------------ f
ax = fig.add_subplot(gs[1, 2]); st.panel_label(ax, "f", dx=-0.26, dy=1.09)
ax.set_title("Overlap at every cut", fontsize=7.5, pad=6, color=INK)
qs = [50, 75, 90, 95]
rows = {}
for line in TXT.splitlines():
    m = re.match(r"\s+([AB])\s+(acquired union|labelled, not engaged)\s+(.*)$", line)
    if m:
        vals = [float(t.strip().rstrip("%")) for t in m.group(3).split()]
        if len(vals) == 4:
            rows[(m.group(1), m.group(2))] = vals
for arm, ls in (("A", "-"), ("B", (0, (4, 2)))):
    for pop, col in (("acquired union", ORANGE), ("labelled, not engaged", GREY)):
        if (arm, pop) in rows:
            ax.plot(qs, rows[(arm, pop)], linestyle=ls, color=col, marker="o",
                    ms=3.2, lw=1.1)
ax.plot(qs, qs, linestyle=":", color=BLUE, lw=1.0)
ax.set_xticks(qs); ax.set_xlabel("Size of the engaged region (%)", fontsize=6.7)
ax.set_ylabel("Sites inside it (%)", fontsize=6.7)
ax.set_ylim(35, 102)
ax.text(50, 99, "acquired", color=ORANGE, fontsize=5.9, ha="left", va="top")
ax.text(50, 93, "labelled, not engaged", color=GREY, fontsize=5.9, ha="left", va="top")
ax.text(95, 96, "engaged", color=BLUE, fontsize=5.9, ha="right", va="bottom")
_F_ROWS = dict(rows)   # panel f, printed so the legend does not have to hold unprinted numbers
for _k in (("A", "acquired union"), ("B", "acquired union"),
           ("A", "labelled, not engaged"), ("B", "labelled, not engaged")):
    assert _k in _F_ROWS, "step28b report no longer carries %r" % (_k,)
fd.panel("figS9", "f",
         [dict(arm=k[0], slug=SLUG[k[0]], population=k[1], engaged_region_pct=q, pct_inside=val)
          for k, vals in sorted(_F_ROWS.items()) for q, val in zip(qs, vals)],
         note="four quantile cuts, both legitimate arms, both comparison populations")

# NO FIGURE-LEVEL TITLE — see the note in fig8_v1.py. 2026-09-05 (M.A.). The gridspec top
# margin was raised from 0.925 to 0.96 to reclaim the band the title occupied.

# ---------------------------------------------------------------- ledger
_D = R.set_index("slug")
_RT = RAT.set_index(RAT.arm.str[0].map(SLUG))
HEAD = {}
for _sl in ("cb", "nz", "cb_symmetric", "mismatched"):
    HEAD["a_d_union_engaged_%s" % _sl] = round(float(_D.at[_sl, "d_union_eng"]), 3)
    HEAD["a_d_control_%s" % _sl] = round(float(_D.at[_sl, "d_control"]), 3)
    HEAD["b_ratio_%s" % _sl] = round(float(_D.at[_sl, "ratio_to_control"]), 2)
for _sl in ("cb", "nz"):
    HEAD["d_ratio_every_site_%s" % _sl] = round(float(_RT.at[_sl, "ratio_site"]), 2)
    HEAD["d_ratio_one_per_protein_%s" % _sl] = round(float(_RT.at[_sl, "ratio_prot"]), 2)
HEAD["e_nz_error_median_angstrom"] = round(med, 2)
HEAD["e_n_measured_lysines"] = int(len(v))
HEAD["e_cb_to_nz_angstrom"] = CB_NZ
for (_arm, _pop), _vals in _F_ROWS.items():
    _who = "acquired" if _pop == "acquired union" else "unengaged"
    for _q, _val in zip(qs, _vals):
        HEAD["f_%s_%s_pct_inside_%d" % (_who, SLUG[_arm], _q)] = _val

st.save(fig, "figS9")
fd.head("figS9", HEAD)
import shutil
for e in ("pdf", "png"):
    shutil.copy("figS9.%s" % e, os.path.join(ROOT, "results", "figS9.%s" % e))
# HEADLINE NUMBERS, added 2026-09-04. This printed only its panel and arm counts, so the
# distances a legend has to state — the three legitimate anchorings against the mismatched one,
# and the median error of the built side chain — would have been copied from prose. They are
# printed now.
print("figS9 headline numbers")
for _k, _v in _F_ROWS.items():
    print("  f arm %s, %-22s inside the 50/75/90/95%% engaged region: %s"
          % (_k[0], _k[1], ", ".join("%.0f%%" % x for x in _v)))
for _, r in R.iterrows():
    print("  %-34s union-vs-engaged %.3f | union-vs-unengaged %.3f | control %.3f | ratio %.2f"
          % (r.arm, r.d_union_eng, r.d_union_not, r.d_control,
             r.d_union_eng / r.d_control if r.d_control else float("nan")))
_leg = R[~R.art]; _art = R[R.art]
if len(_leg) and len(_art):
    print("  legitimate arms   union-vs-engaged %s"
          % " / ".join("%.3f" % v for v in _leg.d_union_eng))
    print("  legitimate arms   control          %s"
          % " / ".join("%.3f" % v for v in _leg.d_control))
    print("  mismatched arm    %.3f against control %.3f  =  %.2fx"
          % (_art.d_union_eng.iloc[0], _art.d_control.iloc[0],
             _art.d_union_eng.iloc[0] / _art.d_control.iloc[0]))
print("  built side chain  median error %.2f A over %s measured lysines (CB->NZ is 6.4 A)"
      % (med, format(len(v), ",")))
print("  quantile series parsed %d | panels 6 | arms %d" % (len(rows), len(R)))
for k in sorted(_F_ROWS):
    print("  f %s  %-22s %s" % (k[0], k[1], _F_ROWS[k]))
