#!/usr/bin/env python3
"""
Figure 5, Ligandability is reproducible, and structurally unpredictable.

THE CENTRAL COMPARISON (panel a)
A structural model of ligandability reaches AUC 0.617 (95% CI 0.594-0.640), cross-validated by
PROTEIN so no test protein is ever seen in training. On its own that number says nothing: an AUC
is only interpretable against how reproducible the labels themselves are. Step 29 measures that
directly, split the compounds tested at each site into halves, and ask how well engagement in
one half predicts engagement in the other. That is the best any site-level classifier could do.

WHICH CEILING IS THE FAIR ONE
Three datasets give three ceilings, and they are not interchangeable:
  NatChem        median 41 compounds/site, 17.5% of sites engaged -> deepest split 20/half
  CyMPSA         median 10 compounds/site,  2.9% engaged          -> deepest split  5/half
  Meroterpenoid  median 10 compounds/site,  1.5% engaged          -> deepest split  5/half
Only NatChem is screened deeply enough to estimate a ceiling at a realistic engaged rate, and its
ceiling RISES with split depth (0.628 at 2/half to 0.723 at 20/half) because a deeper half is a
less noisy label. The shallow screens return higher numbers (0.88-0.94) at engaged rates of 1-3%,
where a handful of positives makes AUC both unstable and easy to inflate. NatChem at 20/half,
0.723 [0.706, 0.744], is therefore the binding comparator, and step 29 is explicit that all of
these are LOWER bounds. figS6 shows all three so the choice can be checked rather than trusted.

A CONTRADICTION IN THE REPO, RESOLVED HERE
step10_report.txt still carries an interpretation guide from an earlier pass: "the Part 1 ceiling
analysis put the maximum achievable AUC at 0.55-0.65 ... an AUC in that band means the model is at
the ceiling, not that it is weak." Read that way, 0.617 sits AT the ceiling and the conclusion
inverts, there would be nothing left to measure. Step 29 supersedes it with an empirical,
per-dataset estimate, and the gap is real. The stale guide should be struck from step10_report.
"""
import os, re, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
from matplotlib.patches import Rectangle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import abbasov_style as st
import figure_data as fd
# One body type size for this figure; every other size is a ratio to it.
BODY = 7.0
st.apply(BODY)
# 300, not the style module's 600. THIS SCRIPT NEVER SET IT, so its PNG came out at twice the
# canvas and every pixel measurement taken off it was out by a factor of two: 4322 px read as
# 14.407 in at 300 when the sheet was really 7.205 in. Figures 1 to 4 all pin it. One display
# pixel is one output pixel.
plt.rcParams["savefig.dpi"] = 300

def find(n):
    for d in (".", "..", "ligand_rebuild", "../ligand_rebuild"):
        p = os.path.join(d, n)
        if os.path.exists(p):
            return p
    raise FileNotFoundError(n)          # a hard failure, never a silent blank panel

R10 = open(find("step10_report.txt")).read()
R29 = open(find("step29_report.txt")).read()
R29B = open(find("step29b_report.txt")).read()
R16 = open(find("step16_report.txt")).read()
R13 = open(find("step13_report.txt")).read()
ROC_S = pd.read_csv(find("step10_roc_splits.csv"))
ROC_P = pd.read_csv(find("step10_roc_pooled.csv"))
PCS = [(int(m.group(1)), float(m.group(2)), float(m.group(3)))
       for m in re.finditer(r"PC(\d)\s+([\d.]+)%\s+cumulative\s+([\d.]+)%", R10)]

# ---------------------------------------------------------------- parse, then assert
def auc_row(txt, name):
    m = re.search(re.escape(name) + r"\s+([\d.]+)\s+\[([\d.]+), ([\d.]+)\]", txt)
    return tuple(float(x) for x in m.groups())

COMP = auc_row(R10, "COMPOSITE (test)")
TRAIN = auc_row(R10, "composite (train)")
BEST = auc_row(R10, "dist_centroid")
SINGLES = [(f, auc_row(R10, f)[0]) for f in
           ("dist_centroid", "plddt_site", "propka_pKa", "hydrophobic_8A", "rsa_tien")]

# ceiling curves: dataset -> [(compounds per half, n sites, % engaged, auc, lo, hi)]
CEIL, cur = {}, None
for line in R29.splitlines():
    m = re.match(r"^(\w+) — (\d+) sites \| compounds per site: median (\d+)", line)
    if m:
        cur = m.group(1); CEIL[cur] = dict(n=int(m.group(2)), med=int(m.group(3)), rows=[])
    m = re.match(r"^\s*(\d+) compounds/half \|\s*(\d+) sites \|\s*([\d.]+)% engaged \| "
                 r"ceiling AUC ([\d.]+) \[([\d.]+), ([\d.]+)\]", line)
    if m and cur:
        CEIL[cur]["rows"].append(tuple(float(x) for x in m.groups()))

# ceiling by structural-confidence band
BANDS, cur = {}, None
for line in R29B.splitlines():
    if re.match(r"^(NatChem|CyMPSA|Meroterpenoid)\s*$", line.strip()):
        cur = line.strip(); BANDS[cur] = []
    m = re.match(r"^\s*(< 50\s+disordered|50-70 low|70-90 confident|>= 90 very high|no structure)"
                 r"\s+n=\s*(\d+) \| deepest split\s+(\d+)/half \((\d+) sites\) \| "
                 r"ceiling ([\d.]+) \[([\d.]+), ([\d.]+)\] \|\s*([\d.]+)% engaged", line)
    if m and cur:
        BANDS[cur].append((m.group(1).strip(), int(m.group(2)), int(m.group(3)), int(m.group(4)),
                           float(m.group(5)), float(m.group(6)), float(m.group(7)),
                           float(m.group(8))))

# breadth: multivariate distance from the 0-chemotype group
BRE = []
for line in R16.splitlines():
    m = re.match(r"^\s*(\d|\d-\d|>=\d)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)", line)
    if m and "chemotypes TESTED" not in line:
        BRE.append((m.group(1), int(m.group(2)), float(m.group(3)), float(m.group(4))))
BRE = [b for b in BRE if b[0] in ("1", "2-3", "4-5", ">=6")]

# enantiomer pairs
ENA = []
for m in re.finditer(r"\((.)\)-([\w ]+?)\s+vs\s+\((.)\)-([\w ]+?)\s+n=\s*(\d+)\n"
                     r"\s+both\s+(\d+) \(expected ([\d.]+)\) \| only\(\+\)\s+(\d+) \| "
                     r"only\(-\)\s+(\d+) \| r = ([\d.]+)", R13):
    ENA.append(dict(name=m.group(2).strip(), n=int(m.group(5)), both=int(m.group(6)),
                    exp=float(m.group(7)), plus=int(m.group(8)), minus=int(m.group(9)),
                    r=float(m.group(10))))

NAT = CEIL["NatChem"]["rows"]
NAT_DEEP = NAT[-1]                      # 20 compounds/half
assert abs(COMP[0] - 0.617) < 1e-9 and abs(NAT_DEEP[3] - 0.723) < 1e-9, (COMP, NAT_DEEP)
assert len(BRE) == 4 and len(ENA) == 2, (BRE, ENA)
GAP = NAT_DEEP[3] - COMP[0]

HEAD = {
    "a_composite_auc": COMP[0], "a_composite_lo": COMP[1], "a_composite_hi": COMP[2],
    "a_train_auc": TRAIN[0], "a_optimism": round(TRAIN[0] - COMP[0], 3),
    "a_ceiling_auc": NAT_DEEP[3], "a_ceiling_lo": NAT_DEEP[4], "a_ceiling_hi": NAT_DEEP[5],
    "a_gap": round(GAP, 3),
    "a_ceiling_shallowest": NAT[0][3], "a_compounds_per_half_deepest": int(NAT_DEEP[0]),
    "a_best_single_auc": BEST[0],
    "a_composite_over_best_single": round(COMP[0] - BEST[0], 3),
}
fd.panel("fig5", "a", [dict(row="composite (test)", auc=COMP[0], lo=COMP[1], hi=COMP[2]),
                       dict(row="composite (train)", auc=TRAIN[0], lo=TRAIN[1], hi=TRAIN[2]),
                       dict(row="best single feature (dist_centroid)",
                            auc=BEST[0], lo=BEST[1], hi=BEST[2])] +
         [dict(row="NatChem ceiling, %d compounds/half" % int(r[0]), auc=r[3], lo=r[4], hi=r[5])
          for r in NAT],
         note="AUC and 2.5-97.5 percentile interval; ceiling series is the split-half ladder")

fd.panel("fig5", "b", [dict(dataset=ds, band=r[0], n_sites=r[1], deepest_split_per_half=r[2],
                            n_sites_at_split=r[3], ceiling_auc=r[4], lo=r[5], hi=r[6],
                            pct_engaged=r[7])
                       for ds, v in BANDS.items() for r in v],
         note="reproducibility ceiling within AlphaFold-confidence bands")
_nb = {r[0]: r for r in BANDS["NatChem"]}
HEAD["b_natchem_ceiling_min"] = round(min(r[4] for r in BANDS["NatChem"]), 3)
HEAD["b_natchem_ceiling_max"] = round(max(r[4] for r in BANDS["NatChem"]), 3)
for _k, _lab in ((">= 90 very high", "plddt90"),):
    if _k in _nb:
        HEAD["b_natchem_ceiling_%s" % _lab] = round(_nb[_k][4], 3)
        HEAD["b_natchem_n_sites_%s" % _lab] = int(_nb[_k][1])
for _ds in BANDS:
    HEAD["b_n_bands_%s" % _ds] = len(BANDS[_ds])

fd.panel("fig5", "d", [dict(component=p_[0], variance_pct=p_[1], cumulative_pct=p_[2])
                       for p_ in PCS], note="principal components of the 14 features")
HEAD["d_pc1_variance_pct"] = PCS[0][1]
HEAD["d_n_components_to_80pct"] = int(next(p_[0] for p_ in PCS if p_[2] >= 80))
HEAD["d_n_features"] = 14

fd.panel("fig5", "e", [dict(chemotypes=b[0], n_sites=b[1], observed_distance=b[2],
                            null_distance=b[3]) for b in BRE],
         note="multivariate distance from a never-engaged site; null over 2,000 label shuffles")
for _b in BRE:
    _t = _b[0].replace(">=", "ge").replace("-", "_")
    HEAD["e_observed_%s" % _t] = _b[2]
    HEAD["e_null_%s" % _t] = _b[3]
    HEAD["e_n_sites_%s" % _t] = _b[1]

fd.panel("fig5", "f", [dict(pair=e["name"], n_sites_measured=e["n"], engaged_by_both=e["both"],
                            expected_by_chance=e["exp"], only_plus=e["plus"],
                            only_minus=e["minus"], correlation=e["r"]) for e in ENA],
         note="enantiomer pairs with enough shared coverage to compare")
for _e in ENA:
    _t = _e["name"].replace(" ", "_")
    HEAD["f_%s_both" % _t] = _e["both"]
    HEAD["f_%s_one_hand" % _t] = _e["plus"] + _e["minus"]
    HEAD["f_%s_expected" % _t] = _e["exp"]
    HEAD["f_%s_n_sites" % _t] = _e["n"]
HEAD["f_n_pairs"] = len(ENA)

# ---------------------------------------------------------------- figure
# 0.63 -> 0.76. Panel a labels five model variants down its y axis, three of them on two
# lines, and at the raised body they ran into each other. Panel a is published as
# Figure 5C, so this is a main-text panel, not one of the three that get cut.
# ALL FOUR MARGINS IN INCHES, AND THE CANVAS SOLVED FROM THE PANEL BLOCK, which is how
# Figures 1 to 4 are built. Fractional margins shrink with the canvas and walk the figure
# into _fit_once, which scales every axes at once and takes the y axes off the family's
# 483.6 px. _PS is the sum of the three panel widths and is HELD: matplotlib spends a
# gridspec's width on panels and gaps together, so changing wspace against a held block
# width would resize all six panels instead of changing the gaps.
_ML, _MR = 2.00, 0.60          # left and right margins, inches. _ML is large because panel
                               # a names its four model variants down its own left side.
_MB, _MT = 0.90, 0.26          # bottom and top margins, inches
_PS = 9.039                    # sum of the three panel widths, held. Chosen so the trimmed
                               # sheet is 3461 px, which is fig1.png and fig4.png, so this figure
                               # finalizes at Figure 1's 12.00 in on the shared raster.
_WS, _HS = 0.30, 0.66          # column and row spacing, fractions of the average panel
_YAX = 1.55020                 # the family y axis, 483.6 px at the 311.96 raster
_PW = _PS * (1 + 2 * _WS / 3)
_PH = (2 + _HS) * _YAX
_FW, _FH = _PW + _ML + _MR, _PH + _MB + _MT
fig = plt.figure(figsize=(_FW, _FH))
# (text artist, bar patch) pairs whose x is settled only after st.save has fitted the
# bar thickness and re-centred the groups. See panel f.
_SNAP = []
gs = fig.add_gridspec(2, 3, hspace=_HS, wspace=_WS,
                      left=_ML / _FW, right=1 - _MR / _FW,
                      top=1 - _MT / _FH, bottom=_MB / _FH,
                      width_ratios=[1.22, 1.0, 1.0])

# ---------------------------------------------------------------- a: the central comparison
ax = fig.add_subplot(gs[0, 0]); st.head(fig, ax, "a")
ax._panel_head_flush = True
# A slides RIGHT to put its letter on D's, which is the move Figure 3 panel a makes in
# the other direction. Flush letters sit on their panel's leftmost ink, and a names its
# four model variants down its own left side while d only carries a y-axis title, so A
# printed 340 px left of D. The trim follows the ink, so the figure loses that 340 px
# off its left edge rather than gaining it as white.
# The divisor is 300, the raster this script renders at, so the numerator is the pixel
# count measured in fig5.png. Figures 1 to 4 write /311.96 here and tune the numerator
# against a 300 dpi render, which makes their numerators not pixels of anything.
ax._axes_dx_in = 150.0 / 300.0
# AND NARROWED BY THE SLIDE PLUS THE OVERLAP. a's x-axis spine ran 74 px past the left
# edge of b's y-axis title: the two clear each other only because they sit at different
# heights, which is not a gap. _axes_dw_in holds x0, so the panel shortens to the right.
ax._axes_dw_in = -188.0 / 300.0
# The two comparator rows are set in the darker neutral, not the light one: their value labels
# sit partly over the light-grey gap band, where light grey on light grey is unreadable.
rows = [("chance", 0.50, None, None, st.GREY),
        ("best single feature", BEST[0], BEST[1], BEST[2], st.GREY),
        ("all features", COMP[0], COMP[1], COMP[2], st.BLUE),
        # "how reproducible the labels are" over two lines becomes one word, and the word is the one
# panel b's y-axis title and the caption both already use.
        ("reproducibility", NAT_DEEP[3], NAT_DEEP[4], NAT_DEEP[5], st.ORANGE)]
y = np.arange(len(rows))[::-1]
# THE BAND STOPS ABOVE THE RAIL'S OWN LABELS. axvspan runs the full height of the axes, so
# the 2, the 20 and "compounds per half" were set on light grey: the rail spans 0.628 to
# 0.723 and the band spans 0.617 to 0.723, so a label at either end of the rail is inside it
# by construction and no horizontal move can get it out. The band covers the four rows and
# stops at -0.30; the labels and the gap arrow sit below it on white.
_A_TOP, _A_BOT = len(rows) - 0.45, -1.45      # the floor came up with the arrow
ax.add_patch(Rectangle((COMP[0], -0.20), NAT_DEEP[3] - COMP[0], _A_TOP + 0.20,
                       facecolor=st.LIGHTGREY, edgecolor="none", zorder=0))
for yy, (lab, v, lo, hi, c) in zip(y, rows):
    if lo is not None:
        # CAPPED, as panel b's are. A bare line through a dot reads as a connector, and these
        # are 95% intervals: the 2.5th to 97.5th percentile across 20 protein-held-out splits
        # for the two model rows and across 50 random half-splits for the ceiling.
        ax.errorbar([v], [yy], xerr=[[v - lo], [hi - v]], fmt="none", ecolor=c,
                    elinewidth=1.4, capsize=3.0, capthick=1.4, zorder=3)
    ax.plot([v], [yy], "o", color=c, ms=5.0, mew=0, zorder=4)
    ax.text(v, yy + 0.26, "%.3f" % v, ha="center", fontsize=st.sz(6.3), color=c)
# The ceiling is not one number: it rises with how deeply each half is screened, from 0.628 at
# 2 compounds per half to 0.723 at 20, and it is still climbing at the deepest split available.
# Putting the whole series on this axis lets the comparator justify itself here rather than only
# in figS6 a, and shows 0.723 for what it is — a lower bound.
_cy = y[-1]
# The rail sits inside the light-grey gap band, so it has to be the darker neutral or it
# disappears into the band.
ax.plot([r[3] for r in NAT], [_cy] * len(NAT), "-", color=st.GREY, lw=1.0, zorder=2)
for r in NAT[:-1]:
    # 2.2 -> 3.2 pt. These are five measurements, not texture; still well under the deepest
    # point's 5.0 pt, which is the one the row is reporting.
    ax.plot([r[3]], [_cy], "o", color=st.ORANGE, ms=3.2, mew=0, zorder=3)
# TICKED, NOT JUST LABELLED. Two numbers floating under a rail do not say which dot each one
# belongs to, and the rail's five dots span only 0.095 of AUC, so proximity decides nothing.
# A tick under each labelled end turns the rail into the small axis it already is.
for _x_, _l_ in ((NAT[0][3], "2"), (NAT_DEEP[3], "20")):
    # The band's own light grey, so the ticks read as the rail's furniture rather than as
    # two more orange marks competing with the five that carry the data.
    ax.plot([_x_, _x_], [_cy - 0.24, _cy - 0.44], color=st.LIGHTGREY, lw=1.1, zorder=3)
    ax.text(_x_, _cy - 0.50, _l_, ha="center", va="top", fontsize=st.sz(6.3), color=st.ORANGE)
# AND WHAT THE 2 AND THE 20 ARE. Two bare numbers under a rail say nothing; the rail is the
# ceiling measured at every screen depth the data support, and the depth is compounds in each
# half of the split. Centred on the rail's own midpoint, one line below the two ends.
ax.text((NAT[0][3] + NAT_DEEP[3]) / 2.0, _cy - 0.95, "compounds per half", ha="center",
        va="top", fontsize=st.sz(6.3), color=st.ORANGE)
# linespacing 1.35 put a two-line label at 90 px against a row pitch of 88.7, so the second
# line of "all 14 features together" and the first of "how reproducible the labels are" met.
# At 1.15 a two-line block is 76 px and the two clear each other by 13.
ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=st.sz(6.3), linespacing=1.15)
ax.set_xlabel("how well a site can be sorted\ninto engaged or not (AUC)",
              fontsize=st.sz(7), linespacing=1.4)
# dist_centroid's CI runs from 0.375 — below chance — because a single feature is unstable
# across protein-held-out splits. Truncating it would flatter the comparison, so the axis
# opens far enough to show it.
# The two-line gap statement below the arrow needs its own room; at -2.00 its second line
# ran along the x-axis spine.
# ROUND LIMITS ON A ROUND STEP. It ran 0.355 to 0.79 and matplotlib answered with ticks at
# 0.4 through 0.75. The lowest thing drawn is the single feature's interval, which starts
# at 0.375 because one feature is unstable across protein-held-out splits; truncating it
# would flatter the comparison. 0.3 to 0.8 on a step of 0.1 holds it and reads round.
ax.set_xlim(0.3, 0.8); ax.set_xticks([0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
ax.set_ylim(_A_BOT, _A_TOP)
# THE GAP ARROW IS GONE. It spanned 0.617 to 0.723, which is exactly what the grey block
# spans, so the panel was drawing one quantity twice, and once the block stopped short of the
# rail's labels the arrow sat 1.4 units below it with nothing connecting the two. Caption (A)
# states the gap as 0.106.
for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
ax.tick_params(axis="y", length=0)

# ---------------------------------------------------------------- b: ceiling vs structural order
ax = fig.add_subplot(gs[0, 1]); st.head(fig, ax, "b")
ax._panel_head_flush = True
ax._axes_dx_in = 132.0 / 300.0     # row 0 evened
# THE KEY CITES, IT DOES NOT NAME. "NatChem", "CyMPSA" and "meroterpenoids" are this lab's
# internal names for three published screens and mean nothing to a reader. Author and year,
# not a reference number: the numbers are assigned at build time by order of first appearance
# in the manuscript, so a number typed in here would go stale silently the first time a
# citation is added ahead of it. The caption carries the full identification.
CITE = {"NatChem": "Abbasov 2021", "CyMPSA": "Ryan 2024", "Meroterpenoid": "Bracken 2024"}
ORDER = ["< 50  disordered", "50-70 low", "70-90 confident", ">= 90 very high", "no structure"]
NICE = ["< 50\ndisordered", "50–70\nlow", "70–90\nconfident", "≥ 90\nvery high", "no\nstructure"]
# (colour, marker). All three lines are SOLID now. They were solid, dashed and dotted, which
# is a second encoding on top of colour and marker and costs the dashed and dotted series
# their legibility at this line width. CyMPSA takes the house dark #46494C and the
# meroterpenoids take teal #1985A1.
# NOTE: teal carries "the composite model" in panels a, c, d and e. In this panel it now
# carries a screen instead, and the model's own line has gone grey. That is a colour meaning
# two things inside one figure, which is worth a look before the figure is finalized.
# (colour, marker, z). The meroterpenoids are only in two bands and both of them are bands
# CyMPSA also occupies, within a hundredth of an AUC, so drawn in dict order the teal series
# went under the dark one and read as two stray marks. +3 puts it on top of both.
STYLE = {"NatChem": (st.ORANGE, "o", 0), "CyMPSA": (st.GREEN, "s", 0),
         "Meroterpenoid": (st.BLUE, "^", 3)}
for ds, rowset in BANDS.items():
    c, mk, _z = STYLE[ds]
    xs = [ORDER.index(r[0]) for r in rowset if r[0] in ORDER]
    ys = [r[4] for r in rowset if r[0] in ORDER]
    los = [r[5] for r in rowset if r[0] in ORDER]
    his = [r[6] for r in rowset if r[0] in ORDER]
    # CAPPED. The intervals were bare lines, so an interval end and a gap in the series read
    # alike; a short perpendicular at each end says where the interval stops.
    ax.errorbar(xs, ys, yerr=[np.array(ys) - np.array(los), np.array(his) - np.array(ys)],
                fmt="none", ecolor=c, elinewidth=0.9, capsize=2.5, capthick=0.9, zorder=2 + _z)
    ax.plot(xs, ys, color=c, lw=1.0, ls="-", zorder=3 + _z)
    ax.plot(xs, ys, mk, color=c, ms=4.5, mew=0, zorder=4 + _z, label=CITE[ds])
ax.axhline(COMP[0], color=st.GREY, lw=1.0, ls=(0, (3, 1.6)))
# UNDER THE RULE AND CENTRED ON IT. Above it and hard left the label sat in the corner where
# the first screen's interval reaches down to 0.583; at x = 2.0 it spans 1.35 to 2.65 and the
# nearest thing drawn below the rule is that same interval, four categories away.
ax.text(2.0, COMP[0] - 0.012, "the model, 0.617", va="top", ha="center", fontsize=st.sz(6.3),
        color=st.GREY)
ax.set_xticks(range(5)); ax.set_xticklabels(NICE, fontsize=st.sz(6.3), linespacing=1.3)
ax.set_xlabel("AlphaFold confidence at the lysine", fontsize=st.sz(7))
ax.set_ylabel("reproducibility\nceiling (AUC)", fontsize=st.sz(7), linespacing=1.4)
ax.set_xlim(-0.5, 4.5)
# ON A TICK AT BOTH ENDS. It ran 0.36 to 1.05 with ticks every 0.1, so the bottom sat 0.04
# below the lowest tick. Everything drawn lies inside 0.5 to 1.0: the lowest interval end is
# NatChem disordered at 0.583, the highest the meroterpenoids' confident band at 0.992, and
# the model's reference line is 0.617.
ax.set_ylim(0.4, 1.0); ax.set_yticks([0.4, 0.6, 0.8, 1.0])
# ONE ROW ALONG THE FOOT. Stacked at centre left the key sat across the error bars of the two
# upper series at the first two bands, which the collision gate cannot see because it compares
# text to text. The band below 0.57 is the only clear space in this panel: the lowest thing
# drawn is Abbasov's interval at 0.583 and the model's rule is at 0.617. Three entries in one
# row measure about 856 px against a panel of 1015.
ax.legend(frameon=False, fontsize=st.sz(6.3), loc="lower center", bbox_to_anchor=(0.5, 0.05),
          ncol=3, handlelength=1.4, handletextpad=0.4, columnspacing=0.8, borderaxespad=0.0)
# MOVED TO THE CAPTION, 7 November. Four lines of prose stood inside the frame here. The
# artwork carries the data, the axis labels and the key; a conclusion is prose and reads as
# prose. Caption (B) ends "Panel a therefore cannot be dismissed as an artefact of disordered
# sites", which is the same sentence, so nothing is lost and the panel gains its lower half.

# ---------------------------------------------------------------- e: breadth across chemotypes
ax = fig.add_subplot(gs[1, 1]); st.head(fig, ax, "e")
ax._panel_head_flush = True
ax._axes_dx_in = -72.0 / 300.0      # row 1 evened against row 0's right edge
# SHORTER, by the same 12% panel f took. _axes_dw_in holds x0, so it shortens to the right.
ax._axes_dw_in = -0.12 * _PS / 3.22
x = np.arange(len(BRE)); w = st.grp(2)
obs = [b[2] for b in BRE]; nul = [b[3] for b in BRE]
_EB = ax.bar(x - w/2, obs, w, color=st.BLUE, label="observed")
_EN = ax.bar(x + w/2, nul, w, color=st.LIGHTGREY, label="2,000 label shuffles")
for i, b in enumerate(BRE):
    # BOTH BARS CARRY THEIR VALUE. The null was drawn and left unlabelled, so the one number
    # the reader has to compare the observed against had to be estimated off the axis. It is
    # set in the darker neutral rather than the bar's own light grey, which is unreadable on
    # white; panel a's comparator rows make the same swap for the same reason.
    # Snapped to the drawn patch: st.save re-centres a grouped container inside its slot
    # after these are placed, so a label written at the x its bar was drawn at ends up
    # beside it. See panel f.
    _SNAP.append((ax.text(i - w/2, b[2] + 0.06, "%.2f" % b[2], ha="center",
                          fontsize=st.sz(6.3), color=st.BLUE), _EB.patches[i]))
    _SNAP.append((ax.text(i + w/2, b[3] + 0.06, "%.2f" % b[3], ha="center",
                          fontsize=st.sz(6.3), color=st.GREY), _EN.patches[i]))
    # The n used to hang at y = -0.20, below the lowest tick, which is what held this axis
    # off 0. It belongs under its own category, which is the house form for a count with
    # its denominator.
ax.set_xticks(x)
ax.set_xticklabels(["%s\nn=%s" % (b[0].replace(">=", "≥"), f"{b[1]:,}") for b in BRE],
                   fontsize=st.sz(6.3), linespacing=1.35)
ax.set_xlabel("chemotypes that engage the site\n(of 34 tested)",
              fontsize=st.sz(7), linespacing=1.4)
ax.set_ylabel("distance from a\nnever-engaged site", fontsize=st.sz(7), linespacing=1.4)
# ON A TICK AT BOTH ENDS. It ran -0.32 to 4.03 to hold the n labels that have just moved
# under the categories. The bars run 0.120 to 2.650 and the value labels sit 0.06 above
# them, so 0 to 3 on a step of 1 holds everything.
ax.set_ylim(0, 3); ax.set_yticks([0, 1, 2, 3])
ax.legend(frameon=False, fontsize=st.sz(6.3), loc="upper left", handlelength=1.1,
          bbox_to_anchor=(0.0, 1.0), borderaxespad=0.0)

# ---------------------------------------------------------------- f: enantiomers
ax = fig.add_subplot(gs[1, 2]); st.head(fig, ax, "f")
ax._panel_head_flush = True
ax._axes_dx_in = 133.0 / 300.0      # row 1 evened against row 0's right edge
# SHORTER. Two categories were spread over a 502 px pitch against a floor of about 400,
# which is the wider of the two tick labels plus clearance, so the panel was carrying
# white rather than information. _axes_dw_in holds x0, so it shortens to the right.
ax._axes_dw_in = -0.12 * _PS / 3.22
lab, both, one, exp = [], [], [], []
for e in sorted(ENA, key=lambda z: -z["both"]):
    lab.append("%s\n(%s sites measured)" % (e["name"], f"{e['n']:,}"))
    both.append(e["both"]); one.append(e["plus"] + e["minus"]); exp.append(e["exp"])
x = np.arange(len(lab)); w = st.grp(2)
_FB = ax.bar(x - w/2, both, w, color=st.BLUE, label="engaged by both enantiomers")
_FO = ax.bar(x + w/2, one, w, color=st.ORANGE, label="engaged by one enantiomer")
for i in range(len(lab)):
    # SNAPPED TO THE DRAWN BAR, not to the x it was drawn at. st.save re-fits bar thickness
    # to the page and RE-CENTRES a grouped container inside its slot, which happens after
    # these labels are placed, so a label written at i - w/2 ends up beside its bar rather
    # than over it: measured here at 19 px for a 75 px bar. _SNAP carries the pair and the
    # label is moved onto the patch's own centre once the fit has run.
    _SNAP.append((ax.text(i - w/2, both[i] + 0.45, "%d" % both[i], ha="center",
                          fontsize=st.sz(6.3), color=st.BLUE), _FB.patches[i]))
    _SNAP.append((ax.text(i + w/2, one[i] + 0.45, "%d" % one[i], ha="center",
                          fontsize=st.sz(6.3), color=st.ORANGE), _FO.patches[i]))
    # A rule at 0.1 sites is invisible next to a bar of 19, so the expectation is stated
    # rather than drawn: the whole point is that it is two orders of magnitude smaller.
    # Set inside the bar it ran wider than the bar itself and the overflow was white on
    # white paper, so both expectations sit above their bar, above the count.
    # ONE LINE, AND SHORTER. Two lines reached 4.1 units up a 30-unit axis from a label that
    # already starts above the count, which put its top inside the key now that the key is at
    # the top of the panel. One line of "vs 0.1 by chance" is 266 px against the panel's 900,
    # so the two groups' labels also clear each other, which a one-line "vs 0.1 expected by
    # chance" at 415 px did not.
    _SNAP.append((ax.text(i - w/2, both[i] + 2.6, "vs %.1f by chance" % exp[i],
                          ha="center", va="bottom", fontsize=st.sz(6.3),
                          color=st.BLUE), _FB.patches[i]))
ax.set_xticks(x); ax.set_xticklabels(lab, fontsize=st.sz(6.3), linespacing=1.35)
ax.set_ylabel("sites engaged at\n≥60% occupancy", fontsize=st.sz(7), linespacing=1.4)
# 0 to 30 on a step of 10. The tallest bar is 19 and the expectation label above it reaches
# 26.2, so 30 holds everything drawn and drops the empty tier 40 was carrying.
ax.set_ylim(0, 30); ax.set_yticks([0, 10, 20, 30])
ax.legend(frameon=False, fontsize=st.sz(6.3), loc="upper right", handlelength=1.1,
          bbox_to_anchor=(1.0, 1.0), borderaxespad=0.0)
# MOVED TO THE CAPTION, 7 November. Four lines of prose stood across the top of this panel,
# over the key and over the two expectation labels; it is caption (F), which already says
# "Each pair overwhelmingly engages the same sites" and "cannot be intrinsic reactivity".


# ---------------------------------------------------------------- c: is it weak everywhere?
# The AUC says the model is weak on average. A reviewer's next hope is that it is sharp at high
# specificity — useful for picking a shortlist even if useless overall. It is not: at 5% false
# positives it recovers 11% of engaged sites, at 10% it recovers 19%. Predictions regenerated
# 2026-09-02 (step10_composite.py saves none); the re-run reproduces 0.617 exactly.
ax = fig.add_subplot(gs[0, 2]); st.head(fig, ax, "c")
ax._panel_head_flush = True
# SHORTER AGAIN, 30% off its column against the 12% of the first cut. At 12% the frame was
# 741 px and the x-axis title, 39 characters at 8.860 pt, is 721, so the title was the floor
# and there was nothing left to give. Breaking it over two lines takes the floor to 444 and
# the frame can come to 589. _axes_dw_in holds x0, so it shortens to the right.
ax._axes_dw_in = -0.30 * _PS / 3.22
ax._axes_dx_in = 252.0 / 300.0     # and back out to row 0's right edge
ax.plot([0, 1], [0, 1], color=st.LIGHTGREY, lw=0.9, ls=(0, (3, 2)), zorder=1)
ax.fill_between(ROC_S.fpr, ROC_S.lo_tpr, ROC_S.hi_tpr, color=st.LIGHTGREY, zorder=2,
                label="95% of 20 splits")
# Key labels shortened. A ROC curve leaves only the lower right free, and at the raised
# body the full labels were wider than that corner, so the key reached back across the
# plot into the 11% and 19% operating-point markers. Legend (C) gives both in full.
ax.plot(ROC_P.fpr, ROC_P.tpr, color=st.BLUE, lw=1.6, zorder=4, label="each site once")
for _fp in (0.05, 0.10):
    _tp = float(np.interp(_fp, ROC_P.fpr, ROC_P.tpr))
    ax.plot([_fp], [_tp], "o", color=st.ORANGE, ms=4.5, mew=0, zorder=6)   # the house size
    # The two operating points are only 0.08 apart in recovery; put one label under its
    # point and one over so they read as two numbers rather than one stack.
    # STRAIGHT UP, CENTRED, AND FAR ENOUGH TO CLEAR THE CURVE. Both operating points sit in
    # the bottom-left corner where the curve climbs steeply out of the origin, so there is no
    # room beside them: down and right puts the label on the x-axis spine and its 0.00 tick,
    # left puts it on the y-axis spine and 0.25, and up and right is where the curve itself
    # runs. The collision gate sees none of this, because it compares text to text.
    # 11 pt clears the curve at both points: the curve reaches 0.155 across the first label's
    # own width and 0.235 across the second's, against label floors of 0.20 and 0.284.
    # A LEADER EACH, so the label no longer has to sit close enough to be read by proximity
    # alone. That constraint is what had both of them wedged into the corner where the curve
    # climbs out of the origin, and it moved them three times. 11% goes right of its dot, into
    # the empty wedge under the chance diagonal, which at x 0.164 to 0.262 runs well above the
    # label's own 0.073 to 0.141. 19% goes straight up over the curve, which across that
    # label's width reaches 0.26 against a floor of 0.363.
    _right = _fp < 0.075
    ax.annotate("%.0f%%" % (100 * _tp), (_fp, _tp),
                xytext=((16, 0) if _right else (0, 20)), textcoords="offset points",
                fontsize=st.sz(6.3), color=st.ORANGE,
                ha="left" if _right else "center",
                va="center" if _right else "bottom",
                arrowprops=dict(arrowstyle="->", color=st.ORANGE, lw=0.8,
                                shrinkA=1, shrinkB=3))
fd.panel("fig5", "c", [dict(false_positive_rate=float(_f),
                            engaged_recovered=float(np.interp(_f, ROC_P.fpr, ROC_P.tpr)))
                       for _f in np.linspace(0, 1, 201)],
         note="pooled ROC of the composite, every site predicted once")
for _fpv in (0.05, 0.10):
    HEAD["c_recovered_at_%dpct_fpr" % int(100 * _fpv)] = round(
        100 * float(np.interp(_fpv, ROC_P.fpr, ROC_P.tpr)), 0)
ax.set_xlabel("false positives\n(of sites never engaged)", fontsize=st.sz(7), linespacing=1.4)
ax.set_ylabel("engaged sites recovered", fontsize=st.sz(7))
ax.set_xlim(0, 1); ax.set_ylim(0, 1.0)
ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0]); ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
ax.legend(frameon=False, fontsize=st.sz(6.3), loc="lower right", handlelength=1.2,
          handletextpad=0.4, borderaxespad=0.4)

# ---------------------------------------------------------------- d: would more features help?
# The reflex objection to a modest AUC is that the feature set is too small. It is not a
# feature-count problem: 14 features carry about five independent properties.
ax = fig.add_subplot(gs[1, 0]); st.head(fig, ax, "d")
ax._panel_head_flush = True
# MUCH SHORTER, 35% off its column against the 12% e and f took. Six bars need little width
# and this panel was the widest in the figure because its column carries panel a's four row
# names above it. What sets the floor is the x-axis title, 28 characters at 8.860 pt, which
# is 465 px under a frame of 668. _axes_dw_in holds x0, so it shortens to the right.
ax._axes_dw_in = -0.35 * _PS * 1.22 / 3.22
_x = [p[0] for p in PCS]; _ind = [p[1] for p in PCS]; _cum = [p[2] for p in PCS]
ax.bar(_x, _ind, color=st.LIGHTGREY, width=st.BARW, label="each property on its own")
ax.plot(_x, _cum, "-o", color=st.BLUE, lw=1.3, ms=4.5, mew=0, label="running total")
ax.axhline(80, color=st.ORANGE, lw=1.0, ls=(0, (3, 1.6)))
# The "80%" label is gone from the artwork. The rule stays, and caption (D) already reads
# "five clear the 80% line", so the number is still printed by the script that draws it and
# still reaches the reader.
ax.set_xticks(_x); ax.set_xticklabels(_x, fontsize=st.sz(6.3))
ax.set_xlabel("independent properties kept", fontsize=st.sz(7))
ax.set_ylabel("how much is captured (%)", fontsize=st.sz(7))
ax.set_xlim(0.4, len(_x) + 0.45)      # the 0.15 that held the 80% label is gone with it
ax.set_ylim(0, 100); ax.set_yticks([0, 20, 40, 60, 80, 100])
ax.legend(frameon=False, fontsize=st.sz(6.3), loc="center right", handlelength=1.2)
# MOVED TO THE CAPTION, 7 November. Two lines of prose ran across the top of this panel and
# over panel e's y-axis title beside it. Caption (D) carries both sentences already.

st.save(fig, "fig5")
# The bar fit has run, so every value label can now be put on its own bar.
for _t, _q in _SNAP:
    _t.set_position((_q.get_x() + _q.get_width() / 2.0, _t.get_position()[1]))
# Trim to the drawing, 18 px all round, the figures-final standard Figure 1 set and Figures 2,
# 3 and 4 follow. Measured off the RENDERED PIXELS rather than get_tightbbox, because an
# axis-off panel still contributes its rectangle to the tight bbox.
_PAD_IN = 18.0 / 300.0
from PIL import Image as _Im
_px = np.asarray(_Im.open("fig5.png").convert("RGB")).astype(int)
_ys, _xs = np.nonzero(_px.sum(2) < 735)
_H_in = fig.get_size_inches()[1]
_trim = Bbox.from_extents(_xs.min() / 300.0 - _PAD_IN,
                          _H_in - _ys.max() / 300.0 - _PAD_IN,
                          _xs.max() / 300.0 + _PAD_IN,
                          _H_in - _ys.min() / 300.0 + _PAD_IN)
for _ext in (".png", ".pdf"):
    fig.savefig("fig5" + _ext, bbox_inches=_trim)
print("trimmed to %.3f x %.3f in (%.0f x %.0f px at 300 dpi), 18 px pads"
      % (_trim.width, _trim.height, _trim.width * 300, _trim.height * 300))
fd.head("fig5", HEAD)

if __import__("os").environ.get("FIG5_FINAL"):
    # THE FAMILY SHARES A RASTER, NOT A WIDTH. Figure 1 at 12.00 in fixes it at
    # 300 x 12.00 / 11.5400 = 311.96 dpi; every figure renders at that raster and is declared
    # 300 dpi, so all of them print 9.213 pt body type whatever their width.
    import PIL.Image as _PI, pymupdf as _mu
    _ref = _mu.open("fig1.pdf")[0].rect.width / 72.0
    _dpi = 300.0 * 12.0 / _ref
    fig.savefig("fig5_final.png", dpi=_dpi, bbox_inches=_trim)
    _im = _PI.open("fig5_final.png")
    _im.save("fig5_final.png", dpi=(300, 300))
    print("final: native %.4f in at the shared %.2f dpi raster -> %d x %d px declared 300 dpi "
          "(%.3f x %.3f in), body type %.3f pt against Figure 1's %.3f"
          % (_trim.width, _dpi, _im.size[0], _im.size[1], _im.size[0] / 300.0,
             _im.size[1] / 300.0, 8.860 * 12.0 / _ref, 8.860 * 12.0 / _ref))
print("composite AUC %.3f [%.3f, %.3f] | NatChem ceiling %.3f [%.3f, %.3f] | gap %.3f"
      % (COMP + NAT_DEEP[3:6] + (GAP,)))
print("ceiling curves parsed:", {k: len(v["rows"]) for k, v in CEIL.items()})
print("bands parsed:", {k: len(v) for k, v in BANDS.items()})
print("breadth bins:", BRE)
print("enantiomer pairs:", [(e["name"], e["both"], e["plus"] + e["minus"], e["exp"]) for e in ENA])
print("c: operating points  " + " | ".join(
      "%.0f%% FPR -> %.1f%% recovered" % (100*f, 100*float(np.interp(f, ROC_P.fpr, ROC_P.tpr)))
      for f in (0.05, 0.10)))
print("d: PCs %d | individual %% %s | cumulative %% %s"
      % (len(PCS), [round(p[1], 1) for p in PCS], [round(p[2], 1) for p in PCS]))
print("d: PCs to clear 80%%: %d" % next(p[0] for p in PCS if p[2] >= 80))
_bb = [(ds, r[0], r[4], r[5], r[6], r[1]) for ds, v in BANDS.items() for r in v]
print("b: ceiling per band %.3f-%.3f over %d band-dataset points" %
      (min(b[2] for b in _bb), max(b[2] for b in _bb), len(_bb)))
for b in _bb:
    print("     %-14s %-18s ceiling %.3f [%.3f, %.3f]  n=%d" % b)
