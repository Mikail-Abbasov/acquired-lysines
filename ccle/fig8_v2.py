#!/usr/bin/env python3
"""
fig8_v2.py, FIGURE 8. Acquired lysines sit in ordinary lysine environments.

THE CLAIM, worded to what the data support
Not "indistinguishable" -- the residual separation is real and survives every check. What the
data support is: acquired lysines occupy the same structural space as engaged lysines; 82% fall
inside the region containing 90% of engaged lysines against 89% for measured lysines that were
labelled but never engaged; they sit closer to engaged lysines than to labelled-but-unengaged
ones; and the residual excess is 1.2x the between-population yardstick once each protein votes
once. The barrier is detection, not chemistry.

NO SENTENCES INSIDE PANELS. Per the 2026-09-02 rule adopted for Figure 6 -- a figure is looked
at, not read. Everything explanatory belongs in the legend.

PANELS
  a  the anchor: why both sides must be measured from the same atom     step21 / step24
  b  the three populations, feature by feature                          step28c
  c  overlap -- do acquired lysines fall inside the engaged region      step28
  d  real crystal structures: does the new lysine move the site         step31
  e  AKT1 E17K -- drugged in vivo, invisible here. THE CLOSING IMAGE.   acquired_K_searched...

ROBUSTNESS (the same comparison under three anchorings, plus the artefact) is figS9.

INPUT   ligand_rebuild/step28c_medians.csv, step28_results.csv, step31_crystal.csv,
        ccle_plex05/acquired_K_searched_vs_detected.csv
OUTPUT  fig8.pdf, fig8.png
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.transforms import Bbox

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import abbasov_style as st
import figure_data as fd
# One body type size for this figure; every other size is a ratio to it.
BODY = 7.0
st.apply(BODY)
# 300, not the style module's 600, and AFTER apply(), which resets the rcParams. fig8.png
# was 4322 x 2593, twice its canvas, for exactly this reason.
plt.rcParams["savefig.dpi"] = 300
# FIVE NAMES, FIVE COLOURS. GREY is the charcoal SERIES colour and INK is the slate
# ANNOTATION colour: they sit 8.5 apart and never carry two series in one panel.
ORANGE, BLUE, LIGHTGREY = st.ORANGE, st.BLUE, st.LIGHTGREY
GREY = st.GREEN                # charcoal, "labeled, not engaged" and "other substitutions"
INK  = st.GREY                 # slate, every rule label, count and annotation

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
LR = os.path.join(ROOT, "ligand_rebuild")
MED = pd.read_csv(os.path.join(LR, "step28c_medians.csv"))
RES = pd.read_csv(os.path.join(LR, "step28_results.csv"))
XT  = pd.read_csv(os.path.join(LR, "step31_crystal.csv"))
AK  = pd.read_csv(os.path.join(HERE, "acquired_K_searched_vs_detected.csv"))
AK  = AK[AK.gene == "AKT1"]

armA = RES[RES.arm.str.startswith("A")].iloc[0]

# P(a random acquired lysine sits closer to the engaged centre than a random engaged one).
# This was typed into panel c as "P = 0.443". It is computed by step28d and is now read from
# that report, so the artwork cannot drift from the analysis.
import re as _re
_r28d = open(os.path.join(LR, "step28d_report.txt")).read()
P_CLOSER = float(_re.search(r"=\s*([\d.]+)\s+\(0\.5 would mean", _r28d).group(1))

# AKT1 E17K: how many peptides did the search actually contain, and in how many plexes?
# The manuscript said "seven peptides in the 42 experiments". Counted from the search
# databases: the mutation is carried by one cell line, so it is searchable in one plex.
import glob as _glob
_akt = {}
for _f in sorted(_glob.glob(os.path.join(HERE, "plex*_varpep.fasta")) +
                 _glob.glob(os.path.join(HERE, "plex*_downpep.fasta"))):
    _plex = os.path.basename(_f).split("_")[0]
    _hdr = None
    for _ln in open(_f):
        if _ln.startswith(">"):
            _p = _ln[1:].strip().split("|")
            _hdr = _p if (len(_p) > 5 and _p[0] in ("VAR", "DOWN")
                          and _p[1] == "AKT1" and _p[5] == "17") else None
        elif _hdr:
            _akt.setdefault(_plex, {})[_ln.strip()] = _hdr[0]
            _hdr = None
AKT_PEPS = sorted({q for v in _akt.values() for q in v})
AKT_PLEXES = sorted(_akt)
AKT_VAR = sum(1 for v in _akt.values() for q, k in v.items() if k == "VAR")
AKT_DOWN = sum(1 for v in _akt.values() for q, k in v.items() if k == "DOWN")
# If the search databases are not present -- a deposit without plex*_{varpep,downpep}.fasta --
# these counts come out as zero, and zero here reads as "the site was never searched", which is
# the opposite of what the paper claims and is the number the caption quotes. Refuse to emit a
# count rather than emit a wrong one.
AKT_OK = bool(_akt)
if not AKT_OK:
    print("  WARNING fig8e: no plex*_varpep.fasta / plex*_downpep.fasta found, so the AKT1 E17K "
          "peptide count cannot be computed.\n"
          "          The deposited value is 3 peptides (2 variant, 1 downstream) in plex36. "
          "No count is written to the ledger.")

# ================================================================ geometry
# THE HOUSE STANDARD, the same one figures 1 and 4 to 7 are built on.
#   every drawn y axis is _YAX = 483.6 px
#   three type sizes only: 13.733 letter, 8.860 = sz(7), 7.974 = sz(6.3)
#   margins in INCHES and the canvas solved from the panel block, never fractions
#   each row's columns are [panel, gap, panel, ..., slack] at wspace 0, so a panel's
#   width and the space beside it are set independently, in inches
#   the gap is set on the INK, which means each panel's y labels hanging left of its
#   frame and its last x tick label hanging right of it are part of the arithmetic
# The figure does not have to span the page: figures 2, 3 and 7 do not either.
_ML, _MR = 1.60, 1.20
_MB, _MT = 1.60, 0.80
_YAX = 1.55020                 # 483.6 px, the family's y axis
_PW = 10.5455                  # the canvas the panels sit on; the trim follows the ink
_HS = 0.62                     # row gap, as a fraction of a panel

# Panel frames, in inches. b's follows its LABEL PITCH: 9 slots across 544 px put the
# angled feature names 56.7 px apart, which is figure 6b's and 7b's pitch exactly.
_AW, _BW_, _CW = 2.25, 1.81, 1.75   # a carries four two-line names under four bars
_DW = 3.20                     # narrowed from 3.9; e takes what is left
_GINK = 100.0                  # ONE ink gap everywhere, px at 300 dpi, as in figure 7
# Overhangs, measured off the render and written here so the gaps are arithmetic rather
# than a nudge. _OL_x is x's ink left of its frame, _OR_x its ink right of it.
_OL_A, _OR_A = 0.393, 0.000
_OL_B, _OR_B = 0.630, 0.000
_OL_C, _OR_C = 0.597, 0.030
_OL_D, _OR_D = 0.473, 0.063
# e's render has white at its own edges, so the ink the gate measures sits INSIDE the frame
# on both sides: these two are negative for that reason, not as a nudge.
_OL_E, _OR_E = -0.0633, -0.16
# d sits 24 px right of its cell so that its letter stands in the same column as a's: a's y
# numbers are two characters and d's are one, so their frames do not start at the same ink.
_DXD = 24.0 / 300.0
_G_AB = (_GINK / 300.0) + _OL_B + _OR_A
_G_BC = (_GINK / 300.0) + _OL_C + _OR_B
_G_DE = (_GINK / 300.0) + _OL_E + _OR_D + _DXD
# ROW 1: d IS SET, e TAKES THE REST, AND ROW 1'S HEIGHT FOLLOWS FROM e. Panel e is a rendered
# structure with a fixed aspect, so once its width is known its height is too, and row 1's
# cell is that height rather than the family's axis height. Panel d's own axes is shortened
# back to 483.6 px inside the taller cell, which is what figure 7's panel b does. The two rows
# still end at the same x because e's width is the remainder.
_ROW0 = _OL_A + _AW + _G_AB + _BW_ + _G_BC + _CW + _OR_C
_EIMG = os.path.join(HERE, "results", "structures", "akt1_8uvy_chimerax.png")
_EPIX0 = plt.imread(_EIMG)
_IH0, _IW = _EPIX0.shape[0], _EPIX0.shape[1]
_EW = _ROW0 - _OL_D - _DXD - _G_DE - _DW - _OR_E
# CROPPED TO PANEL D'S AXIS. The render is 2.165 wide to 1 tall, which at this width stood
# taller than the 0-to-7 axis beside it. The crop takes the difference off the top, where the
# render carries ribbon and nothing else, and leaves the site centred on the zinc.
_IHC = int(round(_YAX * _IW / _EW))
_ETOP = max(0, min(_IH0 - _IHC, 584 - _IHC // 2))
_EPIX = _EPIX0[_ETOP:_ETOP + _IHC]
_IH = _IHC
_Y1 = _YAX + 2.4 * (41.5 / 300.0)              # the image, now d's own axis, plus e's band
_S0 = _PW - _AW - _G_AB - _BW_ - _G_BC - _CW
_S1 = _PW - _DW - _G_DE - _EW

_HG = _HS * _YAX
_PH = _YAX + _HG + _Y1
_FW, _FH = _PW + _ML + _MR, _PH + _MB + _MT
fig = plt.figure(figsize=(_FW, _FH))
gt = fig.add_gridspec(1, 6, left=_ML / _FW, right=1 - _MR / _FW,
                      top=1 - _MT / _FH, bottom=1 - (_MT + _YAX) / _FH,
                      wspace=0, width_ratios=[_AW, _G_AB, _BW_, _G_BC, _CW, _S0])
gb = fig.add_gridspec(1, 4, left=_ML / _FW, right=1 - _MR / _FW,
                      top=(_MB + _Y1) / _FH, bottom=_MB / _FH, wspace=0,
                      width_ratios=[_DW, _G_DE, _EW, _S1])

# ---------------------------------------------------------------- a  the anchor
ax = fig.add_subplot(gt[0, 0]); st.head(fig, ax, "a")
ax._panel_head_flush = True
XS = [0.0, 1.0, 2.2, 3.2]
# READ, NOT TYPED. These four were literals, [12, 6, 12, 11], in a panel whose y axis is a
# count: the one rule this figure set has about numbers is that the script that draws a value
# has to compute it. The acquired side is the union's own crowding_8A median; the two measured
# bars are step 24's NZ-anchored and CB-anchored medians for the same lysines.
_A_UNION = float(MED[(MED.population == "acquired union")
                     & (MED.feature == "crowding_8A")]["median"].iloc[0])
_m24 = _re.search(r"crowding_8A\s+([\d.]+)\s+([\d.]+)\s+[+\-][\d.]+",
                  open(os.path.join(LR, "step24_report.txt")).read())
_A_NZ, _A_CB = float(_m24.group(1)), float(_m24.group(2))
VALS = [_A_UNION, _A_NZ, _A_UNION, _A_CB]
COLS = [ORANGE, LIGHTGREY, ORANGE, BLUE]
# NO BLACK EDGE. Black is matplotlib's text ink, not a palette entry, and an outline on a
# fill is a sixth colour on the sheet.
ax.bar(XS, VALS, 0.62, color=COLS, lw=0, zorder=3)
for x_, v in zip(XS, VALS):
    ax.text(x_, v + 0.35, "%g" % v, ha="center", va="bottom", fontsize=st.sz(6.3), color=INK)
for i0_, i1_, c_ in ((0, 1, ORANGE), (2, 3, BLUE)):
    x0_, x1_ = XS[i0_], XS[i1_]
    nm_ = "gap %g" % (VALS[i0_] - VALS[i1_])
    # Raised from 14.0: the bar's own value label is a 41.5 px box, 1.72 units at this scale,
    # so "12" ran from 12.35 to 14.07 and the bracket started inside it.
    ax.plot([x0_, x0_, x1_, x1_], [14.8, 15.6, 15.6, 14.8], color=c_, lw=0.8, zorder=3)
    ax.text((x0_ + x1_) / 2, 15.9, nm_, ha="center", va="bottom",
            fontsize=st.sz(6.3), color=c_)
ax.text((XS[0] + XS[1]) / 2, 18.0, "mismatched", ha="center", va="bottom",
        fontsize=st.sz(6.3), color=ORANGE)
ax.text((XS[2] + XS[3]) / 2, 18.0, "matched", ha="center", va="bottom",
        fontsize=st.sz(6.3), color=BLUE)
ax.set_xticks(XS)
ax.set_xticklabels(["acquired\nat Cβ", "measured\nat Nζ",
                    "acquired\nat Cβ", "measured\nat Cβ"],
                   fontsize=st.sz(6.3), linespacing=1.3)
ax.set_xlim(-0.70, 3.90)
ax.set_ylim(0, 20); ax.set_yticks([0, 5, 10, 15, 20])
ax.set_ylabel("residues within 8 Å", fontsize=st.sz(7))
for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
HEAD_A = dict(a_union_residues_within_8A_at_CB=_A_UNION,
              a_measured_residues_within_8A_at_NZ=_A_NZ,
              a_measured_residues_within_8A_at_CB=_A_CB)

# ---------------------------------------------------------------- b  three populations
# TRANSPOSED, the move figures 4b, 4c, 6b and 7b make. Nine feature names down a y axis
# cost AXIS LENGTH; angled at 45 under an x axis their length costs ROW HEIGHT instead,
# and the effect-size axis is then free to be exactly the family's 483.6 px.
ax = fig.add_subplot(gt[0, 2]); st.head(fig, ax, "b")
ax._panel_head_flush = True
PRETTY = {"crowding_8A": "neighbors", "atoms_10A": "atom density",
          "hydrophobic_8A": "greasy", "aromatic_8A": "rings", "polar_8A": "polar",
          "anionic_8A": "acidic", "cationic_8A": "basic", "net_charge_8A": "net charge",
          "n_lys_8A": "other lysines"}
piv = MED.pivot(index="feature", columns="population", values="median")
sd = MED[MED.population == "acquired union"].set_index("feature")["sd"]
order = ["crowding_8A", "atoms_10A", "hydrophobic_8A", "aromatic_8A", "polar_8A",
         "anionic_8A", "cationic_8A", "net_charge_8A", "n_lys_8A"][::-1]
x = np.arange(len(order))
z_e = [(piv.loc[f, "engaged"] - piv.loc[f, "acquired union"]) / sd[f] for f in order]
z_n = [(piv.loc[f, "labelled, not engaged"] - piv.loc[f, "acquired union"]) / sd[f]
       for f in order]
# FURNITURE, not a series. The zero line was rust at lw 1.1, which is a colour two of this
# figure's populations use; it is now light grey dotted behind the data, on its own tick.
ax.axhline(0, color=LIGHTGREY, lw=0.8, ls=(0, (1.2, 2.0)), zorder=1)
# DODGED. Four of the nine features have identical medians in both populations, which is the
# panel's point, and undodged the charcoal square sat exactly on top of the teal circle and
# hid it: a reader saw one marker and could not tell the two agreed.
# NO CONNECTOR. The two markers are dodged either side of the tick and that is the whole of
# it: a rule joining them read as a slope between two things that are a pair of medians.
_DG = 0.15
ax.plot(x - _DG, z_e, "o", color=BLUE, ms=4.5, mew=0, ls="none", zorder=4)
ax.plot(x + _DG, z_n, "s", color=GREY, ms=4.5, mew=0, ls="none", zorder=4)
ax.set_xticks(x)
ax.set_xticklabels([PRETTY[f] for f in order], fontsize=st.sz(6.3), rotation=45,
                   ha="right", rotation_mode="anchor")
ax.set_xlim(-0.8, len(order) - 0.2)
# -0.8 to 0.4, both ends on a major tick, and the extremes are -0.730 and +0.330, so
# nothing is clipped off the frame.
ax.set_ylim(-0.8, 0.4); ax.set_yticks([-0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4])
ax.set_ylabel("difference from\nacquired lysines (SD)", fontsize=st.sz(7), linespacing=1.3)
for j_, (nm_, c_, mk_) in enumerate([("engaged", BLUE, "o"),
                                     ("labeled", GREY, "s")]):
    yy_ = 0.96 - j_ * 0.105
    ax.plot([0.13], [yy_], mk_, color=c_, ms=4.5, mew=0, transform=ax.transAxes,
            clip_on=False)
    ax.text(0.165, yy_, nm_, fontsize=st.sz(6.3), color=c_, va="center",
            transform=ax.transAxes)
for s_ in ("top", "right"): ax.spines[s_].set_visible(False)

HEAD = dict(HEAD_A)
fd.panel("fig8", "b", [dict(feature=_f, plain_name=PRETTY[_f],
                            median_acquired=float(piv.loc[_f, "acquired union"]),
                            median_engaged=float(piv.loc[_f, "engaged"]),
                            median_labelled_not_engaged=float(piv.loc[_f,
                                                        "labelled, not engaged"]),
                            sd_acquired=float(sd[_f]),
                            z_engaged=float((piv.loc[_f, "engaged"]
                                             - piv.loc[_f, "acquired union"]) / sd[_f]),
                            z_not_engaged=float((piv.loc[_f, "labelled, not engaged"]
                                                 - piv.loc[_f, "acquired union"]) / sd[_f]))
                       for _f in order[::-1]],
         note="medians of the nine environment features, three populations, and the engaged "
              "and unengaged offsets in SD of the acquired distribution")
HEAD["b_n_features"] = len(order)
HEAD["b_n_identical_medians"] = int(sum(
    1 for _f in order
    if piv.loc[_f, "acquired union"] == piv.loc[_f, "engaged"] == piv.loc[_f,
                                                             "labelled, not engaged"]))
for _f in order:
    HEAD["b_z_engaged_%s" % _f] = round(float((piv.loc[_f, "engaged"]
                                        - piv.loc[_f, "acquired union"]) / sd[_f]), 2)
HEAD["b_largest_abs_z"] = round(max(abs(_z) for _z in z_e), 2)

# ---------------------------------------------------------------- c  overlap, no threshold
ax = fig.add_subplot(gt[0, 4]); st.head(fig, ax, "c")
ax._panel_head_flush = True
CV = pd.read_csv(os.path.join(LR, "step28d_curve.csv"))
ax.plot(CV.radius, 100 * CV.engaged, color=BLUE, lw=1.3, zorder=3)
ax.plot(CV.radius, 100 * CV.not_engaged, color=GREY, lw=1.3, zorder=3)
ax.plot(CV.radius, 100 * CV.union, color=ORANGE, lw=1.3, zorder=3)
# Two lines: one line is 730 px against a 525 px frame and reached back under panel b.
ax.set_xlabel("distance from the engaged\ncenter (engaged SD)", fontsize=st.sz(7),
              linespacing=1.3)
# Two lines, because a 37-character rotated label is 614 px against a 483.6 px axis and
# would print past both ends of its own spine.
ax.set_ylabel("cumulative share\nof each population (%)", fontsize=st.sz(7), linespacing=1.3)
# 0 to 8, not 0 to 6.2: the curve runs to 7.28 and an axis that stopped at 6.2 was cutting
# real data off the frame. Both ends on a major tick.
ax.set_xlim(0, 8); ax.set_xticks([0, 2, 4, 6, 8])
ax.set_ylim(0, 100); ax.set_yticks([0, 25, 50, 75, 100])
for j_, (nm_, c_) in enumerate([("engaged", BLUE), ("labeled", GREY),
                                ("acquired", ORANGE)]):
    ax.text(0.04, 0.97 - j_ * 0.105, nm_, color=c_, fontsize=st.sz(6.3), ha="left",
            va="center", transform=ax.transAxes)
ax.text(0.98, 0.04, "P = %.3f\n(0.5 = interchangeable)" % P_CLOSER, ha="right", va="bottom",
        fontsize=st.sz(6.3), color=INK, linespacing=1.35, transform=ax.transAxes)
for s_ in ("top", "right"): ax.spines[s_].set_visible(False)

fd.panel("fig8", "c", CV, note="cumulative share of each population within a given "
                                "distance of the engaged centre, in engaged SD")
HEAD["c_p_closer"] = P_CLOSER
HEAD["c_n_union_positions"] = int(armA.n_union)
HEAD["c_n_engaged"] = int(armA.n_eng)
HEAD["c_n_labelled_not_engaged"] = int(armA.n_not)
HEAD["c_separation_armA"] = round(float(armA.d_union_eng), 3)
HEAD["c_control_armA"] = round(float(armA.d_control), 3)
HEAD["c_overlap_union_pct"] = round(100 * float(armA.overlap_union), 0)
HEAD["c_overlap_not_engaged_pct"] = round(100 * float(armA.overlap_not), 0)

# ---------------------------------------------------------------- d  crystal validation
ax = fig.add_subplot(gb[0, 0]); st.head(fig, ax, "d")
ax._panel_head_flush = True
ax._axes_dx_in = _DXD              # d's letter into a's column
ax._axes_dh_in = -(_Y1 - _YAX)     # back to the family's 483.6 px inside a taller cell
SITES = [("PTPN11 76", "SHP2\nE76K"), ("NRAS 61", "NRAS\nQ61K"), ("AKT1 17", "AKT1\nE17K"),
         ("PIK3CA 542", "PI3Kα\nE542K"), ("PIK3CA 545", "PI3Kα\nE545K"),
         ("TP53 280", "p53\nR280K")]
rng = np.random.default_rng(3)
floor = {}
for site, _ in SITES:
    floor[site] = XT[(XT.site == site) & (XT.kind == "WT vs WT")].rmsd.median()
# FURNITURE. The wild-type floor was a teal rule at lw 1.2; it is the reference every point
# is divided by, not a fourth series, and at y = 1 it now lands on its own tick.
ax.axhline(1.0, color=LIGHTGREY, lw=0.8, ls=(0, (1.2, 2.0)), zorder=1)
for i, (site, lab) in enumerate(SITES):
    f = floor[site]
    for kind, col, off in (("other sub vs WT", GREY, -0.19), ("K vs WT", ORANGE, 0.19)):
        v = XT[(XT.site == site) & (XT.kind == kind)].rmsd.values / f
        if not len(v):
            continue
        # 4.5 pt, the house marker, and alpha 0.55, which is what lets fifteen points in one
        # column stay countable.
        ax.plot(i + off + rng.uniform(-0.075, 0.075, len(v)), v, "o", color=col, ms=4.5,
                ls="none", zorder=4, mew=0, alpha=0.55)
        ax.plot([i + off - 0.15, i + off + 0.15], [np.median(v)] * 2, color=col, lw=1.6,
                zorder=5)
        # ABOVE THE COLUMN, not beside the median line: at the side it sat inside the
        # neighbouring column's points at three of the six sites.
        ax.text(i + off, v.max() + 0.16, "%.1f" % np.median(v), ha="center",
                va="bottom", fontsize=st.sz(6.3), color=col)
ax.set_xticks(range(len(SITES)))
ax.set_xticklabels([s[1] for s in SITES], fontsize=st.sz(6.3), linespacing=1.3)
# Only NRAS is annotated. The caption explains that one floor rests on a single pair; the
# other five n values are in the panel's deposited table, and six copies of the same
# annotation collided with each other at this width.
# TOP RIGHT. The left column, SHP2, reaches 5.13 of 6, so the only corner this panel leaves
# empty is the one above its two lowest columns.
for j_, (nm_, c_, mk_) in enumerate([("acquired lysine", ORANGE, "o"),
                                     ("other substitutions at the same position", GREY, "o")]):
    yy_ = 0.97 - j_ * 0.105
    # Left-aligned from a fixed column so the dot can sit at a known x; right-aligned text
    # starts wherever its own width puts it, which leaves nowhere to anchor a marker.
    ax.plot([0.36], [yy_], mk_, color=c_, ms=4.5, mew=0, transform=ax.transAxes,
            clip_on=False)
    ax.text(0.395, yy_, nm_, color=c_, fontsize=st.sz(6.3), ha="left", va="center",
            transform=ax.transAxes)
# Named where it is drawn rather than in the key, and in slate: the rule is light grey,
# which is a fill and not a text colour.
# AT THE LEFT END OF THE RULE, not the right. At the right it printed over the p53 column's
# low points and its median label, which the collision gate cannot see: that gate compares
# text to TEXT and has nothing to say about a word landing on a marker.
ax.text(-0.55, 1.12, "wild-type floor", color=INK, fontsize=st.sz(6.3), ha="left",
        va="bottom")
# 0 to 6 in steps of 1, so the floor at 1.0 is a tick rather than a line between two, and
# the largest point, SHP2 at 5.13, is inside the frame.
ax.set_ylim(0, 7); ax.set_yticks([0, 1, 2, 3, 4, 5, 6, 7])
ax.set_xlim(-0.60, 5.60)
# A 39-character line is 647 px against a 483.6 px axis, so the old label printed past both
# ends of its own spine and into the panel letter. Two short lines instead.
ax.set_ylabel("local shift, relative to\nthe wild-type floor",
              fontsize=st.sz(7), linespacing=1.3)
for s_ in ("top", "right"): ax.spines[s_].set_visible(False)

_drows = []
for _site, _lab in SITES:
    _f = floor[_site]
    _k = XT[(XT.site == _site) & (XT.kind == "K vs WT")].rmsd.values / _f
    _o = XT[(XT.site == _site) & (XT.kind == "other sub vs WT")].rmsd.values / _f
    _nfl = int(((XT.site == _site) & (XT.kind == "WT vs WT")).sum())
    _drows.append(dict(site=_site, label=_lab.replace(chr(10), " "),
                       wt_vs_wt_floor_angstrom=float(_f), n_floor_pairs=_nfl,
                       median_K_vs_WT_relative=float(np.median(_k)) if len(_k) else None,
                       n_K_pairs=len(_k),
                       median_other_sub_relative=float(np.median(_o)) if len(_o) else None,
                       n_other_sub_pairs=len(_o)))
    _key = _lab.replace(chr(10), "_")
    if len(_k):
        HEAD["d_%s_relative" % _key] = round(float(np.median(_k)), 2)
    if len(_o):
        HEAD["d_%s_other_sub_relative" % _key] = round(float(np.median(_o)), 2)
fd.panel("fig8", "d", _drows,
         note="local Ca displacement in a co-local shell, normalized to that protein's own "
              "wild-type-vs-wild-type floor")
HEAD["d_n_sites"] = len(SITES)

# ---------------------------------------------------------------- e  the drugged site
# THE PI'S OWN CHIMERAX RENDER, not a projection redrawn here. The panel is the image, drawn
# in its own pixel coordinates so every label and leader is placed in the same frame the
# render was measured in, and the cell was sized from the image's aspect above so it fills
# the panel rather than floating in a box of the wrong shape.
AT = pd.read_csv(os.path.join(LR, "step32_akt1_atoms.csv"))
BD = pd.read_csv(os.path.join(LR, "step32_akt1_bonds.csv"))
ax = fig.add_subplot(gb[0, 2]); st.head(fig, ax, "e")
ax._panel_head_flush = True
ax.imshow(_EPIX, extent=[0, _IW, _IH, 0], aspect="auto", interpolation="antialiased",
          zorder=1)
# OUTLINED. A render is a block of ink on a white page and, unlike a framed panel, has no
# edge of its own; without one its pale ribbon dissolves into the paper at the corners.
ax.add_patch(Rectangle((0, 0), _IW, _IH, fill=False, edgecolor=LIGHTGREY, lw=1.0, zorder=5))
# COMPUTED, NOT TYPED. The five distances were literals in the artwork; step 32 deposits the
# shell's real coordinates, so they are measured here from the zinc to each of its contacts
# and to the carbon that carries the covalent bond.
_XYZ = AT[["x", "y", "z"]].values
_ZN = _XYZ[AT.index[AT.part == "zinc"][0]]
_NZI = AT.index[(AT.part == "lys17") & (AT.name == "NZ")][0]
_DZ = {r["name"].split()[0]: float(np.linalg.norm(np.array([r.x, r.y, r.z]) - _ZN))
       for _, r in AT[AT.part == "znligand"].iterrows()}
_LGD = AT[AT.part == "ligand"]
_DCOV = float(np.min(np.linalg.norm(_LGD[["x", "y", "z"]].values - _XYZ[_NZI], axis=1)))
# A LABEL BAND UNDER THE IMAGE, not labels on it. The render carries ribbon across its whole
# frame, so text placed anywhere inside it printed over secondary structure and could not be
# read, and a white halo behind a word is the trick the PI deleted from figure 6f. The panel
# is the image plus two lines of white below it, and the cell height solved above carries the
# band as well.
_LINEPX = (41.5 / 300.0) * _IW / _EW          # one 7.974 pt line box, in image pixels
_BANDPX = 2.4 * _LINEPX
ax.set_xlim(0, _IW); ax.set_ylim(_IH + _BANDPX, 0); ax.set_axis_off()
# ARROWS, NOT LABELS LAID ON THE STRUCTURE. Each text sits on a patch of the render measured
# at 100% white rather than placed by eye. Positions are in the FULL image's pixels and the
# crop is subtracted at draw time, so they survive a change to the crop.
# FIGURE 6f'S ARROW, an open "->" head at lw 0.8, and ALL THREE THE SAME LENGTH. The arrow
# is a leader, not a measurement: it starts clear of its own word and runs _ALEN px towards
# the thing it names without reaching it, so the three read as one mark repeated rather than
# three different ones. Only the direction comes from the target, which is the atom's own
# position in the render, located by colour.
# 90 px, AND NO ARROW TOUCHES THE STRUCTURE. The three label boxes, their start points and
# this length were solved against the render's own ink: a 4 px corridor along each shaft
# and 8 px around each head, walked pixel by pixel, with every arrow required to stay on
# white for its whole run. 90 is what the zinc allows, which is the binding one, because
# it is the only target approached between two of its own ligands. All three then point
# exactly at their atom, 0 degrees off, rather than at an angle chosen to dodge ink.
_ALEN = 90.0                       # leader length, image px, the same for all three
for _nm, _c, (_tx, _ty), (_sx, _sy), (_px, _py) in (
        # Off the "3", with air between them. The corridor the leader has to run through is
        # narrow, clear only to x 310 at this height, so the gap could not come from moving
        # the arrow further right: the corridor is clear only to x 304 at this height and
        # narrows further as the tail rises, so the gap comes from the label's own position.
        # Its last glyph ends at 279 and the tail starts 25 px right of it and 20 px above
        # the line.
        ("compound 3", ORANGE, (155, 770), (304, 750), (455, 650)),
        ("K17 N\u03b6", INK, (1580, 560), (1479, 541), (1182, 486)),
        ("Zn", BLUE, (1250, 750), (1230, 709), (1171, 584))):
    _vx, _vy = _px - _sx, _py - _sy
    _n = (_vx * _vx + _vy * _vy) ** 0.5
    ax.text(_tx, _ty - _ETOP, _nm, fontsize=st.sz(6.3), color=_c, ha="center",
            va="center", zorder=6)
    ax.annotate("", xy=(_sx + _ALEN * _vx / _n, _sy + _ALEN * _vy / _n - _ETOP),
                xytext=(_sx, _sy - _ETOP), zorder=6,
                arrowprops=dict(arrowstyle="->", color=_c, lw=0.8, shrinkA=0, shrinkB=0))
ax.text(0, _IH + 0.75 * _LINEPX,
        "compound 3 covalently bonded to K17 N\u03b6, %.2f \u00c5" % _DCOV,
        ha="left", va="center", fontsize=st.sz(6.3), color=ORANGE)
ax.text(0, _IH + 1.85 * _LINEPX,
        "Zn contacts: K17 N\u03b6 %.2f \u00c5, C296 %.2f \u00c5, C310 %.2f \u00c5, compound 3 %.2f \u00c5"
        % (_DZ["LYS17"], _DZ["CYS296"], _DZ["CYS310"], _DZ["XOO"]),
        ha="left", va="center", fontsize=st.sz(6.3), color=BLUE)

fd.panel("fig8", "e", AT, note="PDB 8UVY local shell as plotted: atom identities and "
                                "2-D projection coordinates")
HEAD["e_n_atoms_8uvy_shell"] = int(len(AT))
HEAD["e_akt1_e17k_spectra"] = int(AK[AK.change == "E17K"].detected.sum())
if AKT_OK:
    HEAD["e_akt1_e17k_peptides_searched"] = len(AKT_PEPS)
    HEAD["e_akt1_e17k_variant_peptides"] = AKT_VAR
    HEAD["e_akt1_e17k_downstream_peptides"] = AKT_DOWN
    HEAD["e_akt1_e17k_plexes_searchable"] = len(AKT_PLEXES)
    HEAD["e_akt1_e17k_peptide_sequences"] = AKT_PEPS
fd.head("fig8", HEAD)
st.save(fig, "fig8")

# Trim to the drawing, 18 px all round, the standard Figure 1 set. Measured off the RENDERED
# PIXELS rather than get_tightbbox, because an axis-off panel still contributes its rectangle
# to the tight bbox and panel e has its axis off.
_PAD_IN = 18.0 / 300.0
from PIL import Image as _Im
_px = np.asarray(_Im.open("fig8.png").convert("RGB")).astype(int)
_ys, _xs = np.nonzero(_px.sum(2) < 735)
_H_in = fig.get_size_inches()[1]
_trim = Bbox.from_extents(_xs.min() / 300.0 - _PAD_IN,
                          _H_in - _ys.max() / 300.0 - _PAD_IN,
                          _xs.max() / 300.0 + _PAD_IN,
                          _H_in - _ys.min() / 300.0 + _PAD_IN)
for _ext in (".png", ".pdf"):
    fig.savefig("fig8" + _ext, bbox_inches=_trim)
print("trimmed to %.3f x %.3f in (%.0f x %.0f px at 300 dpi), 18 px pads"
      % (_trim.width, _trim.height, _trim.width * 300, _trim.height * 300))

if __import__("os").environ.get("FIG8_FINAL"):
    # THE FAMILY SHARES A RASTER, NOT A WIDTH. Figure 1 at 12.00 in fixes it at
    # 300 x 12.00 / 11.5400 = 311.96 dpi; every figure renders at that raster and is declared
    # 300 dpi, so all of them print 9.213 pt body type whatever their width.
    import PIL.Image as _PI, pymupdf as _mu
    _ref = _mu.open("fig1.pdf")[0].rect.width / 72.0
    _dpi = 300.0 * 12.0 / _ref
    fig.savefig("fig8_final.png", dpi=_dpi, bbox_inches=_trim)
    _im = _PI.open("fig8_final.png")
    _im.save("fig8_final.png", dpi=(300, 300))
    print("final: native %.4f in at the shared %.2f dpi raster -> %d x %d px declared 300 dpi "
          "(%.3f x %.3f in), body type %.3f pt against Figure 1's %.3f"
          % (_trim.width, _dpi, _im.size[0], _im.size[1], _im.size[0] / 300.0,
             _im.size[1] / 300.0, 8.860 * 12.0 / _ref, 8.860 * 12.0 / _ref))
import shutil
for e in ("pdf", "png"):
    shutil.copy("fig8.%s" % e, os.path.join(ROOT, "results", "fig8.%s" % e))
print("fig8: 5 panels | P(closer) %.3f | sites in d: %d | 8UVY atoms %d"
      % (P_CLOSER, len(SITES), len(AT)))
