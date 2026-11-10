#!/usr/bin/env python3
"""
Figure 7, The probe is not a neutral reporter: reagent charge steers labelling.

Answers the exception Figure 6 leaves open. Figure 6d shows no compound property predicts which
lysines are engaged, and figS7 f shows charge was never among the properties tested. Here charge
is tested directly, on the reagent rather than the compound.

WHAT CHANGED FROM THE PLANNED VERSION (2026-09-03)
The architecture planned a second half arguing "proximity to a concavity, not occupancy of a
pocket", with "41.6% of lysines within 5 A of a cavity, only 10.3% large enough to matter".
Three problems: those numbers appear in no step report and do not reproduce (40.4% ungated,
48.3% at pLDDT>=70; no volume threshold gives 10.3%); step 19's `in_pocket` IS
`dist_nearest_pocket <= 5`, so the two halves of the sentence are one variable; and step 36 had
already found fpocket's cavities here are surface depressions, not binding sites. Panel e
replaces it with enclosure (step 43), measured from the AlphaFold structure and independent of
both fpocket and solvent accessibility.

THE DECOMPOSITION (panels c, d; step 44)
A57 said A49's charge effect is "at least partly" reagent. Bootstrap CIs over proteins split it
three ways rather than hedging: cationic attraction is the REAGENT (the neutral probe's interval
brackets zero and excludes the anionic probe's estimate); acidic avoidance, predicted pKa and
salt bridges are the SITE (no probe difference); net local charge is MIXED, and for that one
"at least partly" is still the right words.

NOT USED: the inference that cationic is all reagent because STP-vs-unlabelled is significant
while NASA-vs-unlabelled is not. A difference between significant and non-significant is not a
significant difference. Panels c and d show intervals so the reader can see the actual basis.
"""
import os, re, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import abbasov_style as st
import figure_data as fd
from matplotlib.transforms import Bbox
# One body type size for this figure; every other size is a ratio to it.
BODY = 7.0
st.apply(BODY)
# 300, not the style module's 600, and it has to be set AFTER apply(), which resets the
# rcParams. THIS SCRIPT NEVER SET IT, so its PNG came out at twice the canvas and every
# pixel measurement taken off it was out by a factor of two. Set before apply() it is
# silently undone, which is how the first pass here still trimmed to 6899 px.
plt.rcParams["savefig.dpi"] = 300

def find(n):
    for d in (".", "..", "ligand_rebuild", "../ligand_rebuild"):
        p = os.path.join(d, n)
        if os.path.exists(p): return p
    raise FileNotFoundError(n)

R22 = open(find("step22_report.txt")).read()
T22 = pd.read_csv(find("step22_two_probe.csv"))
D44 = pd.read_csv(find("step44_charge_decomposition.csv"))
P43 = pd.read_csv(find("step43b_probe_contrast.csv"))
S43 = pd.read_csv(find("step43b_site_vs_unlabelled.csv"))

# step 22's own feature families, so the figure cannot drift from the analysis
CHARGE = ["anionic_8A","cationic_8A","net_charge_8A","propka_pKa","salt_bridge",
          "hbond_partners_NZ","n_lys_8A"]
POCKET = ["in_pocket","dist_nearest_pocket","pocket_volume","pocket_score","pocket_drugscore"]
SHAPE  = ["crowding_8A","atoms_10A","hydrophobic_8A","aromatic_8A","polar_8A","rsa_tien",
          "sasa_sidechain","nz_sasa","dist_centroid"]
FAM = {**{f:"charge" for f in CHARGE}, **{f:"pocket" for f in POCKET},
       **{f:"shape" for f in SHAPE}}
FCOL = {"charge": st.ORANGE, "pocket": st.BLUE, "shape": st.GREY, "context": st.LIGHTGREY}

NICE = lambda f: (f.replace("_8A"," 8Å").replace("_10A"," 10Å").replace("_"," ")
                   .replace("rsa tien","exposure").replace("propka pKa","predicted pKa")
                   .replace("nz sasa","NZ exposure").replace("sasa sidechain","side-chain exposure")
                   .replace("dist centroid","distance to centroid").replace("n lys","nearby lysines")
                   .replace("dist nearest pocket","distance to pocket")
                   .replace("hbond partners NZ","H-bonds at NZ")
                   .replace("plddt","pLDDT"))

NSTP, NNASA, NBOTH = [int(x) for x in
    re.search(r"STP-only\s+(\d+)\s+\|\s+NASA-only\s+(\d+)\s+\|\s+both\s+(\d+)", R22).groups()]
CON = T22[T22.group == "STP_ONLY vs NASA_ONLY"].copy()
CON["fam"] = CON.feature.map(lambda f: FAM.get(f, "context"))
NSIG = int((CON.q < 0.05).sum())
NPROT = int(CON.n_proteins.iloc[0])
assert NSIG == 2, NSIG

HEAD = {"a_n_stp_only": NSTP, "a_n_nasa_only": NNASA, "a_n_both": NBOTH,
        "a_n_proteins_with_both": NPROT,
        "b_n_significant": NSIG, "b_n_features": int(len(CON)),
        "b_n_pocket_significant": int((CON[CON.fam == "pocket"].q < 0.05).sum()),
        "b_n_pocket": int((CON.fam == "pocket").sum()),
        "b_n_shape_significant": int((CON[CON.fam == "shape"].q < 0.05).sum()),
        "b_n_shape": int((CON.fam == "shape").sum())}

fd.panel("fig7", "b", CON[["feature", "fam", "d", "q", "n_proteins"]],
         note="all 22 site features, anionic-only vs neutral-only sites, within protein")
for _, _r in CON[CON.q < 0.05].iterrows():
    HEAD["b_%s_d" % _r.feature] = round(float(_r.d), 3)
    HEAD["b_%s_q" % _r.feature] = "%.1e" % float(_r.q)
# the next-largest difference is quoted in 2.7 as the one that does not clear the threshold
_ns = CON[CON.q >= 0.05].sort_values("d", ascending=False)
if len(_ns):
    HEAD["b_next_largest_feature"] = _ns.feature.iloc[0]
    HEAD["b_next_largest_d"] = round(float(_ns.d.iloc[0]), 3)

fd.panel("fig7", "cd", D44, note="bootstrap decomposition: each probe against unlabeled "
                                 "lysines, and the probe-minus-probe difference")

fd.panel("fig7", "e", pd.concat([S43.assign(contrast="labelled vs unlabelled"),
                                 P43.assign(contrast="anionic vs neutral probe")],
                                ignore_index=True),
         note="enclosure and pocket occupancy, site contrast against probe contrast")
for _lab, _df in (("site", S43), ("probe", P43)):
    _r = _df[_df.feature == "enclosure"]
    if len(_r):
        HEAD["e_enclosure_%s_d" % _lab] = round(float(_r.d.iloc[0]), 3)
        HEAD["e_enclosure_%s_q" % _lab] = ("%.1e" % float(_r.q.iloc[0])
                                           if float(_r.q.iloc[0]) < 0.01
                                           else round(float(_r.q.iloc[0]), 2))

# Panel b lists 22 features down its y-axis. At the house body those labels print at 6.0 pt
# and need about 7.5 pt of row pitch, so 22 rows need 2.35 in of axes height. The top row had
# 1.56 in, which is why adjacent feature names were printing through each other: "nearby
# lysines 8A" into "atoms 10A", "net charge 8A" into "aromatic 8A". Type cannot come down to
# fix it, because 6.0 pt is the floor the whole figure set was just raised to meet, so the
# panel has to grow. The figure is therefore taller and the TOP row takes all of the increase:
# the bottom row keeps its old 1.60 in and the gaps and margins keep their old inches too.
# MARGINS IN INCHES, CANVAS SOLVED FROM THE PANEL BLOCK, which is how figures 4, 5 and 6 are
# built. Fractions of a figsize cannot hold a panel width while the figure height changes, and
# they are why this figure's two rows do not line up at either edge.
_ML, _MR = 1.60, 1.20          # left and right margins, inches
_MB, _MT = 1.60, 0.80          # bottom and top margins, inches
_YAX = 1.55020                 # the family y axis, 483.6 px at the 311.96 raster
# Row 0's CELL is 2.20 in, but its drawn y axis is still the family 483.6: panel b is
# shortened inside the cell with _axes_dh_in. The cell is taller than the axis for panel a,
# which stacks two structures, and a structure printed under about 220 px has atom labels
# below the 6 pt floor this figure set is built to.
_Y0 = 2.20
# THE BLOCK WIDTH IS HELD AND THE SPACING IS SOLVED FROM IT, the reverse of Figure 6, because
# every panel in row 1 hangs a column of feature names outside its frame. At Figure 6's 0.42
# those names ate the whole gap: the row measured 75 px and 16 px between panels, and e's
# labels reached 24 px back into d before "accessibility removed" was shortened.
_HS = 0.45                     # row spacing, as a fraction of a panel
# _PW IS NO LONGER THE SHEET WIDTH. Figures 1 and 4 to 6 trim to 3,461 px because their ink
# spans the block; this figure's does not. Both rows are now built from explicit column widths
# and stop where their content stops, so _PW is only the canvas the panels are placed on and
# the trim is set by the ink. The sheet comes out near 3,050 px, which the PI chose knowing
# figure 7 will therefore print at a different scale from the rest of the set.
_PW = 10.5455
# ROW 1, EXPLICIT WIDTHS. Each panel hangs its feature names left of its own frame and its
# last x tick label past the right of it, so an even gridspec cannot hold an even INK gap:
# the old one needed a -29 px nudge on d to even 229 and 228. These are the measured overhangs,
# and the cell gaps are built from them, so the ink gap below is the number that lands.
_YL_D, _YL_E = 238.0, 306.0        # feature names hanging left of d's and e's frames
_OH_C, _OH_D = 33.0, 24.0          # last x tick label hanging right of c's and d's frames
# C AND D ARE ONE WIDTH; E IS SET BY ITS TICK PITCH. c and d both run 5 x ticks across 564 px,
# so their ticks stand 141 px apart. e has 4 ticks, 3 intervals, and at 564 px they stood 188
# apart. 423 px puts e's ticks on the same 141 px pitch as c's and d's, which is what makes the
# three x axes read as one scale rather than three.
_CW = _DW = 564.0 / 300.0          # narrowed from 704 px
_EW = 423.0 / 300.0                # 3 intervals at c and d's 141 px
_GINK = 100.0                      # ONE ink gap, so the three panels are evenly spaced
_GCD, _GDE = (_GINK + _YL_D + _OH_C) / 300.0, (_GINK + _YL_E + _OH_D) / 300.0
_S1_IN = _PW - _CW - _GCD - _DW - _GDE - _EW
# ROW 0 DOES NOT HAVE TO REACH THE RIGHT EDGE. Every earlier version of this row forced it
# to, and that one constraint is what made a's bar and b's label pitch fight each other: the
# only place the width b gave up could go was panel a. Released, each panel is simply as wide
# as its own content needs and the leftover sits at the right of the row as _S_IN.
# A IS ITS CONTENT. The two structures are 662.0 and 639.1 px wide at _H = 224 px tall, and
# _GAPPX is the gap between them, so a's cell is their sum and panel a's 0..100 unit span is
# exactly the pair. The stacked bar spans the same 0..100, which is why it ends where the
# probes do.
# A IS SET BY ROW 1 NOW. Panel a's ink width is whatever makes row 0 span the same length as
# row 1, and _H follows from it: the two structures and the gap between them are scaled
# together by _ASCALE, so the pair still fills 0..100 and the bar still ends where it ends.
_A_W_PX = 1135.0               # row 1's ink (2,614 px) less b's ink (1,379) less the gap (100)
_ASCALE = _A_W_PX / (662.0 + 639.1 + 104.0)
_GAPPX = 104.0 * _ASCALE       # gap between the two structures, px at 300 dpi
_A_IN = _A_W_PX / 300.0
# B IS ITS LABEL PITCH. 22 slots across a 1,232 px axes put the angled names 56.7 px apart,
# which is figure 6b's pitch exactly and inside 4b's 61.7 and 4c's 56.8.
_B_IN = 1232.0 / 300.0
# The gap is set on the INK, not the cells: a's axes is shifted 0.788 in left of its cell to
# put its letter on panel c's, and b hangs its y axis label 0.263 in left of its own frame.
_AFLUSH = 0.812405             # where a's x = 0 lands, which is panel c's letter column
_G_IN = _AFLUSH + 100.0 / 300.0 + 0.490 - _ML      # 100 px, row 1's own gap
_S_IN = _PW - _A_IN - _G_IN - _B_IN        # white at the right of row 0, and that is allowed
_AW = _A_IN
_PXX, _PXY = _AW * 300.0 / 100.0, _Y0 * 300.0 / 10.0
_HG = _HS * _YAX
_PH = _Y0 + _HG + _YAX
_FW, _FH = _PW + _ML + _MR, _PH + _MB + _MT
fig = plt.figure(figsize=(_FW, _FH))
# FOUR COLUMNS, NOT TWO, with the gap and the right-hand slack as columns of their own and
# wspace 0. A two-column gridspec sets its gap as a fraction of the mean panel, so a and b
# could not be sized independently of the space between them; here each width is in inches.
gt = fig.add_gridspec(1, 4, left=_ML / _FW, right=1 - _MR / _FW,
                      top=1 - _MT / _FH, bottom=1 - (_MT + _Y0) / _FH,
                      wspace=0, width_ratios=[_A_IN, _G_IN, _B_IN, _S_IN])
gb = fig.add_gridspec(1, 6, left=_ML / _FW, right=1 - _MR / _FW,
                      top=(_MB + _YAX) / _FH, bottom=_MB / _FH, wspace=0,
                      width_ratios=[_CW, _GCD, _DW, _GDE, _EW, _S1_IN])

# ================================================================ a  the design
# THE PI'S OWN DRAWINGS, cropped from the sheet he supplied, not redrawn here. They already
# carry the point in colour: the desthiobiotin core is teal in both and only the leaving group
# takes the probe's colour. SIDE BY SIDE at a COMMON SCALE on a COMMON BASELINE, so the two
# teal halves read as one shared molecule and each dotted rule marks where its own probe stops
# being that molecule. Probe 1 is dTB-STP, probe 2 is the NASA reagent, the numbering the PI set.
_PDIR = os.path.join("results", "probe_structures")
_P1 = plt.imread(os.path.join(_PDIR, "probe1_stp.png"))     # 854 x 289, teal ends at col 494
_P2 = plt.imread(os.path.join(_PDIR, "probe2_nasa.png"))    # 836 x 293, teal ends at col 490

ax = fig.add_subplot(gt[0, 0]); st.head(fig, ax, "a")
ax._panel_head_flush = True
# a's ink starts at x = 0, so the shift is just the distance from the gridspec's left margin
# to panel c's letter column. Nothing here depends on the cell width.
ax._axes_dx_in = _AFLUSH - _ML
ax.set_xlim(0, 100); ax.set_ylim(0, 10); ax.set_axis_off()
# Both structures are drawn _H px tall and each one's width follows from its own aspect, so the
# two are at the same scale. _H FOLLOWS THE PANEL WIDTH: at 224 px the atom labels measure
# 25 px of cap height, about 8.4 pt, which is the house body size; at this scale they are
# 8.4 * _ASCALE pt. Printed under about 160 px they fall below 6 pt.
_H = 224.0 * _ASCALE
# The gap between the structures is held in PIXELS and converted, so it cannot drift when the
# cell moves. _X0 to 100 is the pair and the bar alike.
_X0 = 0.0
_GAPX = _GAPPX / _PXX
# THE PANEL HANGS FROM ITS TOP LINE, not from the structures' baseline. The panel letter is
# placed on a panel's topmost ink, so if the structures shrink and the block above them follows
# them down, a's letter drops below b's and row 0 stops being flush. "leaving group" is
# therefore the anchor, the rules hang 0.45 below it, the structures 0.26 below those, and
# everything under the structures moves up by _DY, which is the height they gave back.
_LBLY = 9.80                       # the top line, fixed
_RTOP = _LBLY - 0.45               # where the dotted rules end
_PY = _RTOP - 0.26 - _H / _PXY     # the common baseline; the bar goes UNDER both structures
_DY = _PY - 5.70
# Teal ends at 494/854 = 0.578 of probe 1's width and 490/836 = 0.586 of probe 2's. Each rule
# sits _NUDGE left of ITS OWN boundary rather than at a typed fraction, so the two move
# together and each still lands on the acyl bond its leaving group departs from.
_NUDGE = 0.010
# Centred, because the bar is now the full-width element and the pair is narrower than it.
_PAIRW = _GAPX + sum(_H * n[0] / n[1] / _PXX for n in ((854, 289), (836, 293)))
_x, _divs = _X0 + (100.0 - _X0 - _PAIRW) / 2.0, []
for _img, _nat, _frac, _nm, _c in ((_P1, (854, 289), 494.0 / 854.0, "probe 1, anionic", st.ORANGE),
                                   (_P2, (836, 293), 490.0 / 836.0, "probe 2, neutral", st.GREY)):
    _w = _H * _nat[0] / _nat[1] / _PXX
    ax.imshow(_img, extent=[_x, _x + _w, _PY, _PY + _H / _PXY], aspect="auto",
              interpolation="antialiased", zorder=2)
    ax.text(_x + 4.0, _PY - 0.25, _nm, fontsize=st.sz(6.3), color=_c, fontweight="bold",
            va="top")
    _divs.append(_x + (_frac - _NUDGE) * _w)
    _x += _w + _GAPX
_RIGHT = _x - _GAPX                # right edge of probe 2
# IN FRONT of the drawings, at the PI's request, one rule per structure, and each one named.
for _d in _divs:
    ax.plot([_d, _d], [_PY - 0.20, _RTOP], color=st.LIGHTGREY, lw=1.0,
            ls=(0, (1.2, 2.0)), zorder=4)
    ax.text(_d + 0.8, _LBLY, "leaving group", fontsize=st.sz(6.3), color=st.GREY, va="top")

BAR = [("probe 1 only", NSTP, st.ORANGE), ("both", NBOTH, st.LIGHTGREY),
       ("probe 2 only", NNASA, st.GREY)]
# EDGE TO EDGE. The bar is what makes panel a as wide as its cell, so that widening the cell
# to squeeze panel b's names adds bar rather than white.
tot = sum(b[1] for b in BAR); _BW = 100.0 - _X0; x0 = _X0
for nm, n, c in BAR:
    w = _BW * n / tot
    ax.add_patch(Rectangle((x0, 3.60 + _DY), w, 1.10, facecolor=c, lw=0))
    # LIGHT GREY IS A FILL, NOT A TEXT COLOUR. "12,298" and "both" in #C5C3C6 on white are the
    # faintest marks on the sheet; the bar keeps the fill, the words take slate.
    _tc = st.GREY if c is st.LIGHTGREY else c
    ax.text(x0 + w/2, 4.15 + _DY, "%s" % f"{n:,}", fontsize=st.sz(6.3), color=st.PAPER,
            fontweight="bold", ha="center", va="center")
    ax.text(x0 + w/2, 3.40 + _DY, nm, fontsize=st.sz(6.3), color=_tc, ha="center", va="top")
    x0 += w
# Centred on the bar, with the arrow, rather than left-aligned under its left-hand end.
_BC = _X0 + _BW / 2.0
# mutation_scale is what sizes the HEAD. Left unset, annotate takes it from the annotation's
# font size, which is the body size and gives a head too small to read at 300 dpi.
ax.annotate("", xy=(_BC, 1.67 + _DY), xytext=(_BC, 2.70 + _DY),
            arrowprops=dict(arrowstyle="-|>", color=st.GREY, lw=1.05, mutation_scale=10.0,
                            joinstyle="miter", capstyle="butt"))
ax.text(_BC, 1.54 + _DY, "compared within protein", fontsize=st.sz(6.3), color=st.GREY,
        ha="center", va="top")
ax.text(_BC, 0.74 + _DY, "%d proteins carry both" % NPROT, fontsize=st.sz(6.3), color=st.BLUE,
        fontweight="bold", ha="center", va="top")

# ================================================================ b  all 22 features
# TRANSPOSED, the move Figure 4 panels b and c make and for the same reason. Twenty-two feature
# names down a y axis need 930 px of spine before any clearance at all, which is what held this
# row at twice the family height. Angled under an x axis their length costs ROW HEIGHT instead,
# and the effect-size axis is then free to be exactly Figure 1's 483.6 px. No value changes.
ax = fig.add_subplot(gt[0, 2]); st.head(fig, ax, "b")
ax._panel_head_flush = True
# The cell is 2.20 in and every drawn y axis in the family is 1.55020. _axes_dh_in holds
# the top and brings the bottom up, which also puts the angled names inside the cell
# rather than hanging them into the row gap.
ax._axes_dh_in = -(2.20 - _YAX)
C = CON.sort_values("d")
x = np.arange(len(C))
for i, (_, r) in enumerate(C.iterrows()):
    c = FCOL[r.fam] if r.q < 0.05 else st.LIGHTGREY
    ax.plot([i], [r.d], "o", color=c, ms=4.5, mew=0, zorder=3)
ax.axhline(0, color=st.LIGHTGREY, lw=0.8, ls=(0, (1.2, 2.0)), zorder=1)
ax.set_xticks(x)
# ha="right" with rotation_mode="anchor" puts each label's last character under its own tick,
# which is what makes an angled label read as belonging to one position rather than drifting
# between two.
ax.set_xticklabels([NICE(f) for f in C.feature], fontsize=st.sz(6.3), rotation=45,
                   ha="right", rotation_mode="anchor")
for t, f in zip(ax.get_xticklabels(), C.feature):
    if float(C.loc[C.feature == f, "q"].iloc[0]) < 0.05:
        t.set_color(st.ORANGE); t.set_fontweight("bold")
ax.set_ylabel("higher in probe 1's sites (d)", fontsize=st.sz(7))
ax.set_xlim(-0.8, len(C) - 0.2)
# -0.20, NOT -0.10. pLDDT 8Angstrom is d = -0.1098, so at -0.10 its dot was clipped off the
# axes and the leftmost tick label stood over nothing. The PI caught it.
ax.set_ylim(-0.20, 0.30); ax.set_yticks([-0.2, -0.1, 0.0, 0.1, 0.2, 0.3])
# Top left. The features sort ascending, so the left end is the negative one and everything
# above it is empty in this panel.
ax.text(0.015, 0.97, "%d of %d differ\nboth are charge" % (NSIG, len(C)), transform=ax.transAxes,
        fontsize=st.sz(6.3), color=st.ORANGE, fontweight="bold", va="top", linespacing=1.45)
ax.text(0.015, 0.68, "0 of 5 pocket\n0 of 9 shape", transform=ax.transAxes, va="top",
        fontsize=st.sz(6.3), color=st.GREY, linespacing=1.45)
for s_ in ("top", "right"): ax.spines[s_].set_visible(False)

# ================================================================ c  what each probe prefers
ORDER = ["cationic_8A", "net_charge_8A", "anionic_8A", "propka_pKa", "salt_bridge"]
ax = fig.add_subplot(gb[0, 0]); st.head(fig, ax, "c")
ax._panel_head_flush = True
CMPS = [("STP-only vs unlabelled", st.ORANGE, "o"), ("NASA-only vs unlabelled", st.GREY, "s")]
for i, f in enumerate(ORDER):
    for j, (cmp_, c, mk) in enumerate(CMPS):
        r = D44[(D44.feature == f) & (D44.comparison == cmp_)]
        if not len(r): continue
        r = r.iloc[0]; yy = i + (0.17 if j == 0 else -0.17)
        ax.errorbar([r.d], [yy], xerr=[[r.d - r.lo], [r.hi - r.d]], fmt="none", ecolor=c,
                    elinewidth=1.1, capsize=2.8, capthick=1.1, zorder=2)
        ax.plot([r.d], [yy], mk, color=c, ms=4.5, mew=0, zorder=3)
ax.axvline(0, color=st.LIGHTGREY, lw=0.8, ls=(0, (1.2, 2.0)), zorder=1)
ax.set_yticks(range(len(ORDER))); ax.set_yticklabels([NICE(f) for f in ORDER], fontsize=st.sz(6.3))
ax.set_xlabel("preference vs unlabeled lysines (d)", fontsize=st.sz(7))
# net charge's anionic CI reaches +0.45; a tighter limit clipped its marker off the panel.
_lo, _hi = -0.50, 0.50
ax.set_xticks([-0.5, -0.25, 0.0, 0.25, 0.5])
# -0.7, NOT -1.15. The old lower limit existed only to hold the key below the bottom row, and
# it also meant c's five rows sat 34 px off d's five rows, which are the SAME five features.
ax.set_ylim(-0.7, len(ORDER) - 0.4); ax.set_xlim(_lo, _hi)
# TOP RIGHT, in axes fractions. The two top rows carry their data on the left half, so the
# corner is empty; the panel's only right-hand ink is net charge, four rows down.
for j, (nm, c, mk) in enumerate([("probe 1", st.ORANGE, "o"),
                                 ("probe 2", st.GREY, "s")]):
    yy = 0.96 - j * 0.105
    ax.plot([0.80], [yy], mk, color=c, ms=4.5, mew=0, transform=ax.transAxes, clip_on=False)
    ax.text(0.835, yy, nm, fontsize=st.sz(6.3), color=c, va="center", transform=ax.transAxes)
for s in ("top", "right"): ax.spines[s].set_visible(False)

# ================================================================ d  the difference
ax = fig.add_subplot(gb[0, 2]); st.head(fig, ax, "d")
ax._panel_head_flush = True
VERD = {"cationic_8A": ("reagent", st.ORANGE), "net_charge_8A": ("mixed", st.GREEN),
        "anionic_8A": ("site", st.BLUE), "propka_pKa": ("site", st.BLUE),
        "salt_bridge": ("site", st.BLUE)}
for i, f in enumerate(ORDER):
    r = D44[(D44.feature == f) & (D44.comparison == "STP-only vs NASA-only")]
    if not len(r): continue
    r = r.iloc[0]; lab, c = VERD[f]
    ax.errorbar([r.d], [i], xerr=[[r.d - r.lo], [r.hi - r.d]], fmt="none", ecolor=c,
                elinewidth=1.3, capsize=3.0, capthick=1.3, zorder=2)
    ax.plot([r.d], [i], "o", color=c, ms=4.5, mew=0, zorder=3)
    ax.text(r.hi + 0.018, i, lab, fontsize=st.sz(6.3), color=c, va="center")
ax.axvline(0, color=st.LIGHTGREY, lw=0.8, ls=(0, (1.2, 2.0)), zorder=1)
ax.set_yticks(range(len(ORDER))); ax.set_yticklabels([NICE(f) for f in ORDER], fontsize=st.sz(6.3))
ax.set_xlabel("probe 1 minus probe 2 (d)", fontsize=st.sz(7))
ax.set_ylim(-0.7, len(ORDER) - 0.4); ax.set_xlim(-0.20, 0.60)
ax.set_xticks([-0.2, 0.0, 0.2, 0.4, 0.6])
for s in ("top", "right"): ax.spines[s].set_visible(False)

# ================================================================ e  concavity is the site's
ax = fig.add_subplot(gb[0, 4]); st.head(fig, ax, "e")
ax._panel_head_flush = True
ROWS = [("labeled vs\nunlabeled", S43, st.BLUE), ("anionic vs\nneutral probe", P43, st.ORANGE)]
FE = ["enclosure", "enclosure_adj", "in_pocket"]
# "accessibility removed" is 21 characters, 349 px, and it is what pushed e's ink 24 px back
# into panel d. "exposure" is the word this figure already uses for the same quantity, in
# "NZ exposure" and "side-chain exp.", and it is five characters shorter.
FL = {"enclosure": "enclosure", "enclosure_adj": "enclosure,\nexposure removed",
      "in_pocket": "in a pocket"}
for i, f in enumerate(FE):
    for j, (nm, df, c) in enumerate(ROWS):
        r = df[df.feature == f]
        if not len(r): continue
        r = r.iloc[0]; yy = i + (0.18 if j == 0 else -0.18)
        sig = r.q < 0.05
        ax.plot([r.d], [yy], "o" if j == 0 else "s", color=c if sig else st.LIGHTGREY,
                ms=4.5, mew=0, zorder=3)
# BOUNDED, not a full-height axvline. e borrows c and d's limits so its rows share their
# pitch, which leaves two rows of empty axis at the top; the spine and the zero rule stop
# 0.6 above the top row, the same overhang c and d have above theirs, so the empty band
# carries the key and nothing else.
_ETOP = len(FE) - 0.4
ax.plot([0, 0], [-0.7, _ETOP], color=st.LIGHTGREY, lw=0.8, ls=(0, (1.2, 2.0)), zorder=1)
ax.spines["left"].set_bounds(-0.7, _ETOP)
ax.set_yticks(range(len(FE))); ax.set_yticklabels([FL[f] for f in FE], fontsize=st.sz(6.3),
                                                  linespacing=1.35)
ax.set_xlabel("effect size (d)", fontsize=st.sz(7))
ax.set_xlim(-0.10, 0.20); ax.set_xticks([-0.1, 0.0, 0.1, 0.2])
# len(ORDER), NOT len(FE). The limits are c's and d's, so e's three rows fall exactly on their
# bottom three and the row pitch is one number across the whole of row 1. e simply has two
# rows' worth of empty axis above its data, which is where the key goes.
ax.set_ylim(-0.7, len(ORDER) - 0.4)
# THE KEY CANNOT NAME A COLOUR THE PANEL DOES NOT DRAW. Every probe-vs-probe point here is
# non-significant (q = 0.71, 0.18, 0.052), so all three of its squares are drawn light grey,
# and a rust square in the key sent the eye looking for rust ink that does not exist. Each
# entry now takes its series' colour only if that series has a significant point in the panel.
# The swatch goes light grey like the ink; the WORDS go slate, because light grey is a fill.
# THE KEY CARRIES THE VERDICT NOW. "grey = not significant" sat at the bottom right and could
# not be read: all three labeled-vs-unlabeled points are significant and all three
# probe-vs-probe points are not, so every teal mark is a circle and every grey mark is a
# square, and the note was indistinguishable from "grey = probe vs probe". The second entry
# says what grey means in this panel instead.
# ONE TEXT CALL PER LINE, each va="center" on its own baseline, and the marker on the FIRST
# line's baseline. A two-line string with va="top" does not put its first line's centre where
# the arithmetic says it does: matplotlib measures the whole block, so the dot sat off its
# own word. Drawing the lines separately removes the guess.
_KLH = 41.5 / 483.6                # one 7.974 pt line box, as a fraction of the axes height
_ky = 0.96
for j, (nm, df, c, mk) in enumerate([("labeled vs unlabeled", S43, st.BLUE, "o"),
                                     ("probe vs probe\n(none significant)", P43, st.ORANGE, "s")]):
    _any = bool((df[df.feature.isin(FE)].q < 0.05).any())
    # TOP LEFT, not top right as in panel c: "labeled vs unlabeled" is 332 px and the frame
    # is now 564, so right-aligning it would hang the longer entry over the axis.
    ax.plot([0.03], [_ky], mk, color=c if _any else st.LIGHTGREY, ms=4.5, mew=0,
            transform=ax.transAxes, clip_on=False)
    for _li, _line in enumerate(nm.split("\n")):
        ax.text(0.065, _ky - _li * _KLH, _line, fontsize=st.sz(6.3),
                color=c if _any else st.GREY, va="center", transform=ax.transAxes)
    _ky -= (len(nm.split("\n")) - 1) * _KLH + 0.105
for s in ("top", "right"): ax.spines[s].set_visible(False)

st.save(fig, "fig7")
# Trim to the drawing, 18 px all round, the standard Figure 1 set. Measured off the RENDERED
# PIXELS rather than get_tightbbox, because an axis-off panel still contributes its rectangle
# to the tight bbox and panel a has its axis off.
_PAD_IN = 18.0 / 300.0
from PIL import Image as _Im
_px = np.asarray(_Im.open("fig7.png").convert("RGB")).astype(int)
_ys, _xs = np.nonzero(_px.sum(2) < 735)
_H_in = fig.get_size_inches()[1]
_trim = Bbox.from_extents(_xs.min() / 300.0 - _PAD_IN,
                          _H_in - _ys.max() / 300.0 - _PAD_IN,
                          _xs.max() / 300.0 + _PAD_IN,
                          _H_in - _ys.min() / 300.0 + _PAD_IN)
for _ext in (".png", ".pdf"):
    fig.savefig("fig7" + _ext, bbox_inches=_trim)
print("trimmed to %.3f x %.3f in (%.0f x %.0f px at 300 dpi), 18 px pads"
      % (_trim.width, _trim.height, _trim.width * 300, _trim.height * 300))

if __import__("os").environ.get("FIG7_FINAL"):
    # THE FAMILY SHARES A RASTER, NOT A WIDTH. Figure 1 at 12.00 in fixes it at
    # 300 x 12.00 / 11.5400 = 311.96 dpi; every figure renders at that raster and is declared
    # 300 dpi, so all of them print 9.213 pt body type whatever their width.
    import PIL.Image as _PI, pymupdf as _mu
    _ref = _mu.open("fig1.pdf")[0].rect.width / 72.0
    _dpi = 300.0 * 12.0 / _ref
    fig.savefig("fig7_final.png", dpi=_dpi, bbox_inches=_trim)
    _im = _PI.open("fig7_final.png")
    _im.save("fig7_final.png", dpi=(300, 300))
    print("final: native %.4f in at the shared %.2f dpi raster -> %d x %d px declared 300 dpi "
          "(%.3f x %.3f in), body type %.3f pt against Figure 1's %.3f"
          % (_trim.width, _dpi, _im.size[0], _im.size[1], _im.size[0] / 300.0,
             _im.size[1] / 300.0, 8.860 * 12.0 / _ref, 8.860 * 12.0 / _ref))
fd.head("fig7", HEAD)
print("a: STP-only %d | both %d | NASA-only %d | %d proteins" % (NSTP, NBOTH, NNASA, NPROT))
print("b: %d of %d features differ; families %s"
      % (NSIG, len(CON), CON[CON.q < 0.05].fam.tolist()))
print("c/d: %s" % {f: (VERD[f][0]) for f in ORDER})
print("e: enclosure labelled-vs-unlab d %.3f (q %.3g) | probe-vs-probe d %.3f (q %.2f)"
      % (S43[S43.feature=="enclosure"].d.iloc[0], S43[S43.feature=="enclosure"].q.iloc[0],
         P43[P43.feature=="enclosure"].d.iloc[0], P43[P43.feature=="enclosure"].q.iloc[0]))
