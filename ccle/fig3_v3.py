#!/usr/bin/env python3
"""
Figure 3, how the false-discovery rate was controlled, and what it cost.

THE PROBLEM THIS FIGURE ANSWERS
Standard target-decoy FDR assumes the candidate population has a reasonable prior of being
present. A variant database does not: most of the 8,470 acquired-lysine peptides are absent from
any given cell line, so a 1% FDR computed against the whole search is meaningless for the
variant subset. Two things were done about it, a two-stage search that establishes the
reference population first, and a paired entrapment set that measures the false-discovery
proportion directly rather than estimating it.

a  the design: two-stage search, and what each stage is for
b  the score margin dose-response, how FDP falls as the margin requirement rises. This is what
   justifies the >= 10 threshold; without it the threshold is arbitrary
c  per-plex consistency: the two populations side by side, one point per plex, so the
   question of whether the entrapment rate is stable across 42 experiments or driven by a
   few outliers is answered by looking at the two columns
d  spectrum quality, scored by SpecQuality independently of the search engine, the database
   and the FDR, so it is not circular with b or c
e  mass accuracy, observed minus calculated, zoomed to the monoisotopic population. The +1 Da
   cluster is what withdrew YBX1 E4K
f  the pairing used directly: where a variant peptide's score distribution overlaps that of
   its own entrapment partner is where an identification cannot be distinguished from a
   peptide known to be absent

ENTRAPMENT DESIGN. Each variant peptide is paired with an entrapment peptide of matched length
and composition that cannot be present. An entrapment identification is therefore a direct
observation of a false positive, not an estimate.
"""
import os, sys, re, csv, glob
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Polygon
from matplotlib.transforms import Bbox
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import abbasov_style as st
import figure_data as fd

# --- portability -------------------------------------------------------------------------
# The data search paths below were written as ~/Desktop/DepMap/... , which only resolves on
# the machine the analysis was built on. They are now anchored to this script's own location,
# so the build works wherever the folder lives (a collaborator's checkout, a container, a
# different user account). On the original machine the two are identical.
_TOP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# One body type size for this figure; every other size is a ratio to it.
BODY = 7.0
st.apply(BODY)
# 300, not the house default of 600. Figures 1 and 2 both set this, and this one never did,
# so its PNG came out at twice the canvas and every pixel measurement taken off it, by me or
# by the gate tools, was out by a factor of two. One display pixel is now one output pixel.
plt.rcParams["savefig.dpi"] = 300

def find(n):
    for d in (".", "..", os.path.join(_TOP, 'ccle_plex05')):
        p = os.path.join(d, n)
        if os.path.exists(p): return p
    return None

# ---------------------------------------------------------------- gather margins and labels
CACHE = "fig3_margins.csv"
if os.path.exists(CACHE):
    D = pd.read_csv(CACHE)
    print("margins from cache: %s rows" % f"{len(D):,}")
else:
    rows = []
    for pf in sorted(glob.glob("plex*_varpair.txt")):
        plex = pf.split("_")[0]
        g = glob.glob("out_%s_pep/*/psm.tsv" % plex)
        if not g:
            continue
        lab = {}
        with open(pf) as fh:
            fh.readline()
            for ln in fh:
                f = ln.split("\t", 5)
                if len(f) > 3 and f[3] in ("target", "p_target"):
                    lab[f[0]] = f[3]
        best = {}
        for r in csv.DictReader(open(g[0]), delimiter="\t"):
            p = r.get("Peptide", "").strip()
            if p not in lab:
                continue
            try:
                d = float(r.get("Hyperscore") or 0) - float(r.get("Nextscore") or 0)
            except ValueError:
                continue
            if p not in best or d > best[p][0]:
                best[p] = (d, lab[p])
        for p, (d, l) in best.items():
            rows.append(dict(plex=plex, peptide=p, margin=d, label=l))
    D = pd.DataFrame(rows)
    # Never write an empty cache. In a deposit without the raw search output this branch finds no
    # psm.tsv, and an empty fig3_margins.csv here is worse than none: it is written once, then
    # every later run reads it and dies on "No columns to parse from file" with no clue why.
    if D.empty:
        raise SystemExit(
            "fig3: no PSM tables found, so no margins could be computed, and an empty cache was "
            "NOT written.\n       This figure needs either fig3_margins.csv (deposited) or the "
            "out_plex*_pep/*/psm.tsv search output.")
    D.to_csv(CACHE, index=False)
    print("margins computed: %s rows across %d plexes" % (f"{len(D):,}", D.plex.nunique()))

S = pd.read_csv(find("downstream_summary.csv"))

# 6 Nov: THREE ROWS TO TWO, and the sheet onto the family's width. Figures 1 and 2 are
# landscape two-row sheets at an 11.5333 in canvas; this was a 14.41 x 16.12 in portrait
# sheet in three rows, a|b / c / d|e|f, which could not sit beside them. Same six panels,
# same letters, same reading order, re-laid as a|b|c over d|e|f.
#
# The outer grid is 2 x 1 and each row carries its own subgridspec, the Figure 2 pattern,
# because the two rows need different width ratios: c has 42 rotated plex labels and needs
# the widest frame in row 0, while row 1 is three plots of similar appetite.
# Height set so the drawn y axes match FIGURE 1's, which is the reference for this family
# the same way it is for width, raster and type. Figure 1's panels measure 483.6 px.
#
# THE FIRST ATTEMPT AT THIS WAS WRONG AND THE ERROR IS WORTH KEEPING. It read 488 px off
# the published PNG, found 501 here, and scaled the canvas by 488/501 to 4.2966. That put
# the canvas BELOW the height at which the bottom row's tick labels fit, so _fit_once saw
# content past the edge and scaled EVERY axes by the same factor, which is all-or-nothing.
# The panels landed at 468.3, not 488, and panel b at 460.3 because it was shrunk twice.
# Ratio arithmetic cannot cross that cliff: the shrink is a step, not a slope.
#
# The margins are therefore held in INCHES, not in figure fractions. Fractional margins
# move with the canvas, so lowering the height also lowers the room the x labels have and
# walks the figure into the shrink. With _MB and _MT fixed, the height solves in closed
# form: h = (H - _MB - _MT) / (rows + hspace) and H = 2.42 x 1.55020 + _MB + _MT. That is
# 4.36958, and it puts all five drawn y axes on 483.6 px with 0.0 px spread and no fit.
# ALL FOUR MARGINS IN INCHES, AND THE CANVAS SOLVED FROM THE PANEL BLOCK.
# The left margin is 0.780 rather than the 0.634 that left=0.055 used to give. Panel a has
# to slide 0.743 in left to put its letter on panel d's, because d's leftmost ink is a tick
# label and a's is a box edge, and 0.743 in did not fit inside 0.634: a's boxes ran off the
# canvas, _fit_once fired, and the whole sheet came down from 483.6 px to 476.6 with panel b
# at 469.8. The canvas is 0.15 in wider to hold that slide. Nothing in the figure moves or
# resizes because of it: the panel block keeps its own 10.726 x 3.75148 in, and the extra
# margin is outside the trim, so the published file is unchanged except for panel a.
_ML, _MR = 0.780, 0.173        # left and right margins, inches
_MB, _MT = 0.3753, 0.2428      # bottom and top margins, inches
_YAX = 1.55020                 # the family y axis, 483.6 px at the 311.96 raster
_HS = 0.52                     # row spacing, a fraction of the panel height
# _HS went 0.42 -> 0.52 on 7 November to drop row 1 away from row 0; the same 0.10, which is
# 48.4 px, goes into Figure 4. The gap and the block height are tied here rather than written
# as two numbers: matplotlib spends the block on panels and gaps together, so raising hspace
# against a held _PH would have taken the 48.4 px out of all five y axes instead of the canvas.
_PW, _PH = 10.726, (2 + _HS) * _YAX      # the panel block itself
_FW, _FH = _PW + _ML + _MR, _PH + _MB + _MT
fig = plt.figure(figsize=(_FW, _FH))
gs = fig.add_gridspec(2, 1, hspace=_HS,
                      left=_ML / _FW, right=1 - _MR / _FW,
                      top=1 - _MT / _FH, bottom=_MB / _FH)
gs_top = gs[0, 0].subgridspec(1, 3, wspace=0.34, width_ratios=[1000, 820, 1180])
gs_bot = gs[1, 0].subgridspec(1, 3, wspace=0.42, width_ratios=[720, 900, 900])

# ============================================== six panels
# Panel d of the previous version — the Wen et al. estimator over the whole plex-2 search —
# stated in its own annotation that the populations differ and the figure is not comparable
# with the 0.58% next to it. A main-text panel showing a number that cannot be compared to
# our own belongs in the supplement, and is now figS_fdr h.


# ---------------------------------------------------------------- a: the design
ax = fig.add_subplot(gs_top[0, 0]); st.head(fig, ax, "a")
ax._panel_head_flush = True
# a's leftmost ink is its box edge; d's is a y tick label, which reaches further left. With
# the letters flush to their own ink, a's landed 166 px right of d's and the two row openers
# stopped sharing the figure's left edge, which is the one alignment the figures-final README
# names. The frame moves left to put them back. Converged against the render.
# RE-CONVERGED 6 Nov after the canvas height changed. This number cannot survive a geometry
# change: moving a's frame left also moves a's INK left, _fit_once then pushes the whole page
# right to keep it inside the canvas, and both letter columns travel. At -175 the pair had
# drifted to A 62.6 against D 17.6. The window where headrow groups them as one column is
# -222 to -224, so -223 sits in the middle of it. Any future change to this figure's height
# or margins has to be followed by a headrow read of the letter columns.
ax._axes_dx_in = -220.0 / 300.0
# And 184 px wider. At the house body size "MSFragger 4.4.1 . Percolator . Philosopher"
# ran past the white card and the last word was clipped by the border. The width is there
# to take: the a|b gap measured 266 px against the 100 px this paper packs rows at, so
# widening a both fixes the overflow and closes the gap.
# 184 px wider for the text, then 127.3 px back off for the frame slack, so +56.7 net.
# ONE assignment: these are attributes, not accumulators, and writing the second line
# under the first silently discarded the 184 and left the frame 127 px NARROWER than it
# started, which took 80 px off the boxes before the measurement caught it. Panel a's axis is OFF, so its frame is invisible, but
# intrude.py reads frames and panel a's ran 190 px past its own boxes. The moment b came
# left to close the gap, b's y-axis label was sitting 100% inside that empty right half.
# Exactly the defect panel a of Figure 2 had, and the same fix: take the slack off the
# frame, and rewrite the panel's own constants in PIXELS over PXDU_A so the drawing does
# not shrink with it. Net width is 1052.6 - 127.3 = 925.3 px.
# RESTORED TO THE PRE-RESIZE DRAWING. Panel a's contents are laid out in data units over a
# fixed 0..10 in both directions, so the drawing scales with the FRAME while its type does
# not. Taking the figure from 4.2966 to 4.36958 in grew the frame from 962.9 x 468.3 px to
# 968.2 x 483.6 and the boxes grew 3.3% taller around text that had not changed, which is
# the air the PI saw. The frame is put back to 962.9 x 468.3 exactly, so panel a renders
# pixel for pixel as it did. _axes_dh_in holds the TOP edge, so the 15.3 px comes off the
# bottom of a cell whose axis is off and which nothing is drawn in.
ax._axes_dw_in = 51.60 / 300.0
# AND THE PANEL GROWS DOWNWARD. The shrink that was here, -15.3 px to put the frame back to
# its pre-resize 468.3, made the panel WORSE and the PI caught it: the connectors clipped the
# descenders of "acquired-lysine" and the sub-captions sat on both borders of their white
# bodies. Measuring showed why, and it was never a resize problem. A two-line sub-caption is
# 84.7 px tall and the white body it sits in is 1.73 data units, which is 80.9 px. The text
# has always been taller than the box that holds it; the resize changed that by 3 px and the
# restore changed it back.
# Three boxes with 10 px of air around an 84.7 px caption, a 46 px header band and two
# connectors need 608.6 px and this cell is 483.6. Panel a's axis is off and the row gap
# below it is 203 px of white, of which the D letter needs about 60, so the frame takes 125
# of it. Nothing else in the figure moves: _axes_dh_in holds the top edge.
A_GROW_PX = 125.0
ax._axes_dh_in = A_GROW_PX / 311.96
# ylim is 10.8 rather than 10 because the three boxes grew. Each holds a bold title over a
# three-line sub-caption, and at the raised body the sub-caption's top line was printing
# through the title (42% on all three boxes). The box is taller and the two are further
# apart; the extra 0.8 units keeps the top box inside the panel.
# Each box is a bold title over TWO lines of specification. It was three lines, and the
# third of each was commentary rather than specification: "target-decoy FDR behaves here",
# "8,470 sites, DepMap 26Q1", "every hit is an observed false positive". Those went to the
# caption, with the sentence that used to sit under the boxes, on the same rule that moved
# panel f of Figure 2. At the house body size the three-line form did not fit any frame
# this row can give panel a, and the fourth line was prose in a schematic either way.
ax.set_xlim(0, 10); ax.set_ylim(0, 10.0); ax.set_axis_off()
# 6 Nov, rebuilt. Three changes, all asked for.
#
# COLOUR. Each box was a white card with a coloured hairline, so the only colour in the
# panel was a 1.1 pt outline and the stage 1 title printed in C5C3C6, the palette's
# faintest value, which is a rule colour and not a type colour. Each box now carries a
# SOLID COLOURED HEADER with the title in paper over a white body, and the colour runs
# round the body as a hairline border. No tints: the header is the palette colour at full
# strength, the body is paper, and the contrast comes from which of the two carries type.
# Stage 1 takes 46494C rather than C5C3C6, because white on light grey cannot be read, and
# the three then read as the plain reference search, our variant search and the control.
#
# PADDING. Text sat 0.4 units inside a box edge. It is 0.55 now, and the header band gives
# the title its own 0.95 units, so neither line is against a border.
#
# ARROWS. "-|>" at the default mutation_scale draws a head about 42 px long, and the whole
# connector was 0.6 units, about 42 px, so the head WAS the arrow and it printed as a bare
# triangle. mutation_scale 9 puts the head at about 12 px and leaves a 33 px shaft at
# lw 1.6, which reads as an arrow.
# WIDTH IS MEASURED, NOT SET, AND IT IS ONE WIDTH FOR ALL THREE. They were all 9.0 units,
# a number, which fitted the longest line in box 1 and left box 2 ending at 72% of its card
# and box 3 at 63%. Sizing each box to its OWN longest line removed that slack and gave
# three ragged widths, which the PI did not want: the boxes are a set and read as one.
# So the width is the widest line in the PANEL plus the padding, measured off FT2Font
# ADVANCE widths at the drawn size, the same way panel a of Figure 2 lays out its peptide.
# Advances, not the ink box, because a text drawn ha="left" is anchored on its layout box.
# The panel still loses the slack 9.0 was carrying; what it does not do is make the loss
# different for each box.
#
# HEIGHT. The body is 1.73 units, 182 px, against 110 px of type, so 36 px of air above
# and below the two lines instead of the 23 px they had. Taken as a 2.7 unit box with a
# 0.9 unit header rather than 2.8 and 1.0, because the first attempt at 2.8 left only 0.6
# units between boxes and the arrows went back to being bare triangles. The gap is 0.95
# now, 100 px, which is what an arrow needs to read as one.
PXDU_A = 925.3 / 10.0                      # px per data unit across panel a's frame
def _aw(txt, bold=False):
    """Advance width of a string at the drawn body size, in data units."""
    from matplotlib import font_manager as _fm
    from matplotlib.ft2font import FT2Font as _FT
    fp = _fm.FontProperties(family=plt.rcParams["font.sans-serif"][0], size=st.sz(7),
                            weight="bold" if bold else "normal")
    f = _FT(_fm.findfont(fp)); f.set_size(st.sz(7), 300)
    tot = 0.0
    for ch in txt:
        f.set_text(ch); tot += f.get_width_height()[0] / 64.0
    return tot / PXDU_A
# BOX_X and PADX are PIXELS over PXDU_A, not data-unit literals, so narrowing the frame
# above moved the frame and not the drawing. The vertical three are unaffected by a width
# change and stay in data units.
BOX_X, PADX = 52.6 / PXDU_A, 57.9 / PXDU_A
# VERTICAL GEOMETRY IN PIXELS, for the same reason the horizontal geometry is in pixels.
# BOX_H and HEAD_H were data-unit literals, so the boxes scaled with the frame while the type
# inside them did not. Everything vertical is now derived from the MEASURED text instead.
PXDU_AY = (483.6 + A_GROW_PX) / 10.0   # px per y data unit across panel a's frame
# Two-line caption, air each side of it, header band, inner inset, and the clear space the
# stack keeps above the D and E letters that sit in the row gap below this panel. The first
# build of this left 3.2 px there, which reads as touching.
SUB_PX, AIR_PX, HEAD_PX, INSET_PX, BOTTOM_PX = 84.7, 8.0, 46.0, 4.0, 24.0
INSET = INSET_PX / PXDU_AY
BODY_H = (SUB_PX + 2 * AIR_PX) / PXDU_AY
HEAD_H = HEAD_PX / PXDU_AY
BOX_H = BODY_H + HEAD_H + INSET
GAP_H = (10.0 - BOTTOM_PX / PXDU_AY - 3 * BOX_H) / 2.0
_Y0 = 10.0
_YS = [_Y0 - BOX_H, _Y0 - 2 * BOX_H - GAP_H, _Y0 - 3 * BOX_H - 2 * GAP_H]
print("  a box %.1f px (white body %.1f for %.1f px of caption, header %.1f), gap %.1f px, "
      "%.1f px clear below the stack"
      % (BOX_H * PXDU_AY, BODY_H * PXDU_AY, SUB_PX, HEAD_H * PXDU_AY, GAP_H * PXDU_AY,
         _YS[2] * PXDU_AY))
boxes = [(_YS[0], "stage 1, reference space", st.GREEN,
          "MSFragger 4.4.1 \u00b7 Percolator \u00b7 Philosopher\n"
          "MANE v1.2 + contaminants, reversed decoys"),
         (_YS[1], "stage 2, variant space", st.BLUE,
          "same engine, database restricted to\nacquired-lysine peptides not in stage 1"),
         (_YS[2], "paired entrapment", st.ORANGE,
          "FDRBench 1.1.1 \u00b7 shuffled partners of\nmatched length and composition")]
BOX_W = max(max([_aw(_t, bold=True)] + [_aw(_l) for _l in _s.split("\n")])
            for _, _t, _, _s in boxes) + 2 * PADX
for y, title, col, sub in boxes:
    # outer card, solid colour; inner card, paper, inset so a header band and a hairline
    # border of the colour remain
    ax.add_patch(FancyBboxPatch((BOX_X, y), BOX_W, BOX_H, boxstyle="round,pad=0.06",
                                facecolor=col, edgecolor="none"))
    ax.add_patch(FancyBboxPatch((BOX_X + INSET, y + INSET), BOX_W - 2 * INSET,
                                BOX_H - HEAD_H - INSET, boxstyle="round,pad=0.04",
                                facecolor=st.PAPER, edgecolor="none"))
    ax.text(BOX_X + PADX, y + BOX_H - HEAD_H / 2.0, title, fontsize=st.sz(7),
            fontweight="bold", color=st.PAPER, va="center")
    ax.text(BOX_X + PADX, y + INSET + (BOX_H - HEAD_H - INSET) / 2.0, sub,
            fontsize=st.sz(7), color=st.GREY, va="center", linespacing=1.5)
# Sharp heads, and each arrow takes the colour of the box it LEAVES, so the eye follows
# the stage down rather than reading a neutral connector. Two things had to be traded
# against each other on an 83 px connector: a head long enough to come to a point, and a
# shaft long enough to be seen. head_length 1.4 at mutation_scale 7 gave the point and ate
# 41 px, half the arrow, so it read as a bare wedge. 1.2 at 6.5 is a 32 px head 28 px
# across over a 51 px shaft at lw 2.0, which is a head wider than its shaft and shorter
# than it, the proportions that read as an arrow. The x is fixed for both, under the
# titles, so the ragged box widths do not make the spine wander.
ARROW_X = BOX_X + 1.4
# TRIANGLES, NOT ARROWS, and each one hangs from the box it leaves. An arrow on a 66 px
# connector has to be a head long enough to come to a point and a shaft long enough to be
# seen, and at this size one of the two always loses. A triangle has no such problem: it
# sits on the bottom edge of its own box, it points at the next one, and it takes the colour
# of the box it leaves, so the eye still follows the stage down.
# Sized in PIXELS over the panel's two scales, which differ, so the shape is the triangle
# it is meant to be rather than whatever the frame's aspect makes of it.
# The triangle SPANS the gap: its flat top on the drawn bottom edge of the box it leaves and
# its tip on the drawn top edge of the box it points at. The height is therefore not a chosen
# number but the gap itself, less the 0.06 round pad each box adds, so it stays correct if
# the stack is ever resized again.
# EQUILATERAL. The height is fixed by the gap, so the width follows from it rather than being
# chosen: side = height x 2 / sqrt(3). It is converted through the frame's TRUE x scale, 96.29
# px per unit, not through PXDU_A, which is the 92.53 reference the box widths were written
# against. Using PXDU_A would have drawn the triangle 4% wide of equilateral.
PXDU_AX = 962.9 / 10.0
# SMALLER. Spanning the whole 58.9 px gap made the connector the loudest thing in the panel,
# so the height is set at 40 px and the width follows from equilateral. The tip no longer
# reaches the box below: with the flat top held on the box it leaves, a triangle cannot be
# both smaller and still touching, and the PI asked for smaller last.
TRI_H_PX = 40.0
_TH = TRI_H_PX / PXDU_AY
TRI_W_PX = TRI_H_PX * 2.0 / 3.0 ** 0.5
_TW = TRI_W_PX / PXDU_AX
for _i, col in ((0, st.GREEN), (1, st.BLUE)):
    _yt = _YS[_i] - 0.06                       # the drawn bottom edge of the box above
    ax.add_patch(Polygon([[ARROW_X - _TW / 2.0, _yt], [ARROW_X + _TW / 2.0, _yt],
                          [ARROW_X, _yt - _TH]], closed=True,
                         facecolor=col, edgecolor="none", zorder=3))
print("  a connector triangles %.1f wide x %.1f tall px, equilateral, hung on the box edge "
      "with %.1f px clear below in a %.1f px gap"
      % (TRI_W_PX, TRI_H_PX, (GAP_H - 0.12) * PXDU_AY - TRI_H_PX, GAP_H * PXDU_AY))

# ---------------------------------------------------------------- b: margin dose-response
ax = fig.add_subplot(gs_top[0, 1]); st.head(fig, ax, "b")
ax._panel_head_flush = True
# 80 px narrower, taken off the RIGHT so the a|b gap keeps the 100 px it is converged at
# and the b|c gap, which measured 67 px, gets the width back. This figure's offsets arrive
# whole: it trims rather than overflowing, so there is none of the 0.795 damping Figure 2
# applies, and -80 px of dw is -80 px on the page.
ax._axes_dw_in = -80.0 / 300.0
# And 233 px left, which takes the a|b gap from 332 to 100. Panel c moves the same 233 so
# the b|c gap it is converged at travels with it instead of reopening one place along,
# the lesson row 0 of Figure 2 taught.
ax._axes_dx_in = -233.0 / 300.0
th = np.arange(0, 31, 1)
fdp, n_id = [], []
for t in th:
    sub = D[D.margin >= t]
    nt = int((sub.label == "target").sum()); ne = int((sub.label == "p_target").sum())
    fdp.append(100 * ne / max(nt + ne, 1)); n_id.append(nt)
HEAD = {}
_ne_lad = []
for _t in th:
    _s = D[D.margin >= _t]
    _ne_lad.append(int((_s.label == "p_target").sum()))
fd.panel("fig3", "b", [dict(margin_threshold=int(_t), fdp_pct=_f,
                            n_variant_peptides=int(_n), n_entrapment_peptides=int(_e))
                       for _t, _f, _n, _e in zip(th, fdp, n_id, _ne_lad)],
         note="distinct peptides, best PSM each")
HEAD["b_n_rows"] = len(D)
for _t in (0, 2, 5, 10, 15, 20):
    _i = list(th).index(_t)
    HEAD["b_fdp_margin%d_pct" % _t] = round(fdp[_i], 2)
    HEAD["b_n_variant_margin%d" % _t] = int(n_id[_i])
    HEAD["b_n_entrapment_margin%d" % _t] = int(_ne_lad[_i])
# COLOUR FOLLOWS THE SERIES, and the threshold marker is neutral. The left axis and its
# solid curve are orange and always were. The right axis, its label, its ticks and its
# dotted curve were all 4C5C68, the same slate the figure uses for annotation text, so the
# second series had no colour of its own; they are teal now. The vertical marker at
# margin 10 and the block beside it were teal, which read as a third series rather than as
# a threshold; they are grey. Each axis is now the colour of the line it measures, and the
# only thing that is not a series is the only thing drawn in a neutral.
ax.plot(th, fdp, color=st.ORANGE, lw=1.4)
ax.axvline(10, color=st.GREY, lw=1.0, ls="--")
# THE TICKS ARE SET, THE LIMITS KEEP THEIR AIR. The only thing wrong with the panel as it
# stood was the left y axis, which ran 0 to 5 by 1 and spent six labels on a curve with one
# feature; it is 0 to 6 by 2, which is four. Everything else about the framing was right
# and clamping the limits onto the ticks was my addition, not an instruction. Flush limits
# put both curves on the y spine and laid the flat orange tail along the x spine from
# margin 21 to 30.
#
# So: ticks at 0, 10, 20, 30 and 0, 2, 4, 6, with 1.5 of margin at each x end and 5% of
# the range below each y zero. Both curves still begin at x = 0 and therefore still begin
# together, which is the property to protect; what the margin changes is that they begin
# just inside the frame rather than drawn over it.
ax.set_xlim(-1.5, 31.5); ax.set_xticks([0, 10, 20, 30])
# THE TWO CURVES ARE MADE TO START AND END AT THE SAME HEIGHT. They share an x ladder and
# they are read against each other, so the eye wants to see one pair of endpoints, not two.
# With the right axis pinned at a round 0 to 1500, the teal curve's own first and last
# values fix the two page fractions the orange curve has to hit, and solving for the left
# limits is then arithmetic rather than a guess:
#
#     f0 = n_id[0]  / 1500          fraction of the frame the pair starts at
#     f1 = n_id[-1] / 1500          fraction it ends at
#     span = (fdp[0] - fdp[-1]) / (f0 - f1)
#     L0   = fdp[-1] - f1 * span
#
# Computed, not typed, so it survives the ladder being recomputed. It lands on a top of
# 5.00 to two decimals, which is why 0 to 5 by 1 still labels the axis exactly; the assert
# is there because nothing guarantees that and a drifting top would put the last tick off
# the end of the spine.
_R0, _R1 = 0.0, 1500.0
_f0 = (n_id[0] - _R0) / (_R1 - _R0)
_f1 = (n_id[-1] - _R0) / (_R1 - _R0)
_span = (fdp[0] - fdp[-1]) / (_f0 - _f1)
_L0 = fdp[-1] - _f1 * _span
_L1 = _L0 + _span
assert abs(_L1 - 5.0) < 0.02, "panel b: left axis tops out at %.3f, not 5; reset the ticks" % _L1
ax.set_ylim(_L0, 5.0); ax.set_yticks([0, 1, 2, 3, 4, 5])
print("b: curves matched at both ends, left axis %.4f to 5.00, right %g to %g"
      % (_L0, _R0, _R1))
ax.set_xlabel("score margin threshold", fontsize=st.sz(7))
ax.set_ylabel("false-discovery proportion (%)", fontsize=st.sz(7), color=st.ORANGE)
ax.tick_params(axis="y", colors=st.ORANGE)
# THE SPINE TOO, not just the ticks and the label. tick_params colours the tick marks and
# their numbers and stops there, so the left spine was still drawing in the default near
# black while everything attached to it was orange. The x spine is left alone: it is the
# one piece of furniture both series share, so it is the one piece that should not take a
# series colour.
ax.spines["left"].set_color(st.ORANGE)
ax.spines["right"].set_visible(False)
ax2 = ax.twinx()
ax2.plot(th, n_id, color=st.BLUE, lw=1.2, ls=":")
ax2.set_ylabel("identifications retained", fontsize=st.sz(7), color=st.BLUE)
# The right axis is the round one and the left is solved against it, so this stays 0 to
# 1500. The two zeros no longer sit at the same height, and they cannot: the orange curve
# ends at 0 and the teal ends at 155, so matching the ends and aligning the zeros are
# different requests. Matching the ends is the one that was asked for.
ax2.set_ylim(_R0, _R1); ax2.set_yticks([0, 500, 1000, 1500])
for _sp in ("left", "top", "bottom"):
    ax2.spines[_sp].set_visible(False)
# 4 labels, not the 7 the automatic 250 step gave. Rule across this paper: round limits,
# round step, 4 to 6 ticks, and the top tick IS the top of the axis.
ax2.tick_params(axis="y", colors=st.BLUE); ax2.spines["right"].set_visible(True)
ax2.spines["right"].set_color(st.BLUE)
i10 = list(th).index(10)
# This curve counts DISTINCT PEPTIDES, best PSM each, which is what the entrapment pairing is
# defined over. downstream_summary counts identifications per plex and gives 1,039 at the same
# threshold, 1,033 variant and 6 entrapment. Both are correct; they are not the same
# denominator, and the panel says so rather than letting the two numbers look like a
# contradiction.
# Anchored by its TOP at 4.85 rather than by its baseline at a fixed y. At the old 4.3
# with the axis topping out at 6 the block sat inside the frame; taking the axis to 5 left
# its first line printing above the top spine, outside the panel. va="top" ties it to the
# frame instead of to a number, so it cannot leave again when the limit moves.
ax.text(11, 4.85, "margin \u2265 10\nFDP %.2f%%\n%s distinct peptides"
        % (fdp[i10], f"{n_id[i10]:,}"), fontsize=st.sz(7), color=st.GREY,
        va="top", linespacing=1.45)
# "distinct peptides, best PSM each" moved to the caption, where it belongs: it is the
# panel's unit statement, which is prose, and it was the widest thing in panel b.

# ---------------------------------------------------------------- c: per-plex consistency
ax = fig.add_subplot(gs_top[0, 2]); st.head(fig, ax, "c")
ax._panel_head_flush = True
ax._axes_dx_in = -233.0 / 300.0          # follows b, see the note there
# TWO COLUMNS, NOT A RANKED INDEX. The panel used to put the 42 plexes on an x-axis ordered
# by yield. That ordering is a rank index: the claim is 6 entrapment identifications against
# 1,033 variant in 4 of 42 plexes, and which plex is 7th does not enter it. The index cost
# 3.439 in of sheet, the widest panel in the figure, plus a legend inside the axes and four
# "plex NN" labels. Drawn as two columns the same facts need 1.300 in, so the frame is
# narrowed by the difference and nothing else in the figure moves.
ax._axes_dw_in = -641.7 / 300.0          # 3.439 in allocated, 1.300 in drawn
Sp = S.sort_values("var", ascending=False)
mask = (Sp.var_entrap > 0).values
fd.panel("fig3", "c", Sp[["plex", "var", "var_entrap"]].rename(
    columns={"var": "n_variant_identifications", "var_entrap": "n_entrapment_identifications"}),
    note="identifications at margin >= 10, per plex")
HEAD["c_n_entrapment_identifications"] = int(S.var_entrap.sum())
HEAD["c_n_variant_identifications"] = int(S["var"].sum())
HEAD["c_n_total_identifications"] = int(S["var"].sum()) + int(S.var_entrap.sum())
HEAD["c_fdp_pooled_pct"] = round(
    100 * S.var_entrap.sum() / max(S["var"].sum() + S.var_entrap.sum(), 1), 2)
HEAD["c_n_plexes_with_entrapment"] = int((S.var_entrap > 0).sum())
HEAD["c_n_plexes"] = int(len(S))
HEAD["c_var_per_plex_min"] = int(Sp["var"].min())
HEAD["c_var_per_plex_max"] = int(Sp["var"].max())

def _fan(vals, span):
    """Deterministic symmetric spread for points sharing a y, so a column does not stack.
    Five plexes sit at 13 identifications and three at 27; drawn on one x they would be one
    dot and the column would understate how many plexes it holds."""
    out = np.zeros(len(vals), float)
    for v in np.unique(vals):
        idx = np.where(vals == v)[0]
        out[idx] = 0.0 if len(idx) == 1 else np.linspace(-span, span, len(idx))
    return out

# Teal, not 4C5C68. Teal means "variant" everywhere else in this figure: it is the
# identifications-retained curve in panel b, "variant (1,934)" in panel e and
# "variant targets (1,419)" in panel f.
_v = Sp["var"].values
# Bigger and semi-transparent, on the PI's instruction of 6 November, the same call he made
# for Figure 2 panel d. At s=22 the dots in the variant column do overlap, which is what
# makes the opacity report density rather than invent a colour, and both collections carry
# _alpha_ok so the translucency check records an approved exception instead of a defect.
_pc_v = ax.scatter(0 + _fan(_v, 0.26), _v, s=22, color=st.BLUE, alpha=0.65,
                   linewidths=0, zorder=3)
_e = Sp.var_entrap[mask].values
_pc_e = ax.scatter(1 + _fan(_e, 0.12), _e, s=30, color=st.ORANGE, alpha=0.65,
                   linewidths=0, zorder=4)
_pc_v._alpha_ok = "panel c overplotting, 42 plexes in one column at s=22"
_pc_e._alpha_ok = "panel c, matched to the variant column so opacity is not an emphasis channel"
ax.set_yscale("log"); st.plain_log(ax)
# Top at 100, a decade tick. The PI chose 0.5 for the bottom over the matched-clearance
# 0.8 I had proposed, so the two entrapment dots at 1 clear the bottom spine by 87 px
# against 30 px at 0.8, and the data occupies 82.5% of the panel rather than 90.5%.
# Candidates measured in figtools/bottom_3c.py; at a bottom of exactly 1 the dots sit ON
# the spine at 0.0 px, which is what ruled out the tidier two-decade framing.
ax.set_ylim(0.5, 100)
# The series are named by their own tick labels, so the legend that used to sit inside the
# axes is gone, and so are the four "plex NN" labels. The pooled 0.58% went with them: at a
# top of 100 the tallest point reaches 95% of the panel and there is no honest room for a
# text block. 0.58% is already in the body text ("Six in 1,039 is 0.58%").
# OPEN: the plex ids 30, 09, 38, 33 are in fig3_head.json and in the deposited fig3_c.csv,
# and they are cited in NO prose anywhere. The comment this replaces claimed they were "named
# in the caption", which was never true. Either the caption names them or the claim goes.
ax.set_xlim(-0.62, 1.62); ax.set_xticks([0, 1])
# Rung 6.3, the tick size, not rung 7. These are tick labels sitting in the tick role, and
# they are also what replaced the legend, which was rung 6.3 too. At rung 7 typescan showed
# xtick carrying two sizes, 7.974 x26 against 8.860 x2, which breaks one size per role.
ax.set_xticklabels(["variant\n%d plexes" % len(Sp), "entrapment\n%d plexes" % int(mask.sum())],
                   fontsize=st.sz(6.3))
ax.tick_params(axis="x", length=0)
# labelpad 0 rather than the matplotlib default of 4. The gap this closes is small,
# about 6 px: what actually separates the title from the spine is the width of the
# "100" tick-label column, not the pad.
ax.set_ylabel("identifications at margin \u2265 10", fontsize=st.sz(7), labelpad=0)
HEAD["c_plexes_with_entrapment"] = ", ".join("%02d" % p for p in Sp.plex[mask])
# Identifications per plex, from downstream_summary, on the estimator e/(t+e) the Methods
# state. Panel b counts DISTINCT PEPTIDES and gives 7 in five plexes against 6 in four here.
# They are NOT the same events: the extra one is a peptide at margin 10.53 in plex 34, which
# clears the threshold in the recomputed margins and is absent from this summary. Calling them
# one set of events under two denominators was wrong, and the legend now states the difference.
# The "6 in 1,039 identifications, FDP 0.58%, falling in 4 of 42 plexes" block moved to
# the caption. It is a statistic, and a statistic is prose: the same call made for panel f
# of Figure 2 and for panel b here. Every number in it is already in the ledger, as
# c_n_entrapment_identifications, c_n_total_identifications, c_fdp_pooled_pct,
# c_n_plexes_with_entrapment and c_n_plexes, so the caption reads them rather than the
# artwork. The key stays: it decodes the two colours and is not prose.


# ---------------------------------------------------------------- d: SpecQuality
ax = fig.add_subplot(gs_bot[0, 0]); st.head(fig, ax, "d")
ax._panel_head_flush = True
# SpecQuality scores a spectrum on its own properties, with no reference to the search engine,
# the database or the FDR.
#
# SQS_GM CANNOT BE USED HERE. It is a geometric mean across features, and the hit spectra were
# scored from an MGF export in which the Intensity feature is zero (median 0.000 against 6.818
# for the mzML background). A single zero collapses a geometric mean, which is the entire
# 1.31-vs-6.96 difference an earlier version of this panel showed — a file-format artifact, not
# a quality difference. PeakCount and SNR came through the export intact and are compared here.
import glob as _g
SQ = None
for c in ("../SpecQuality", os.path.join(_TOP, 'SpecQuality')):
    if os.path.isdir(c):
        SQ = c; break
if SQ:
    H = pd.read_csv(os.path.join(SQ, "hit_spectra_mgf_features.tsv"), sep="\t")
    bg = []
    # [:20] is an ALPHABETICAL slice of 77 available runs, 4 "c" and 16 "d", excluding every
    # "f" and "g" run. It is not a chosen subset. Reading all 77 takes 1.2 s and MOVES BOTH
    # NUMBERS: background peak-count median 695 -> 629 and SNR 0.1855 -> 0.2275, so the
    # ratios go 0.91x -> 1.00x and 1.19x -> 0.97x. The PI reviewed this on 6 November 2026
    # and chose to keep the published figure rather than reopen the text. DO NOT change the
    # slice without him. The key and the captions now say 20 of the 77 runs.
    for f in sorted(_g.glob(os.path.join(SQ, "sq_out", "*_features.tsv")))[:20]:
        try:
            bg.append(pd.read_csv(f, sep="\t", usecols=["PeakCount", "SNR", "Intensity"]))
        except Exception:
            continue
    B = pd.concat(bg, ignore_index=True) if bg else pd.DataFrame()
    feats = [("PeakCount", "peak count"), ("SNR", "signal-to-noise")]
    xs = np.arange(len(feats)); w = st.grp(2)
    hv = [pd.to_numeric(H[f], errors="coerce").dropna() for f, _ in feats]
    bv = [pd.to_numeric(B[f], errors="coerce").dropna() for f, _ in feats]
    # each feature on its own scale: peak count is in the hundreds, SNR near 0.2
    for k, ((f, lab), a, b) in enumerate(zip(feats, hv, bv)):
        norm = b.median() if b.median() else 1
        ax.bar(k - w/2, a.median() / norm, w, color=st.BLUE)
        ax.bar(k + w/2, 1.0, w, color=st.LIGHTGREY)
        # the bar is a ratio, so the label on it is the ratio; the raw medians, which are on a
        # different scale entirely (hundreds vs 0.2), go underneath as text
        ax.text(k - w/2, a.median() / norm + 0.042, "%.2f\u00d7" % (a.median() / norm),
                ha="center", va="bottom", fontsize=st.sz(7), color=st.BLUE)
        ax.text(k + w/2, 1.0 + 0.042, "1.00\u00d7", ha="center", va="bottom",
                fontsize=st.sz(7), color=st.GREY)
    ax.set_xticks(xs); ax.set_xticklabels([l for _, l in feats], fontsize=st.sz(6.3))
    for _i, (_t, _c) in enumerate(((("acquired-lysine hits (%d)" % len(H)), st.BLUE),
                                   (("background (%s)" % f"{len(B):,}"), st.GREY))):
        # 0.985, not 0.96. With the counts back the key is 1.203 in wide in a 1.400 in
        # frame, so it runs the full width and the only separation available from the
        # "1.19x" label is vertical. At 0.96 the two boxes overlapped by 20%; lifting the
        # block 0.025 of the axes clears it with 17 px of paper to spare.
        # The two lines were set on a 0.072 pitch, which at 483.6 px is 34.8 px against a
        # 35 px line box: no leading at all, so they read as one block of type rather than
        # two entries. The first line alone lifts to 1.030, which opens 56 px between the
        # two. The second line then came up to 0.945, which leaves 41 px between the line
        # tops, 6 px of leading on a 35 px line box, and lifts it clear of "1.19x".
        ax.text(0.02, (1.030 if _i == 0 else 0.945), _t, transform=ax.transAxes,
                fontsize=st.sz(7), color=_c, va="top")
    # CONDENSED, on the Figure 1 panel g recipe. With the two text blocks gone to the caption
    # the widest ink left is "acquired-lysine hits (38)" at 1.305 in, and the binding
    # constraint is the two tick labels, "peak count" at 0.543 in and "signal-to-noise" at
    # 0.721 in. Floor between the category centres = (0.543 + 0.721)/2 + 0.10 clearance =
    # 0.732 in, which is an x span of 1.913 units centred on 0.5.
    # The key went too, for the same reason the two bottom lines did. It read
    # "acquired-lysine hits (38)" and "background (802,067)", and the caption already says
    # "the 38 acquired-lysine identifications ... 802,067 background spectra", word for word.
    # It was also the floor on the width: at a fixed x span the "1.19x" value label runs into
    # "background (802,067)" at every frame below 2.00 in, and even at 2.00 the two are
    # visually touching while scoring under the 12% gate. Teal is the hits and grey the
    # background throughout this figure, the grey bars carry "1.00x" and sit on the dashed
    # reference, and the y-axis label says "relative to background", so nothing is lost.
    ax.set_xlim(0.5 - 1.913 / 2, 0.5 + 1.913 / 2)
    ax._axes_dw_in = -(2.394 - 1.400)
    # 1.683 in of text on a 1.564 in axis overflowed the panel and ran into the letter
    # D at 17% once the canvas came down to Figure 1's axis height. Two lines, no words lost.
    ax.set_ylabel("median, relative\nto background", fontsize=st.sz(7), linespacing=1.4)
    ax.set_ylim(0, 1.6); ax.set_yticks([0, 0.4, 0.8, 1.2, 1.6])
    # The two lines that used to sit under this panel are gone to the caption. The SQS_GM note
    # was already in both captions word for word, so it was a duplicate. The raw medians were
    # NOT in prose anywhere and were added to both captions in the same pass, because the
    # standing rule is that moving a number to the caption never means dropping it.
    _drows = []
    for (f, lab), a, b in zip(feats, hv, bv):
        _drows.append(dict(feature=lab, hit_median=float(a.median()),
                           background_median=float(b.median()),
                           ratio_hit_over_background=float(a.median() / (b.median() or 1))))
        _k = lab.replace("-", "_").replace(" ", "_")
        HEAD["d_%s_ratio" % _k] = round(float(a.median() / (b.median() or 1)), 2)
        HEAD["d_%s_hit_median" % _k] = float(a.median())
        HEAD["d_%s_background_median" % _k] = float(b.median())
    fd.panel("fig3", "d", _drows, note="SpecQuality medians, hits vs background (20 of 77 runs, 802,067 spectra)")
    HEAD["d_n_hit_spectra"] = int(len(H))
    HEAD["d_n_background_spectra"] = int(len(B))

# ---------------------------------------------------------------- e: precursor mass error
ax = fig.add_subplot(gs_bot[0, 1]); st.head(fig, ax, "e")
ax._panel_head_flush = True
# Observed minus calculated, from the PSM tables. An earlier version plotted the fractional part
# of the observed mass alone, which is a property of the peptide's own composition and says
# nothing about error. A cluster at +1 Da is the signature that withdrew YBX1 E4K: one neutron,
# a deamidation or isotope error rather than the peptide claimed.
CACHE_DM = "fig3_deltamass.csv"
if os.path.exists(CACHE_DM):
    DM = pd.read_csv(CACHE_DM)
else:
    rows = []
    for f in sorted(_g.glob("out_plex*_pep/*/psm.tsv")):
        try:
            q = pd.read_csv(f, sep="\t", low_memory=False,
                            usecols=lambda c: c in ("Protein", "Delta Mass"))
        except Exception:
            continue
        if "Delta Mass" not in q.columns:
            continue
        pr = q.get("Protein", pd.Series([""] * len(q))).astype(str)
        rows.append(pd.DataFrame(dict(
            dm=pd.to_numeric(q["Delta Mass"], errors="coerce"),
            grp=np.where(pr.str.startswith("VAR|"), "variant", "background"))))
    DM = pd.concat(rows, ignore_index=True).dropna()
    DM.to_csv(CACHE_DM, index=False)
v = DM.loc[DM.grp == "variant", "dm"]
b = DM.loc[DM.grp == "background", "dm"].sample(min(200000, int((DM.grp == "background").sum())),
                                                random_state=0)
# Zoomed to the monoisotopic population. Over the full +/-1.5 Da range 91% of the mass sits in
# one bin at zero and the plot shows nothing; here the actual accuracy is visible, and variant
# errors sit on top of background rather than merely falling in the same bins.
_v1 = 100 * ((v - 1.00336).abs() < 0.05).mean()
_b1 = 100 * ((b - 1.00336).abs() < 0.05).mean()
_vz, _bz = v[v.abs() < 0.06], b[b.abs() < 0.06]
HEAD["e_pct_variant_within_0.05_da"] = round(100 * float((v.abs() < 0.05).mean()), 1)
HEAD["e_median_variant_da"] = round(float(_vz.median()), 4)
HEAD["e_median_background_da"] = round(float(_bz.median()), 4)
HEAD["e_plus1_variant_pct"] = round(float(_v1), 1)
HEAD["e_plus1_background_pct"] = round(float(_b1), 1)
HEAD["e_n_variant_psms"] = int(len(v))
bins = np.linspace(-0.06, 0.06, 49)
# PERCENT OF THE PLOTTED WINDOW, not density. matplotlib's density=True returns a
# probability density in units of 1/Da, so the old axis topped out at "160" meaning 160 per
# dalton, which is not a quantity a reader can check. The bins are fixed at 0.0025 Da, so
# the same histogram in percent is exact, and it matches panel 2e's "sites (%)" and the
# house forms "regions covered (%)" and "mutant recovered (%)".
_wz = np.full(len(_vz), 100.0 / len(_vz))
_wb = np.full(len(_bz), 100.0 / len(_bz))
ax.hist(_bz, bins=bins, weights=_wb, color=st.LIGHTGREY,
        label="background (%s)" % f"{len(b):,}")
ax.hist(_vz, bins=bins, weights=_wz, histtype="step", lw=1.4, color=st.BLUE,
        label="variant (%s)" % f"{len(v):,}")
# the dashed vertical at 0 is gone, on the same reasoning that removed 3d's horizontal one:
# there is a labelled x tick at 0, so the line drew the reference a second time.
_eb = np.linspace(-0.06, 0.06, 49)
_hv_, _ = np.histogram(_vz, bins=_eb); _hb_, _ = np.histogram(_bz, bins=_eb)
fd.panel("fig3", "e", [dict(bin_left_da=float(_eb[_k]), bin_right_da=float(_eb[_k + 1]),
                            n_variant=int(_hv_[_k]), n_background=int(_hb_[_k]))
                       for _k in range(len(_eb) - 1)],
         note="precursor mass error, monoisotopic window; background is a 200k sample")
ax.set_xlabel("observed \u2212 calculated mass (Da)", fontsize=st.sz(7))
ax.set_ylabel("PSMs (%)", fontsize=st.sz(7))
# 0 to 50 by 10, not 0 to 40. The tallest bin is 35.2% background and 32.2% variant, so 40
# is the tighter round ceiling, but the legend is 6.6 units tall and the peak sits at x = 0,
# the centre of the panel, directly under it. At a 40 ceiling the legend bottom landed at
# 33.1 and overlapped the 35.2 bar by 2.1, measured by figtools/legend_clear.py; tightening
# borderaxespad and labelspacing was not enough. 50 is the next round ceiling and buys the
# clearance. Six ticks, ending on a major tick.
ax.set_ylim(0, 50); ax.set_yticks([0, 10, 20, 30, 40, 50])
# starts and ends on a labelled tick, 5 of them rather than the 7 matplotlib was giving
ax.set_xlim(-0.06, 0.06); ax.set_xticks([-0.06, -0.03, 0, 0.03, 0.06])
# CONDENSED, but not to 2.000 in. The x-axis label at 1.802 in is the ink floor, and at
# 2.000 in the house bar-fit pass reported "slot allows only 2.76 pt for 48 bars, below the
# 3.0 pt floor". 48 bars at 3.0 pt need 144 pt of bar, which is 2.000 in, and the bars take
# 92% of the frame, so the frame has to be 2.000 / 0.92 = 2.174 in. 2.180 clears it.
ax._axes_dw_in = -(2.993 - 2.200)   # 2.180 gave 2.99 pt per bar after the height change
# Packed left to the house gap. Condensing d from 2.394 to 1.400 in left 525 px of ink-to-ink
# white between d and e, measured off the render by figtools/pxgaps.py against the 100 px
# house gap that Figure 1's top row sets (103, 100, 97). 425 px of that comes out here.
# ONE assignment: a second _axes_dx_in would silently discard the first.
ax._axes_dx_in = -440.0 / 311.96
# Pulled up against the top spine. The tallest bin is 35.2% of 40, so there are only 4.8
# units of headroom, and at matplotlib's default borderaxespad and labelspacing the second
# legend row sat at about 33.9 and the peak ran through it. Tightening the block is the
# only lever that does not move the ceiling: the peak is at x = 0, the centre of the panel,
# and "background (200,000)" is 1.086 in in a 2.000 in frame, so the legend reaches the
# centre from either side and no horizontal placement avoids it.
ax.legend(frameon=False, fontsize=st.sz(6.3), loc="upper left", handlelength=1.1,
          borderaxespad=0.12, labelspacing=0.22, handletextpad=0.5)
# Both text blocks moved to the caption, where drafted.py already carried them word for
# word. The "91% within 0.05 Da" block also sat ON the histogram's right tail, which the
# collision gate compares text to text and cannot see. Every number is still printed, in
# HEAD above, and every one reaches the caption.

# ---------------------------------------------------------------- f: paired score separation
ax = fig.add_subplot(gs_bot[0, 2]); st.head(fig, ax, "f")
ax._panel_head_flush = True
# The empirical version of the argument. An earlier panel here plotted false positives computed
# FROM the nominal FDR and then expressed them as a fraction of that same total, which is
# algebraically just the FDR — two flat lines proving nothing.
#
# This instead uses the pairing directly: every variant peptide has exactly one entrapment
# partner (peptide_pair_index in the varpair files). Where the two score distributions overlap
# is where an identification cannot be distinguished from a peptide known to be absent.
tg = D.loc[D.label == "target", "margin"]
en = D.loc[D.label == "p_target", "margin"]
bins = np.linspace(0, 45, 34)
# Two translucent fills overlap into a third colour. Outlines overlap into nothing.
# PERCENT OF EACH GROUP, not density, for the same reason as panel e. density=True returns
# 1/margin-unit, so the old ceiling of 0.15 meant 0.15 per margin unit. Bins are fixed at
# 1.3636 units, so percent is exact. Each curve is normalised to ITS OWN total, which is the
# point of the panel (1,419 against 70 cannot be compared on counts) and is said in the caption.
_wt = np.full(len(tg), 100.0 / len(tg))
_we = np.full(len(en), 100.0 / len(en))
# FILLED teal, OUTLINED orange, which is the construction panels 2e and 3e already use:
# fill the larger reference population, outline the one being judged against it. The old
# comment here warned that two translucent fills overlap into a third colour, and that is
# still true, which is why only ONE of the two is filled and the fill is fully opaque.
# The entrapment outline is solid now rather than dashed. Against a fill the two series are
# already told apart twice over, by colour and by weight, so the dash was a third encoding
# of the same distinction, and on 70 peptides in 33 bins it broke the curve into fragments.
# It also frees dashing to mean one thing in this figure, the threshold.
# stepfilled, so the 33 bins are ONE polygon rather than 33 abutting rectangles. As
# separate bars the white seams between them read as structure that is not in the data.
ax.hist(tg, bins=bins, weights=_wt, histtype="stepfilled", color=st.BLUE, lw=0, zorder=2,
        label="variant targets (%s)" % f"{len(tg):,}")
ax.hist(en, bins=bins, weights=_we, histtype="step", color=st.ORANGE, lw=1.5, zorder=4,
        label="entrapment (%d)" % len(en))
# The threshold marker STAYS. In panel e the dashed vertical sat on zero, a labelled tick,
# so it drew the axis twice and went. Here the x ticks are 0, 15, 30, 45 and the threshold is
# 10, which falls between two of them and cannot be read off the axis at all. It is also the
# quantity this panel exists to show, and panel B marks the same 10.
ax.axvline(10, color=st.LIGHTGREY, lw=1.2, ls="--", zorder=5)
ax.set_xlabel("score margin", fontsize=st.sz(7))
ax.set_ylabel("peptides (%)", fontsize=st.sz(7))
# tallest bin is 20.0% entrapment and 7.8% variant. 0 to 20 by 5 is the roundest pair that
# holds them, and the precedent is panel 2e, where a 19.96% bar was set against a 20 ceiling
# on the PI's instruction: a round ceiling outranks clearance above the tallest bar.
ax.set_ylim(0, 20); ax.set_yticks([0, 5, 10, 15, 20])
# starts and ends on a labelled tick
ax.set_xlim(0, 45); ax.set_xticks([0, 15, 30, 45])
# CONDENSED. Widest ink is the legend's "variant targets (1,419)" at 1.085 in, so 1.700 in
# carries it with room. a, b, c, d and e do not move.
# Stretched to 2.854 in. At 1.700 the row ended 360 px short of row 1 and the two score
# curves were cramped. The family tolerates some ragged right, Figure 2's rows differ by
# 283 px, but Figure 1's differ by 16, and Figure 1 is the reference the raster comes from,
# so this goes flush. The sheet grows with it, which is the wanted direction: Figure 3's
# native is well under the family's 11.540 in, not over.
ax._axes_dw_in = -(2.993 - 2.854)
# Packed left to the house gap, the same move condensing d forced on e. Taking f's prose
# off and narrowing the frame left 951 px of ink-to-ink white between e and f against the
# 100 px house gap, so 851 px comes out here. ONE assignment of _axes_dx_in.
ax._axes_dx_in = -825.0 / 311.96   # 56 px less, because e grew 0.180 in to its right
ax.legend(frameon=False, fontsize=st.sz(6.3), loc="upper right", handlelength=1.1)
above = (en >= 10).sum()
fd.panel("fig3", "f", D[["plex", "peptide", "margin", "label"]],
         note="every distinct peptide with its best score margin; label target = variant, "
              "p_target = entrapment")
HEAD["f_n_entrapment_above_margin10"] = int(above)
HEAD["f_n_entrapment_peptides"] = int(len(en))
HEAD["f_n_variant_peptides"] = int(len(tg))
HEAD["f_pct_entrapment_below_margin10"] = round(100 * float((en < 10).mean()), 1)
# The two-line block that used to sit under this panel is in the captions. Its second
# sentence, that panel c counts identifications rather than peptides, was already there in
# the italic note after (F); its first, 7 of 70 clearing the threshold and 90% sitting below
# it, was in no prose anywhere and was added to both captions in the same pass.
st.save(fig, "fig3")
if __import__("os").environ.get("FIG3_PROBE"):
    _r = fig.canvas.get_renderer()
    _L = {getattr(_a, "_panel_letter", "?"): _a for _a in fig.axes}
    for _ln in "abcdef":
        _a = _L.get(_ln)
        if _a is None: continue
        _bb = [_a.get_window_extent()]
        for _ch in (list(_a.texts) + list(_a.collections) + list(_a.patches) + list(_a.lines)
                    + [_a.xaxis.label, _a.yaxis.label]
                    + _a.get_xticklabels() + _a.get_yticklabels()):
            try:
                _e = _ch.get_window_extent(_r)
                if _e.width > 0 and _e.height > 0: _bb.append(_e)
            except Exception: pass
        _lg = _a.get_legend()
        if _lg is not None:
            try: _bb.append(_lg.get_window_extent(_r))
            except Exception: pass
        print("probe %s  ink %8.1f .. %8.1f  (w %7.1f)   frame %8.1f .. %8.1f"
              % (_ln, min(e.x0 for e in _bb), max(e.x1 for e in _bb),
                 max(e.x1 for e in _bb) - min(e.x0 for e in _bb),
                 _a.get_window_extent().x0, _a.get_window_extent().x1))
# Trim to the drawing, 18 px all round, the figures-final standard that Figure 1 set and
# Figure 2 follows. Measured off the RENDERED PIXELS rather than get_tightbbox, because an
# axis-off panel still contributes its rectangle to the tight bbox and panel a has its axis
# off. This figure had no trim at all and carried whatever margin the gridspec left.
_PAD_IN = 18.0 / 300.0
from PIL import Image as _Im
_px = np.asarray(_Im.open("fig3.png").convert("RGB")).astype(int)
_ys, _xs = np.nonzero(_px.sum(2) < 735)
_H_in = fig.get_size_inches()[1]
_trim = Bbox.from_extents(_xs.min() / 300.0 - _PAD_IN,
                          _H_in - _ys.max() / 300.0 - _PAD_IN,
                          _xs.max() / 300.0 + _PAD_IN,
                          _H_in - _ys.min() / 300.0 + _PAD_IN)
for _ext in (".png", ".pdf"):
    fig.savefig("fig3" + _ext, bbox_inches=_trim)
print("trimmed to %.3f x %.3f in (%.0f x %.0f px at 300 dpi), 18 px pads"
      % (_trim.width, _trim.height, _trim.width * 300, _trim.height * 300))

if __import__("os").environ.get("FIG3_FINAL"):
    # The finalized file for CCB-submission/figures-final/, on Figure 2's block exactly.
    # THE FAMILY SHARES A RASTER, NOT A WIDTH. Figure 1 at 12.00 in fixes the raster at
    # 300 x 12.00 / 11.5400 = 311.96 dpi. Figure 3 renders at that same raster and is
    # declared 300 dpi, so its body type prints at the same 9.213 pt as Figures 1 and 2.
    # Scaling this figure to 12.00 in instead would scale it by 12.00 / 8.37 = 1.434 against
    # Figure 1's 1.0399 and print its type 38% larger, which is the one thing the three
    # figures must not do. A narrower figure is the honest outcome of having less in it.
    # Saved WITH _trim, because the trim here is a bbox_inches on save and not a resize, so
    # fig.get_size_inches() still reports the untrimmed canvas.
    import PIL.Image as _PI, pymupdf as _mu
    _ref = _mu.open("fig1.pdf")[0].rect.width / 72.0      # Figure 1's native width, measured
    _dpi = 300.0 * 12.0 / _ref
    fig.savefig("fig3_final.png", dpi=_dpi, bbox_inches=_trim)
    _im = _PI.open("fig3_final.png")
    _im.save("fig3_final.png", dpi=(300, 300))
    print("final: native %.4f in at the shared %.2f dpi raster -> %d x %d px declared 300 dpi "
          "(%.3f x %.3f in), body type %.3f pt against Figure 1's %.3f"
          % (_trim.width, _dpi, _im.size[0], _im.size[1], _im.size[0] / 300.0,
             _im.size[1] / 300.0, 8.860 * 12.0 / _ref, 8.860 * 12.0 / _ref))
fd.head("fig3", HEAD)
for _e in (".png", ".pdf"):          # review copies alongside fig4_v2 / figS_react_v2
    fig.savefig("../results/fig3_v3" + _e, dpi=400, bbox_inches="tight")
print("FDP at margin>=10: %.2f%% | variant %d, entrapment %d"
      % (fdp[i10], int(S["var"].sum()), int(S.var_entrap.sum())))
try:    print("  e precursor mass: %.0f%% within 0.05 Da | medians %.4f vs %.4f Da | +1 Da tail %.1f%% variant vs %.1f%% background"
              % (100 * (v.abs() < 0.05).mean(), _vz.median(), _bz.median(), _v1, _b1))
except Exception as _e: print("  e: %s" % _e)
print("  b margin ladder (distinct peptides, best PSM each):")
for _t in (0, 2, 5, 10, 15, 20):
    _sub = D[D.margin >= _t]
    _nt = int((_sub.label == "target").sum()); _ne = int((_sub.label == "p_target").sum())
    print("     margin >= %-2d  FDP %5.2f%%   variant %5d   entrapment %3d"
          % (_t, 100 * _ne / max(_nt + _ne, 1), _nt, _ne))
print("  c per-plex (identifications, downstream_summary): %d entrapment in %s identifications "
      "= %.2f%%; entrapment appears in %d of %d plexes"
      % (int(S.var_entrap.sum()),
         f"{int(S['var'].sum()) + int(S.var_entrap.sum()):,}",
         100 * S.var_entrap.sum() / max(S["var"].sum() + S.var_entrap.sum(), 1),
         int((S.var_entrap > 0).sum()), len(S)))
print("  NOTE two denominators: panel b counts DISTINCT PEPTIDES and finds 7 in five plexes, "
      "panel c counts IDENTIFICATIONS and finds 6 in four. Not interchangeable, and not the "
      "same events: one peptide at margin 10.53 in plex 34 clears the threshold in the "
      "recomputed margins and is absent from downstream_summary.")
