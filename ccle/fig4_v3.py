"""
Figure 4, proposed rebuild, six panels instead of nine.

WHAT CHANGED AND WHY

Old b and c were two ranked bar charts of the same features asked two different ways, but with
different feature lists in different orders, so the comparison the figure is built on could not
be made by eye. They are now one chart drawn twice (b, c) on one shared feature order: before
the gate and after it. The argument is now visible rather than asserted.

Old e, f and g were three panels on feature redundancy, two scree plots and an adjusted-vs-alone
chart. That is scaffolding, not result, and it moves to the SI.

Old i reported four correlation coefficients over 5 and 9 hand-classified features with no null.
It is replaced by e, an exact permutation test over all C(15,6) = 5005 ways of splitting these
features 5-against-9, with the ungated version of the same test as its control (p = 0.008 gated,
p = 0.487 ungated).

Old d and h are kept, cleaned: d loses its unlabelled twin-axis traces (moved to SI) and gains
distinguishable line styles; h moves its tier counts out of the plotting area.

HONEST LABELLING: under the gate only the LABELLING effect is recomputed on well-modelled
lysines. The ligandability effect is as measured. Panel c says so.
"""
import os, re, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.lines import Line2D
from matplotlib.transforms import Bbox, ScaledTranslation
from itertools import combinations
import abbasov_style as st
import figure_data as fd
# One body type size for this figure; every other size is a ratio to it.
BODY = 7.0
st.apply(BODY)
# 300, not the house default of 600. Figures 1, 2 and 3 all set this and this one never did,
# so its PNG came out at twice the canvas and every pixel measurement taken off it, by me or
# by the gate tools, was out by a factor of two. 4322 px read as 14.407 in at 300 when the
# sheet is really 7.205 in. One display pixel is now one output pixel.
plt.rcParams["savefig.dpi"] = 300   # these two scripts imported the style module for its colours but
             # never applied it, so their typeface, sizes, spines and DPI were
             # matplotlib defaults while every other figure used the BMC style
st.slate()         # palette from the supplied ramp; see abbasov_style.SLATE

def find(n):
    for d in [".", "..", "ligand_rebuild", "../ligand_rebuild"]:
        p = os.path.join(d, n)
        if os.path.exists(p): return p
    raise FileNotFoundError(n)

S9  = pd.read_csv(find("step9_results.csv"))
G   = pd.read_csv(find("step11_gates.csv"))           # reactivity, gate ladder
G2  = pd.read_csv(find("step11b_ligandability_gates.csv"))   # ligandability, same ladder
G2  = G2[~G2.dropped]
R11 = open(find("step11_report.txt")).read()
R11B = open(find("step11b_report.txt")).read()    # carries the engaged count panel a needs
R38 = open(find("step38_report.txt")).read()

NICE = lambda f: (f.replace("_8A", " 8Å").replace("_10A", " 10Å").replace("_", " ")
                   .replace("rsa tien", "exposure").replace("propka pKa", "predicted pKa")
                   .replace("plddt site", "pLDDT").replace("plddt 8Å", "pLDDT 8Å")
                   .replace("nz sasa", "NZ exposure").replace("sasa sidechain", "side-chain exposure")
                   .replace("dist centroid", "distance to centroid").replace("n lys", "nearby lysines")
                   .replace("salt bridge", "salt bridge"))
SHAPE = {"crowding_8A", "atoms_10A", "hydrophobic_8A", "aromatic_8A", "polar_8A",
         "rsa_tien", "sasa_sidechain", "nz_sasa", "dist_centroid"}
CCOL  = {"charge": st.BLUE, "shape": st.ORANGE}

# ---------------------------------------------------------------- shared feature table
# Both contrasts are now measured on the same pLDDT ladder (step 11 for reactivity, step 11b
# for ligandability), so b and c differ only in the gate — not in which side was gated.
# pLDDT features are excluded throughout: they are the gating variable, not a site property.
def ladder(df, col):
    return {g: df[df.gate == g].set_index("feature")[col] for g in sorted(df.gate.unique())}
RE, LI = ladder(G, "d"), ladder(G2, "d")
TIER = 60          # the occupancy threshold that defines "engaged" in comparison B
GATE = max(RE)
F = pd.DataFrame({"react": RE[0], "lig": LI[0],
                  "react_g": RE[GATE], "lig_g": LI[GATE]}).dropna().reset_index()
F = F.rename(columns={"index": "feature"})
F["kind"] = np.where(F.feature.isin(SHAPE), "shape", "charge")
F = F.sort_values("react", ascending=False).reset_index(drop=True)
NF, KC = len(F), int((F.kind == "charge").sum())
F["gap_0"]  = (F.react   - F.lig).abs()
F["gap_90"] = (F.react_g - F.lig_g).abs()

SPLITS = list(combinations(range(NF), KC))
def split_test(rc, lc):
    """Exact permutation over every way of splitting these features KC-against-the-rest."""
    R = lambda idx: np.corrcoef(F.loc[idx, rc], F.loc[idx, lc])[0, 1]
    ci, si = F.index[F.kind == "charge"], F.index[F.kind == "shape"]
    obs = R(ci) - R(si)
    g = np.array([R(pd.Index(c)) - R(F.index.difference(pd.Index(c))) for c in SPLITS])
    return obs, g, float((g >= obs).mean()), R(ci), R(si)

NULL = {"gated": split_test("react_g", "lig_g"), "ungated": split_test("react", "lig")}
LADDER = []
for g in sorted(RE):
    if g not in LI: continue
    T = pd.DataFrame({"react": RE[g], "lig": LI[g]}).dropna()
    T = T[T.index.isin(F.feature)]
    T["kind"] = np.where(T.index.isin(SHAPE), "shape", "charge")
    Rr = lambda k: np.corrcoef(T[T.kind == k].react, T[T.kind == k].lig)[0, 1]
    idx = {k: F.index[F.kind == k] for k in ("charge", "shape")}
    obs = Rr("charge") - Rr("shape")
    Fg = F.assign(rr=F.feature.map(T.react), ll=F.feature.map(T.lig))
    Rg = lambda i: np.corrcoef(Fg.loc[i, "rr"], Fg.loc[i, "ll"])[0, 1]
    perm = np.array([Rg(pd.Index(c)) - Rg(Fg.index.difference(pd.Index(c))) for c in SPLITS])
    LADDER.append(dict(gate=g, charge=Rr("charge"), shape=Rr("shape"),
                       p=float((perm >= obs).mean())))
LADDER = pd.DataFrame(LADDER)

# ---------------------------------------------------------------- the setting, from the report
# These four numbers used to be typed into panel a as a string. They are now read from the same
# report the effects come from, so the artwork cannot drift from the analysis.
_m = re.search(r"(\d+) lysines on (\d+) proteins \| labelled (\d+) \(([\d.]+)%\)", R11)
SET_NK, SET_NP, SET_NL, SET_PCT = (int(_m.group(1)), int(_m.group(2)),
                                   int(_m.group(3)), float(_m.group(4)))
# the third level of the nest, read from the ligandability report for the same reason: panel a
# draws it to scale, so it cannot be a number typed into a string
_m2 = re.search(r"engaged at >=60%: (\d+) \(([\d.]+)%\)", R11B)
SET_NE, SET_PCT_E = int(_m2.group(1)), float(_m2.group(2))
SET_NU, SET_NUE = SET_NK - SET_NL, SET_NL - SET_NE

HEAD = {"a_n_lysines": SET_NK, "a_n_proteins": SET_NP, "a_n_labelled": SET_NL,
        "a_pct_labelled": SET_PCT, "a_n_features": NF, "a_n_charge_features": KC,
        "a_n_unlabelled": SET_NU, "a_n_engaged": SET_NE, "a_pct_engaged_of_labelled": SET_PCT_E,
        "a_n_labelled_not_engaged": SET_NUE,
        "a_n_shape_features": NF - KC, "a_occupancy_tier_pct": TIER,
        "a_n_regroupings": len(SPLITS)}

# spread under the gate -- the range-restriction control quoted in 2.4
_SD = {g: G[G.gate == g].set_index("feature")["sd_all"] for g in sorted(G.gate.unique())}
for _f in ("crowding_8A", "polar_8A", "aromatic_8A", "rsa_tien", "dist_centroid"):
    if _f in _SD[0].index and _f in _SD[GATE].index:
        HEAD["bc_sd_%s_ungated" % _f] = round(float(_SD[0][_f]), 2)
        HEAD["bc_sd_%s_gated" % _f] = round(float(_SD[GATE][_f]), 2)
        _ch = 100 * (float(_SD[GATE][_f]) / float(_SD[0][_f]) - 1)
        HEAD["bc_sd_%s_pct_change" % _f] = round(_ch, 0)
        if _ch < 0:                       # 2.4 quotes these as spread LOST, positive
            HEAD["bc_sd_%s_pct_spread_lost" % _f] = round(-_ch, 0)

# ================================================================== canvas
# 11.5400 in, Figure 1's native width, not st.DOUBLE at 7.2047. Panel a's cell was 1.65 in
# and the nest it now draws needs about 3.0, and a figure of six content panels has no
# business being 40% narrower than every other figure in the set. wspace comes down from
# 0.62 to 0.40 at the same time: 62% of the average axes width in gutters was most of why
# the cells were starved. Figure 3 runs 0.34 and 0.42.
# HEIGHT AND MARGINS. Panels b, c, d and f fill their gridspec cells exactly, so the cell
# height IS the y-axis length, and the family's y axes are 483.6 px. Two things had to change
# before this figure could reach that number.
#
# THE MARGINS ARE HELD IN INCHES, NOT IN FIGURE FRACTIONS. _fit_once is all-or-nothing: when
# any drawn label runs past the canvas edge it scales EVERY axes by the same factor, so the
# relation between canvas height and panel height is a step, not a slope. Fractional margins
# shrink with the canvas, so lowering the height also lowers the room the bottom labels have
# and walks the figure straight into that step. Figure 3 had been set by ratio arithmetic
# across exactly this cliff and came out at 468.3 px instead of 488. With _MB and _MT fixed
# the height solves in closed form and lands on the target at the first build.
#
# THE SECOND WAS PANEL B, and that is why b and c are transposed. See their block below.
#
#   h = (H - _MB - _MT) / (rows + hspace)  ->  H = 1.55020 x 2.66 + 0.700 + 0.252 = 5.07553
#
# hspace is 0.66 rather than 0.40 because the gap between the rows is no longer empty: it
# carries b and c's 15 rotated feature names. 0.64 is the first value at which they clear
# panel letters E and F; 0.66 keeps a margin. _MB is 0.700 rather than 0.588 for the same
# kind of reason, panel d's sublabel now hanging a fixed 0.56 in below its axes.
# ALL FOUR MARGINS IN INCHES, AND THE CANVAS SOLVED FROM THE PANEL BLOCK. _PS is the sum of
# the three panel widths and it is HELD: matplotlib spends a gridspec's width on panels and
# gaps together, so dropping wspace with the block width fixed would have widened all six
# panels instead of closing the gaps. Solving the block from _PS and wspace keeps every panel
# exactly the width it had and takes the saving out of the canvas, which the trim removes.
_ML, _MR = 0.9000, 0.8000      # left and right margins, inches
# _MR went 0.2500 -> 0.6000 when panel f moved 85 px right to open the e|f gap. At
# 0.4500, and again at 0.6000 when f moved further, f's own ink ran past the canvas and
# _fit_once fired: every panel came down
# from 483.6 px to 477.9, and panel d to 472.2 because its twin axis took a second
# shrink from the tick fit on top of that. The margin is outside the trim, so it
# costs the published figure nothing.
_MB, _MT = 1.000, 0.252        # bottom and top margins, inches
_PS = 9.0754                   # sum of the three panel widths, held
# 8.5003 -> 9.1538 on 7 November, which puts this sheet's trimmed render at 3461 px, the width
# fig1.png has. THE TARGET IS THE FINALIZED WIDTH, NOT THIS FILE'S. Every finalized figure is
# the same trimmed canvas re-rendered at 311.96 dpi and declared 300, so fig1.png's 3461 px
# become Figure1.png's 3600 px, 12.00 in. Matching Figure1.png's 3600 px here instead would
# have finalized at 3744 px, 12.48 in, and printed this figure 4% wider than Figure 1 with its
# body type 4% larger. Nothing vertical changes; _PS only buys width.
_WS, _HS = 0.18, 0.72          # column and row spacing, fractions of the average panel
# _HS went 0.62 -> 0.72 on 7 November, the same 48.4 px Figure 3's row 1 drops by.
# _WS came down from 0.40. 0.16 is the last value at which panel f's y-axis title clears
# panel e's series labels and panel b's title clears panel a's lower block; 0.18 keeps a
# margin. _HS is 0.62 rather than 0.66 because the feature names are angled now instead
# of upright, which costs them their length times cos 45 in row height; 0.60 is the first
# value at which 'hydrophobic 8Å' clears the E and F letters under it.
_YAX = 1.55020                 # the family y axis, 483.6 px at the 311.96 raster
_PW = _PS * (1 + 2 * _WS / 3)
_PH = (2 + _HS) * _YAX
_FW, _FH = _PW + _ML + _MR, _PH + _MB + _MT
fig = plt.figure(figsize=(_FW, _FH))
# wspace 0.34 -> 0.62: the lower sub-panel of column 1 carries x tick numbers along its
# foot and panel f carries y tick numbers up its side, and at the raised body the two
# met in the gutter between the columns (−0.8 through 1.5 at 84%).
gs  = fig.add_gridspec(2, 3, hspace=_HS, wspace=_WS,
                       left=_ML / _FW, right=1 - _MR / _FW,
                       top=1 - _MT / _FH, bottom=_MB / _FH,
                       width_ratios=[1.02, 1.0, 0.92])

# ---------------------------------------------------------------- a: the design
# WAS two rounded boxes of prose plus the sentence "one contrast cannot answer both
# questions", carrying no data at all and asserting the nesting rather than showing it.
# The three populations were sitting unused in the same reports the effects come from.
#
# NOW the nest drawn to scale. 245,684 lysines, of which 30,805 labelled, of which 3,047
# engaged. The second bar IS the dark slice of the first, opened out, so the containment is
# visible and the old sentence is no longer needed.
#
# ORDERING. A first version put each contrast's heading between the two bars and the zoom
# connector then ran a dashed diagonal straight through the word REACTIVITY. The top
# contrast's heading sits ABOVE its bar and the bottom one's BELOW its bar, which leaves the
# band between them empty for the wedge.
#
# COLOUR. Neutrals only, and specifically the figure's OWN neutral. A first version used
# 46494C for the minority blocks and the headings, and an inventory showed 46494C appears
# in NO other panel of this figure: panel a was the only one carrying it, and the only one
# with neither teal nor rust. That is a coherence problem in its own right. Everything here
# is now C5C3C6 and 4C5C68, both of which b, c, d, e and f already use.
# What panel a must NOT do is borrow teal and rust. Teal and rust are locked to charge and shape across b, c, d, e and
# f, and the old panel set REACTIVITY in teal, which gave one colour two meanings inside one
# figure. The marker key is no longer a separate list either: a filled circle IS the
# reactivity contrast and an open diamond IS the ligandability contrast, so each sits on its
# own heading. The charge-against-shape colour key moves to panel b, where those colours are
# actually used.
ax = fig.add_subplot(gs[0, 0]); st.head(fig, ax, "a")
ax._panel_head_flush = True
# Panel a's own leftmost ink is a box edge and panel d's is a y-axis title, which reaches
# much further left, so flush alone leaves the two row openers on different columns. The
# frame slides left to put them back, the same move Figure 3 panel a makes. _ML is 0.9000
# rather than 0.6347 to hold the slide inside the canvas; the panel block is unchanged and
# the extra margin is outside the trim. Converged against the render.
ax._axes_dx_in = -186.0 / 311.96
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
_BH, _X0, _X1 = 0.115, 0.015, 0.985
_W = _X1 - _X0
# The right-hand count sat flush at _X1, which is also where the right zoom connector
# runs, so the dashed line grazed the last digit and had to stop 0.048 short of the
# lower block to stay off the text. Pulling the label left by _RPAD clears that corridor
# and lets both connectors run all the way onto the block they point at.
_RPAD = 0.022
_yA, _yB = 0.620, 0.250

def _nest_bar(y, frac, nleft, nright, lleft, lright):
    xs = _X0 + _W * (1 - frac)
    ax.add_patch(Rectangle((_X0, y), xs - _X0, _BH, fc=st.LIGHTGREY, ec="none", zorder=2))
    ax.add_patch(Rectangle((xs, y), _X1 - xs, _BH, fc=st.GREY, ec="none", zorder=2))
    ax.text(_X0 + 0.014, y + _BH / 2, "%s  %s" % (lleft, f"{nleft:,}"), va="center",
            ha="left", fontsize=st.sz(6.3), color=st.GREY)
    xm = (xs + _X1) / 2
    ax.plot([xm, xm], [y + _BH, y + _BH + 0.045], color=st.GREY, lw=0.8, zorder=3)
    ax.text(_X1 - _RPAD, y + _BH + 0.052, "%s  %s" % (lright, f"{nright:,}"), va="bottom",
            ha="right", fontsize=st.sz(6.3), color=st.GREY)
    return xs

# The two "what makes a ..." questions are gone from the artwork. They are word for word in
# the body text already, which also carries the "same protein" qualifier on the reactivity
# contrast, and the bar under each heading now shows what the contrast splits, so the
# question restated what the picture says. Removing them loses no fact.
def _contrast(y, marker, mfc, name):
    ax.plot([_X0 + 0.010], [y], marker=marker, ms=4.4, mfc=mfc, mec=st.GREY, mew=1.1,
            ls="none", zorder=4, clip_on=False)
    ax.text(_X0 + 0.040, y, name, fontsize=st.sz(7), fontweight="bold", color=st.GREY,
            va="center", ha="left")

_contrast(_yA + _BH + 0.105, "o", st.GREY, "REACTIVITY")
_sA = _nest_bar(_yA, SET_NL / SET_NK, SET_NU, SET_NL, "not labeled", "labeled")
_sB = _nest_bar(_yB, SET_NE / SET_NL, SET_NUE, SET_NE,
                "labeled, not engaged", "engaged at \u2265%d%%" % TIER)
_contrast(_yB - 0.080, "D", st.PAPER, "LIGANDABILITY")
for _xa, _xb in ((_sA, _X0), (_X1, _X1)):
    ax.plot([_xa, _xb], [_yA, _yB + _BH], color=st.LIGHTGREY, lw=1.0,
            ls=(0, (2.4, 1.8)), zorder=1, clip_on=False)
ax.text(0.5, 0.995, "%s lysines on %s proteins, every block to scale"
        % (f"{SET_NK:,}", f"{SET_NP:,}"), transform=ax.transAxes, ha="center", va="top",
        fontsize=st.sz(6.3), color=st.GREY)
# THE COLOUR HALF OF THE KEY BELONGS HERE TOO. I had moved it to panel b on the reasoning
# that a key belongs where its colours are used, but they are used in b, c, d, e AND f, so
# that argument picks one of five panels arbitrarily. Panel a is the design panel and the
# original figure defined both channels here. Shape is the contrast, carried by the two
# markers on the headings above; colour is the feature kind, carried by these two squares.
# One place to look either of them up.
for _dx, _c, _t in ((0.010, st.BLUE, "charge feature"), (0.330, st.ORANGE, "shape feature")):
    ax.plot([_X0 + _dx], [0.055], marker="s", ms=4.0, color=_c, ls="none", clip_on=False)
    ax.text(_X0 + _dx + 0.028, 0.055, _t, fontsize=st.sz(6.3), color=_c,
            va="center", ha="left")

# ---------------------------------------------------------------- b, c: the paired display
# TRANSPOSED, and the reason is the y axis. These were dumbbells with the 15 feature names
# running down the y axis, and that is what stopped Figure 4 reaching the family's 483.6 px:
# 15 tick labels with 37.4 px glyph boxes need 561.5 px of spine before any clearance at all,
# so the panel could not be shorter than about 600 px without colliding with itself. Built at
# 483.6 px the figure went from 1 text collision to 16.
# Turning the panel through 90 degrees moves the names OUT of the axis and under it, where
# their length costs ROW HEIGHT instead of AXIS LENGTH, and the effect-size axis is then free
# to be exactly as long as Figure 1's. No value changes and no feature is dropped. The pair
# is now vertical and the zero rule horizontal; everything else about the display is as it was.
for j, (rc, lc, gapcol) in enumerate([("react", "lig", "gap_0"),
                                     ("react_g", "lig_g", "gap_90")]):
    ax = fig.add_subplot(gs[0, 1 + j]); st.head(fig, ax, "bc"[j])
    if j == 0:
        # ROW 0 FILLS THE SHEET, the way both of Figure 1's rows do. Row 1 reaches 3441 px
        # and row 0 was stopping at 3281, 160 px short, which read as a short top row.
        # Figure 1 does not hold its gaps equal BETWEEN rows for this: it runs 61, 97 and
        # 98 px across its top row and 423 and 397 across its bottom, because each row
        # spreads its own panels over the same width. Row 0's three panels hold 2,996 px
        # of ink in a 3,346 px span, so the two gaps are 175 px each.
        # EQUAL GAPS WITHIN ROW 0. Measured off the rendered page, panel a's rightmost ink
        # ended at 885 px and panel b's leftmost at 1043, a 158 px gap, while b's rightmost
        # ink ended at 2073 and c's leftmost began at 2093, a gap of 20. The three panels
        # hold 1030 px of ink between a's right edge and c's left, leaving 178 px of white
        # to share, so panel b moves left to put it there, converged against the render.
        # RE-converged at -84 after the type ladder was fixed and panel a slid a further
        # 25 px left to keep its letter on D's: a's right ink travelled with it, which
        # opened a|b to 107 against b|c's 83. Now 96 px and 95 px. a and c do
        # not move; only the white between them is redistributed.
        ax._axes_dx_in = 6.5 / 311.96
    if j == 1:
        # Panel c moves right to put row 0 on the full sheet width. Closing _WS to 0.18 pulled all three columns
        # together evenly, and c had the least to give: its y-axis title and tick labels
        # ended up nearer panel b than a|b are to each other. The frame shifts without
        # resizing, and _MR goes from 0.1385 to 0.2500 in so the move stays inside the
        # canvas; the extra margin is outside the trim.
        ax._axes_dx_in = 201.5 / 311.96
    # Flush, the way every panel of Figures 1, 2 and 3 is placed. Those three set
    # _panel_head_flush on 6, 8 and 6 panels; this figure set it on none, so its letters
    # sat on the axes FRAME while each panel's y-axis title and tick labels stuck out to
    # the left of them. Flush puts the letter over the panel's own leftmost ink.
    ax._panel_head_flush = True
    x = np.arange(NF)
    for k, r_ in F.iterrows():
        c = CCOL[r_.kind]
        ax.plot([x[k]] * 2, [r_[rc], r_[lc]], color=st.LIGHTGREY, lw=2.0,
                solid_capstyle="round")
        ax.scatter(x[k], r_[rc], s=16, color=c, zorder=3, lw=0.0)
        ax.scatter(x[k], r_[lc], s=15, marker="D", facecolor="none", edgecolor=c, lw=1.0, zorder=4)
    # Dotted, not solid. The zero rule is a reference the eye crosses 15 times in this
    # panel and it is not data; a dotted rule reads as the reference it is.
    ax.axhline(0, color=st.LIGHTGREY, lw=0.8, ls=(0, (1.2, 2.0)))
    ax.set_xticks(x)
    # rung 6.3, not 5.0. Panels a, d and f already tick at 7.974 and b and e were at 6.645,
    # so "tick" carried two sizes. One size per role.
    # BOTH panels carry the names now. They shared one set when the features were on the y
    # axis and c's rows lined up with b's labels across the gap; side by side on x there is
    # no such alignment to borrow, so c labels itself.
    # Angled rather than upright. At 90 degrees the longest name, "nearby lysines 8Å",
    # costs its full 17 characters in row height; at 45 it costs them times cos 45 and
    # the row gap comes down with it. ha="right" anchors each label's last character
    # under its own tick, which is what makes an angled label readable as belonging to
    # one position rather than drifting between two.
    ax.set_xticklabels([NICE(f) for f in F.feature], fontsize=st.sz(6.3),
                       rotation=45, ha="right", rotation_mode="anchor")
    for t, kk in zip(ax.get_xticklabels(), F.kind): t.set_color(CCOL[kk])
    ax.set_xlim(-0.8, NF - 0.2)
    # Ends on a major tick at both ends, the house rule for a y axis. Data run -0.454 to
    # +0.352, so -0.5 to 0.5 holds every point with the axis symmetric about the zero rule,
    # which a reader comparing a positive effect with a negative one needs it to be.
    ax.set_ylim(-0.5, 0.5); ax.set_yticks([-0.5, -0.25, 0.0, 0.25, 0.5])
    ax.set_ylabel("effect size, d", fontsize=st.sz(7))
    g = F.groupby("kind")[gapcol].mean()
    # Bottom left. The features that reach below -0.25 are all at the right-hand end, so the
    # lower left corner of the frame is empty in both panels.
    ax.text(0.015, 0.02,
            "average distance between the pair\ncharge %.2f     shape %.2f" % (g["charge"], g["shape"]),
            transform=ax.transAxes, ha="left", va="bottom", fontsize=st.sz(6.3), color=st.GREY,
            linespacing=1.55)
    for s in ("top", "right"): ax.spines[s].set_visible(False)

fd.panel("fig4", "bc", F[["feature", "kind", "react", "lig", "react_g", "lig_g",
                          "gap_0", "gap_90"]],
         note="react/lig = ungated effect size d; _g = at pLDDT >= %d; gap = |react - lig|" % GATE)
_g0, _g9 = F.groupby("kind").gap_0.mean(), F.groupby("kind").gap_90.mean()
HEAD["b_mean_gap_charge"] = round(float(_g0["charge"]), 3)
HEAD["b_mean_gap_shape"] = round(float(_g0["shape"]), 3)
HEAD["c_mean_gap_charge"] = round(float(_g9["charge"]), 3)
HEAD["c_mean_gap_shape"] = round(float(_g9["shape"]), 3)
for _f in ("crowding_8A", "atoms_10A", "polar_8A", "anionic_8A", "salt_bridge"):
    _r = F[F.feature == _f]
    if len(_r):
        _r = _r.iloc[0]
        for _c, _n in (("react", "react_ungated"), ("react_g", "react_gated"),
                       ("lig", "lig_ungated"), ("lig_g", "lig_gated")):
            HEAD["bc_%s_%s" % (_f, _n)] = round(float(_r[_c]), 3)

# ---------------------------------------------------------------- e: the gate, step by step
ax = fig.add_subplot(gs[1, 0]); st.head(fig, ax, "d")
ax._panel_head_flush = True
# Row 1 re-evened after "dist. to centroid" and "side-chain exp." were spelled out in
# panels b and c: the longer names reach further left and down, which moved both rows.
ax._axes_dx_in = -13.5 / 311.96
gates = sorted(G.gate.unique()); xs = np.arange(len(gates))
SERIES = [("net_charge_8A", st.BLUE,   "-",  "o"), ("anionic_8A", st.BLUE,   "--", "s"),
          ("crowding_8A",   st.ORANGE, "-",  "o"), ("aromatic_8A", st.ORANGE, "--", "s")]
for f, c, ls, mk in SERIES:
    sub = G[G.feature == f].set_index("gate").reindex(gates)
    # clip_on=False: the frame now ends ON the >= 90 rung, so a 4.5 pt mark centred there
    # would lose its right half to the axes box. The lines stay inside it either way.
    ax.plot(xs, sub.d, ls=ls, marker=mk, color=c, lw=1.25, ms=4.5, mew=0,
            clip_on=False, label=NICE(f))
# Dotted, as in b and c. The zero rule is the same reference in all four panels and it
# should not read as a drawn series in one of them and a reference in the others.
ax.axhline(0, color=st.LIGHTGREY, lw=0.8, ls=(0, (1.2, 2.0)))
GATE_N, cur = {}, None
for ln in R11.splitlines():
    m = re.match(r"GATE\s+(\S+)", ln.strip())
    if m: cur = 0 if m.group(1) == "none" else int(re.sub(r"\D", "", m.group(1)) or 0)
    m2 = re.search(r"labelled ([\d,]+) \((\d+)% of all labelled\) \| unlabelled ([\d,]+)", ln)
    if m2 and cur is not None:
        GATE_N[cur] = (int(m2.group(1).replace(",", "")), int(m2.group(3).replace(",", "")))
fd.panel("fig4", "d", [dict(gate=int(_g_), feature=_f_,
                           effect_on_being_labelled_d=float(
                               G[(G.feature == _f_) & (G.gate == _g_)].d.iloc[0]),
                           n_labelled_kept=GATE_N.get(_g_, (None, None))[0],
                           n_unlabelled_kept=GATE_N.get(_g_, (None, None))[1])
                      for _g_ in gates for _f_, _c_, _l_, _m_ in SERIES
                      if len(G[(G.feature == _f_) & (G.gate == _g_)])],
         note="reactivity effect at each AlphaFold-confidence rung, with lysines retained")
for _g_ in gates:
    if _g_ in GATE_N:
        HEAD["d_n_labelled_gate%d" % _g_] = int(GATE_N[_g_][0])
        HEAD["d_n_unlabelled_gate%d" % _g_] = int(GATE_N[_g_][1])
        # The artwork no longer prints the pair. It prints the TOTAL kept under each bar,
        # rounded to thousands, and shows the labelled share as the bar's dark segment, so
        # the rounded total is the form that has to be checkable against the caption.
        HEAD["d_kept_gate%d_as_printed" % _g_] = "%dk" % round(sum(GATE_N[_g_]) / 1000.0)
        HEAD["d_pct_kept_gate%d" % _g_] = round(
            100.0 * sum(GATE_N[_g_]) / float(sum(GATE_N[gates[0]])))
ax.set_xticks(xs)
ax.set_xticklabels(["none" if g == 0 else "\u2265%d" % g for g in gates], fontsize=st.sz(6.3))
ax.set_xlabel("AlphaFold confidence gate", fontsize=st.sz(7))
# Anchored a fixed 0.56 in below the axes BOTTOM, not at -0.30 of the axes height. In axes
# fraction this sublabel tracked the panel: at 608.8 px it sat 183 px down, clear of the
# xlabel, and when the panel came to the family's 483.6 px the same fraction became 145 px
# and landed on "AlphaFold confidence gate" at 21%. The xlabel's own offset is in points, so
# the two only agree if this one is physical too.
# A STRIP, NOT A RATIO. The second tick row read "30k/214k", "28k/161k", "25k/141k",
# "17k/80k", which asks the reader to divide two numbers in their head four times over to see
# what the gate costs. One bar per rung instead, drawn to scale against the ungated total, so
# the 60% of lysines the tightest gate throws away is a LENGTH. The totals are printed under
# the bars because a bar with no number cannot be read off the page.
# NEUTRALS ONLY, and the same two panel a uses. Teal is the charge feature and rust the shape
# feature in every panel of this figure, so a teal bar here would give one colour two meanings
# inside panel d, which is the defect panel a was rebuilt to remove. Slate is the labelled
# minority and light grey the unlabelled majority, exactly as panel a draws the same two
# populations, so the strip and the nested bars read as one idea.
_TOT0 = float(sum(GATE_N[gates[0]]))
# THE X LIMIT IS SET HERE, not left to the line further down. The bars are positioned by
# converting each tick's data x into an axes fraction, and the first build read get_xlim()
# while the axis was still on matplotlib's autoscale range: the bars came out on a 405 px
# pitch against the ticks' 270 and the last one sat 380 px right of the rung it belongs to.
# A conversion from data to axes fraction is only valid after the limits it divides by exist.
# THE AXIS STOPS AT THE LAST RUNG. It ran to 4.62 against a last tick at 3.00, so 1.62 of
# its 4.97 units, a third of the spine, was empty rail to the right of >= 90. The width it
# was holding was for the four direct labels, which sit at 3.18; they are outside the frame
# now, beside the panel rather than inside it, which is where a series label belongs when the
# series has stopped.
# AND THE FRAME COMES IN WITH IT. Shortening the axis alone would have spread the four rungs
# from 185 px apart to 267, so the panel is narrowed by the same ratio: 3.45 units at the
# same 185.11 px per unit is 638.6 px, against 920.0, and the rungs do not move relative to
# each other. _axes_dw_in holds x0, so the panel shortens to the right.
ax.set_xlim(-0.35, 3.00)      # ON the last rung, not 0.10 past it
# The letter goes on the twin's own top line, not on the frame. d is the only panel in the
# family carrying apparatus above its frame, and with the letter on the frame D printed
# beside the four counts, in the middle of the header, with the line that says what they are
# running on above it. The opt-in is per panel: it lifts D off row 1's baseline and leaves E
# and F where they are.
ax._panel_head_twin_top = True
ax._axes_dw_in = (3.35 / 4.97 - 1.0) * _PS * 1.02 / 2.94
# Written as a fraction of its own column, not as 638.6 px against 920.0, so the rung pitch
# stays whatever the full column would have given at the original 4.97-unit span. At _PS
# 8.5003 the two forms are the same number; past that the literal one would have let the
# rungs spread while every other panel kept its scale.
_X0D, _X1D = ax.get_xlim()
_BARW = 0.80 / (_X1D - _X0D)                  # axes fraction: 80% of one tick spacing
_BARH = 0.034
# THE STRIP IS A TOP AXIS, not an annotation hanging under the panel. It reads as a second
# scale on the same four rungs, so it gets the furniture of one: its own tick labels and its
# own axis title, above the frame, with the bars sitting between the labels and the spine.
# Hung underneath it was four strata deep below the x axis and pushed the bottom margin to
# 1.000 in; on top it shares the row gap the D letter already occupies.
_tr = ax.transAxes + ScaledTranslation(0, 0.055, fig.dpi_scale_trans)
for _g_, _x_ in zip(gates, xs):
    if _g_ not in GATE_N:
        continue
    _nl, _nu = GATE_N[_g_]
    _xc = (_x_ - _X0D) / (_X1D - _X0D)
    _w = _BARW * (_nl + _nu) / _TOT0
    _wl = _w * _nl / float(_nl + _nu)
    ax.add_patch(Rectangle((_xc - _w / 2.0, 1.0), _wl, _BARH, transform=_tr,
                           facecolor=st.GREY, edgecolor="none", clip_on=False, zorder=4))
    ax.add_patch(Rectangle((_xc - _w / 2.0 + _wl, 1.0), _w - _wl, _BARH, transform=_tr,
                           facecolor=st.LIGHTGREY, edgecolor="none", clip_on=False, zorder=4))
_axt = ax.twiny()
_axt.set_xlim(_X0D, _X1D)
_axt.set_xticks(list(xs))
_axt.set_xticklabels(["%dk" % round(sum(GATE_N[_g_]) / 1000.0) if _g_ in GATE_N else ""
                      for _g_ in gates], fontsize=st.sz(6.3), color=st.GREY)
# pad clears the bars, which occupy 0.055 to 0.107 in above the spine. 10 pt is 0.139 in.
_axt.tick_params(axis="x", length=0, pad=10, colors=st.GREY)
_axt.set_xlabel("lysines kept at each gate, labeled part in dark",
                fontsize=st.sz(7), color=st.GREY)
for _sp_ in ("top", "right", "bottom", "left"):
    _axt.spines[_sp_].set_visible(False)
_axt.set_yticks([])
print("  d strip, lysines kept against the ungated %s: " % f"{int(_TOT0):,}"
      + ", ".join("%s %d%% (%s labelled)"
                  % ("none" if _g_ == 0 else ">=%d" % _g_,
                     round(100.0 * sum(GATE_N[_g_]) / _TOT0),
                     "%d%%" % round(100.0 * GATE_N[_g_][0] / sum(GATE_N[_g_])))
                  for _g_ in gates if _g_ in GATE_N))
# The key line that was here is the top axis's title now, set on _axt above.
# "effect on reactivity, d", not "effect on being labeled". Panel a defines REACTIVITY as
# the labelled-against-unlabelled contrast, so the panel that tracks that contrast across
# the gate ladder should use the word the figure already defines. It also stops the axis
# naming the assay step, the probe label, where what is plotted is the property.
ax.set_ylabel("effect on reactivity, d", fontsize=st.sz(7))
# the x limit is set above, with the strip that is positioned against it
# ON b AND c's SCALE, AND ON A TICK AT BOTH ENDS. The axis ran -0.56 to 0.52 with ticks every
# 0.2, so end_on_tick raised the top to 0.6 and the BOTTOM sat 0.16 below the lowest tick,
# which the loose-axis report now catches. The series run -0.464 to +0.401, so -0.5 to 0.5
# holds them, and that is the same range and the same 0.25 step b and c use. All three panels
# plot Cohen's d, so a reader comparing an effect in d against one in b is now reading one
# scale instead of three.
ax.set_ylim(-0.5, 0.5); ax.set_yticks([-0.5, -0.25, 0.0, 0.25, 0.5])
for f, c, ls, mk in SERIES:                      # direct labels, no legend over the data
    v = G[(G.feature == f) & (G.gate == gates[-1])].d.iloc[0]
    ax.text(len(gates) - 0.82, v, "  " + NICE(f), color=c, fontsize=st.sz(6.3), va="center", ha="left")
for s in ("top", "right"): ax.spines[s].set_visible(False)

# ---------------------------------------------------------------- f: the test
# Now that both contrasts run on the same gate ladder, the separation can be shown as a
# dose-response rather than a single before/after: charge stays coupled at every gate while
# shape comes apart progressively, and the exact permutation p falls with it.
# Panel e's width as a fraction of its own column, one number for both halves. 3.45/4.05 is
# the cut that took the rail off the right of the ladder without letting its four rungs spread;
# 0.95 of that is a further condense, asked for by eye, and it moves the rungs.
_E_W = 0.95 * 3.45 / 4.05
sub_gs = gs[1, 1].subgridspec(2, 1, height_ratios=[1.35, 1.0], hspace=0.62)

ax = fig.add_subplot(sub_gs[0, 0]); st.head(fig, ax, "e")
# ROW 1 HAD NO GAPS. Measured panel to panel, including every label each one hangs
# outside its frame, row 0 runs 87.9 px and 100.9 px while row 1 ran 11.4 px and -18.3,
# the second of those an overlap of e's series labels with f's y-axis title that misses
# the collision gate only because they sit at different heights. wspace is one number
# for the whole gridspec and row 0 is already even, so the two row-1 panels move
# individually instead. Converged against the render, and e re-converged at 37 after its
# y tick labels went from "1.5" to "0.25" and took its left ink 19 px further out.
ax._axes_dx_in = 74.8 / 311.96
# Flush, like the other five. It was the only panel here that was not: the y-axis title
# ran two lines on a 212.1 px axis, so it was taller than the panel and overran its top,
# and moving the letter onto that ink put E on the title's first line at 23%. The title
# is one line now and sits inside the panel, so the letter can go where it belongs.
ax._panel_head_flush = True
xs = np.arange(len(LADDER))
for k, c, mk in (("charge", st.BLUE, "o"), ("shape", st.ORANGE, "s")):
    # clip_on=False: charge sits at 0.98 to 0.99 against a top rung of 1.00, and 0.01 of this
    # axis is 2.8 px against a 4.5 pt marker's 9 px radius, so the axes box was cutting the
    # top off every charge mark. The line itself never leaves the frame.
    ax.plot(xs, LADDER[k], marker=mk, color=c, lw=1.3, ms=4.5, mew=0, clip_on=False)
    ax.text(xs[-1] + 0.12, LADDER[k].iloc[-1], "  " + k, color=c, fontsize=st.sz(6.3), va="center")
# ONE COLOUR AND ONE WEIGHT FOR ALL FOUR. The row used to print p >= 0.05 grey and light and
# p < 0.05 orange and bold, which drew the 0.05 line for the reader instead of letting the
# four values do it. The values are printed; a reader reads 0.412 and 0.003 and needs no help.
for x, pv in zip(xs, LADDER.p):
    ax.text(x, 1.055, "p %s" % (("%.3f" % pv).lstrip("0") if pv >= 0.001 else "<.001"),
            ha="center", fontsize=st.sz(6.3), color=st.ORANGE, fontweight="bold")
ax.set_xticks(xs)
ax.set_xticklabels(["none" if g == 0 else "≥%d" % g for g in LADDER.gate], fontsize=st.sz(6.3))
# ON A TICK AT BOTH ENDS. It ran 0.20 to 1.10 with ticks at 0.5, 1.0 and 1.5, so end_on_tick
# lifted the top to 1.5 and left the BOTTOM at 0.20, between nothing. 0.25 to 1.00 on a 0.25
# step ends on a tick at both ends and holds the series, which run 0.35 to 0.99. The p row is
# at 1.055, outside the frame, where it reads as the annotation it is rather than as data.
ax.set_ylim(0.25, 1.00); ax.set_yticks([0.25, 0.50, 0.75, 1.00])
# THE RAIL PAST >= 90 IS GONE, the same cut panel d took. The axis ran to 3.75 against a last
# rung at 3, and the 0.75 units it was holding were for the two direct labels, which sit at
# 3.12; they are outside the frame now, beside the panel rather than inside it. AND THE FRAME
# COMES IN WITH IT: shortening the axis alone would have spread the four rungs by 17%, so both
# halves of e are narrowed by the same 3.45/4.05, which keeps the rung pitch and keeps the
# histogram below exactly as wide as the ladder above it.
ax.set_xlim(-0.35, 3.10)
ax._axes_dw_in = (_E_W - 1.0) * _PS / 2.94
ax.set_ylabel("agreement", fontsize=st.sz(7))
ax.tick_params(axis="y", labelsize=st.sz(6.3))
for sp in ("top", "right"): ax.spines[sp].set_visible(False)

ax = fig.add_subplot(sub_gs[1, 0])
ax._axes_dx_in = 74.8 / 311.96        # the other half of e moves with it
ax._axes_dw_in = (_E_W - 1.0) * _PS / 2.94             # and narrows with it
obs, g, p, _, _ = NULL["gated"]
ax.hist(g, bins=58, color=st.LIGHTGREY, edgecolor="none")
ax.axvline(obs, color=st.ORANGE, lw=1.5, zorder=3)
ax.axvline(NULL["ungated"][0], color=st.GREY, lw=1.2, ls=(0, (2, 1.5)), zorder=3)
# ma="center" so the two lines of each marker sit centred on one another. They were
# right-aligned, which left the shorter second line ragged down its left side. The value goes
# in parentheses as a value under its name rather than a second bare number.
# ONE TOP FOR BOTH. They were anchored at 0.62 and 0.90 of the axis with the default vertical
# alignment, so the two markers started a line apart and read as a stagger rather than as two
# labels on one row. Same y and va="top" puts their first lines on one baseline. "no gate"
# takes the deeper horizontal offset because its rule sits at 0.01, in the body of the
# histogram, where 7 pt left it against the bars.
# 0.98, not 0.72. With va="top" the block hangs DOWN from its anchor, so 0.72 printed
# both markers across the body of the histogram: each one reaches left from its own
# rule into bars 20 to 30% of full height. Hung from 0.98 the pair clears them.
_E_LBL_Y = ax.get_ylim()[1] * 0.98
ax.annotate("no gate\n(%.2f)" % NULL["ungated"][0], (NULL["ungated"][0], _E_LBL_Y),
            xytext=(-14, 0), textcoords="offset points", ha="right", va="top", ma="center",
            fontsize=st.sz(6.3), color=st.GREY, linespacing=1.5)
# RIGHT OF ITS OWN RULE. Left of it the two lines reached back across the body of the
# histogram, over the bars the rule is being compared against; right of it there is nothing
# but the tail, and the rule reads as the edge of the label rather than as something behind it.
ax.annotate("pLDDT ≥ 90\n(%.2f)" % obs, (obs, _E_LBL_Y),
            xytext=(4, 0), textcoords="offset points", ha="left", va="top", ma="center",
            fontsize=st.sz(6.3), color=st.ORANGE, linespacing=1.5)
# -1.0 to 1.0 ON A 0.5 STEP, ends on a tick at both ends. It clips 32 of the 5,005
# permutation values, which run to -1.5709, and that is the LEFT tail. The test reads off the
# right: the observed statistic is +0.639 and the p-value counts permutations at or above it,
# where the null tops out at 0.9153, well inside the frame. Nothing the panel is used for is
# cut off, and the bars being lost are one and two counts tall.
ax.set_xlim(-1.0, 1.0); ax.set_xticks([-1.0, -0.5, 0.0, 0.5, 1.0]); ax.set_yticks([])
ax.tick_params(axis="x", labelsize=st.sz(6.3))
for sp in ("top", "right", "left"): ax.spines[sp].set_visible(False)
ax.set_xlabel("%s regroupings of the same features" % f"{len(g):,}", fontsize=st.sz(7))

fd.panel("fig4", "e", LADDER.rename(columns={
    "charge": "charge_coupling_r", "shape": "shape_coupling_r", "p": "exact_permutation_p"}),
    note="how closely reactivity and ligandability track, by gate; p over all %d regroupings"
         % len(SPLITS))
for _r in LADDER.itertuples():
    _k = "none" if _r.gate == 0 else "%d" % _r.gate
    HEAD["e_charge_coupling_gate%s" % _k] = round(float(_r.charge), 2)
    HEAD["e_shape_coupling_gate%s" % _k] = round(float(_r.shape), 3)
    HEAD["e_p_gate%s" % _k] = round(float(_r.p), 3)
for _k in ("ungated", "gated"):
    HEAD["e_perm_obs_%s" % _k] = round(float(NULL[_k][0]), 3)
    HEAD["e_perm_p_%s" % _k] = round(float(NULL[_k][2]), 3)
    HEAD["e_corr_charge_%s" % _k] = round(float(NULL[_k][3]), 2)
    HEAD["e_corr_shape_%s" % _k] = round(float(NULL[_k][4]), 2)

# ---------------------------------------------------------------- f: occupancy dose-response
ax = fig.add_subplot(gs[1, 2]); st.head(fig, ax, "f")
ax._axes_dx_in = 161.0 / 311.96       # opens e|f as well as carrying e's own move
ax._panel_head_flush = True
rows, hdr = [], None
for ln in R38[R38.find("Q1a"):R38.find("Q1b")].splitlines():
    p_ = ln.split()
    if len(p_) > 3 and p_[0] == "tier": hdr = p_
    elif hdr and len(p_) == len(hdr): rows.append(p_)
T = pd.DataFrame(rows, columns=hdr)
T = T[T.tier.str.match(r"^[\d>=\u2265\-]+$", na=False)].reset_index(drop=True)
for c in T.columns:
    if c != "tier": T[c] = pd.to_numeric(T[c], errors="coerce")
T["tier"] = T.tier.str.replace(">=", "\u2265", regex=False)
xs = np.arange(len(T))
# All four of these are SHAPE features (see SHAPE above, the panel b legend and panel d), so
# all four carry the shape colour. The old version gave crowding the CHARGE colour and put
# aromatic on a third hue, which contradicted the class code the rest of the figure uses.
# Identity comes from the marker and the direct label; exposure keeps its dashed rule because
# it is the one that runs the other way.
SER_F = (("hydrophob", st.ORANGE, "-",  "o", "hydrophobic 8Å"),
         ("crowding_", st.ORANGE, "-",  "s", "crowding 8Å"),
         ("aromatic_", st.ORANGE, "-",  "^", "aromatic 8Å"),
         ("rsa_tien",  st.ORANGE, "--", "D", "exposure"))
for f, c, ls, mk, lab in SER_F:
    if f in T.columns:
        ax.plot(xs, T[f], ls=ls, marker=mk, color=c, lw=1.25, ms=4.5, mew=0,
                clip_on=False)                      # d and e's marker size
ax.axhline(0, color=st.LIGHTGREY, lw=0.8, ls=(0, (1.2, 2.0)))   # dotted, as in b, c and d
# ON A TICK AT BOTH ENDS. The axis was left to matplotlib and came out -0.668 to 1.0 with
# ticks every 0.5, so the bottom sat 0.168 below the lowest tick. The series run -0.598
# to +0.807 and the direct labels reach 0.836, so -0.75 to 1.0 holds all of it on the
# same 0.25 step b, c and d use.
ax.set_ylim(-1.0, 1.0)
ax.set_yticks([-1.0, -0.5, 0.0, 0.5, 1.0])
# ANGLED, not upright, which is what panels b and c do with their feature names. Same 45
# degrees and the same anchor, so the three panels carrying a long category name read alike.
ax.set_xticks(xs)
ax.set_xticklabels(T.tier, fontsize=st.sz(6.3), rotation=45, ha="right", rotation_mode="anchor")
ax.set_xlabel("occupancy of the lysine (%)", fontsize=st.sz(7))
ax.set_ylabel("shift from the 0–10% tier, d", fontsize=st.sz(7))
# DIRECT LABELS CANNOT ALL SIT EXACTLY ON THEIR LINE. hydrophobic ends at +0.807 and
# crowding at +0.740, 0.067 apart, and one label's glyph box measures 0.111 in these data
# units. They overlapped by 24% at the old panel height and by 40% once the panel came to
# the family's 483.6 px. The ends are pushed apart to 0.125, symmetrically about their own
# midpoint, so each label stays as near its line as the type allows and neither series is
# privileged over the other. Identity is still the marker and the word, not the colour.
# COMPUTED, NOT SET. 0.125 was measured against the old 6.645 pt label on the old -0.668 to
# 1.0 axis. The type ladder took the label to 7.974 and the tick rule took the axis to
# -0.75 to 1.0, and the pair overlapped again at 13%. The separation is now derived from the
# drawn type and the axis it sits on, so neither change can stale it: one line box is
# fontsize/72 x 1.25 inches, over a _YAX-inch panel, times the span in data units.
_MINSEP = (st.sz(6.3) / 72.0 * 1.25) / _YAX * (ax.get_ylim()[1] - ax.get_ylim()[0]) * 1.10
_ends = sorted([[float(T[f].iloc[-1]), c, lab] for f, c, ls, mk, lab in SER_F
                if f in T.columns], key=lambda t: t[0])
for _i in range(len(_ends) - 1):
    _d = _ends[_i + 1][0] - _ends[_i][0]
    if _d < _MINSEP:
        _sh = (_MINSEP - _d) / 2.0
        _ends[_i][0] -= _sh
        _ends[_i + 1][0] += _sh
for _y, _c, _lab in _ends:                       # identity is not colour alone
    ax.text(len(T) - 0.75, _y, "  " + _lab, color=_c, fontsize=st.sz(6.3),
            va="center", ha="left")
print("  f direct labels, drawn y against series end: "
      + ", ".join("%s %+.3f (end %+.3f)" % (_l, _y, float(T[_f].iloc[-1]))
                  for (_y, _c, _l), _f in zip(_ends, [dict((l, f) for f, c, ls, mk, l in SER_F)[_l]
                                                      for _y, _c, _l in _ends])))
# THE RAIL PAST >= 90 IS GONE, the third panel in this row to take that cut. The axis ran to
# 12.4 against a last tier at 9, and the 2.4 units it was holding were for the four direct
# labels at 9.25; they are outside the frame now. The frame comes in by the same 9.5/12.8, so
# the ten tiers keep their pitch AND the labels land on the page exactly where they were,
# which is why this costs row 1 nothing: f's ink still ends where it ended.
ax.set_xlim(-0.4, len(T) - 1.0)      # ON the >= 90 tier
ax._axes_dw_in = (9.4 / 12.8 - 1.0) * _PS * 0.92 / 2.94
fd.panel("fig4", "f", T, note="shift from the 0-10%% occupancy tier, d, by tier")
HEAD["f_n_tiers"] = int(len(T))
# MOVED TO THE CAPTION. "monotonic across all ten tiers, no discontinuity at the 60% cut"
# was a conclusion printed on the artwork. It now reads, in the (F) entry of CAP_2_4 in
# manuscript/drafted.py: "Monotonic across all ten tiers, with no discontinuity at the
# 60% engagement cut." Same claim, same two facts, and the panel is left showing the
# curves the claim is read off. Figures 1, 2 and 3 moved their equivalents the same way.
for s in ("top", "right"): ax.spines[s].set_visible(False)


# NOTE: writes to results/ only. It does not overwrite figures/main/fig4 — adopting this
# version is a decision for the PI, not a side effect of running the script.
# via st.save so the blank-panel check covers this figure as well
st.save(fig, "fig4")
# Trim to the drawing, 18 px all round, the figures-final standard Figure 1 set and
# Figures 2 and 3 follow. This figure had no trim and carried the gridspec margin, which
# put its panel letters 198.8 px from the left edge where the other three sit at 18 to
# 22 px: side by side, Figure 4 was the one indented sheet. A bbox crop moves no drawn
# element relative to any other and resizes no panel, so the letters travel with their
# content and nothing inside the figure changes.
# Measured off the RENDERED PIXELS rather than get_tightbbox, because an axis-off panel
# still contributes its rectangle to the tight bbox and panel a has its axis off.
_PAD_IN = 18.0 / 300.0
from PIL import Image as _Im
_px = np.asarray(_Im.open("fig4.png").convert("RGB")).astype(int)
_ys, _xs = np.nonzero(_px.sum(2) < 735)
_H_in = fig.get_size_inches()[1]
_trim = Bbox.from_extents(_xs.min() / 300.0 - _PAD_IN,
                          _H_in - _ys.max() / 300.0 - _PAD_IN,
                          _xs.max() / 300.0 + _PAD_IN,
                          _H_in - _ys.min() / 300.0 + _PAD_IN)
for _ext in (".png", ".pdf"):
    fig.savefig("fig4" + _ext, bbox_inches=_trim)
print("trimmed to %.3f x %.3f in (%.0f x %.0f px at 300 dpi), 18 px pads"
      % (_trim.width, _trim.height, _trim.width * 300, _trim.height * 300))
fd.head("fig4", HEAD)

if __import__("os").environ.get("FIG4_FINAL"):
    # The finalized file for CCB-submission/figures-final/, on the family's raster.
    # THE FAMILY SHARES A RASTER, NOT A WIDTH. Figure 1 at 12.00 in fixes it at
    # 300 x 12.00 / 11.5400 = 311.96 dpi; every figure renders at that raster and is declared
    # 300 dpi, so all four print 9.213 pt body type whatever their width. This figure's native
    # trimmed sheet is deliberately fig1.png's 3461 px, so it finalizes at Figure 1's 12.00 in.
    # Saved WITH _trim, because the trim here is a bbox_inches on save and not a resize, so
    # fig.get_size_inches() still reports the untrimmed canvas.
    import PIL.Image as _PI, pymupdf as _mu
    _ref = _mu.open("fig1.pdf")[0].rect.width / 72.0      # Figure 1's native width, measured
    _dpi = 300.0 * 12.0 / _ref
    fig.savefig("fig4_final.png", dpi=_dpi, bbox_inches=_trim)
    _im = _PI.open("fig4_final.png")
    _im.save("fig4_final.png", dpi=(300, 300))
    print("final: native %.4f in at the shared %.2f dpi raster -> %d x %d px declared 300 dpi "
          "(%.3f x %.3f in), body type %.3f pt against Figure 1's %.3f"
          % (_trim.width, _dpi, _im.size[0], _im.size[1], _im.size[0] / 300.0,
             _im.size[1] / 300.0, 8.860 * 12.0 / _ref, 8.860 * 12.0 / _ref))
fig.savefig("../results/fig4_v3.png", dpi=400, bbox_inches="tight")
fig.savefig("../results/fig4_v3.pdf", bbox_inches="tight")
print("features used: %d (%d charge / %d shape); ligandability tier = %d%%" % (NF, KC, NF-KC, TIER))
for c in NULL: print("  %-12s obs %+.3f   p %.4f" % (c, NULL[c][0], NULL[c][2]))
print("  mean pair gap  ungated:", F.groupby("kind").gap_0.mean().round(3).to_dict())
print("  mean pair gap    gated:", F.groupby("kind").gap_90.mean().round(3).to_dict())
for _k in ("ungated", "gated"):
    _o, _g, _p = NULL[_k][0], NULL[_k][1], NULL[_k][2]
    print("  e permutation %-8s observed %+.3f, p = %.3f, over %d regroupings" % (_k, _o, _p, len(_g)))
print("  e gate ladder:", ", ".join("%s charge %.2f shape %.2f p %.3f" %
      ("none" if r.gate == 0 else ">=%d" % r.gate, r.charge, r.shape, r.p) for r in LADDER.itertuples()))
print("  panels b,c display gate: pLDDT >= %d (GATE = max(RE)); ladder rungs %s" % (GATE, sorted(RE)))
print("  features (%d): charge = %s | shape = %s"
      % (NF, ", ".join(sorted(F.feature[F.kind == "charge"])),
         ", ".join(sorted(F.feature[F.kind == "shape"]))))
print("  per-feature, reactivity then ligandability, ungated -> gated:")
for _, r in F.iterrows():
    print("     %-16s %-6s react %+0.3f -> %+0.3f   lig %+0.3f -> %+0.3f   gap %.3f -> %.3f"
          % (r.feature, r.kind, r.react, r.react_g, r.lig, r.lig_g, r.gap_0, r.gap_90))
print("  within-kind correlation react vs lig: ungated charge %.2f shape %.2f | gated charge %.2f shape %.2f"
      % (NULL["ungated"][3], NULL["ungated"][4], NULL["gated"][3], NULL["gated"][4]))