#!/usr/bin/env python3
"""
ga_v1.py, the GRAPHICAL ABSTRACT for the Cell Chemical Biology Resource.

Cell Press requires one for this journal and specifies a SINGLE-PANEL, SQUARE image
(cell.com/information-for-authors/figure-guidelines). It is not a figure of the family:
it has no panel letters, no caption to lean on and no axis furniture, it is displayed at
about a third of a page, and it has to carry the paper's one claim on its own.

THE CLAIM, in the abstract's own words: "Acquired lysines are structurally ordinary, so the
barrier is detection rather than chemistry."

Three acts, read top to bottom, then the verdict:
  1  a missense mutation makes a new lysine, and the class is large and sits on drivers
  2  those lysines sit in the same structural space as the ones compounds already engage
  3  the same mutation makes a peptide too short to identify, so they are not seen,
     and the error rate on the ones that are seen is not the one the search reports

EVERY NUMBER IS READ, not typed: each comes from the deposited ledger of the figure that
computed it, so the abstract, the figures and this image cannot drift apart.

THE DETECTION CENSUS IS GIVEN AT BOTH THRESHOLDS. The Figure 3 caption binds the whole
paper to it, that every error rate is quoted at margin ≥ 10 and that Figure 2's census is
given at both thresholds wherever it appears. This square is a place it appears, so 0.63%
at the search's own 1% filter and 0.35% at margin ≥ 10 are printed together.

INPUT   ../results/panel_data/fig1_head.json, fig2_head.json, fig3_head.json,
        figS3_head.json, fig8_head.json
        ../ligand_rebuild/step28d_curve.csv
OUTPUT  graphical_abstract.png, .pdf, and a copy in ../results/
"""
import os, sys, json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import abbasov_style as st

BODY = 7.0
st.apply(BODY)
plt.rcParams["savefig.dpi"] = 300          # AFTER apply(), which resets the rcParams

ORANGE, BLUE, LIGHTGREY = st.ORANGE, st.BLUE, st.LIGHTGREY
GREY = st.GREEN                            # charcoal, the second series colour
INK  = st.GREY                             # slate, every label

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
# The working tree writes ledgers to results/panel_data; the published deposit ships them at
# the root as panel_data/. Take whichever exists, so this runs from a clone unchanged.
LED  = next(d for d in (os.path.join(ROOT, "results", "panel_data"),
                        os.path.join(ROOT, "panel_data"))
            if os.path.isdir(d) and os.path.exists(os.path.join(d, "fig1_head.json")))
def led(fig):
    d = json.load(open(os.path.join(LED, "%s_head.json" % fig)))
    return d.get("numbers", d)
L1, L2, L3, LS3, L8 = led("fig1"), led("fig2"), led("fig3"), led("figS3"), led("fig8")

N_SITES  = int(L1["a_union_K"])            # 295,411 acquired lysines
N_GENES  = int(L1["a_genes"])              # in 20,365 genes
# d_OR, NOT the ratio of the two percentages. They differ: 1.675 against 1.668 here, and
# the paper names an odds ratio, so the square has to print the statistic the paper names
# rather than one of its own that happens to round the same way.
OR_DRIV  = float(L1["d_OR"])
N_OBS    = int(L2["b_pep_total"])          # 8,470 site-in-cell-line observations
PCT_REC  = float(L2["c_total_pct"])        # 0.63% recover the mutant peptide
PCT_REC_M = float(L2["c_total_pct_margin"])  # 0.35% of them at margin ≥ 10
FDP_RAW  = float(L3["b_fdp_margin0_pct"])  # 4.70% realized, before any score filtering
# The search's own filter SETTING, not a measurement, so no ledger holds it. It is named
# here rather than typed into a string, so the two places it is printed cannot disagree.
FDR_REPORTED = 1
LEN_WT   = int(LS3["f_wt_median_aa"])      # 13 residues
LEN_MUT  = int(LS3["f_mut_median_aa"])     # 6
P_CLOSER = float(L8["c_p_closer"])         # 0.443

# THREE TYPE SIZES, as every figure in the set has, but larger: a graphical abstract prints
# at about a third of a page, so the family's 7.974 pt body would land near 5 pt on the page.
SM, MD, LG = st.sz(8.0), st.sz(9.5), st.sz(12.0)

# SQUARE, because Cell Press asks for a square single panel. The axes is the whole canvas and
# the drawing is in 0..100 units both ways, so every position below is a percentage of the
# side and nothing has to be re-derived if the square changes size.
SIDE = 4.0
fig = plt.figure(figsize=(SIDE, SIDE))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.set_axis_off()

# THREE BANDS OF EQUAL HEIGHT. The acts hold different amounts of ink, so left to
# themselves they take different amounts of the square and the sheet reads as three
# unequal things. Each act is given the SAME band and anchored the same way inside it:
# its title sits _LABEL_DY below the band top and its last line _CLOSE_DY above the band
# bottom, so the three read as one rhythm whatever sits between. The verdict is not an
# act and gets the strip below them.
_VBAND = 7.0                                   # the verdict's strip, 0 to here
_BAND  = (100.0 - _VBAND) / 3.0                # 31.00 units each, by construction
_T = [100.0 - k * _BAND for k in (0, 1, 2)]    # each band's top: 100.00, 69.00, 38.00
_LABEL_DY = 5.0                                # title, below its band top
_CLOSE_DY = 2.8                                # last line, above its band bottom
_PITCH    = 4.6                                # between two closing lines

# The thread down the left is the one object the three acts share: it is the lysine, and it
# changes colour where the story does. It runs the whole height of the three bands.
ax.plot([6.0, 6.0], [_VBAND, _T[0] - _LABEL_DY + 2.2], color=LIGHTGREY, lw=1.2, zorder=1)

def act(y, n, text):
    """A numbered stop on the thread."""
    ax.add_patch(Circle((6.0, y), 2.2, facecolor=ORANGE, lw=0, zorder=3))
    # va="center" centres the LAYOUT box, which runs from ascent to descent, so a digit,
    # having no descender, is set 4 px high in the disc at this size. Measured off the
    # render, not assumed. va="center_baseline" puts the ink bbox dead on the disc centre.
    ax.text(6.0, y, "%d" % n, color=st.PAPER, fontsize=SM, fontweight="bold",
            ha="center", va="center_baseline", zorder=4)
    ax.text(12.0, y, text, color=INK, fontsize=MD, fontweight="bold",
            ha="left", va="center", zorder=3)

# ---------------------------------------------------------------- 1  the class
act(_T[0] - _LABEL_DY, 1, "a missense mutation makes a new lysine")
_X0, _X1, _CY = 13.0, 87.0, _T[0] - 16.0
_N = 11
_xs = np.linspace(_X0, _X1, _N)
for i, x in enumerate(_xs[:-1]):
    ax.plot([x, _xs[i + 1]], [_CY, _CY], color=LIGHTGREY, lw=1.6, zorder=2)
for i, x in enumerate(_xs):
    _is_k = (i == 6)
    ax.add_patch(Circle((x, _CY), 2.6 if _is_k else 1.9,
                        facecolor=ORANGE if _is_k else LIGHTGREY, lw=0, zorder=3))
ax.text(_xs[6], _CY, "K", color=st.PAPER, fontsize=SM, fontweight="bold",
        ha="center", va="center", zorder=4)
# The arrow only. The word "missense" was printed above it and ran into act 1's own line,
# which says the same thing.
ax.annotate("", xy=(_xs[6], _CY + 3.4), xytext=(_xs[6], _CY + 7.6),
            arrowprops=dict(arrowstyle="->", color=ORANGE, lw=0.9))
# CENTRED ON THE CONTENT COLUMN, 54, not on the page. The thread sits at x 6 and every
# act's own line starts at 12, so a line centred at 50 hangs its left end out past the
# column the rest of the sheet is built on.
_CX = 54.0
ax.text(_CX, _T[0] - _BAND + _CLOSE_DY + _PITCH, "%s sites in %s genes" % (f"{N_SITES:,}", f"{N_GENES:,}"),
        color=INK, fontsize=SM, ha="center", va="center")
ax.text(_CX, _T[0] - _BAND + _CLOSE_DY, "%.1f× more often on driver hotspots than other substitutions"
        % OR_DRIV, color=INK, fontsize=SM, ha="center", va="center")

# ---------------------------------------------------------------- 2  ordinary
act(_T[1] - _LABEL_DY, 2, "it sits in the same environments as engaged lysines")
# The two cumulative curves of Figure 8c, drawn without axes: the point is that they lie on
# top of each other, and a frame would invite the eye to read values instead of overlap.
CV = pd.read_csv(os.path.join(ROOT, "ligand_rebuild", "step28d_curve.csv"))
_cx0, _cx1, _cy0, _cy1 = 16.0, 62.0, _T[1] - 23.5, _T[1] - 11.0
_r = CV.radius.values / CV.radius.max()
for _col, _c, _lw in (("engaged", BLUE, 1.8), ("union", ORANGE, 1.8)):
    ax.plot(_cx0 + _r * (_cx1 - _cx0), _cy0 + CV[_col].values * (_cy1 - _cy0),
            color=_c, lw=_lw, zorder=3, solid_capstyle="round")
ax.plot([_cx0, _cx1], [_cy0, _cy0], color=LIGHTGREY, lw=1.0, zorder=2)
ax.text(_cx1 + 3.0, _T[1] - 13.5, "engaged", color=BLUE, fontsize=SM, ha="left", va="center")
ax.text(_cx1 + 3.0, _T[1] - 19.0, "acquired", color=ORANGE, fontsize=SM, ha="left", va="center")
# ONE line under the curve, not two. P alone means nothing to someone browsing, and the
# number plus what 0.5 would have meant is the whole of it.
ax.text(_CX, _T[1] - _BAND + _CLOSE_DY, "P = %.3f, where 0.5 would mean the two are interchangeable"
        % P_CLOSER, color=INK, fontsize=SM, ha="center", va="center")

# ---------------------------------------------------------------- 3  unseen
act(_T[2] - _LABEL_DY, 3, "and it cuts the peptide too short to identify")
_bx0, _bw, _by = 24.0, 46.0, 0.0
_scale = _bw / float(max(LEN_WT, LEN_MUT))
# The two counts sit at ONE x, past the longer bar, so the shorter bar's number cannot land
# on the search floor it is being compared against.
_NX = _bx0 + _bw + 2.5
for _i, (_nm, _n, _c) in enumerate((("wild type", LEN_WT, GREY),
                                    ("mutant", LEN_MUT, ORANGE))):
    _y = _T[2] - 10.4 - _i * 4.8
    ax.add_patch(Rectangle((_bx0, _y - 1.9), _n * _scale, 3.8, facecolor=_c, lw=0, zorder=3))
    ax.text(_bx0 - 2.0, _y, _nm, color=_c, fontsize=SM, ha="right", va="center")
    ax.text(_NX, _y, "%d aa" % _n, color=_c, fontsize=SM, ha="left", va="center")
# The search floor is furniture: light grey, dotted, behind the bars, exactly at 7 residues,
# which is where the mutant peptide falls out of the search.
_fx = _bx0 + 7.0 * _scale
ax.plot([_fx, _fx], [_T[2] - 17.1, _T[2] - 8.0], color=LIGHTGREY, lw=1.0, ls=(0, (1.2, 2.0)), zorder=2)
ax.text(_fx, _T[2] - 18.0, "7 aa, the search floor", color=INK, fontsize=SM, ha="center", va="top")
# TWO closing lines, because the detection claim has two halves and the square carried one.
# The first is the census at both thresholds. The second is the half that was missing: the
# peptides that ARE recovered do not carry the error rate the search reports for them.
ax.text(_CX, _T[2] - _BAND + _CLOSE_DY + _PITCH, "recovered in %.2f%% of %s observations, %.2f%% at margin %s 10"
        % (PCT_REC, f"{N_OBS:,}", PCT_REC_M, chr(8805)),
        color=INK, fontsize=SM, ha="center", va="center")
# "where the search reports 1%" is not decoration, it PINS THE THRESHOLD. 4.70% is the
# proportion before any score filtering, which is the search's own 1% output, and it is
# 0.83% at margin >= 10. Stating 4.70% flat, one line under "0.35% at margin >= 10", would
# read as if the two sat at the same cut. They do not.
ax.text(_CX, _T[2] - _BAND + _CLOSE_DY, "entrapment measures %.2f%% error where the search reports %d%%"
        % (FDP_RAW, FDR_REPORTED), color=INK, fontsize=SM, ha="center", va="center")

# ---------------------------------------------------------------- the verdict
ax.text(_CX, 3.4, "the barrier is detection, not chemistry",
        color=ORANGE, fontsize=LG, fontweight="bold", ha="center", va="center")

# NO TRIM. The square IS the deliverable, so the canvas is saved whole; trimming to the ink
# would make it a rectangle and Cell Press asks for a square.
for _e in (".png", ".pdf"):
    fig.savefig("graphical_abstract" + _e)
import shutil
for _e in ("png", "pdf"):
    shutil.copy("graphical_abstract.%s" % _e,
                os.path.join(ROOT, "results", "graphical_abstract.%s" % _e))
# THE COLLISION GATE, run here rather than after the fact. This sheet has no panels for
# figtools/headrow.py to group, but the same check applies: no two words may overlap.
_coll = st.text_collisions(fig, 0.12)
print("text collisions at 12%%: %d%s"
      % (len(_coll), "" if not _coll else "  " + "; ".join(
          "%r over %r" % (a, b) for a, b, *_ in _coll[:4])))
from PIL import Image as _Im
_im = _Im.open("graphical_abstract.png")
# THE BANDS ARE MEASURED, not asserted. Each act's ink is found in the rendered square and
# its extent reported, so "the same dimensions" is a reading off the page.
import numpy as _np
_a = _np.asarray(_im.convert("L")).astype(int)
_ink = (_a < 245)
_ink[:, : int(0.085 * _im.size[0])] = False     # the thread runs the full height, so drop it
print("  bands %.2f units each, tops %s" % (_BAND, ", ".join("%.2f" % t for t in _T)))
for _k, _t in enumerate(_T):
    _y0 = int((1 - _t / 100.0) * _im.size[1]); _y1 = int((1 - (_t - _BAND) / 100.0) * _im.size[1])
    _rows = _np.nonzero(_ink[max(0, _y0):_y1].any(axis=1))[0]
    _h = (_rows.max() - _rows.min() + 1) / float(_im.size[1]) * 100.0
    print("    act %d  band %.2f to %.2f   ink spans %.2f units, starts %.2f below the band top"
          % (_k + 1, _t, _t - _BAND, _h, _rows.min() / float(_im.size[1]) * 100.0))
print("graphical abstract %d x %d px, %.2f x %.2f in at 300 dpi, square=%s"
      % (_im.size[0], _im.size[1], _im.size[0] / 300.0, _im.size[1] / 300.0,
         _im.size[0] == _im.size[1]))
print("  read: %s sites | %s genes | OR %.2f | %s obs | %.2f%% (%.2f%% at margin %s 10)"
      " | %d vs %d aa | P %.3f | FDP %.2f%% vs %d%% reported"
      % (f"{N_SITES:,}", f"{N_GENES:,}", OR_DRIV, f"{N_OBS:,}", PCT_REC, PCT_REC_M,
         chr(8805), LEN_WT, LEN_MUT, P_CLOSER, FDP_RAW, FDR_REPORTED))
