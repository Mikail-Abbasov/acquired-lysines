#!/usr/bin/env python3
"""
PUBLISHED FIGURE 5, the sheet the paper actually prints, which until now did not exist.

THE PAPER HAS SEVEN FIGURES AND THE REPO HAS EIGHT SHEETS. `manuscript/panels.py` sets
MERGE_INTO, MERGE_FROM = 5, 6, so source figures 5 and 6 become one published Figure 5 and
six of their twelve panels go to Extended Results. `results/panel_data/figure_numbering.json`
is the authority and this script implements it:

    published   source      what it is
    a           fig6 a      the assay, two rails
    b           fig6 b      the fourteen site features, within compound
    c           fig5 a      the central comparison, model against reproducibility ceiling
    d           fig5 c      is it weak everywhere, the ROC
    e           fig6 d      does any compound descriptor predict site preference
    f           fig5 e      breadth across chemotypes

    cut to Extended Results: 5b, 5d, 5f, 6c, 6e, 6f

THE PANELS ARE NOT REDRAWN HERE. Each one is read out of its own source script at run time
and executed against this figure's axes, so the six blocks keep every decision taken in their
own review rounds and cannot drift from the sources. What this script owns is the FRAME:
which cell each panel sits in, the column widths, the gaps, the panel letters and the trim.
A crop-and-paste composite was the other option, and `build.py` calls that a reading copy
for good reason: it inherits the old layout's proportions and loses anything drawn outside a
panel's own rectangle.

INPUT   fig5_v2.py, fig6_v2.py and everything their own heads read
OUTPUT  fig5pub.png, .pdf, and with FIG5PUB_FINAL=1 the figures-final sheet
"""
import os, re, sys, io
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
from matplotlib.patches import Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import abbasov_style as st
import figure_data as fd

BODY = 7.0
st.apply(BODY)
plt.rcParams["savefig.dpi"] = 300          # AFTER apply(), which resets the rcParams

# ---------------------------------------------------------------- read the two sources
_FIGMARK = "---------------------------------------------------------------- figure"


def source_head(path):
    """Everything the source computes before it draws: the data, the parse and the ledger."""
    src = io.open(os.path.join(HERE, path), encoding="utf-8").read()
    assert src.count(_FIGMARK) == 1, path
    ns = {"__name__": "__head__", "__file__": os.path.join(HERE, path)}
    exec(compile(src[:src.index(_FIGMARK)], path, "exec"), ns)
    return ns


def source_block(path, letter):
    """One panel's drawing code, from its own section marker to the next one."""
    lines = io.open(os.path.join(HERE, path), encoding="utf-8").read().splitlines()
    pat = re.compile(r"^# ={10,} ([a-f])\s*$|^# -{10,} ([a-f]):")
    marks = [(i, (m.group(1) or m.group(2)))
             for i, l in enumerate(lines) for m in [pat.match(l)] if m]
    for k, (i, let) in enumerate(marks):
        if let == letter:
            j = marks[k + 1][0] if k + 1 < len(marks) else len(lines)
            return lines[i:j]
    raise KeyError("%s has no panel %s" % (path, letter))


NS5, NS6 = source_head("fig5_v2.py"), source_head("fig6_v2.py")

# ---------------------------------------------------------------- figure
# THE FRAME IS THIS SCRIPT'S, the panels are not. Margins in inches and the canvas solved
# from the panel block, as Figures 1 to 6 are built. _PS is the sum of the three panel widths
# and is HELD, so changing wspace closes the gaps instead of widening all six panels.
_ML, _MR = 1.60, 1.20          # left and right margins, inches
_MB, _MT = 1.60, 0.80          # bottom and top margins, inches
_PS = 8.475                    # sum of the three panel widths, held. Chosen so the trimmed
                               # sheet is 11.540 in native, which is Figure 1's width on the shared
                               # raster, so this sheet finalizes at 12.00 in like Figures 1 and 5
_WS, _HS = 0.42, 0.72          # column and row spacing, fractions of the average panel
_YAX = 1.55020                 # the family y axis, 483.6 px at the 311.96 raster
# COLUMN 0 CARRIES THE WIDE PANEL OF ROW 0 and the narrow one of row 1, which is what
# _axes_dw_in exists for: the column is sized for a, and d shortens inside it.
_RATIOS = [1.55, 1.10, 1.35]
_PW = _PS * (1 + 2 * _WS / 3)
_PH = (2 + _HS) * _YAX
_FW, _FH = _PW + _ML + _MR, _PH + _MB + _MT
fig = plt.figure(figsize=(_FW, _FH))
_SNAP = []                     # (text, patch) pairs snapped after st.save fits the bars
_GS = fig.add_gridspec(2, 3, hspace=_HS, wspace=_WS,
                       left=_ML / _FW, right=1 - _MR / _FW,
                       top=1 - _MT / _FH, bottom=_MB / _FH,
                       width_ratios=_RATIOS)

# published letter -> (source script, source letter, row, col, dx px, dw px)
# dx and dw are THIS layout's, measured off this render. The ported blocks carry the source
# layout's values and those lines are stripped, because a number tuned against a different
# set of columns is worse than no number at all.
PLAN = [("a", "fig6_v2.py", "a", 0, 0, -146.0,    0.0, -59.3),
        ("b", "fig6_v2.py", "b", 0, 1, -113.0,  105.0,   0.0),
        ("c", "fig5_v2.py", "a", 0, 2,  120.0, -120.0,   0.0),
        ("d", "fig5_v2.py", "c", 1, 0,    0.0, -260.0,   0.0),
        ("e", "fig6_v2.py", "d", 1, 1,  -78.0,    0.0,   0.0),
        ("f", "fig5_v2.py", "e", 1, 2,   24.0,    0.0,   0.0)]
# dw: b takes the 15% of its column its fourteen angled names need, as it does in fig6, and
# without it the tick gate reports 14 labels overlapping by 7%. c gives back 120 px so row 0's
# gaps open from 137 to 197. d is the ROC and shortens by 260 inside a column sized for a.
# dx: solved, not estimated. figtools/solve5.py renders, reads the ink extents panelgaps
# measures and moves each panel until both rows span the sheet with equal gaps, which lands
# at 197, 197 in row 0 and 314, 313 in row 1.
# ldx moves the LETTER only. Panel a has its axis off, so st.head puts A on the axes frame
# while the panel's leftmost ink is its rail and key 59.3 px further left, which is the
# offset fig6 carried as a panel shift. Here the panel is where the gaps want it, so the
# letter moves instead and lands on D.

# TUNING IS DATA, NOT CODE, while the frame is being solved. FIG5PUB_TUNE names a JSON of
# {letter: [dx_px, dw_px]} that overrides PLAN, so the solver can iterate without rewriting
# this file. The settled values are written back into PLAN above and the file stands alone.
if os.environ.get("FIG5PUB_TUNE"):
    import json as _json
    _ov = _json.load(open(os.environ["FIG5PUB_TUNE"]))
    PLAN = [(n, p_, o, r, c, _ov.get(n, [dx, dw])[0], _ov.get(n, [dx, dw])[1], ld)
            for n, p_, o, r, c, dx, dw, ld in PLAN]

_HEADSET = re.compile(r'^ax = fig\.add_subplot\(\w+\[[^\]]*\]\); st\.head\(fig, ax, "[a-f]"\)')
_TUNE = re.compile(r"^ax\._axes_d[xw]_in\s*=")
AXES = {}
for _new, _path, _old, _r, _c, _dx, _dw, _ldx in PLAN:
    _ns = NS5 if _path.startswith("fig5") else NS6
    _ns.update(fig=fig, _GS=_GS, gs=_GS, _PS=_PS, _SNAP=_SNAP, _YAX=_YAX)
    _body = []
    for _l in source_block(_path, _old):
        if _HEADSET.match(_l):
            _body.append('ax = fig.add_subplot(_GS[%d, %d]); '
                         'st.head(fig, ax, "%s", dx_in=%.6f)'
                         % (_r, _c, _new, _ldx / 300.0))
        elif _TUNE.match(_l):
            _body.append("# [frame] " + _l)          # stripped, see PLAN
        else:
            _body.append(_l)
    exec(compile("\n".join(_body), "%s:%s" % (_path, _old), "exec"), _ns)
    _ax = _ns["ax"]
    _ax._axes_dx_in = _dx / 300.0
    _ax._axes_dw_in = _dw / 300.0
    AXES[_new] = _ax

# ---------------------------------------------------------------- save
st.save(fig, "fig5pub")
for _t, _q in _SNAP:
    _t.set_position((_q.get_x() + _q.get_width() / 2.0, _t.get_position()[1]))
_PAD_IN = 18.0 / 300.0
from PIL import Image as _Im
_px = np.asarray(_Im.open("fig5pub.png").convert("RGB")).astype(int)
_ys, _xs = np.nonzero(_px.sum(2) < 735)
_H_in = fig.get_size_inches()[1]
_trim = Bbox.from_extents(_xs.min() / 300.0 - _PAD_IN,
                          _H_in - _ys.max() / 300.0 - _PAD_IN,
                          _xs.max() / 300.0 + _PAD_IN,
                          _H_in - _ys.min() / 300.0 + _PAD_IN)
for _ext in (".png", ".pdf"):
    fig.savefig("fig5pub" + _ext, bbox_inches=_trim)
print("trimmed to %.3f x %.3f in (%.0f x %.0f px at 300 dpi), 18 px pads"
      % (_trim.width, _trim.height, _trim.width * 300, _trim.height * 300))

if os.environ.get("FIG5PUB_FINAL"):
    import PIL.Image as _PI, pymupdf as _mu
    _ref = _mu.open("fig1.pdf")[0].rect.width / 72.0
    _dpi = 300.0 * 12.0 / _ref
    fig.savefig("fig5pub_final.png", dpi=_dpi, bbox_inches=_trim)
    _im = _PI.open("fig5pub_final.png")
    _im.save("fig5pub_final.png", dpi=(300, 300))
    print("final: native %.4f in at the shared %.2f dpi raster -> %d x %d px declared 300 dpi "
          "(%.3f x %.3f in), body type %.3f pt"
          % (_trim.width, _dpi, _im.size[0], _im.size[1], _im.size[0] / 300.0,
             _im.size[1] / 300.0, 8.860 * 12.0 / _ref))

print("panels: " + " | ".join("%s<-%s %s" % (n, p.split("_")[0], o)
                              for n, p, o, _, _, _, _, _ in PLAN))
