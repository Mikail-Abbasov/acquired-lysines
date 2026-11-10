"""
abbasov_style.py — figure defaults for Bioorganic & Medicinal Chemistry submission.

Elsevier expects vector output, sans-serif labels, and figures legible at column width.
Liberation Sans is metrically identical to Arial, so text sizes and line breaks match what
Arial would produce at the printer.

Sizes: BMC single column ~85 mm, double column ~180 mm.
"""
import os as _os, sys as _sys, math as _math
# The extra figure libraries (matplotlib-venn, adjustText, rdkit, spectrum_utils, colorspacious…)
# live in figtools/py rather than site-packages: this Mac's VM has no root and is wiped between
# sessions, so they sit on real disk in the project folder instead. Put them on the path here so
# every figure script picks them up without anyone having to remember to source figtools/env.sh.
_FIGTOOLS = _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))),
                          "figtools", "py")
if _os.path.isdir(_FIGTOOLS) and _FIGTOOLS not in _sys.path:
    _sys.path.append(_FIGTOOLS)

import matplotlib as mpl
import matplotlib.text as mtext
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

MM = 1/25.4
# Nature's "Guide to preparing final artwork", read 23 September 2026:
#   "89 mm width (1 column)", "1.5 columns where necessary - 120 mm or 136 mm",
#   "183 mm width (2 column)", "The full depth of the page is 247 mm."
#   Panel labels "8 pt bold, upright (not italic) a, b, c's"; other text max 7 pt, min 5 pt.
# Elsevier's own grid, which this submission is going to, is 90 / 140 / 190 with a 7 pt
# finished lettering target and a 6 pt floor on sub/superscripts. At 183 mm the figures sit
# 4% narrow for Elsevier and production would scale them up by 1.038.
SINGLE = 89*MM      # 3.50 in
ONEHALF = 120*MM    # 4.72 in   (136 mm is the other permitted 1.5-column width)
DOUBLE = 183*MM     # 7.20 in

# THE FAMILY SHEET. Every main figure is built on an 11.5400 in canvas and rendered at
# 300 x 12.00 / 11.5400 = 311.96 dpi, so its body type prints 9.213 pt. A sheet built on the
# 7.20 in DOUBLE canvas carries the same source type on a canvas 1.6 times smaller, so once
# both are placed at the text width that type prints 1.6 times LARGER. The supplementary
# sheets were left on DOUBLE when the main set moved to the family canvas on 4 November, and
# the SI has been printing 8.0 pt body type against the main text's 5.0 pt ever since.
# sheet() brings a DOUBLE-based figsize onto the family canvas without touching its aspect.
FAMILY = 11.5400
FAMK   = FAMILY / DOUBLE          # 1.60176


YAX = 1.55020   # the family y axis: 483.6 px at the 311.96 dpi raster


def block_h(rows, hspace, ratios=None, yax=YAX):
    """Height of the panel block, in inches, for a gridspec whose unit row is one y axis.

    A figure's height follows from its y axis, never the other way round. matplotlib spends
    the block on rows in proportion to height_ratios and on gaps of hspace x the MEAN row, so
    a row of ratio 1.0 lands on yax when the block is this tall. With no ratios this is the
    (rows + hspace x gaps) x yax that every main figure uses.
    """
    r = [1.0] * int(rows) if ratios is None else [float(x) for x in ratios]
    mean = sum(r) / len(r)
    return yax * (sum(r) + hspace * (len(r) - 1) * mean)


def sheet(w, h):
    """A DOUBLE-based figsize brought to the family canvas, aspect unchanged."""
    return (w * FAMK, h * FAMK)
MAX_DEPTH = 247*MM  # the full depth of a Nature page

# THE FIVE COLOURS. Supplied by M.A. on 2026-09-01 and used as given, with nothing added:
#   CC3F0C rust  ·  C5C3C6 light grey  ·  46494C charcoal  ·  4C5C68 slate  ·  1985A1 teal
# Earlier builds also carried a hand-made tint (#80BCCB) and shade (#125D71) of the teal for
# ordered ramps. They were not in the brief and they were not needed: the five already form a
# monotone lightness sequence, so the ramps are now built from the five themselves. Both
# derived values are deleted from every palette below so no script can reach for them again.
#
# L* (0 black, 100 white):  light grey 79.0 · teal 51.3 · rust 47.5 · slate 38.2 · charcoal 30.9
# dE2000, normal vision:    teal/rust 55.5 · light grey/charcoal 45.8 · light grey/rust 39.2
#                           slate/rust 38.1 · light grey/slate 37.7 · rust/charcoal 33.8
#                           light grey/teal 33.0 · teal/charcoal 25.7 · teal/slate 18.0
#                           charcoal/slate 8.5  <- the one pair that must never share a panel
#
# ROLES, one each:
#   teal      BLUE       primary series (charge, reactivity, "found")
#   rust      ORANGE     the mutation itself, and anything flagged
#   charcoal  GREEN      second series (shape, ligandability, "missed")
#   slate     GREY       annotation ink: rules, tick labels, footnotes, n= counts
#   light grey LIGHTGREY background fills, "not covered", reference bars
# Charcoal and slate sit 8.5 apart, so they never carry two series in the same panel: charcoal
# is a series colour, slate is ink. Everything else clears a dE of 18.
#
# PAPER and INK are the page and matplotlib's default text colour, not palette entries.
PAPER = "#FFFFFF"
INK   = "#000000"
VIVID = dict(BLUE="#1F6FBF", ORANGE="#D85A30", GREEN="#2E8B57", GREY="#888780", LIGHTGREY="#CFCEC9")
CALM  = dict(BLUE="#3A6BA5", ORANGE="#C2764B", GREEN="#0E7D6B", GREY="#8A8A84", LIGHTGREY="#DAD8D3")
WARM  = dict(BLUE="#4A5BA0", ORANGE="#BC6B40", GREEN="#3F7A3C", GREY="#8A8078", LIGHTGREY="#DCD3C7")
EARTH = dict(BLUE="#0A9396", ORANGE="#BB3E03", GREEN="#EE9B00", GREY="#7A736A", LIGHTGREY="#E9D8A6")
SLATE = dict(BLUE="#1985A1", ORANGE="#CC3F0C", GREEN="#46494C", GREY="#4C5C68", LIGHTGREY="#C5C3C6")
BLUE, ORANGE = SLATE["BLUE"], SLATE["ORANGE"]
GREEN, GREY, LIGHTGREY = SLATE["GREEN"], SLATE["GREY"], SLATE["LIGHTGREY"]

def slate(): globals().update(SLATE)

# ---------------------------------------------------------------- the R61 house standard
# NCI-R61/scripts/housestyle.py, spec in NCI-R61/results/panel_style_spec_2026-10-28.md.
# One body type size per figure and every other size a fixed ratio to it. The palette is
# already the house five (SLATE is this module's default), so the delta is type and
# furniture, not color.
LETTER_K    = 1.55      # panel letter, multiple of body. Was 1.96 against body 7.0;
                        # lowered so the letter still prints at 12.4 pt after BODY_K
                        # raised the body. The letters are unchanged on the page.
LINE_K      = 0.0935    # spine and tick width, pt per pt of body
TICKMAJ_K   = 0.44
TICKMIN_K   = 0.234
VALUE_OFF_K = 0.55      # value label offset above its mark, pt per pt of body
TEXT = "#000000"
FURN = SLATE["GREEN"]   # #46494C, the house dark, for edges and ticks

_BODY = 7.0             # the body the manuscript figures were drawn against
BODY  = 7.0             # the active body, set by apply(body)

# Measured 4 November against Cell Press, which asks for "about 6-8 pt at the desired print
# size": all 1,514 body strings in the 18 figures printed BELOW 6 pt, at 4.1 to 5.4. Nothing
# in the set was in range. The figures are drawn 7.20 in wide and placed at 6.5 in, so every
# size rides a 0.903 scale down, and the ladder was never raised to pay for it.
#
# Each figure sets its own BODY = 7.0 and calls apply(BODY), so the knob has to live here
# rather than in eighteen scripts. apply() multiplies whatever body it is given by BODY_K.
BODY_K = float(_os.environ.get('BODYK', '1.2657'))         # 7.0 -> 8.86, chosen so the lowest rung in use prints at exactly 6.0

# The house type ladder. Measured 3 November: the 18 figures were drawing 25 distinct type
# sizes between 3.5 and 13.7, fig1 alone using 16, which is what stopped them reading as one
# set. Four roles cover every call site:
#
#   label  7.0  axis labels, legends, panel titles. This is the body itself.
#   note   6.0  in-panel annotations and value labels
#   dense  5.0  ordinary tick labels
#   micro  4.2  tick labels on an axis with more than about 25 categories
#
# The rungs step by about 1.17, and the top rung is the body, which is what housestyle sets a
# panel title in. Anchoring the ladder on the body rather than below it also stops the set
# drifting smaller: 5 and 6 were the two most used sizes, and a ladder that snapped 6 down to
# 5.6 would have shrunk most of the annotation type in the paper.
#
# The panel letter is not on this ladder: it is LETTER_K x body, 13.7 at body 7.
#
# sz() snaps to the nearest rung rather than making every caller name a role, so the 513 call
# sites keep the intent they were written with and the ladder itself stays adjustable in one
# place. Change a number here and every figure follows.
# The rungs are given at body 7 and scaled by BODY_K, so what each one PRINTS at 0.903 is:
#
#   label  7.00 -> 8.86 drawn -> 8.00 pt   axis labels, legends, panel titles
#   note   6.30 -> 7.97 drawn -> 7.20 pt   in-panel annotations and value labels
#   dense  5.70 -> 7.21 drawn -> 6.51 pt   ordinary tick labels
#   micro  5.25 -> 6.64 drawn -> 6.00 pt   tick labels on a crowded axis
#
# The old ladder (4.2, 5.0, 6.0, 7.0) spanned 1.67 between its ends, and 6-8 pt is a span of
# only 1.33, so no uniform lift could put both ends in the band. The ladder is compressed as
# well as raised. Checked against all 491 st.sz() call sites: every one re-snaps in the same
# order it snapped before, so no panel's type hierarchy inverts.
LADDER = (5.25, 5.7, 6.3, 7.0)

# Bar geometry. Every bar group occupies the same slot footprint, so a single bar,
# a pair and a triplet all read as the same width of ink on the page.
BARW = 0.70


def grp(n):
    """Width for one member of an n-member bar group sharing one category slot."""
    return BARW / float(n)


# BARW above is a fraction of a CATEGORY SLOT, which is not what a reader sees. A two-category
# panel renders a 0.70 slot as a 50 pt slab; a forty-five-category panel renders the same 0.70
# as a 1 pt hairline. Measured across the paper before this pass: 0.9 pt to 56 pt, a 60-fold
# spread, with every bar at the same slot fraction. Page thickness has to be set on the page.
#
# A single constant thickness is wrong in the other direction: two bars at 6 pt in a 1.5 in
# panel read as two needles marooned in white space. So the rendered thickness is CLAMPED into
# a band, and the slot fraction decides where inside the band each panel lands.
BAR_MIN_PT = 3.0
# 10.0 until 2026-11-05. What a reader sees is the FILL FRACTION, bar against slot, and an
# absolute cap of 10 pt was deciding that fraction for most panels instead of BARW. Measured
# across the paper at 10 pt: fig1 d 14%, fig1 g 20%, fig2 c 21%, fig1 e 35%, fig2 b 36%,
# fig1 c 46%, fig8 a 47%, fig1 f 70%, fig2 e 92%, fig3 e 92%. Only the panels too crowded to
# reach the cap kept the 0.70 they were authored with. At 16 pt every panel that has the room
# lands on 0.70, the crowded ones are untouched, and the widest bar in the paper is 0.22 in,
# which is nobody's idea of a slab.
BAR_MAX_PT = 16.0


def _bar_axis(cont, patches):
    """'v' for ax.bar, 'h' for ax.barh. Read it from the container, measure only as fallback."""
    o = getattr(cont, "orientation", None)
    if o in ("vertical", "horizontal"):
        return "v" if o == "vertical" else "h"
    # Fallback: a vertical bar chart varies its patches along x and shares a baseline.
    xs = {round(q.get_x(), 9) for q in patches}
    ys = {round(q.get_y(), 9) for q in patches}
    return "v" if len(xs) >= len(ys) else "h"


def fit_bars(fig, verbose=False):
    """Set every bar's thickness ON THE PAGE, after layout. Never touches a bar's LENGTH.

    Runs last, because points per data unit depends on the final axes position and limits.
    Three cases are handled separately:

      uniform        every bar in the container is the same thickness -> clamp to the band
      encoded        the thickness varies by design (figS4 e scales width by identification
                     count) -> scale by ONE factor so the widest bar lands on BAR_MAX_PT and
                     the encoding survives
      grouped        two or more containers offset within one slot -> share the slot between
                     them and re-centre, so the group keeps its geometry
    """
    from matplotlib.container import BarContainer
    changed = []; crowded = []
    for ax in fig.axes:
        conts = [c for c in getattr(ax, "containers", []) if isinstance(c, BarContainer)]
        if not conts:
            continue
        for axis in ("v", "h"):
            sel = []
            for c in conts:
                pts = [q for q in c.patches
                       if q.get_width() > 0 or q.get_height() > 0]
                if pts and _bar_axis(c, pts) == axis:
                    sel.append((c, pts))
            if not sel:
                continue
            # points per data unit along the CATEGORY axis
            t = ax.transData
            if axis == "v":
                ppd = abs(t.transform((1, 0))[0] - t.transform((0, 0))[0])
                thick = lambda q: q.get_width()
                centre = lambda q: q.get_x() + q.get_width() / 2.0
                def place(q, c_, w_): q.set_x(c_ - w_ / 2.0); q.set_width(w_)
            else:
                ppd = abs(t.transform((0, 1))[1] - t.transform((0, 0))[1])
                thick = lambda q: q.get_height()
                centre = lambda q: q.get_y() + q.get_height() / 2.0
                def place(q, c_, w_): q.set_y(c_ - w_ / 2.0); q.set_height(w_)
            ppd = ppd / fig.dpi * 72.0                      # pixels -> points
            if ppd <= 0:
                continue

            # Stacks: containers whose bar centres coincide are stacked or overlaid, not
            # grouped, and must keep one centre between them.
            stacks = []
            for c, pts in sel:
                key = tuple(round(centre(q), 6) for q in pts)
                for s in stacks:
                    if s["key"] == key:
                        s["members"].append((c, pts)); break
                else:
                    stacks.append({"key": key, "members": [(c, pts)], "n": len(pts)})
            # A set of stacks is a GROUP sharing one slot only if each stack carries the
            # same number of bars at parallel positions. Six single-bar barh() calls in a
            # loop are six bars in six different slots, NOT a six-way group: treating them
            # as one would re-centre every bar away from the row it labels. figS7 panel a
            # is exactly that shape and was being scrambled.
            counts = {len(s["key"]) for s in stacks}
            grouped = len(stacks) > 1 and counts == {max(counts)} and max(counts) >= 2
            if not grouped:
                stacks = [{"key": tuple(sorted(c_ for s in stacks for c_ in s["key"])),
                           "members": [m for s in stacks for m in s["members"]],
                           "n": sum(s["n"] for s in stacks)}]
            n = len(stacks)

            # Slot pitch. For a group, one stack's own spacing IS the slot pitch. Otherwise
            # take it from the union of every bar centre on this axis.
            pitch = None
            for s in stacks:
                cs = sorted(set(s["key"]))
                if len(cs) > 1:
                    d = min(b - a for a, b in zip(cs, cs[1:]))
                    pitch = d if pitch is None else min(pitch, d)
            slot_pt = pitch * ppd if pitch else None

            for k, s in enumerate(stacks):
                pts = [q for _, ps in s["members"] for q in ps]
                cur = [thick(q) * ppd for q in pts]
                if not cur or max(cur) <= 0:
                    continue
                uniform = (max(cur) - min(cur)) <= 1e-6 * max(cur)
                # The slot is a HARD cap. Raising a hairline to BAR_MIN_PT past the room
                # its slot allows makes neighbouring bars touch, which trades one defect for
                # a worse one. A panel whose room is below the floor is too crowded for the
                # band and is reported rather than forced.
                room = (slot_pt * 0.92 / n) if slot_pt else BAR_MAX_PT
                if room < BAR_MIN_PT:
                    crowded.append((getattr(ax, "_panel_letter", "?"), axis,
                                    round(room, 2), len(pts)))
                if uniform:
                    want = min(max(max(cur), BAR_MIN_PT), BAR_MAX_PT, room)
                    scale = want / max(cur)
                else:
                    # keep the encoding, land the widest bar on the band's top
                    scale = min(BAR_MAX_PT, room) / max(cur)
                if abs(scale - 1.0) < 0.02:
                    continue
                for c_, ps in s["members"]:
                    for q in ps:
                        w_new = thick(q) * scale
                        c_old = centre(q)
                        c_orig, w_orig = c_old, thick(q)
                        if n > 1 and slot_pt:
                            # re-centre the group: the slot centre is this bar's centre with
                            # its old offset removed, and the new offset is set from w_new
                            off_old = (k - (n - 1) / 2.0) * (thick(q))
                            slot_c = c_old - off_old
                            c_old = slot_c + (k - (n - 1) / 2.0) * w_new
                        place(q, c_old, w_new)
                        # THE VALUE LABEL GOES WITH ITS BAR. Re-centring a group moves the
                        # bars and nothing else, so "7.3%" stayed where the wide bar used to
                        # be and ended up over the gap beside the narrow one. figS2 a and c
                        # had all four labels adrift that way.
                        _dx = c_old - c_orig
                        if abs(_dx) > 1e-9:
                            # ONCE PER LABEL. The tolerance has to be wide enough to catch a
                            # label sitting over its bar and narrow enough to miss the next
                            # group's, and with two groups a quarter of a slot apart there is
                            # no such width: figS2 c moved 15.0% onto its bar on the first
                            # stack's pass and the second stack's pass moved it straight back.
                            _tol = max(abs(w_orig), abs(w_new)) * 0.75
                            for _t in ax.texts:
                                if getattr(_t, "_bar_label_moved", False):
                                    continue
                                if _t.get_transform() is not ax.transData:
                                    continue
                                _x, _y = _t.get_position()
                                _v = _x if axis == "v" else _y
                                if abs(_v - c_orig) <= _tol:
                                    _t.set_position((_x + _dx, _y) if axis == "v"
                                                    else (_x, _y + _dx))
                                    _t._bar_label_moved = True
                changed.append((getattr(ax, "_panel_letter", "?"), axis,
                                round(max(cur), 2), round(max(cur) * scale, 2),
                                "grouped/%d" % n if n > 1 else
                                ("encoded" if not uniform else "uniform")))
    if verbose:
        for L, a, was, now, kind in changed:
            print("   bar fit  panel %-3s %s  %6.2f -> %5.2f pt  %s" % (L, a, was, now, kind))
        for L, a, room, nb in crowded:
            print("   CROWDED  panel %-3s %s  slot allows only %.2f pt for %d bars, below the "
                  "%.1f pt floor" % (L, a, room, nb, BAR_MIN_PT))
    return changed


def sz(at7):
    """A type size, snapped to the house ladder and scaled to the active body.

    The argument is the size the panel asked for at body 7. It is snapped to the nearest rung
    so that the figures share one set of sizes, then scaled, so at body 7 a rung is itself.
    """
    x = float(at7)
    rung = min(LADDER, key=lambda t: (abs(t - x), t))
    return rung / _BODY * BODY
def earth(): globals().update(EARTH)
def warm():  globals().update(WARM)
def calm():  globals().update(CALM)
def vivid(): globals().update(VIVID)

_SUP = {"-": "\u207b", "0": "\u2070", "1": "\u00b9", "2": "\u00b2", "3": "\u00b3",
        "4": "\u2074", "5": "\u2075", "6": "\u2076", "7": "\u2077", "8": "\u2078",
        "9": "\u2079"}


def sup(v):
    """An exponent as real glyphs instead of mathtext.

    Matplotlib sets a mathtext superscript at 0.7 of the surrounding size, so a 7 pt label
    carrying 10$^{-10}$ puts 4.9 pt on the page, under Elsevier's 6 pt floor for a
    superscript. A Unicode superscript is an ordinary glyph at full size, and Liberation
    Sans carries the whole set."""
    return "".join(_SUP[c] for c in str(v))


def plain_log(ax, axis="y", fmt=None, sci=False):
    """Label a log axis with plain numerals, for the same reason as sup().

    The default used to be "{:,.0f}", zero decimal places, which is correct only while every
    tick is at least 1. Below that it rounds the label to nothing: fig6's occupancy inset put
    ticks at 1, 10^-2 and 10^-4 and the axis read "1, 0, 0". Decimals are now chosen from the
    value, so a tick below 1 prints at the precision it needs. Pass fmt to override, or
    sci=True for an exponent on an axis with no room for the plain form.
    """
    from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter
    import math

    def _f(v, _):
        if fmt is not None:
            return fmt.format(v)
        if v == 0:
            return "0"
        if abs(v) >= 1:
            return "{:,.0f}".format(v)
        # sci is for an axis with no room for the plain form. fig6's occupancy inset is
        # 15 mm wide and "0.0001" reached out of it and across the schematic beside it.
        if sci:
            return "$10^{%d}$" % int(round(math.log10(abs(v))))
        d = int(math.ceil(-math.log10(abs(v))))
        return ("{:,.%df}" % d).format(v)

    a = ax.yaxis if axis == "y" else ax.xaxis
    a.set_major_locator(LogLocator(base=10))
    a.set_major_formatter(FuncFormatter(_f))
    a.set_minor_formatter(NullFormatter())
    return ax


def missing_glyphs(fig):
    """Characters the figure draws that the chosen face does not carry.

    This is not a matplotlib error. A missing glyph falls back to another font silently, or
    prints nothing at all. Substituting Unicode superscripts for mathtext looked like an
    improvement until the render showed Liberation Sans carries only the Latin-1 set
    (U+00B9, U+00B2, U+00B3) and has no U+2070, U+2074 to U+2079, U+207B or U+2082, so the
    exponents on Figure 1c and 1d simply disappeared. Check before believing the file."""
    import matplotlib.font_manager as fm
    try:
        from fontTools.ttLib import TTFont
    except ImportError:
        return []
    name = mpl.rcParams["font.sans-serif"][0]
    try:
        cmap = TTFont(fm.findfont(fm.FontProperties(family=name))).getBestCmap()
    except Exception:
        return []
    bad = set()
    for ax in list(fig.axes) + [fig]:
        for t in list(getattr(ax, "texts", [])) + list(
                getattr(ax, "get_xticklabels", list)() if hasattr(ax, "get_xticklabels") else []) + list(
                getattr(ax, "get_yticklabels", list)() if hasattr(ax, "get_yticklabels") else []):
            for c in t.get_text():
                if c.strip() and ord(c) > 127 and ord(c) not in cmap and c != "$":
                    bad.add(c)
    return sorted(bad)


# The page a sheet has to sit on, in inches: the width it is placed at and the depth of the
# text block it shares with its legend.
TEXT_W, TEXT_H = 6.5, 9.0


def check_depth(fig, stem="", width_in=None, cap_in=None):
    """Depth as the READER gets it, not as the canvas is built.

    This compared the canvas height against Nature's 247 mm page, which only tested anything
    while a sheet was built at its printed size. The main set moved to the 11.54 in family
    canvas on 4 November and has been silently exempt ever since: fig1 is 5.01 in tall at
    face value and passes a 9.72 in cap without ever being measured. What matters is the
    depth the sheet prints at once it is placed at the text width, legend still to come.
    """
    width_in = TEXT_W if width_in is None else width_in
    cap_in = TEXT_H if cap_in is None else cap_in
    w, h = fig.get_size_inches()
    placed = width_in * h / w
    if placed > cap_in + 1e-6:
        print("OVER DEPTH: %s prints %.2f in deep at the %.1f in text width, against a "
              "%.1f in text block, and its legend has to fit under it"
              % (stem or "figure", placed, width_in, cap_in))
        return False
    return True


def fit_inside(fig, pad_mm=0.8, passes=4):
    """Repeat the inward fit until it converges.

    Moving the axes changes where the tick labels sit, so one pass leaves a fraction of a
    millimetre hanging over. Three or four passes settle it."""
    for _ in range(passes):
        before = [tuple(a.get_position().bounds) for a in fig.axes]
        _fit_once(fig, pad_mm)
        if [tuple(a.get_position().bounds) for a in fig.axes] == before:
            break
    return fig


def _fit_once(fig, pad_mm=0.8):
    """Bring everything the figure draws inside the figure rectangle.

    bbox="tight" solved overflow by growing the canvas, which is why no two figures came out
    the same width. The canvas is now fixed at the column width, so overflow has to be
    solved the other way: measure how far the drawn content runs past each edge and move
    every axes inward by that much, keeping their relative positions. The plotting area
    shrinks a little, which is the correct trade: the labels genuinely need the room."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    W, H = fig.get_size_inches()
    pad = pad_mm * MM
    try:
        bb = fig.get_tightbbox(r)
    except TypeError:
        bb = fig.get_tightbbox(renderer=r)
    if bb is None:
        return fig
    over_l = max(0.0, -bb.x0) + pad
    over_r = max(0.0, bb.x1 - W) + pad
    over_b = max(0.0, -bb.y0) + pad
    over_t = max(0.0, bb.y1 - H) + pad
    if max(over_l, over_r, over_b, over_t) <= pad * 1.001:
        return fig
    fl, fr = over_l / W, over_r / W
    fb, ft = over_b / H, over_t / H
    sx, sy = 1.0 - fl - fr, 1.0 - fb - ft
    if sx <= 0.4 or sy <= 0.4:                 # refuse a squeeze that would be a redesign
        print("FIT REFUSED: content runs %.1f mm past the %.0f mm trim; the layout needs "
              "changing, not scaling" % (max(over_l, over_r) / MM, W / MM))
        return fig
    for ax in fig.axes:
        p = ax.get_position()
        ax.set_position([fl + p.x0 * sx, fb + p.y0 * sy, p.width * sx, p.height * sy])
    return fig


def seq_cmap(name="abbasov_seq"):
    """Continuous low-to-high map anchored on the palette: paper -> light grey -> teal
    (L* 100 -> 79 -> 51, saturation rising with value). A continuous map necessarily contains
    values between its anchors; the anchors are the constraint.

    Charcoal was tried as a fourth anchor and removed: it put the only chromatic step in the
    MIDDLE of the scale, so mid-range cells outshone the maximum. The most saturated colour
    belongs at the top of a sequential scale, not in its middle. The same mistake, made with
    five discrete steps, is why Fig. 1c is not a ramp either."""
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list(name, [PAPER, LIGHTGREY, BLUE])


def apply(body=None):
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"],
        "font.size": 7,
        "axes.labelsize": 7,
        "axes.titlesize": 8,
        "xtick.labelsize": 6.5,
        "ytick.labelsize": 6.5,
        "legend.fontsize": 6.5,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "xtick.direction": "out", "ytick.direction": "out",
        "axes.spines.top": False, "axes.spines.right": False,
        "lines.linewidth": 1.0,
        "lines.markersize": 3,
        "legend.frameon": False,
        "figure.dpi": 300,
        "savefig.dpi": 600,
        # bbox="tight" saved whatever the labels happened to need rather than the width the
        # figure declares, so the eighteen figures came out between 176.9 and 216.9 mm
        # against a declared 180. Production scales each one to the column, and that scaling
        # multiplies every type size on the page. save() now fits the content inside the
        # declared rectangle instead, so the trim is exact and nothing is cropped.
        "savefig.bbox": None,
        "savefig.pad_inches": 0.0,
        # Matplotlib's mathtext default is DejaVu Sans, so every $...$ expression was being
        # set in a DIFFERENT typeface from the body text — PyMuPDF found whole alphabets of
        # DejaVu embedded in fig4 and figS5. Point mathtext at the same face as everything else.
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold",
        "mathtext.sf": "Arial",
        # cal/tt/bfit default to "cursive"/"monospace"/"sans:italic:bold". Any figure with a
        # log axis builds the whole custom fontset, resolves "cursive", finds nothing, and
        # silently falls back to DejaVu Sans - the same leak fixed above, by another door.
        "mathtext.cal": "Arial:italic",
        "mathtext.tt": "Arial",
        "mathtext.bfit": "Arial:italic:bold",
        "mathtext.default": "regular",
        "pdf.fonttype": 42,     # embed as TrueType so text stays editable in Illustrator
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })
    if body is None:
        # The other scripts in this directory were drawn against the fixed sizes above.
        # Leaving them alone is deliberate: only the manuscript figures opt in.
        return
    global BODY
    BODY = float(body) * BODY_K
    mpl.rcParams.update({
        "font.size": BODY,
        "axes.labelsize": BODY, "axes.titlesize": BODY,
        "xtick.labelsize": sz(6.5), "ytick.labelsize": sz(6.5),
        "legend.fontsize": sz(6.5),
        "text.color": TEXT, "axes.labelcolor": TEXT, "axes.edgecolor": FURN,
        "xtick.color": FURN, "ytick.color": FURN,
        "xtick.labelcolor": TEXT, "ytick.labelcolor": TEXT,
        "axes.linewidth": LINE_K * BODY,
        "xtick.major.width": LINE_K * BODY, "ytick.major.width": LINE_K * BODY,
        "xtick.minor.width": LINE_K * BODY, "ytick.minor.width": LINE_K * BODY,
        "xtick.major.size": TICKMAJ_K * BODY, "ytick.major.size": TICKMAJ_K * BODY,
        "xtick.minor.size": TICKMIN_K * BODY, "ytick.minor.size": TICKMIN_K * BODY,
    })

def blank_panels(fig):
    """Axes containing nothing but their panel letter.

    Every figure script guards its panels with `if INPUT:` and most of those guards have no
    else branch, so a missing input skips the panel, the script still exits 0, and
    build_figures prints "ok". That is how fig1 e and fig2 d/f/g/h shipped blank for an hour
    on 2026-09-01 with a clean build log. Nothing here prevents the guard from firing; it
    makes the result impossible to miss.
    """
    def _content(ax):
        return (len(ax.lines) + len(ax.patches) + len(ax.collections) + len(ax.images)
                + len(ax.tables) + len(getattr(ax, "containers", []))
                + len([t for t in ax.texts if len(t.get_text().strip()) > 2]))

    boxes = [(a, a.get_position()) for a in fig.axes]
    out = []
    for ax, bb in boxes:
        if _content(ax):
            continue
        # A LABEL-ONLY CONTAINER is not a blank panel. figS2 g is an axis-off axes that exists
        # only to carry its letter while the spectrum is drawn into a subgridspec inside it;
        # flagging that as empty would be a false alarm. An axes with other axes nested inside
        # its own rectangle is a container, so judge it by what those hold instead.
        nested = [a for a, b in boxes
                  if a is not ax and b.x0 >= bb.x0 - 1e-6 and b.x1 <= bb.x1 + 1e-6
                  and b.y0 >= bb.y0 - 1e-6 and b.y1 <= bb.y1 + 1e-6]
        if nested and any(_content(a) for a in nested):
            continue
        out.append(getattr(ax, "_panel_letter", "?"))
    return out


def translucent_artists(fig):
    """Artists drawn at less than full opacity.

    A palette colour at alpha 0.5 on white paper renders as a colour that is not in the
    palette - a tint by another name - and neither the pixel audit nor the eye reliably
    catches it, because a blend toward the paper is indistinguishable from antialiasing.
    The brief was five colours and no derived values, so opacity is not available as an
    emphasis channel; use a different palette colour, line weight, or dash pattern instead.

    Returns (panel letter, artist type, alpha) for every offender. Legend and title artists
    are exempt only when fully opaque, so nothing is silently skipped.

    An artist carrying ``_alpha_ok`` is a case the PI has looked at and approved, and is left
    out of the report with its reason kept on the artist. Opacity is still not an emphasis
    channel; the one accepted use is overplotting, where alpha reports density rather than
    inventing a colour. Without this, an approved exception prints on every run until it is
    noise and the check stops being read.
    """
    out = []
    for ax in fig.axes:
        letter = getattr(ax, "_panel_letter", "?")
        groups = (list(ax.lines) + list(ax.patches) + list(ax.collections)
                  + list(ax.images) + list(ax.texts))
        for art in groups:
            if getattr(art, "_alpha_ok", None):
                continue
            a = art.get_alpha()
            if a is None:
                continue
            try:
                a = float(a)
            except (TypeError, ValueError):
                continue                      # per-element alpha array: reported below
            if a < 0.999:
                out.append((letter, type(art).__name__, round(a, 3)))
    return out


def panel_label(ax, letter, dx=-0.18, dy=1.06):
    ax._panel_letter = letter
    ax.text(dx, dy, letter, transform=ax.transAxes,
            fontsize=8, fontweight="bold", style="normal", va="top", ha="left")

def head(fig, ax, letter, title="", dx_in=0.0):
    """Record the house panel header for this axes. Nothing is drawn yet.

    housestyle draws the letter straight away because its figures are one panel per column
    with a fixed layout. These are gridspec, and save() runs fit_inside(), which moves every
    axes to bring drawn content inside the trim. A letter positioned at panel-creation time
    would be left behind by that move, so the header is recorded here and placed by save()
    from the final axes boxes, with fit_inside() and placement interleaved until both settle.
    """
    ax._panel_letter = letter
    ax._panel_title = title
    ax._panel_dx_in = dx_in
    if not hasattr(fig, "_heads"):
        fig._heads = []
    if ax not in fig._heads:
        fig._heads.append(ax)
    return ax


def text_collisions(fig, frac=0.12):
    """Every pair of text artists whose drawn boxes overlap by more than `frac`.

    This exists because the check it replaces could not fail. The ink-margin test measured
    the figure's OUTER bounding box, so it reported "clear" for a figure whose tick labels
    were printing on top of each other in the middle of a panel: 30 such pairs were sitting
    in eight figures, found only by reading the rendered PDFs. A whole-figure box can never
    see an internal collision, so the guard has to compare the boxes to each other.

    Rotated labels are included: get_window_extent returns the drawn box, not the glyph run.
    """
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    # ax.twinx() and ax.twiny() hide the duplicated AXIS, not the individual label artists:
    # each tick label still answers get_visible() == True while never being drawn. Counting
    # them reported 23 "stacked duplicates" at 100% across fig3 and figS6 - the twin's labels
    # sitting exactly on the originals they are hidden behind. They are not on the page, so
    # they are not collisions. Drop anything belonging to an invisible axis or axes.
    ghost = set()
    for a_ in fig.axes:
        for nm in ("xaxis", "yaxis"):
            axis = getattr(a_, nm, None)
            # ax.axison is the flag axis("off") sets; the individual x/y axis objects keep
            # answering get_visible() == True under it, and their tick labels are rebuilt
            # with DEFAULT styling the moment anything reads majorTicks. Figure 1b is drawn
            # with axis("off") and reported a phantom 18% overlap between its heading and a
            # "10" that is nowhere on the page.
            if axis is None or (getattr(a_, "axison", True) and axis.get_visible()
                                and a_.get_visible()):
                continue
            ghost.add(id(axis.get_label()))
            for tk in getattr(axis, "majorTicks", []) + getattr(axis, "minorTicks", []):
                for lab in (getattr(tk, "label1", None), getattr(tk, "label2", None)):
                    if lab is not None:
                        ghost.add(id(lab))
        if not a_.get_visible():
            for t in list(a_.texts) + [a_.title]:
                ghost.add(id(t))
        # A TICK LABEL OUTSIDE ITS OWN AXIS LIMITS IS NOT ON THE PAGE. A log axis keeps a
        # full decade ladder whatever the view is: figure 3 panel c runs 0.75 to 205 and
        # still carries live, get_visible() == True labels for 0.01, 0.1, 1,000 and 10,000,
        # whose reported boxes span 1,450 px and land on two other panels. That produced
        # one "collision" and one "intrusion" against text a reader cannot see. The same
        # happens on a linear axis whose locator overshoots the limits.
        for axis, lim in ((getattr(a_, "xaxis", None), a_.get_xlim()),
                          (getattr(a_, "yaxis", None), a_.get_ylim())):
            if axis is None:
                continue
            lo, hi = (lim[0], lim[1]) if lim[0] <= lim[1] else (lim[1], lim[0])
            span = hi - lo
            for tk in getattr(axis, "majorTicks", []) + getattr(axis, "minorTicks", []):
                try:
                    v = float(tk.get_loc())
                except Exception:
                    continue
                if lo - span * 1e-6 <= v <= hi + span * 1e-6:
                    continue
                for lab in (getattr(tk, "label1", None), getattr(tk, "label2", None)):
                    if lab is not None:
                        ghost.add(id(lab))
    items = []
    for t in fig.findobj(mtext.Text):
        try:
            if id(t) in ghost or not t.get_visible() or not str(t.get_text()).strip():
                continue
            # An Annotation's window extent is the union of its TEXT and its LEADER LINE, so
            # two callouts whose leaders merely cross each other measure as a text collision.
            # figS10's panel b reported 'MDAMB231' at 111 px tall for a 30 px label and gave
            # five false overlaps at 40-98% between labels that do not touch. Measure the text.
            if isinstance(t, mtext.Annotation):
                bb = mtext.Text.get_window_extent(t, renderer=r)
            else:
                bb = t.get_window_extent(renderer=r)
        except Exception:
            continue
        if bb.width <= 0 or bb.height <= 0:
            continue
        items.append((bb, str(t.get_text()).strip(), t))
    out = []
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            a, b = items[i][0], items[j][0]
            if (min(a.x1, b.x1) - max(a.x0, b.x0) <= 0
                    or min(a.y1, b.y1) - max(a.y0, b.y0) <= 0):
                continue
            f = _overlap_frac(items[i][2], a, items[j][2], b, r)
            if f > frac:
                out.append((f, items[i][1], items[j][1]))
    out.sort(reverse=True)
    return out


def _twin_extra(fig, ax, bb, r, FH, lx1=None):
    """Height a twin axis has stacked on this panel's top edge, in inches.

    A panel carrying an ax.twiny() puts that twin's tick numbers and axis label ALONG
    ITS TOP EDGE, which is exactly where the letter goes. figS6 panel c printed its
    letter C on top of the twin's "0.4" at 72%. Panels without a twin return 0.0.

    ONLY WHAT THE LETTER WOULD ACTUALLY HIT. lx1 is the letter's right edge in inches.
    A letter placed flush to its panel's leftmost ink sits out beside the y-axis title,
    well left of the plot frame, so the twin's furniture is nowhere near it and lifting
    the letter clear of something it cannot touch is pure damage: it put fig4's D 116.3 px
    above its own row-mate F for a twin whose nearest tick number starts 295 px to its
    right. Artists entirely right of the letter are skipped. With lx1 unset the old
    behaviour stands, so nothing that does not pass it changes.
    """
    out = 0.0
    for a2 in fig.axes:
        if a2 is ax or not a2.get_visible():
            continue
        p2 = a2.get_position()
        if abs(p2.x0 - bb.x0) > 1e-6 or abs(p2.y1 - bb.y1) > 1e-6:
            continue
        for t2 in _live(a2, "x") + [a2.xaxis.get_label()]:
            try:
                if not t2.get_visible() or not str(t2.get_text()).strip():
                    continue
                e2 = t2.get_window_extent(r)
                if lx1 is not None and e2.x0 / fig.dpi > lx1:
                    continue
                out = max(out, e2.y1 / fig.dpi - bb.y1 * FH)
            except Exception:
                pass
    return out


def _ink_x0(ax, r):
    """Leftmost display-x of everything this axes actually draws, or None.

    The axes RECTANGLE is not ink. A venn panel has no spines and no ticks, so its
    frame sits 33 px left of the first thing a reader can see, while a panel whose
    heading overhangs its own cell has ink 29 px LEFT of the frame. A letter placed
    on the frame is therefore ragged against the panel either way.
    """
    best = None
    # AN ARTIST CANNOT BE INK FURTHER LEFT THAN THE AXES' OWN TIGHT BBOX. matplotlib computes
    # that bbox from a real draw; get_window_extent on a tick label or on a blended-transform
    # line can answer from a stale or uninitialised position, and on fig6 panel d it put two
    # lines at -5.7 px and every tick label at 0 px, which dragged the panel letter 0.98 in
    # left of the panel and off the canvas, where _fit_once then shrank every axes in the
    # figure by 7% trying to bring it back. Figure-level text is exempt: a caption anchored to
    # the gridspec cell is that panel's ink and legitimately sits outside the axes bbox.
    try:
        _floor = ax.get_tightbbox(r).x0 - 1.0
    except Exception:
        _floor = None

    def _note(art, text=False, guard=True):
        nonlocal best
        try:
            if not art.get_visible():
                return
            if text and not str(art.get_text()).strip():
                return
            if text and isinstance(art, mtext.Annotation):
                e = mtext.Text.get_window_extent(art, renderer=r)
            elif text:
                e = art.get_window_extent(renderer=r)
            else:
                e = art.get_window_extent()
            if e.width <= 0 or e.height <= 0:
                return
            if guard and _floor is not None and e.x0 < _floor:
                return
            best = e.x0 if best is None else min(best, e.x0)
        except Exception:
            pass
    for art in list(ax.texts):
        _note(art, text=True)
    # A caption pinned to the gridspec CELL rather than to the axes is still that panel's
    # ink, and a reader sees no difference. Figure 1a's union line is one: leaving it out
    # moved the panel letter 10 px right, onto the venn circle, the moment the caption was
    # re-anchored. Count figure-level text that falls inside this panel's own cell.
    try:
        _c = ax.get_subplotspec().get_position(ax.figure)
        _fw, _fh = ax.figure.get_size_inches()
        _dpi = ax.figure.dpi
        _cx0, _cx1 = _c.x0 * _fw * _dpi, _c.x1 * _fw * _dpi
        _cy0, _cy1 = _c.y0 * _fh * _dpi, _c.y1 * _fh * _dpi
        for art in list(ax.figure.texts):
            if not art.get_visible() or not str(art.get_text()).strip():
                continue
            e = art.get_window_extent(renderer=r)
            mx, my = (e.x0 + e.x1) / 2.0, (e.y0 + e.y1) / 2.0
            if _cx0 <= mx <= _cx1 and _cy0 <= my <= _cy1:
                _note(art, text=True, guard=False)
    except Exception:
        pass
    for seq in (ax.collections, ax.patches, ax.lines, ax.images):
        for art in list(seq):
            _note(art)
    # A panel drawn with axis("off") keeps tick label artists that still answer
    # get_visible() == True while never being drawn. Counting those ghosts put a mark
    # 38.7 px left of anything a reader can see in Figure 1b and would have dragged the
    # letter out there with it. Same rule as text_collisions: an invisible axis or axes
    # contributes no ink.
    if ax.get_visible() and getattr(ax, "axison", True):
        for nm in ("xaxis", "yaxis"):
            axis = getattr(ax, nm, None)
            if axis is None or not axis.get_visible():
                continue
            for tk in getattr(axis, "majorTicks", []):
                for lab in (getattr(tk, "label1", None), getattr(tk, "label2", None)):
                    if lab is not None:
                        _note(lab, text=True)
            _note(axis.get_label(), text=True)
    return best


def place_heads(fig, rise_in=0.055, gap_in=0.075):
    """Create or re-anchor the recorded panel headers against the current axes boxes.

    Oversized bold letter at the panel's left edge, title beside it on the same baseline.
    Nothing else in the figure is bold.
    """
    if not getattr(fig, "_heads", None):
        return fig
    FW, FH = fig.get_size_inches()
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    # Letters in one gridspec row share a baseline. fit_inside only ever SHRINKS an
    # axes, so a panel whose labels once overflowed keeps a lower top for the rest of
    # the run and its letter followed it down: Figure 1 printed A 120 px below B, in
    # the same row. Take the row's highest top as the base and add each panel's own
    # clearance to that, so a twin-axis panel still rises above its row-mates where
    # it has to and a row whose tops already agree does not move at all.
    _key, _base, _xtra, _ttop = {}, {}, {}, {}
    for ax in fig._heads:
        if not ax.get_visible():
            continue
        ss = ax.get_subplotspec() if hasattr(ax, "get_subplotspec") else None
        try:
            k = (id(ss.get_gridspec()), ss.rowspan.start, ss.rowspan.stop)
        except Exception:
            k = id(ax)
        _key[ax] = k
        bb = ax.get_position()
        # The base is the row's highest AXES top. Twin-axis clearance is deliberately NOT
        # pooled into it: carrying C's 95.5 px of clearance across figS6's row put E's
        # letter on top of E's own y-axis label at 41%, trading a ragged row for a
        # collision. A row-mate of a twin panel therefore stays on the row baseline and
        # the twin panel alone rises above it.
        # The letter's own x, computed the same way the placement loop below computes it,
        # so the twin check knows where the letter will actually be.
        _lx0 = bb.x0 * FW + getattr(ax, "_panel_dx_in", 0.0)
        if getattr(ax, "_panel_head_flush", False):
            _ix0 = _ink_x0(ax, r)
            if _ix0 is not None:
                _lx0 = _ix0 / fig.dpi + getattr(ax, "_panel_dx_in", 0.0)
        _xtra[ax] = _twin_extra(fig, ax, bb, r, FH, _lx0 + LETTER_K * BODY / 72.0)
        # OPT-IN: put this letter's TOP on the twin apparatus's own top line. The filtered
        # clearance above is the right default, because a letter out beside the y-axis title
        # cannot hit a tick number 295 px to its right. It also leaves the letter beside the
        # twin's tick numbers rather than beside its label, which on fig4 d read as a letter
        # sitting in the middle of the panel's header rather than at the head of it. Panels
        # that want the second behaviour ask for it by name, one at a time, because it lifts
        # the letter off its row's baseline and every other panel in that row stays put.
        if getattr(ax, "_panel_head_twin_top", False):
            _ttop[ax] = _twin_extra(fig, ax, bb, r, FH)
        _base[k] = bb.y1 * FH if k not in _base else max(_base[k], bb.y1 * FH)
    for ax in fig._heads:
        if not ax.get_visible():
            continue
        bb = ax.get_position()
        x0 = bb.x0 * FW + getattr(ax, "_panel_dx_in", 0.0)
        if getattr(ax, "_panel_head_flush", False):
            _ix = _ink_x0(ax, r)
            if _ix is not None:
                x0 = _ix / fig.dpi + getattr(ax, "_panel_dx_in", 0.0)
                # Recorded so an audit can report what the letter was aligned TO. Measuring
                # the ink again after save() does not reproduce it: reading majorTicks on a
                # decoration-off panel rebuilds a default tick set that was never drawn.
                ax._head_ink_x0 = _ix
        _extra = _xtra.get(ax, 0.0)
        if _os.environ.get("NOROWBASE") == "1":
            y = (bb.y1 * FH + _extra + rise_in) / FH
        else:
            y = (_base.get(_key.get(ax), bb.y1 * FH) + _extra + rise_in) / FH
        t = getattr(ax, "_head_letter_artist", None)
        if t is None:
            # Letters are DISPLAYED capitalised. _panel_letter stays lowercase because it
            # is also the identifier the ledger, the crop filenames and the published
            # panel_map are keyed on, and those must not move.
            t = fig.text(x0 / FW, y, str(ax._panel_letter).upper(), ha="left", va="bottom",
                         fontsize=LETTER_K * BODY, fontweight="bold", color=TEXT)
            ax._head_letter_artist = t
        else:
            t.set_position((x0 / FW, y))
        if ax in _ttop:
            # Measured, not computed: the letter's drawn height is what has to land on the
            # twin's top, and 13.733 pt of bold Arial is not 13.733 pt of ink.
            _hh = t.get_window_extent(renderer=r).height / fig.dpi
            y = (bb.y1 * FH + _ttop[ax] - _hh) / FH
            t.set_position((x0 / FW, y))
        if ax._panel_title:
            wl = t.get_window_extent(renderer=r).width / fig.dpi
            tt = getattr(ax, "_head_title_artist", None)
            if tt is None:
                tt = fig.text((x0 + wl + gap_in) / FW, y, ax._panel_title, ha="left",
                              va="bottom", fontsize=BODY, fontweight="bold", color=TEXT)
                ax._head_title_artist = tt
            else:
                tt.set_position(((x0 + wl + gap_in) / FW, y))
    return fig


def value(ax, x, y, text, color, body=None, italic=False):
    """A value label above its mark, in that mark's own color, never bold."""
    body = BODY if body is None else float(body)
    ax.annotate(text, xy=(x, y), xytext=(0, VALUE_OFF_K * body),
                textcoords="offset points", ha="center", va="baseline",
                fontsize=body, color=color, style="italic" if italic else "normal")

# Raising the body 27% to reach Cell Press's 6-8 pt put 335 text pairs into collision across
# the 18 figures. Classified by what each side was, 189 of them - 49% - were tick labels
# hitting OTHER TICK LABELS ON THEIR OWN AXIS. That is not a design problem, it is the same
# problem fit_bars solves for bars: how much room a label needs is a PAGE fact, and the
# script that wrote the label cannot know it. So it is settled here, once, after layout.
#
# Three moves, in order of how little they cost the reader:
#   numeric axis      thin the ticks. Fewer gridlines on a continuous scale costs nothing.
#   categorical x     rotate further. A category label cannot be dropped without orphaning
#                     its bar, but it can be stood up.
#   last resort       step the type down the ladder, and never below the micro rung, since
#                     the whole point of the exercise is that nothing prints under 6 pt.
# An axis that still collides after all three is reported, not silently mangled.
TICK_ROT_STEPS = (42.0, 55.0, 70.0, 90.0)
TICK_FRAC = 0.02        # any measurable overlap between two labels on one axis is too much


def _is_numeric_axis(labels):
    """Do all these tick labels read as plain numbers? Mathtext and names do not."""
    ok = 0
    for t in labels:
        x = str(t.get_text()).strip().replace("\u2212", "-").replace("\u2013", "-")
        if not x:
            continue
        try:
            float(x.replace(",", "").replace("%", ""))
            ok += 1
        except ValueError:
            return False
    return ok > 0


def _oriented_rect(t, r):
    """The four display-space corners of a text's OWN box, or None if it is upright.

    An axis-aligned bounding box is the wrong shape for rotated type. Two tick labels set
    at 45 degrees and stepping down a diagonal have boxes that overlap by two thirds while
    not a pixel of ink touches, and that phantom is what drove fit_ticks to stand labels up
    to 70 degrees across this figure set. The AABB of a rotated rectangle shares its centre,
    so the unrotated size plus that centre reconstructs the real box.
    """
    try:
        orig = t.get_rotation()
        rot = float(orig) % 180.0
    except Exception:
        return None
    if min(rot, 180.0 - rot) < 1.0:
        return None                      # upright: the AABB already IS the box
    try:
        t.set_rotation(0)
        up = t.get_window_extent(renderer=r)
        w, h = up.width, up.height
    except Exception:
        return None
    finally:
        t.set_rotation(orig)             # exactly what it was, always
    return (w, h)


def _corners(t, bb, wh):
    """Corners of the oriented box, given its unrotated size and its rotated AABB."""
    w, h = wh
    cx, cy = (bb.x0 + bb.x1) / 2.0, (bb.y0 + bb.y1) / 2.0
    th = _math.radians(float(t.get_rotation()))
    c, s = _math.cos(th), _math.sin(th)
    return [(cx + dx * c - dy * s, cy + dx * s + dy * c)
            for dx, dy in ((-w / 2, -h / 2), (w / 2, -h / 2), (w / 2, h / 2), (-w / 2, h / 2))]


def _poly_area(poly):
    a = 0.0
    for i in range(len(poly)):
        x0, y0 = poly[i]
        x1, y1 = poly[(i + 1) % len(poly)]
        a += x0 * y1 - x1 * y0
    return abs(a) / 2.0


def _poly_clip(subject, clip):
    """Sutherland-Hodgman: intersection of two CONVEX polygons, counter-clockwise."""
    out = list(subject)
    for i in range(len(clip)):
        if not out:
            return []
        ax_, ay_ = clip[i]
        bx_, by_ = clip[(i + 1) % len(clip)]
        ex, ey = bx_ - ax_, by_ - ay_
        inside = lambda p: (ex * (p[1] - ay_) - ey * (p[0] - ax_)) >= -1e-9
        new, prev = [], out[-1]
        for cur in out:
            if inside(cur):
                if not inside(prev):
                    new.append(_seg_cross(prev, cur, (ax_, ay_), (bx_, by_)))
                new.append(cur)
            elif inside(prev):
                new.append(_seg_cross(prev, cur, (ax_, ay_), (bx_, by_)))
            prev = cur
        out = [q for q in new if q is not None]
    return out


def _seg_cross(p0, p1, q0, q1):
    x1, y1 = p0; x2, y2 = p1; x3, y3 = q0; x4, y4 = q1
    d = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    if abs(d) < 1e-12:
        return None
    t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / d
    return (x1 + t * (x2 - x1), y1 + t * (y2 - y1))


def _ccw(poly):
    a = 0.0
    for i in range(len(poly)):
        x0, y0 = poly[i]
        x1, y1 = poly[(i + 1) % len(poly)]
        a += x0 * y1 - x1 * y0
    return poly if a >= 0 else poly[::-1]


def _overlap_frac(ta, bba, tb, bbb, r):
    """Overlap of two labels as a fraction of the smaller, honest about rotation.

    Falls back to the axis-aligned ratio when both are upright, which is exact there.
    """
    qa = _oriented_rect(ta, r) if ta is not None else None
    qb = _oriented_rect(tb, r) if tb is not None else None
    if qa is None and qb is None:
        ov = (max(0.0, min(bba.x1, bbb.x1) - max(bba.x0, bbb.x0))
              * max(0.0, min(bba.y1, bbb.y1) - max(bba.y0, bbb.y0)))
        sm = min(bba.width * bba.height, bbb.width * bbb.height)
        return (ov / sm) if sm > 0 else 0.0
    box = lambda bb: [(bb.x0, bb.y0), (bb.x1, bb.y0), (bb.x1, bb.y1), (bb.x0, bb.y1)]
    pa = _corners(ta, bba, qa) if qa else box(bba)
    pb = _corners(tb, bbb, qb) if qb else box(bbb)
    aa, ab = _poly_area(pa), _poly_area(pb)
    sm = min(aa, ab)
    if sm <= 0:
        return 0.0
    inter = _poly_clip(_ccw(pa), _ccw(pb))
    return (_poly_area(inter) / sm) if inter else 0.0


def _worst_overlap(labels, r):
    """Largest pairwise overlap among these labels, as a fraction of the smaller box."""
    bs = []
    for t in labels:
        try:
            bb = t.get_window_extent(renderer=r)
        except Exception:
            continue
        if bb.width > 0 and bb.height > 0:
            bs.append((bb, t))
    worst = 0.0
    for i in range(len(bs)):
        for j in range(i + 1, len(bs)):
            (a, ta), (b, tb) = bs[i], bs[j]
            if (min(a.x1, b.x1) - max(a.x0, b.x0) <= 0
                    or min(a.y1, b.y1) - max(a.y0, b.y0) <= 0):
                continue
            worst = max(worst, _overlap_frac(ta, a, tb, b, r))
    return worst


def _live(ax, which):
    """The tick label artists ALREADY on the page, read without disturbing them.

    ax.get_xticklabels() is not safe here. It routes to Axis.get_major_ticks(), which
    materialises Tick objects against the current locator, and doing that before the layout
    loop changes what the figure ends up drawing: measured on fig3, simply calling it took
    the figure from 6 text collisions to 15 without this function altering anything. Reading
    axis.majorTicks after a draw touches only what the draw already built.
    """
    axis = ax.xaxis if which == "x" else ax.yaxis
    out = []
    for tk in getattr(axis, "majorTicks", []):
        for lab in (getattr(tk, "label1", None), getattr(tk, "label2", None)):
            if lab is not None and lab.get_visible() and str(lab.get_text()).strip():
                out.append(lab)
    return out


def _tick_state(fig):
    """Snapshot everything fit_ticks may touch, so it can be put back exactly.

    Stored BY POSITION, not by artist identity: changing a locator makes matplotlib build a
    different set of Tick objects, so holding references to the old label artists restores
    attributes onto objects the figure no longer draws. That is what left fig5 at 36
    collisions after a revert that should have returned it to 29.
    """
    snap = []
    for ax in fig.axes:
        for which in ("x", "y"):
            axis = ax.xaxis if which == "x" else ax.yaxis
            labs = _live(ax, which)
            snap.append((ax, which, axis.get_major_locator(),
                         [(t.get_rotation(), t.get_ha(), t.get_fontsize(),
                           t.get_rotation_mode()) for t in labs]))
    # Axes POSITIONS as well. fit_inside() moves panels and is not reversible, so restoring
    # only the labels left figS10 at 74 collisions after a revert meant to return it to 61.
    return (snap, [(a_, tuple(a_.get_position().bounds)) for a_ in fig.axes])


def _tick_restore(fig, state):
    snap, boxes = state
    for ax, which, loc, _ in snap:
        axis = ax.xaxis if which == "x" else ax.yaxis
        axis.set_major_locator(loc)
    for a_, bounds in boxes:
        a_.set_position(bounds)
    fig.canvas.draw()
    for ax, which, _, labs in snap:
        live = _live(ax, which)
        for t, (rot, ha, size, rmode) in zip(live, labs):
            t.set_rotation(rot); t.set_ha(ha); t.set_fontsize(size)
            t.set_rotation_mode(rmode)
    for a_, bounds in boxes:
        a_.set_position(bounds)
    fig.canvas.draw()


def fit_ticks(fig, verbose=True):
    """Make tick labels stop overlapping, but only where that helps the whole page.

    The first version of this optimised each axis on its own and made two figures worse:
    standing fig5's panel b categories up from 0 to 70 degrees cured that axis and drove the
    taller labels into the panel's annotation block, 29 collisions to 37. An axis does not
    get to improve itself at the figure's expense, so every candidate change is scored on the
    FIGURE's total collision count and reverted unless the total falls.

    Three moves, in order of how little they cost the reader:
      numeric axis    thin the ticks. Fewer gridlines on a continuous scale costs nothing.
      categorical x   rotate further. A category label cannot be dropped without orphaning
                      its bar, but it can be stood up.
      last resort     step the type down the ladder, never below the micro rung, since the
                      point of the exercise is that nothing prints under 6 pt.
    """
    MICRO = min(LADDER) / _BODY * BODY
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    total = lambda: len(text_collisions(fig, frac=0.12))
    best = total()
    report = []

    for ax in fig.axes:
        for which in ("x", "y"):
            labs = _live(ax, which)
            if len(labs) < 2 or _worst_overlap(labs, r) <= TICK_FRAC:
                continue
            tag = "%s%s" % (getattr(ax, "_panel_letter", "?"), which)
            axis = ax.xaxis if which == "x" else ax.yaxis
            loc0 = axis.get_major_locator()
            rot0 = labs[0].get_rotation()
            ha0 = labs[0].get_ha()
            sz0 = labs[0].get_fontsize()

            def restore():
                axis.set_major_locator(loc0)
                fig.canvas.draw()
                for t in _live(ax, which):
                    t.set_rotation(rot0); t.set_ha(ha0); t.set_fontsize(sz0)
                fig.canvas.draw()

            cands = []
            if _is_numeric_axis(labs):
                vals = list(loc0())
                for k in (2, 3, 4, 5, 6):
                    keep = vals[::k]
                    if vals and vals[-1] not in keep:
                        keep = keep + [vals[-1]]
                    cands.append(("thinned %d -> %d labels" % (len(vals), len(keep)),
                                  ("loc", keep)))
            elif which == "x":
                for rot in TICK_ROT_STEPS:
                    if rot > rot0:
                        cands.append(("rotated %.0f -> %.0f deg" % (rot0, rot), ("rot", rot)))
            for rung in sorted(LADDER, reverse=True):
                cand = rung / _BODY * BODY
                if cand < sz0 - 1e-9 and cand >= MICRO - 1e-9:
                    cands.append(("type %.2f -> %.2f pt" % (sz0, cand), ("size", cand)))

            won = None
            for desc, (kind, val) in cands:
                if kind == "loc":
                    axis.set_major_locator(mticker.FixedLocator(val))
                elif kind == "rot":
                    for t in _live(ax, which):
                        t.set_rotation(val); t.set_ha("right"); t.set_rotation_mode("anchor")
                else:
                    for t in _live(ax, which):
                        t.set_fontsize(val)
                fig.canvas.draw()
                now = total()
                if now < best and _worst_overlap(_live(ax, which), r) <= TICK_FRAC:
                    best = now; won = desc
                    break
                restore()
            if won:
                if verbose:
                    print("   tick fit %-4s %s" % (tag, won))
            else:
                w = _worst_overlap(_live(ax, which), r)
                report.append((tag, len(labs), w))
    if report and verbose:
        for tag, n, w in report:
            print("   TICKS CROWDED  %-4s %d labels overlap by %.0f%%, and no rotation, "
                  "thinning or 6 pt-floor type change improved the figure"
                  % (tag, n, 100 * w))
    return report


def nudge_axes(fig):
    """Apply each axes' declared _axes_dx_in / _axes_dy_in / _axes_dw_in change, once.

    _axes_dh_in changes the HEIGHT and holds the TOP edge, for a short panel sharing a row
    with a tall one. Figure 1's row of a, b, c and d is as tall as b, and c and d are bar
    panels with no aspect lock, so without this they stretch to b's height instead of
    keeping their own. Holding the top is what keeps all four letters on one baseline.

    _axes_dw_in changes the WIDTH and holds x0, for a panel whose spacing is set by how
    far apart its categories can sit rather than by the gridspec. Figure 1g is two bars
    whose tick labels are wide, so a full-cell axes forces the bars apart or the end
    margins open. Narrowing the axes costs nothing on the page, because repanel_extract.py
    crops every panel to its own tight ink box and discards the white around it.

    A panel whose content does not reach its own frame cannot be aligned with its
    neighbours by the gridspec: Figure 1a's venn is centred in its axes and aspect-locked,
    so its leftmost ink sat 33 px right of panel c's frame directly below it. Shifting the
    whole axes is the only move that takes the diagram, its caption and its letter together.
    Applied before the layout loop so place_heads reads the moved ink, and guarded so the
    loop's repeated calls do not stack the offset.
    """
    FW, FH = fig.get_size_inches()
    for ax in fig.axes:
        if getattr(ax, "_axes_nudged", False):
            continue
        dx = getattr(ax, "_axes_dx_in", 0.0)
        dy = getattr(ax, "_axes_dy_in", 0.0)
        dw = getattr(ax, "_axes_dw_in", 0.0)
        dh = getattr(ax, "_axes_dh_in", 0.0)
        ax._axes_nudged = True
        if not dx and not dy and not dw and not dh:
            continue
        b = ax.get_position()
        ax.set_position([b.x0 + dx / FW, b.y0 + dy / FH - dh / FH,
                         b.width + dw / FW, b.height + dh / FH])
    return fig


def _categorical_y(ax):
    """True when the y tick labels are names rather than numbers.

    A categorical axis has no values to read off it, so "end on a major tick" does not apply:
    fig2 d plots two named rows, missed and found, and rounding its top up to the next step
    only added two unlabelled ticks to the spine.

    The formatter is ASKED for each label rather than inspected. set_yticklabels usually
    installs a FixedFormatter, but this module's own tick handling leaves a FuncFormatter in
    place, so an isinstance check on FixedFormatter missed fig2 d entirely. Calling the
    formatter works whatever it is, and a numeric axis that has not been drawn yet returns
    empty strings, which are skipped, so nothing is wrongly exempted.
    """
    fmt = ax.yaxis.get_major_formatter()
    try:
        locs = list(ax.get_yticks())
    except Exception:
        return False
    for i, loc in enumerate(locs):
        try:
            t = str(fmt(loc, i)).strip()
        except Exception:
            continue
        if not t:
            continue
        try:
            float(t.replace("\u2212", "-").replace(",", ""))
        except ValueError:
            return True
    return False


def end_on_tick(fig, verbose=False):
    """Extend each linear y axis so it stops on a major tick, never between two.

    An axis whose top is 1.30 x the tallest bar lands wherever the data put it, so the
    spine runs past the last labelled tick into unlabelled space: Figure 1e ended at 36.7
    with its last tick at 30. Round the top up to the next whole step of the tick spacing
    the locator already chose and relabel to match. Axes that already end on a tick, and
    anything not linear, are left alone.
    """
    out, bot = [], []
    for ax in fig.axes:
        if not ax.get_visible() or not getattr(ax, "axison", True):
            continue
        if ax.get_yscale() != "linear" or not ax.yaxis.get_visible():
            continue
        if _categorical_y(ax):
            continue
        try:
            y0, y1 = ax.get_ylim()
            locs = [float(t) for t in ax.get_yticks()]
        except Exception:
            continue
        inside = [t for t in locs if y0 - 1e-9 <= t <= y1 + 1e-9]
        if len(inside) < 2 or y1 <= y0:
            continue
        gaps = [b - a for a, b in zip(inside, inside[1:])]
        step = min(gaps)
        if step <= 0 or max(gaps) - step > step * 1e-6:
            continue                                  # not a uniform grid, leave it alone
        if abs(y1 - inside[-1]) <= step * 1e-6:
            continue                                  # already ends on a tick
        # Anchor on the EXISTING ticks, never on the axis bottom. Anchoring on y0 shifted
        # the whole grid whenever the axis had headroom below zero, and turned fig2 c's
        # clean 0, 2, 4 into -0.9, 0.1, 1.1, 2.1 ... which is worse than the defect.
        import math as _m
        base = inside[0]
        k = _m.ceil((y1 - base) / step - 1e-9)
        top = base + k * step
        added = k - (len(inside) - 1)
        if added < 1 or added > 2:
            # More than two new ticks means the author's tick set was deliberately sparse
            # and this would be redesigning their axis, not finishing it.
            continue
        ticks = [base + i * step for i in range(k + 1)]
        ax.set_ylim(y0, top)
        ax.set_yticks(ticks)
        out.append((getattr(ax, "_panel_letter", "?"), round(y1, 4), round(top, 4)))

    # THE BOTTOM IS AN END OF THE AXIS TOO. The pass above finishes the top and leaves the
    # bottom wherever the data margin put it, so thirteen supplementary panels had a spine
    # running past its own lowest label: figS5 d started at 35 under a first tick of 40.
    # Two ways to finish it, and which one is right depends on the data, not on the number:
    # raise the bottom onto the tick when nothing is drawn below that tick, and only extend
    # downward when something is.
    for ax in fig.axes:
        if not ax.get_visible() or not getattr(ax, "axison", True):
            continue
        if ax.get_yscale() != "linear" or not ax.yaxis.get_visible() or _categorical_y(ax):
            continue
        try:
            y0, y1 = ax.get_ylim()
            locs = [float(t) for t in ax.get_yticks()]
            dlo = float(ax.dataLim.y0)
        except Exception:
            continue
        inside = [t for t in locs if y0 - 1e-9 <= t <= y1 + 1e-9]
        if len(inside) < 2 or y1 <= y0:
            continue
        gaps = [b - a for a, b in zip(inside, inside[1:])]
        step = min(gaps)
        if step <= 0 or max(gaps) - step > step * 1e-6:
            continue
        low = inside[0]
        if abs(y0 - low) <= step * 1e-6:
            continue                                  # already ends on a tick
        if not (dlo == dlo) or abs(dlo) == float("inf"):
            continue                                  # no usable data limit, leave it alone
        if dlo >= low - step * 1e-6:
            ax.set_ylim(low, y1)                      # nothing below the tick: close the gap
            new_lo = low
        else:
            import math as _m
            m = _m.ceil((low - dlo) / step - 1e-9)
            if m > 2:
                continue                              # redesigning the axis, not finishing it
            new_lo = low - m * step
            ax.set_ylim(new_lo, y1)
            ax.set_yticks([new_lo + i * step
                           for i in range(int(round((y1 - new_lo) / step)) + 1)])
        bot.append((getattr(ax, "_panel_letter", "?"), round(y0, 4), round(new_lo, 4)))

    if verbose and out:
        for L, was, now in out:
            print("   y axis  panel %s   top %g -> %g, now a labelled tick" % (L, was, now))
    if verbose and bot:
        for L, was, now in bot:
            print("   y axis  panel %s   bottom %g -> %g, now a labelled tick" % (L, was, now))
    return fig


def final_sheet(fig, stem, trim=None, ref="fig1.pdf"):
    """The finalized file, written at the raster the whole family shares.

    THE FAMILY SHARES A RASTER, NOT A WIDTH. Figure 1 at 12.00 in fixes it at
    300 x 12.00 / 11.5400 = 311.96 dpi; every sheet renders at that raster and is declared
    300 dpi, so all of them print 9.213 pt body type whatever their width. Each main figure
    carries its own copy of this block. The ten supplementary sheets share this one, so the
    raster cannot drift between them.
    """
    import PIL.Image as _PI, pymupdf as _mu
    _ref = _mu.open(ref)[0].rect.width / 72.0
    dpi = 300.0 * 12.0 / _ref
    out = "%s_final.png" % stem
    fig.savefig(out, dpi=dpi, **({"bbox_inches": trim} if trim is not None else {}))
    im = _PI.open(out)
    im.save(out, dpi=(300, 300))
    print("final: %s at the shared %.2f dpi raster -> %d x %d px declared 300 dpi "
          "(%.3f x %.3f in), body type %.3f pt"
          % (out, dpi, im.size[0], im.size[1], im.size[0] / 300.0, im.size[1] / 300.0,
             8.860 * dpi / 300.0))
    return out


def save(fig, stem):
    """Write both a vector file for the journal and a PNG for review.

    Also refuses to be quiet about an empty panel: see blank_panels(). The build treats the
    marker line below as a failure, so a missing input can no longer produce a clean log.
    """
    # The headers and the fit depend on each other: placing a letter above a panel can push
    # the content past the trim, and fitting moves the panel the letter is anchored to. A
    # fixed number of rounds is not enough, because the last placement can reintroduce the
    # overflow the last fit removed: fig7's letters a and b sat 2.2 pt above the page that
    # way. Iterate until a fit changes nothing, then place once more against those boxes.
    nudge_axes(fig)
    for _ in range(6):
        place_heads(fig)
        before = [tuple(a_.get_position().bounds) for a_ in fig.axes]
        fit_inside(fig, passes=1)
        if [tuple(a_.get_position().bounds) for a_ in fig.axes] == before:
            break
    place_heads(fig)
    # Tick labels are fitted against the SETTLED layout, not before it. Scoring a rotation
    # against a layout that has not happened yet is what made fig5 worse: the change looked
    # like an improvement pre-layout and cost 8 collisions once the panels moved.
    if _os.environ.get('NOFITTICKS') != '1':
        # fit_ticks must not be able to make a figure worse. Its candidates are scored one at
        # a time, but the fit_inside/place_heads that follow shift the page again afterwards,
        # so the only honest verdict is the finished count. Snapshot, fit, re-settle, compare,
        # and put everything back if the page did not improve.
        _snap = _tick_state(fig)
        _base = len(text_collisions(fig, frac=0.12))
        fit_ticks(fig, verbose=True)
        fit_inside(fig, passes=1)
        place_heads(fig)
        if len(text_collisions(fig, frac=0.12)) > _base:
            _tick_restore(fig, _snap)
            place_heads(fig)
            print("   tick fit reverted: the settled page was no better than before")
    # Finishing the y axes comes AFTER the tick fit, not before. Run first, it installs a
    # FixedLocator that fit_ticks can no longer thin, and fig2 c went from three labels to
    # seven. Run here, it only extends a top that the thinning has already settled.
    if end_on_tick(fig, verbose=True) is not None:
        fit_inside(fig, passes=1)
        place_heads(fig)
    # Bar thickness is set on the PAGE, so it has to run after the layout has settled:
    # points per data unit depends on the final axes position and limits. It changes no
    # bar's length, so nothing here can move a plotted value.
    fit_bars(fig, verbose=True)
    _coll = text_collisions(fig)
    if _coll:
        print("TEXT COLLISION: %s has %d overlapping text pair(s); the first few are %s"
              % (stem, len(_coll),
                 "; ".join("%.0f%% %r over %r" % (f * 100, a[:30], b[:30])
                           for f, a, b in _coll[:4])))
    check_depth(fig, stem)
    gone = missing_glyphs(fig)
    if gone:
        print("MISSING GLYPHS: %s draws %s, which %s does not carry; matplotlib "
              "will drop or substitute them silently"
              % (stem, ", ".join("U+%04X %r" % (ord(c), c) for c in gone),
                 mpl.rcParams["font.sans-serif"][0]))
    blank = blank_panels(fig)
    if blank:
        print("BLANK PANEL: %s has %d empty panel(s) — %s — an input guard skipped the body"
              % (stem, len(blank), ", ".join(str(b) for b in blank)))
    trans = translucent_artists(fig)
    if trans:
        from collections import Counter
        tally = Counter("%s:%s@%s" % t for t in trans)
        print("TRANSLUCENT: %s draws %d artist(s) below full opacity — %s — alpha over the "
              "paper is a tint, which the palette does not allow"
              % (stem, len(trans), ", ".join("%s x%d" % (k, v) for k, v in
                                             sorted(tally.items()))))
    fig.savefig(f"{stem}.pdf")
    fig.savefig(f"{stem}.png")
    # The finalized sheet comes off the SAME settled figure, so it cannot differ from the
    # reviewed one by a layout pass. Main figures run their own block after their trim.
    if _os.environ.get("FINAL_SHEET"):
        final_sheet(fig, stem)
    return f"{stem}.pdf", f"{stem}.png"


# Superseded scripts on disk (fig1.py, fig2_v2.py, figS_class_v3.py, the _*.py sketches) still
# reference MIDBLUE and DARKBLUE. Those two values were a tint and a shade of the teal that the
# five-colour brief did not include, and they are gone. Anything reaching for them should say so
# out loud rather than fail with a bare AttributeError.
def __getattr__(name):
    if name in ("MIDBLUE", "DARKBLUE"):
        raise AttributeError(
            "abbasov_style.%s was a hand-made %s of the teal and is not part of the five-colour "
            "palette. For an ordered sequence use seq_cmap(); for a second series use GREEN."
            % (name, "tint" if name == "MIDBLUE" else "shade"))
    raise AttributeError("module 'abbasov_style' has no attribute %r" % name)
