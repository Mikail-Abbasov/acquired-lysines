"""
abbasov_style.py — figure defaults for Bioorganic & Medicinal Chemistry submission.

Elsevier expects vector output, sans-serif labels, and figures legible at column width.
Liberation Sans is metrically identical to Arial, so text sizes and line breaks match what
Arial would produce at the printer.

Sizes: BMC single column ~85 mm, double column ~180 mm.
"""
import os as _os, sys as _sys
# The extra figure libraries (matplotlib-venn, adjustText, rdkit, spectrum_utils, colorspacious…)
# live in figtools/py rather than site-packages: this Mac's VM has no root and is wiped between
# sessions, so they sit on real disk in the project folder instead. Put them on the path here so
# every figure script picks them up without anyone having to remember to source figtools/env.sh.
_FIGTOOLS = _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))),
                          "figtools", "py")
if _os.path.isdir(_FIGTOOLS) and _FIGTOOLS not in _sys.path:
    _sys.path.append(_FIGTOOLS)

import matplotlib as mpl
import matplotlib.pyplot as plt

MM = 1/25.4
SINGLE = 85*MM      # 3.35 in
ONEHALF = 130*MM
DOUBLE = 180*MM     # 7.09 in

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
def earth(): globals().update(EARTH)
def warm():  globals().update(WARM)
def calm():  globals().update(CALM)
def vivid(): globals().update(VIVID)

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


def apply():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Liberation Sans", "Arial", "Helvetica", "DejaVu Sans"],
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
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
        # Matplotlib's mathtext default is DejaVu Sans, so every $...$ expression was being
        # set in a DIFFERENT typeface from the body text — PyMuPDF found whole alphabets of
        # DejaVu embedded in fig4 and figS5. Point mathtext at the same face as everything else.
        "mathtext.fontset": "custom",
        "mathtext.rm": "Liberation Sans",
        "mathtext.it": "Liberation Sans:italic",
        "mathtext.bf": "Liberation Sans:bold",
        "mathtext.sf": "Liberation Sans",
        # cal/tt/bfit default to "cursive"/"monospace"/"sans:italic:bold". Any figure with a
        # log axis builds the whole custom fontset, resolves "cursive", finds nothing, and
        # silently falls back to DejaVu Sans - the same leak fixed above, by another door.
        "mathtext.cal": "Liberation Sans:italic",
        "mathtext.tt": "Liberation Sans",
        "mathtext.bfit": "Liberation Sans:italic:bold",
        "mathtext.default": "regular",
        "pdf.fonttype": 42,     # embed as TrueType so text stays editable in Illustrator
        "ps.fonttype": 42,
        "svg.fonttype": "none",
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
    """
    out = []
    for ax in fig.axes:
        letter = getattr(ax, "_panel_letter", "?")
        groups = (list(ax.lines) + list(ax.patches) + list(ax.collections)
                  + list(ax.images) + list(ax.texts))
        for art in groups:
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
            fontsize=9, fontweight="bold", va="top", ha="left")

def save(fig, stem):
    """Write both a vector file for the journal and a PNG for review.

    Also refuses to be quiet about an empty panel: see blank_panels(). The build treats the
    marker line below as a failure, so a missing input can no longer produce a clean log.
    """
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
