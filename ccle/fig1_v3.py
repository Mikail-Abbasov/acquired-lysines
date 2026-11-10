#!/usr/bin/env python3
"""
fig1.py, Figure 1: cancer mutation creates a large, driver-enriched class of new lysines.

  a  union assembly: three variant sources -> 295,411 acquired-lysine sites
  b  substitution ranking, E->K first of 363 observed pairs
  c  gnomAD allele-frequency bands, with "not assessed" shown separately
  d  driver enrichment among Hess-assessed positions
  e  mutational-signature composition of driver hotspots
  f  acquired-lysine burden per tumour by cancer type, with the uveal control

SCOPE, which the legend must state: panels a-c are union-wide (DepMap + TCGA + COSMIC) and
panel f is TCGA; panels d and e are DepMap only, because the Hess driver and signature
annotations exist only there.

Panels d and e reproduce the logic of driver_figure.py exactly rather than re-deriving it, a
driver enrichment depends entirely on how the denominator is defined.
"""
import os, re, sys
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch, Rectangle
from scipy import stats

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
# NOT scaled. 5 Nov: multiplying the body by 1.04882 to make the INK rather than the
# canvas measure 12 inches looked like a uniform enlargement and was not. The panel
# nudges _axes_dx_in, _axes_dw_in and _axes_dh_in are absolute inches and BAR_MAX_PT is
# absolute points, so none of them scaled with the type and the panels. Letters a, b, c
# and e came out 119, 180, 81 and 123 px left of where they belong while d, f and g
# stayed put. A real scale needs those four to scale too; the figure is 12 inches as it
# stands and the trimmed image is 11.47 only because the margins are not ink.
BODY = 7.0
st.apply(BODY)

# Panel a carries ONE type size for all its wording: region counts, source names, the union
# line and the footnote. 6.8 lands on the ladder's 7 pt rung, which is the largest that still
# fits the Venn crescents, and prints at 5.4 pt where 5.0 pt printed at 3.8.
PANEL_A_PT = 6.8
# Panel b carries ONE type size too. It was running four: the route names and percentages at
# 5.0, the column heads at 6.0, LYSINE and the two tallies at 7.0. 36 of its 42 strings sat
# at 5.0, which prints at 4.1 pt in this figure, the smallest type in the paper.
PANEL_B_PT = 6.8
# Panel a styling variant, for review only. 0 is what the manuscript has shipped.
# APPROVED 2026-11-04: 7, the diagonal wash across the shared territory. The other
# treatments stay reachable with PANELA=<n> for review and are not what ships.
_VA = int(os.environ.get("PANELA", "7"))
# ARROWS: 1 straight, which is what ships. 0 restores the old arc, 2 is a horizontal run
# into a sharp bend, which merged the two rust routes into one mass at this scale and sent
# the top outbound arrow past its own label. Straight gives the sharp reading without that.
_ARR = int(os.environ.get("ARROWS", "0"))


def _sharp(ann, z=None):
    """Square the ends of an arrow and miter its joins.

    FancyArrowPatch uses round caps and joins by default, which turns the tail of a thick
    arrow into a blob and softens every edge along the curve.
    """
    pa = getattr(ann, "arrow_patch", None)
    if pa is None:
        return ann
    try:
        pa.set_capstyle("butt"); pa.set_joinstyle("miter")
        if z is not None:
            pa.set_zorder(z)
    except Exception as e:
        print("arrow cap/join not settable: %s" % e)
    return ann


def _ramp(ax, patch, c0, c1, radial=True, n=600, angle=45.0):
    """Fill one Venn region with a fade, clipped to that region's own outline.

    Both ends of every ramp are palette colours, so the only values that are not in the five
    are the ones the fade passes THROUGH. A region is its own patch, so this is a clip, not a
    blend: no two fills are ever composited on top of each other.
    """
    import numpy as np
    from matplotlib.colors import LinearSegmentedColormap, to_rgb
    cm = LinearSegmentedColormap.from_list("_r", [to_rgb(c0), to_rgb(c1)], N=n)
    ext = patch.get_extents()
    inv = ax.transData.inverted()
    x0, y0 = inv.transform((ext.x0, ext.y0))
    x1, y1 = inv.transform((ext.x1, ext.y1))
    if radial:
        # Centre the glow on the REGION, not on its bounding box. A crescent's box centre
        # can sit outside the crescent, which is what put the hotspot off-centre in the
        # lens and turned the two slivers into bright dots.
        try:
            vtx = patch.get_path().vertices
            cx, cy = float(np.nanmean(vtx[:, 0])), float(np.nanmean(vtx[:, 1]))
            px, py = ax.transData.transform((cx, cy))
            fx = (px - ext.x0) / max(ext.width, 1e-9)
            fy = (py - ext.y0) / max(ext.height, 1e-9)
            fx = min(max(fx, 0.0), 1.0); fy = min(max(fy, 0.0), 1.0)
        except Exception:
            fx = fy = 0.5
        yy, xx = np.mgrid[0:n, 0:n]
        g = np.sqrt(((xx / float(n) - fx) / 0.72) ** 2 + ((yy / float(n) - fy) / 0.72) ** 2)
        g = np.clip(g, 0, 1)
    else:
        th = np.radians(angle)
        yy, xx = np.mgrid[0:n, 0:n]
        u = (xx / float(n)) * np.cos(th) + (yy / float(n)) * np.sin(th)
        g = (u - u.min()) / max(u.max() - u.min(), 1e-9)
    im = ax.imshow(g, extent=[x0, x1, y0, y1], origin="lower", cmap=cm, vmin=0, vmax=1,
                   zorder=patch.get_zorder() + 0.1, interpolation="bilinear", aspect="auto")
    im.set_clip_path(patch)
    return im


def find(name):
    for d in (".", "..", _TOP,
              os.path.join(_TOP, 'ccle_plex05'),
              os.path.join(_TOP, 'Datasets'),
              os.path.join(_TOP, 'Datasets/pandepmap_out')):
        p = os.path.join(d, name)
        if os.path.exists(p): return p
    return None

ALLM   = find("union_all_missense.tsv")
SUBM   = find("union_subst_matrix.csv")
UNIONG = find("union_sites_gnomad_v2.tsv")
PARQ   = find("depmap_variants_annotated.parquet")
MAF    = find("OmicsSomaticMutations.csv")
BURD   = find("tcga_burden_by_lineage.csv")
notes  = []

# GRID REBUILT 2026-09-03. The old 2x3 gave every panel 1.73 in. Two of the panels now here
# cannot live in that: the codon route map needs ~3.7 in for its route labels, and the
# 32-lineage burden needs ~4.2 in for its tick labels. Both were measured, not guessed — an
# earlier attempt to squeeze the route map into a 1.73 in cell put its labels at ~2.5 pt.
# So: 3 rows x 12 columns, with the two wide panels spanning 8 and 9 columns.
#
#   row 0   a  union assembly (0:4)          b  codon routes into and out of lysine (4:12)
#   row 1   c  gnomAD bands (0:4)            d  driver enrichment (4:8)   e  signatures (8:12)
#   row 2   f  burden across every lineage (0:9)        g  germline vs somatic (9:12)
# HEADLINE NUMBERS, added 2026-09-04. Until now this script printed only "wrote fig1.pdf",
# so every number in its legend would have had to be copied out of FIGURE_ARCHITECTURE.md.
# That is the failure mode that produced A36's 8,530 sites, Figure 8's non-existent 0.633/0.384
# and Figure 6's 331/479/372 Da. Each panel now files its headline quantities here and the
# script prints them, so the legend is written from build output and a changed input shows up.
HEAD = {}

# Row 0 is taller, and the figure with it. Panels A and B set their line gaps in data or
# axes units while the type is fixed in points, so a taller row widens every gap in both
# panels at once and leaves every coordinate alone. Raising the gaps instead would have
# worked on the left of B and failed on the right: that column holds seven entries at a
# 1.10 pitch, so a 0.56 within-pair gap leaves 0.54 between pairs and the fourteen lines
# read as one block. Rows 1 and 2 keep their old inches, which is what 1.05 -> 1.13
# against height_ratios 1.38 -> 1.63 buys.
# 5 Nov: two rows on a 12 inch canvas at 300 dpi, a b c d over e f g. The panels are NOT
# resized. Each cell is set to the frame that panel already had in the three-row figure,
# measured in display px at dpi 300, so every width, every x limit and every margin inside
# a panel is the one that was signed off:
#   a 602.46   b 1263.57   c 559.39   d 410.02   e 559.39   f 1313.36   g 298.34
# Heights 665.78 (a, aspect-locked), 907.79 (b) and 465.53 (the five bar panels).
# savefig.dpi drops from 600 to 300 so the PNG is 3900 px wide and one display pixel is
# one output pixel, which is also what every measurement in this file is counted in.
plt.rcParams["savefig.dpi"] = 300
# Canvas trimmed to the ink, 5 Nov. The pads were L 27 R 147 T 80 B 33, so the right and
# top carried 120 and 47 px of nothing and the figure sat off-centre on its own sheet.
# 18 px all round now. The CONTENT box keeps its pixel size, 3300 x 1250.2, and each
# panel's overhang past it is unchanged, so no panel moves relative to any other.
fig = plt.figure(figsize=(11.54, 5.01))
# hspace and the two margins are solved together so the cells land on 907.79 and
# 465.53 exactly. c and d keep their own height at the top of a row sized for b, so
# 442 px of their cell is empty below them; a large hspace on top of that opened a
# white band across the middle of the figure.
gs = fig.add_gridspec(2, 1, hspace=0.2122,
                      left=0.0407279, right=0.9939342, top=0.9415169, bottom=0.1097139,
                      height_ratios=[665.78, 465.53])

# Rows 1 and 2 get their own column spacing. At wspace=0.55 the gutter between two
# side-by-side panels is 61 px, while a panel's y-axis label plus its tick numbers need
# 107 px (measured: c|d 107, d|e 100, f|g 108). Every one of those labels was therefore
# being drawn 46 px INSIDE its left neighbour's plot frame. text_collisions() cannot see
# this — the labels land in whitespace, so no two text boxes intersect — which is why it
# survived the consistency pass. figtools/intrude.py is the check that does see it.
#
# Row 0 keeps the original spacing on purpose. Panels a and b are settled, and panel b's
# arrowhead clearances were tuned at its current width; a global wspace change would
# narrow it by 5% and put those back in play.
GUT = 2.6    # gutter 129 px against the 108 px needed, so ~3 pt of printed clearance.
             # 1.65 (gutter 110) cleared the collision but left under 1 pt, thin enough
             # that the panel crop still caught a sliver of the neighbour's label.
             # Costs panels c-e 7% of their width and panel g 11%; both have the room.
# width_ratios are the panels' own frame widths, so the cells reproduce them exactly.
# wspace is a fraction of the AVERAGE axes width within its own row, so the two rows need
# different numbers. Both rows run the full 3300 px of content, which puts a and e on the
# left edge and d and g on the right. The bottom row carries three panels against the top
# row's four, so its gutters are wider; that is where its spare width goes.
gs_top = gs[0, 0].subgridspec(1, 4, wspace=0.19876,
                              width_ratios=[602.46, 1300.00, 559.39, 410.02])
gs_bot = gs[1, 0].subgridspec(1, 3, wspace=0.77989,
                              width_ratios=[559.39, 1313.36, 298.34])

# ------------------------------------------------------------------ a
ax = fig.add_subplot(gs_top[0, 0]); st.head(fig, ax, "a")
# Letter flush to this panel's leftmost ink rather than to its axes frame: the venn
# has no frame to read, and panel b's heading overhangs its cell by 29 px.
ax._panel_head_flush = True
# Two different left edges are in play and they are not the same number. Panel c's BLOCK
# starts at 11.6, where its rotated y-axis label sits; its FRAME starts at 148.4. Aligning
# panel a on the frame left its block indented against c's, so the axes moves out to the
# block edge instead: 170.2 px left of where the venn was drawn, putting a's leftmost ink
# on 11.6 with c and f. The letter must not follow it out there, so it is pushed back to
# the block edge too, so the letter marks the panel's left edge rather than floating 137 px
# inside it. figtools/headrow.py prints the letter column, so if this constant goes stale
# the misalignment is reported rather than silently shipped.
# The same alignment this offset always did, re-measured for the new grid. a leads the top
# row and e the bottom one, so their letters carry the figure's left edge. The venn is
# centred in its cell and aspect-locked, so its ink stops 26 px inside the frame and the
# "A" sat that far right of the "E". Shifting the axes moves the diagram, its caption and
# its letter together; the panel is not resized.
# Re-measured 5 Nov after "union 295,411 acquired lysines" went to the caption. That line
# was centred under the venn and reached further left than the diagram does, so taking it
# off moved this panel's leftmost ink 16.4 px right and the "A" with it.
# Hand positioning from the PI's own edit of 5 Nov, measured off that image rather than
# guessed. His top row has an even ~95 px of air between ink blocks; the script's had
# 200 px after a and b's labels OVERLAPPING c's by 42, which is why that row renders as
# one merged block. Shifts: a -79, b -179, c -75, d 0. The a and e letters go to the
# figure's left edge, which is where the two row-leading letters mark the margin.
ax._axes_dx_in = -121.4 / 300.0
ax._panel_dx_in = -44.1 / 300.0
# THIS PANEL USED TO LIE, QUIETLY. Until 2026-09-02 the three circles had areas proportional
# to each source's count but POSITIONS CHOSEN BY HAND: the overlaps a reader saw were drawn,
# not computed. They are computed now — matplotlib-venn solves the layout from the seven real
# region counts — and the true structure is not what the old picture implied: DepMap sits
# almost entirely inside COSMIC, and the DepMap-with-TCGA-but-not-COSMIC region is 147 sites,
# not the wedge the old diagram showed.
if ALLM:
    u = pd.read_csv(ALLM, sep="\t", usecols=["gene", "alt_aa", "in_depmap", "in_tcga", "in_cosmic"])
    k = u[u.alt_aa == "K"]
    n_src = k[["in_depmap", "in_tcga", "in_cosmic"]].sum(axis=1)
    multi = int((n_src > 1).sum())
    HEAD["a_union_K"] = len(k); HEAD["a_all_missense"] = len(u)
    HEAD["a_multi_source"] = multi
    HEAD["a_genes"] = int(k.gene.nunique()) if "gene" in k.columns else None
    from matplotlib_venn import venn3, venn3_circles
    D, T, C = k.in_depmap.astype(bool), k.in_tcga.astype(bool), k.in_cosmic.astype(bool)
    reg = (int((D & ~T & ~C).sum()), int((~D & T & ~C).sum()), int((D & T & ~C).sum()),
           int((~D & ~T & C).sum()), int((D & ~T & C).sum()), int((~D & T & C).sum()),
           int((D & T & C).sum()))
    # Three translucent circles blend into six secondary colours that are not in the palette,
    # and a 0.42 blend of teal and rust is not a value anyone chose. Fill every region with the
    # one neutral instead and let the OUTLINES carry the set identity, at full opacity.
    # The seven region counts belong in the ledger: four of them are drawn in the panel and
    # three are too small to label, and a caption that stands alone needs all seven.
    for _k, _n in zip(("depmap_only", "tcga_only", "depmap_tcga", "cosmic_only",
                       "depmap_cosmic", "tcga_cosmic", "all_three"), reg):
        HEAD["a_region_" + _k] = int(_n)
    v = venn3(subsets=reg, set_labels=("DepMap", "TCGA", "COSMIC"), ax=ax,
              set_colors=(st.LIGHTGREY, st.LIGHTGREY, st.LIGHTGREY), alpha=1.0)
    # Every region is its OWN patch at alpha 1.0, so there are no blends to avoid and each
    # can take a solid palette colour. The three exclusive territories carry their source's
    # colour; everything shared stays neutral, which is the panel's whole point.
    if _VA == 2:        # each exclusive territory in its own source colour
        _FILL = {"100": st.BLUE, "010": st.GREEN, "001": st.ORANGE}
    elif _VA == 3:      # the SHARED territory, which is what the 32% claim is about
        _FILL = {"110": st.BLUE, "101": st.BLUE, "011": st.BLUE, "111": st.BLUE}
    else:
        _FILL = {}
    _SHARED = ("110", "101", "011", "111")
    _SRC = {"100": st.BLUE, "010": st.GREEN, "001": st.ORANGE}
    for key in ("100", "010", "001", "110", "101", "011", "111"):
        pt = v.get_patch_by_id(key)
        if pt is None:
            continue
        pt.set_facecolor(_FILL.get(key, st.LIGHTGREY)); pt.set_alpha(1.0)
        pt.set_edgecolor("none")
        _e = pt.get_extents()
        _tiny = (_e.width * _e.height) < 0.004 * (ax.bbox.width * ax.bbox.height)
        if _VA in (4, 6, 7, 8, 9, 10, 11):
            # the shared territory carries the treatment, everything else lies flat
            if key in _SHARED:
                if _VA != 4 and _tiny:
                    pt.set_facecolor(st.BLUE)
                elif _VA == 7:                       # diagonal wash
                    pt.set_facecolor(st.LIGHTGREY)
                    _ramp(ax, pt, st.BLUE, st.LIGHTGREY, radial=False, angle=35.0)
                elif _VA == 8:                       # teal to rust, both parents named
                    pt.set_facecolor(st.LIGHTGREY)
                    _ramp(ax, pt, st.BLUE, st.ORANGE, radial=False, angle=55.0)
                elif _VA == 9:                       # halo: pale core, saturated rim
                    pt.set_facecolor(st.LIGHTGREY); _ramp(ax, pt, st.LIGHTGREY, st.BLUE)
                else:                                # 4, 6, 10, 11: centred glow
                    pt.set_facecolor(st.LIGHTGREY); _ramp(ax, pt, st.BLUE, st.LIGHTGREY)
        elif _VA == 5:
            # every territory fades from its own source outward
            pt.set_facecolor(st.LIGHTGREY)
            c0 = st.BLUE if key in _SHARED else _SRC.get(key, st.LIGHTGREY)
            _ramp(ax, pt, c0, st.LIGHTGREY)
    _lw = 0.0 if _VA == 10 else (0.6 if _VA == 11 else 1.0)
    for circ, col in zip(venn3_circles(subsets=reg, ax=ax, lw=max(_lw, 0.01)),
                         (st.BLUE, st.GREEN, st.ORANGE)):
        circ.set_edgecolor(col)
        if _VA == 10:
            circ.set_visible(False)
    for t, col in zip(v.set_labels or [], (st.BLUE, st.GREEN, st.ORANGE)):
        if t: t.set_fontsize(st.sz(PANEL_A_PT)); t.set_color(col); t.set_fontweight("bold")
    # venn3 drops the three set names wherever its own layout lands them, so DepMap, TCGA and
    # COSMIC each sat a different distance from their circle. Put every name the same gap
    # outside its own circle's edge, on the bearing that points away from the other two.
    _NAMEPOS = []
    # Name-to-count gap, in the venn's own data units. Measured, not guessed: at 0.072
    # the two glyph boxes OVERLAPPED by 2.7 px (name 2408.6..2442.6, count 2377.3..2411.3),
    # a centre distance of 31.3 px for 34 px type, which is tighter than single spacing.
    # 434.7 px per unit here, so 0.106 buys 46 px, about 1.35 lines.
    _LINE = 0.106
    try:
        _cen, _rad = v.centers, v.radii
        _GAP = 0.055                       # of the axes, clear of the circle's edge
        # Each name is centred on ITS OWN circle and sits directly above it, or below it for
        # COSMIC, whose top is buried under the other two. Bearings put each name at a
        # different angle, so none of the three sat over the circle it names.
        _ABOVE = (True, True, False)       # DepMap, TCGA, COSMIC
        for _t, _c, _r, _up in zip(v.set_labels or [], _cen, _rad, _ABOVE):
            if _t is None:
                continue
            _cx = float(_c[0]) if not hasattr(_c, "x") else float(_c.x)
            _cy = float(_c[1]) if not hasattr(_c, "y") else float(_c.y)
            _r = float(_r)
            # the count hangs one line under the name, so a block above a circle starts a
            # line higher or the count lands on the ring
            _py = (_cy + _r + _GAP + _LINE) if _up else (_cy - _r - _GAP)
            _t.set_position((_cx, _py))
            _t.set_ha("center"); _t.set_va("center")
            _NAMEPOS.append((_cx, _py))
    except Exception as _e:
        print("set-label placement skipped: %s" % _e)
    import os
    if os.environ.get("VENN_PROBE"):
        for i, key in enumerate(("100", "010", "001", "110", "101", "011", "111")):
            t = v.get_label_by_id(key)
            if t is not None: print("region", key, t.get_position(), repr(t.get_text()))
        for t in (v.set_labels or []):
            if t: print("setlabel", repr(t.get_text()), t.get_position())
        print("xlim", ax.get_xlim(), "ylim", ax.get_ylim())
    # Region counts span three orders of magnitude. Label the four a reader can act on and
    # drop the rest rather than printing numbers into lens regions too small to hold them.
    KEEP = {"100": True, "010": True, "001": True, "011": True}
    # DepMap-only and TCGA-only are thin crescents: their counts are wider than the region
    # that holds them, so those two sit outside on a leader instead of straddling an edge.
    # (text x, text y, alignment, leader-start dx, leader-start dy in text-relative units)
    OUT = {"100": (-0.645, 0.275, "right", 0.012, 0.0),
           "010": (-0.050, 0.610, "center", 0.0, -0.030)}
    # 167,321 comes off the grey fill onto paper. Rust on C5C3C6 measured 2.8:1 against the
    # pixels actually behind the glyphs, below the 3:1 floor for large text, and this prints
    # at 5.7 pt. Placed from COSMIC's own circle rather than a hardcoded pair, so it stays
    # put if the geometry changes. No leader: the rust is what ties it to COSMIC.
    try:
        _cc, _cr = v.centers[2], v.radii[2]
        _ccx = float(_cc[0]) if not hasattr(_cc, "x") else float(_cc.x)
        _ccy = float(_cc[1]) if not hasattr(_cc, "y") else float(_cc.y)
        _bb = np.radians(318.0)
        OUT["001"] = (_ccx + (float(_cr) + 0.050) * np.cos(_bb),
                      _ccy + (float(_cr) + 0.050) * np.sin(_bb), "left", 0.0, 0.0)
    except Exception as _e:
        print("001 placement skipped, count stays on the fill: %s" % _e)
    for key in ("100", "010", "001", "110", "101", "011", "111"):
        t = v.get_label_by_id(key)
        if t is None: continue
        if key in KEEP:
            t.set_fontsize(st.sz(PANEL_A_PT))
            # Every count takes its own source's colour, which is the house rule for value
            # labels and is also what lets the leader lines go: a count and its database
            # name now carry the same colour, so the colour says which region it belongs to.
            # The shared count sits ON the wash, where its sources are teal and rust both,
            # so it stays white for contrast rather than claiming one of them.
            _CNT = {"100": st.BLUE, "010": st.GREEN, "001": st.ORANGE, "011": st.PAPER}
            t.set_color(_CNT.get(key, st.GREY))
            t.set_text(f"{int(t.get_text()):,}")
            # Each exclusive count sits directly under its own database name, in that
            # database's colour, so the pair reads as one block and nothing needs a leader.
            _UNDER = {"100": 0, "010": 1, "001": 2}
            if key in _UNDER and len(_NAMEPOS) == 3:
                _nx, _ny = _NAMEPOS[_UNDER[key]]
                t.set_position((_nx, _ny - _LINE)); t.set_ha("center"); t.set_va("center")
            elif key in OUT:
                x0, y0 = t.get_position()
                dx, dy, ha, ox, oy = OUT[key]
                t.set_position((dx, dy)); t.set_ha(ha)
                # No leader. The count is coloured to its source instead, so the line that
                # used to point at the crescent is doing nothing the colour does not.
        else:
            t.set_visible(False)
    # venn3 fits its layout to the axes, so the COSMIC set label lands hard against the
    # bottom edge and collided with the caption. Open space below the diagram instead of
    # moving the caption further away from the panel it belongs to.
    _y0, _y1 = ax.get_ylim()
    # 0.057 aligns this panel's topmost ink with panel b's. The set labels sit OUTSIDE the
    # circles, so they overflow whatever headroom the data range leaves: at 0.03 the TCGA
    # label stood 25.9 px above the axes top and so 15.9 px above panel b's heading, at
    # 0.08 it fell 13.3 px below it. Interpolated and checked: TCGA now tops out at 2562.0
    # against b's 2561.9, and the two letter-to-ink gaps read 6.3 and 6.5 px.
    _top = 0.13 if _VA == 0 else 0.057
    _bot = 0.11          # was 0.20, which made room for a footnote that now lives in the caption
    ax.set_ylim(_y0 - _bot * (_y1 - _y0), _y1 + _top * (_y1 - _y0))
    # "union  295,411 acquired lysines" moved to the caption on 2026-11-05, the same move
    # the footnote below made on 11-04. The count is still filed in the ledger as a_union_K
    # and the caption now names it as the UNION of the three sources, which is the word that
    # was doing the work here: 295,411 is not the sum of the three catalogues.
    HEAD["a_union_label"] = len(k)
    # The footnote moved to the caption on 2026-11-04. In the artwork it was 2.56 in wide
    # against a 1.45 in diagram, so it set panel a's width, reached into panel b's column,
    # and widened the composite, which shrinks every label in Figures 1 and 5 at print.
# venn3 sets an equal aspect, so matplotlib shrinks this axes to the diagram's own shape
# and, anchored "C", centres the shrunken box in the cell. That put panel a's content 111 px
# below its letter while panel b's sat 6.5 px below the same letter baseline. Anchoring north
# spends the slack at the bottom, where the caption lives, instead of above the diagram.
ax.set_anchor("N")
ax.set_axis_off()

# ------------------------------------------------------------------ b
# The full codon route map, promoted from figS2 panel b. Spans 8 columns (~4.3 in), which is
# what its route labels need. The compact variant built for the old 1.73 in cell is retired.
ax = fig.add_subplot(gs_top[0, 1]); st.head(fig, ax, "b")
ax._axes_dx_in = -179.0 / 300.0
# Letter flush to this panel's leftmost ink rather than to its axes frame: the venn
# has no frame to read, and panel b's heading overhangs its cell by 29 px.
ax._panel_head_flush = True

ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")

# C/T>A/G, not "C>A/G, T>A/G". Lysine has two codons, so Asn either way runs through
# C>A, C>G, T>A and T>G, which is exactly {C,T} to {A,G}. The collapsed form is the same
# four routes and says what they have in common, pyrimidine to purine. It is also 140 px
# wide against 233, which is what lets panel b come back down from 1560 px.
# Three-letter codes, 5 Nov. "asparagine  N" is 229.6 px and "Asn" is about 75, and that
# is what buys the codon type back onto the line. With the full names a one-line entry
# wanted 598.5 px of left column against 361.7, which is 636 px of extra panel, half the
# figure. With three-letter codes it is about 300, which the top row's gutters can pay.
_IN = [("Glu",      74.3, "C>T",          True),
       ("Arg",      11.2, "C>T",          True),
       ("Gln",       5.5, "C>A",          False),
       ("Asn",       5.4, "C/T>A/G"     , False),
       ("Thr",       2.4, "C>A",          False),
       ("Met/Ile",   1.2, "T>A",          False)]
_OUT = [("Asn", 45.5, "C/T>A/G"     ), ("Thr", 16.0, "T>G"),
        ("Glu", 11.2, "T>C"),          ("Arg", 10.8, "T>C"),
        ("Gln",  6.3, "T>G"),          ("stop", 4.8, "T>A"),
        ("Met/Ile", 5.4, "T>A")]
# Recentred on the one-line stack, 8.75 to 2.60, and 1.60 units tall rather than the 1.20
# the two strings need. The six inbound heads land spread down the left edge, and at 1.20
# that edge was 61 px for heads wanting 126, so they merged into one blob.
KX0, KX1, KY0, KY1 = 4.17, 5.73, 4.87, 6.47
# The box sits IN FRONT of every arrow: arrowheads used to ride over its outline and the
# outbound tails over its right edge. zorder does the work, and the white fill under the
# wash means an arrow that strays beneath it is hidden rather than showing through.
_BOX = FancyBboxPatch((KX0, KY0), KX1-KX0, KY1-KY0,
                      boxstyle="round,pad=0.05,rounding_size=0.06",
                      facecolor=st.PAPER, edgecolor=st.GREEN, lw=1.3, zorder=6)
ax.add_patch(_BOX)
# A wash across the box, the same move as panel a, running light grey to paper so the two
# dark strings on it keep their contrast. Both ends are palette values.
_ramp(ax, _BOX, st.LIGHTGREY, st.PAPER, radial=False, angle=35.0).set_zorder(6.1)
ax.text((KX0+KX1)/2, KY0+(KY1-KY0)*0.65, "Lys", ha="center", va="center", fontsize=st.sz(PANEL_B_PT),
        fontweight="bold", color=st.GREEN, zorder=7)
ax.text((KX0+KX1)/2, KY0+(KY1-KY0)*0.24, "AAA   AAG", ha="center", va="center", fontsize=st.sz(PANEL_B_PT),
        color=st.INK, zorder=7)
_lw = lambda sh: 0.5 + 0.042 * sh      # 74.3% -> 3.6 pt, was 5.9
# Every inbound arrow used to end on the SAME point of the box's left edge, so all six
# arrowheads piled up on each other, the widest one burying the rest. They now land spread
# down that edge, in the order their routes are listed, and still read as converging on it.
# The box's left edge is 1.70 data units, about 25 pt on the page. Six heads at
# mutation_scale 7 want roughly 42 pt of it, so they cannot all fit however they are spread:
# the heads come down to 5 and use the whole edge, which gives about 5 pt each.
# Measured: clipping the tails to a vertical half-plane at 2.98 removes nothing, because the
# arrows already start on that line and leave it close to horizontal, so the butt cut is
# already vertical. The ragged edge was at the OTHER end.
_AY = np.linspace(KY1 - 0.14, KY0 + 0.14, len(_IN))
_HEAD_IN = 4      # 88 px of edge for six heads; 5 wanted 126
# ONE LINE PER ENTRY, 5 Nov. Two lines each was what made this panel 907.8 px tall and so
# set the height of the whole top row. The diagram is unchanged: same two columns, same box,
# same arrows with thickness carrying the share. Only the entry folds onto one line, name in
# its own right-aligned column and the value in another, so nothing shifts as values change.
# The codon classes come off the artwork. Keeping them on one line needs 598.5 px of left
# column against the 361.7 px there is, and 619.1 px on the right against 364.5, which is
# 636 px of extra panel width, half the figure. The C>T routes are still orange, and the two
# footers still count them, 4 of 14 in and none of 16 out.
# Three columns, not two. Right-aligning "5.4%  C>A/G, T>A/G" as one string put that
# row's number hard against its name while every other row had a gap, because the
# string is the widest and right alignment pushes its START left. Name right-aligned,
# number right-aligned so the decimals line up, class left-aligned.
# One run-on line per entry with bullet separators, not three aligned columns, and both
# sides justified to the SAME vertical extent: six entries on the left and seven on the
# right both run 8.75 down to 2.60, so the two stacks start and finish level and the box
# sits centred between them. The left lines are right-aligned and the right lines left-
# aligned, which puts the flush edge against the arrows on each side and leaves the rag
# facing outwards.
_LX, _RX = 2.87, 7.02
for _i, ((lab, sh, cls, ct), yy) in enumerate(zip(_IN, np.linspace(8.75, 2.60, len(_IN)))):
    col = st.ORANGE if ct else st.LIGHTGREY
    # The box carries boxstyle pad=0.05, so its drawn edge IS at KX0-0.05: an arrow ending
    # there has its tip on the outline, and with the box in front at zorder 6 the tip is cut
    # by it. Clearance measured on the page rather than guessed: at 0.34 the heads stood
    # 8.8 to 9.7 pt printed off the box, wider than the 7.2 pt type beside them. 1 data unit
    # is 28.0 pt printed here, so 0.11 leaves a nominal 0.06 units and, once shrinkB's 1 pt
    # is added at this end, about 2.5 pt on the page.
    _ann = ax.annotate("", xy=(KX0-0.11, _AY[_i]), xytext=(2.99, yy),
                arrowprops=dict(arrowstyle="-|>", color=col, lw=_lw(sh), shrinkA=0, shrinkB=1,
                                mutation_scale=_HEAD_IN,
                                connectionstyle=(
                                    "arc3,rad=0" if _ARR == 1 else
                                    "angle,angleA=0,angleB=%.0f,rad=0"
                                    % (-58.0 if yy > 5.4 else 58.0) if _ARR == 2 else
                                    "arc3,rad=%.2f" % (-0.11 if yy > _AY[_i] else 0.11))))
    _sharp(_ann, z=3)
    ax.text(_LX, yy, "%s \u2022 %.1f%% \u2022 %s" % (lab, sh, cls.strip()),
            ha="right", va="center", fontsize=st.sz(PANEL_B_PT),
            color=st.ORANGE if ct else st.GREY)
_TY = np.linspace(KY1 - 0.14, KY0 + 0.14, len(_OUT))
for _j, ((lab, sh, cls), yy) in enumerate(zip(_OUT, np.linspace(8.75, 2.60, len(_OUT)))):
    # Tails leave the box with the same 2.5 pt of printed clearance the inbound heads stop
    # short of it. The constant differs because shrinkB bites at xy, which is the FAR end on
    # this side: no shrink at the box, so 0.14 here matches 0.11 there on the page.
    _ann = ax.annotate("", xy=(6.90, yy), xytext=(KX1+0.14, _TY[_j]),
                arrowprops=dict(arrowstyle="-|>", color=st.BLUE, lw=_lw(sh), shrinkA=0,
                                shrinkB=1, mutation_scale=7,
                                connectionstyle=(
                                    "arc3,rad=0" if _ARR == 1 else
                                    "angle,angleA=%.0f,angleB=0,rad=0"
                                    % (58.0 if yy > 5.4 else -58.0) if _ARR == 2 else
                                    "arc3,rad=%.2f" % (0.11 if yy > _TY[_j] else -0.11))))
    _sharp(_ann, z=3)
    ax.text(_RX, yy, "%s \u2022 %.1f%% \u2022 %s" % (lab, sh, cls.strip()),
            ha="left", va="center", fontsize=st.sz(PANEL_B_PT), color=st.BLUE)
ax.text(2.87, 9.55, "what becomes a lysine", ha="right", va="center", fontsize=st.sz(PANEL_B_PT),
        fontweight="bold", color=st.INK)
ax.text(7.02, 9.55, "what a lysine becomes", ha="left", va="center", fontsize=st.sz(PANEL_B_PT),
        fontweight="bold", color=st.INK)
# Both summaries hang from a shared BASELINE at y = _SUMY, anchored va="bottom". They used
# to be va="top" at 1.70 and 1.45, which put their bottoms 0.25 units apart. Anchoring the
# bottom rather than matching the tops keeps them aligned if either ever gains a line.
# Justifying both stacks to 2.60 brought the last entry down onto these, 18 px apart.
_SUMY = 0.42
ax.text(2.87, _SUMY, "14 codon routes\n4 of them C>T", ha="right", va="bottom",
        fontsize=st.sz(PANEL_B_PT), color=st.ORANGE, fontweight="bold", linespacing=1.5)
ax.text(7.02, _SUMY, "16 codon routes\nnone of them C>T", ha="left", va="bottom",
        fontsize=st.sz(PANEL_B_PT), color=st.BLUE, fontweight="bold", linespacing=1.5)
# The AAA/AAG explanation and the rust key moved to the caption on 2026-11-04. Both were
# prose, and the caption is where this lab's rule puts prose; removing them also frees the
# centre of the panel, which is what let every label go up to one size.


# ------------------------------------------------------------------ c
ax = fig.add_subplot(gs_top[0, 2]); st.head(fig, ax, "c")
ax._axes_dx_in = -79.0 / 300.0
# Keeps its own 465.53 px height in a row sized for b, top held so the
# four letters of this row stay on one baseline.
ax._axes_dh_in = -200.25 / 300.0
ax._panel_head_flush = True   # letter on the block edge, with panel a above it
if UNIONG:
    g = pd.read_csv(UNIONG, sep="\t", usecols=["band"])
    order = ["absent", "lt_1e-5", "1e-5_1e-4", "1e-4_1e-3", "gt_1e-3", "not_assessed"]
    # One line per band, set at an angle. The two-line form stacked "10^-5-" over "10^-4"
    # and "not" over "assessed", which read as eight labels for six bars.
    # "to", not a hyphen. Closing a mathtext group leaves trailing advance and the ASCII
    # hyphen that followed sat at superscript height, so "10^-5-10^-4" read as
    # "10^-5  ^-10^-4": the separator looked like the second number's minus sign, and it
    # was a different glyph from the mathtext minus in the exponents.
    show = ["absent", "<10$^{-5}$", "10$^{-5}$ to 10$^{-4}$", "10$^{-4}$ to 10$^{-3}$",
            ">10$^{-3}$", "not assessed"]
    v = np.array([(g.band == b).sum() for b in order])
    # A five-step lightness ramp was tried here and rejected on sight: the one chromatic step
    # in the middle of it pulled the eye to the 10^-5-10^-4 bin, which is not the point. The
    # point is the annotation, so colour carries the annotation: teal = the two rare bins that
    # make up the 91%, background grey = the three common bins, rust = never assessed.
    cols = [st.BLUE, st.BLUE, st.LIGHTGREY, st.LIGHTGREY, st.LIGHTGREY, st.ORANGE]
    ax.bar(range(len(v)), v / 1000, color=cols, width=st.BARW,
           edgecolor=st.GREY, linewidth=0.5)
    ax.set_xticks(range(len(v)))
    # rotation_mode="anchor" with ha="right" pins each label's right end to its own tick, so
    # the angled strings stay on their bars instead of drifting left as they lengthen.
    # No explicit size. sz(5.4) snapped to the bottom rung and printed these at 5.99 pt
    # while every other tick label in the figure printed at 7.19, which is the house
    # xtick.labelsize. The size was chosen when these were two-line horizontal labels.
    ax.set_xticklabels(show, rotation=45, ha="right", rotation_mode="anchor")
    ax.set_ylabel("sites (thousands)")
    # Ticks every 40 rather than matplotlib's 25: five labels instead of seven, and the top
    # one is a round 160 just above the 154,188 of the tallest bar.
    ax.set_ylim(0, 160); ax.set_yticks([0, 40, 80, 120, 160])
    assessed = v.sum() - v[-1]; som = v[0] + v[1]
    HEAD["c_bands"] = dict(zip(order, [int(x) for x in v]))
    HEAD["c_pct_somatic_of_assessed"] = 100 * som / assessed
    HEAD["c_pct_not_assessed"] = 100 * v[5] / v.sum()
    # Body size, matching the in-panel text of panels a and b. sz(6) snapped to 6.51.
    ax.text(0.60, 0.97, "%.1f%% of assessed sites\nabsent or <10$^{-5}$"
            % (100 * som / assessed), transform=ax.transAxes, ha="center", va="top",
            fontsize=st.sz(7), color=st.GREY, linespacing=1.4)
    ax.text(5, v[5]/1000 + v.max()/1000*0.05, "%.1f%%" % (100*v[5]/v.sum()),
            ha="center", fontsize=st.sz(7), color=st.ORANGE)
else:
    notes.append("1c: union_sites_gnomad_v2.tsv not found")

# ------------------------------------------------------------------ d  (DepMap scope)
# Four columns. Two bars cannot fill a 134 pt axes: squeezing the side margins to 8 pt
# drives the bars 115 pt apart, and holding them 42 pt apart leaves 44 pt empty at each
# edge. Three columns was tried and gives a 93.3 pt axes, which the 96 pt OR line
# overflows, so the annotation's width is the floor on how narrow this cell can go.
ax = fig.add_subplot(gs_top[0, 3]); st.head(fig, ax, "d")
# Keeps its own 465.53 px height in a row sized for b, top held so the
# four letters of this row stay on one baseline.
ax._axes_dh_in = -200.25 / 300.0
# Narrowed and recentred in its cell. Widening the x range was what brought the two bars
# together, because fit_bars holds the bar thickness on the page, but it bought that with
# 146.5 px of dead air inside each end of the frame, the largest end margins in the figure
# and the opposite of what was asked for. Taking the width out of the axes instead keeps
# the bars exactly where they were approved and brings both ends to 30 px, the same margin
# panel g carries. The shift keeps the narrowed panel centred between c and e rather than
# leaving a hole on one side. What stops it reaching 30 px is the OR annotation, which is
# 399.8 px wide on one line. At a 326 px frame it hung 37 px past each side and landed on
# the panel letter at 20% and on the 1.2 tick at 13%, so the frame is held at 410 px, just
# wide enough to carry it. That leaves 71.8 px at each end, down from 146.5. Putting the
# statistic on two lines is what would buy the rest.
# [cell is now the frame; offset retired with the 3-row grid] ax._axes_dw_in = -156.0 / 300.0
# [cell is now the frame; offset retired with the 3-row grid] ax._axes_dx_in = 74.7 / 300.0
# Letter to the block edge, as a, c, f and g already do. It used to sit at the frame's left
# corner, and narrowing the frame brought the OR annotation's left end to within 35 px of
# it, 17% over the "D". Flush to the panel's own leftmost ink, the y-axis label, it clears.
ax._panel_head_flush = True
if PARQ:
    v = pd.read_parquet(PARQ)
    m = v[v.vclass == "missense"].dropna(subset=["alt_aa"]).copy()
    m["hess"] = m.hess_driver.astype(str).str.lower().isin(["true", "yes", "y", "1"])
    ass = m[m.hess_driver.notna()]
    K = ass[ass.alt_aa == "K"]; O = ass[ass.alt_aa != "K"]
    kr = 100 * K.hess.sum() / len(K); orr = 100 * O.hess.sum() / len(O)
    OR, P = stats.fisher_exact([[K.hess.sum(), len(K) - K.hess.sum()],
                                [O.hess.sum(), len(O) - O.hess.sum()]])
    # 6 Nov: the second bar moves in from 1.0 to 0.851, bringing the centres from 199.8 px
    # to 170.0 px. The pair is then re-centred in the frame through xlim below, so the two
    # end margins stay equal and the frame itself does not change width. Positions come from
    # _DX everywhere, the same way panels e and g of figure 2 and panel g here already do it,
    # so there is one place to change and the ticks cannot drift off the bars.
    _DX = [0, 0.851]
    ax.bar(_DX, [orr, kr], color=[st.LIGHTGREY, st.BLUE], width=st.BARW)
    for _x, val in zip(_DX, [orr, kr]):
        # Body size, matching the in-panel text of panels a to c. sz(6.2) snapped to 7.19.
        ax.text(_x, val * 1.03, "%.3f%%" % val, ha="center", fontsize=st.sz(7))
    # No explicit size: the house xtick.labelsize is already the 7.19 these asked for.
    ax.set_xticks(_DX); ax.set_xticklabels(["other\nsubstitutions", "acquired\nlysines"])
    ax.set_ylabel("driver hotspots (%)")
    # 0 to 1.5 by 0.5, the PI's call, replacing 0 to 1.2 by 0.4. Four labels either way and
    # a rounder step, and the extra headroom is what lets the OR line sit clear of the
    # 1.024% bar label.
    ax.set_ylim(0, 1.5); ax.set_yticks([0, 0.5, 1.0, 1.5])
    # Span held at 2.0521 units so the 410 px frame still draws 199.8 px per unit and the
    # bars keep their drawn thickness. Only the centre moves: the pair now spans 0 to 0.851
    # with its centre at 0.4255, so the window is that centre plus and minus half the span.
    # The 29.8 px taken out from between the bars reappears as 14.9 px at each end, which
    # takes the end margins from 71.8 px to 86.7 px.
    _DC = 0.5 * (_DX[0] + _DX[1])
    ax.set_xlim(_DC - 1.02605, _DC + 1.02605)
    HEAD["d_pct_driver_K"] = kr; HEAD["d_pct_driver_other"] = orr
    HEAD["d_OR"] = OR; HEAD["d_p"] = P
    HEAD["d_n_K"] = len(K); HEAD["d_n_other"] = len(O)
    mant = P / (10 ** np.floor(np.log10(P))); expo = int(np.floor(np.log10(P)))
    # Stays above the frame, where it was approved. Hanging it inside instead put it 89%
    # over the 1.024% bar label, because a 1.024 bar on a 1.2 axis leaves no headroom.
    # Held at 93.75% of the axis height, 1.40625 of 1.5 where it was 1.125 of 1.2, and
    # centred on the PAIR rather than on x = 0.5, which stopped being the pair's centre
    # when the second bar came in. On a 1.5 axis the 1.024% bar label tops out near 1.09,
    # so there is now 0.3 of clear air under this line where there was 0.03.
    ax.text(_DC, 1.40625, "OR %.2f, p = %.1f \u00d7 10$^{%d}$" % (OR, mant, expo),
            ha="center", va="bottom", fontsize=st.sz(7), color="k")
else:
    notes.append("1d: depmap_variants_annotated.parquet not found")

# ------------------------------------------------------------------ e  (DepMap scope)
ax = fig.add_subplot(gs_bot[0, 0]); st.head(fig, ax, "e")
ax._panel_dx_in = -123.0 / 300.0
if MAF:
    d = pd.read_csv(MAF, usecols=["HugoSymbol", "ProteinChange", "HessDriver", "HessSignature"],
                    low_memory=False).dropna(subset=["ProteinChange"])
    d["alt"] = d.ProteinChange.astype(str).str.extract(r"([A-Z])$")
    h = d[d.HessDriver.astype(str).str.lower().isin(["true", "yes", "y", "1"])].copy()
    h["sig"] = h.HessSignature.astype(str).str.split(":").str[0]
    hK = h[h.alt == "K"]; hO = h[h.alt != "K"]
    order = ["UV", "MSI", "POLE", "APOBEC", "smoking", "CpG"]
    pk = [100 * (hK.sig == s).sum() / len(hK) for s in order]
    po = [100 * (hO.sig == s).sum() / len(hO) for s in order]
    HEAD["e_sig_K"] = dict(zip(order, [round(x, 2) for x in pk]))
    HEAD["e_sig_other"] = dict(zip(order, [round(x, 2) for x in po]))
    HEAD["e_n_driver_K"] = len(hK); HEAD["e_n_driver_other"] = len(hO)
    x = np.arange(len(order)); w = st.grp(2)
    ax.bar(x - w/2, pk, w, color=st.BLUE, label="acquired lysine")
    ax.bar(x + w/2, po, w, color=st.LIGHTGREY, label="other")
    # 30 degrees let POLE and APOBEC overlap by 31% of the smaller box (found by
    # figtools/overlap_check.py, not by eye). 42 degrees and a hair smaller clears it.
    # No explicit size: the house xtick.labelsize is 7.19, which is what c and d print.
    ax.set_xticks(x); ax.set_xticklabels(order, rotation=45, ha="right",
                                         rotation_mode="anchor")
    ax.set_ylabel("driver hotspots (%)")
    ax.set_ylim(0, max(max(pk), max(po)) * 1.30)
    ax.legend(frameon=False, loc="upper left", handlelength=1.1,
              bbox_to_anchor=(0.0, 1.02))
    for i in range(len(order)):
        r = pk[i] / po[i] if po[i] else 0
        # label every ratio that departs from parity, in both directions. Depletion in
        # smoking and CpG is part of the argument — those processes do not generate the
        # transitions that create lysine — so labelling only enrichments would be selective.
        if r and (r >= 1.25 or r <= 0.8):
            # Body size, matching the in-panel text of panels a to d. sz(5.6) gave 6.51.
            ax.text(i, max(pk[i], po[i]) * 1.04, "%.1f\u00d7" % r, ha="center",
                    fontsize=st.sz(7),
                    color=st.BLUE if r > 1 else st.GREY)
else:
    notes.append("1e: OmicsSomaticMutations.csv not found")

# ------------------------------------------------------------------ f
# Promoted from figS1 panel e: every TCGA lineage, not the top eight. The eight-lineage
# version that used to sit here is retired — the whole point of this panel is the UVM
# CONTROL, and a truncated list cannot carry it. Cutaneous melanoma is ultraviolet-driven
# and uveal melanoma is not; if the acquired-lysine excess were intrinsic to melanocytes
# both would be high. SKCM is 13x pan-TCGA and UVM is zero. Spans 9 columns because 32
# rotated lineage labels need ~4.2 in.
# 0:9 -> 0:8, with panel g left at 9:12 so the gap between them widens rather than g
# growing to absorb it. This is what sets Figure 1's PRINT SCALE. The composite packs
# [a b] [c d] [e f], and the row holding e and f was 1261 px against the 1173 it can
# have: wider than COMMON_W, so Figure 1 alone printed at 0.84 while every other figure
# printed at 0.90, and Figure 1 alone still had type below 6 pt. Narrowing f and NOT
# widening g brings that row to 1160 and the composite under COMMON_W.
ax = fig.add_subplot(gs_bot[0, 1]); st.head(fig, ax, "f")
ax._panel_head_flush = True   # letter on the block edge, with panel a above it
# Bringing this panel up to the common type scale widened its two-line y-axis label, which
# pushed its block out to 25.5 and left it 21 px adrift of the a/c edge at 46.4. Shift the
# whole axes right by that much so the block, not just the letter, joins the column.
# [cell is now the frame; offset retired with the 3-row grid] ax._axes_dx_in = 20.9 / 300.0
# Wider by 6%, which is what buys air between the 32 rotated lineage codes. They step down
# a 45 degree diagonal, so what separates two of them is the centre spacing times sin 45,
# and at the gridspec's own 1242 px the worst neighbouring pair, ACC and PRAD, overlapped
# by 10.6%, under the 12% gate but touching on the page. Each 2% of width buys about 1.7
# points of that overlap back. Clearing it outright wants 13%, which would carry this
# panel's crop past the point where the composed figure still fits 1225 px, and a figure
# wider than COMMON_W is placed at a smaller scale than the rest of the set, so its type
# prints smaller. 6% is the most this panel can take without paying that. The bars widen
# with the axes, 27.2 px to 28.8 px, still short of the 16 pt cap.
# [cell is now the frame; offset retired with the 3-row grid] ax._axes_dw_in = 74.5 / 300.0
if BURD:
    T = pd.read_csv(BURD).sort_values("median", ascending=False)
    HEAD["f_n_lineages"] = len(T)
    HEAD["f_SKCM"] = float(T.loc[T.cancer == "SKCM", "median"].iloc[0]) if (T.cancer == "SKCM").any() else None
    HEAD["f_UVM"] = float(T.loc[T.cancer == "UVM", "median"].iloc[0]) if (T.cancer == "UVM").any() else None
    HEAD["f_SKCM_fold"] = float(T.loc[T.cancer == "SKCM", "fold_vs_pan"].iloc[0]) if "fold_vs_pan" in T.columns and (T.cancer == "SKCM").any() else None
    x = np.arange(len(T))
    cols = [st.ORANGE if c == "SKCM" else (st.GREEN if c == "UVM" else st.BLUE)
            for c in T.cancer]
    ax.bar(x, T["median"], color=cols, width=st.BARW)
        # No explicit size: the house xtick.labelsize is 7.19, as in c, d and e.
    # 45 degrees, the same angle c and e use. Upright was the safe choice when overlap was
    # judged on axis-aligned boxes; measured properly these four-letter codes step down the
    # diagonal with no pair over 12%, and the label band is 19.7 pt rather than 21.2.
    ax.set_xticks(x); ax.set_xticklabels(T.cancer, rotation=45, ha="right",
                                         rotation_mode="anchor")
    ax.set_ylabel("median acquired lysines\nper tumor")
    ax.set_ylim(0, 40); ax.set_yticks([0, 10, 20, 30, 40])
    # No margin at either end: the first bar starts on the y axis and the last tick is the
    # right edge. Half a slot each side is what a categorical axis needs, and BARW is 0.70,
    # so the bars clear the spines by 0.15 of a slot.
    ax.set_xlim(-0.5, len(T) - 0.5)
    # Each callout is anchored on the side it grows AWAY from its own bar. The UVM one sits
    # at the right-hand end of 32 lineages, so left-anchored it grew rightward with the type
    # and ran 15% into panel g's axis label the moment this panel came up to body size.
    # Right-anchored, it can only grow back into panel f's own empty right-hand side.
    # Both callout numbers are READ from the same frame the bars are drawn from. They were
    # literal strings, "13x pan-TCGA" and "median 0", which happened to agree with the data
    # today and would have gone stale silently the first time the burden table was rebuilt.
    # Both callout numbers are READ from the frame the bars are drawn from. They were
    # literal strings, "13x pan-TCGA" and "median 0", which happened to agree with the data
    # today and would have gone stale silently the first time the burden table was rebuilt.
    _fold = HEAD["f_SKCM_fold"]; _uvm = HEAD["f_UVM"]
    # SKCM is the tallest bar and sits against the y axis, so its label needs no leader:
    # there is only one bar it can be naming.
    if "SKCM" in set(T.cancer):
        _i = list(T.cancer).index("SKCM")
        ax.text(_i + 1.1, 31.0,
                ("melanoma, %g\u00d7 pan-TCGA" % round(_fold)) if _fold is not None
                else "melanoma",
                ha="left", va="center", fontsize=st.sz(7), color=st.ORANGE)
    # UVM has NO bar to point at, its median being 0, so the leader carries the whole job of
    # saying which lineage is meant. A vertical arrow dropped onto its tick does that
    # unambiguously where a diagonal line across 30 neighbours did not.
    if "UVM" in set(T.cancer):
        _i = list(T.cancer).index("UVM")
        ax.text(_i + 0.35, 15.2,
                ("uveal melanoma, not UV-driven,\nmedian %g" % _uvm) if _uvm is not None
                else "uveal melanoma, not UV-driven",
                ha="right", va="bottom", fontsize=st.sz(7), color=st.GREEN, linespacing=1.45)
        ax.annotate("", xy=(_i, 0.4), xytext=(_i, 14.4),
                    arrowprops=dict(arrowstyle="-|>", color=st.GREEN, lw=0.8,
                                    shrinkA=0, shrinkB=0, mutation_scale=6))
    ax.text(0.99, 0.97, "%d lineages, %s tumors" % (len(T), f"{int(T.n.sum()):,}"),
            transform=ax.transAxes, ha="right", va="top", fontsize=st.sz(7), color=st.GREY)
else:
    notes.append("1f: tcga_burden_by_lineage.csv not found")

# ------------------------------------------------------------------ g
# Promoted from figS1 panel d. This is the control that excludes the cheapest alternative to
# the whole premise — that acquired lysines are simply what mutation does everywhere, and the
# same excess would appear in healthy people. It does not: somatic mutation creates them at
# twice the germline rate. Nothing else in this figure addresses that.
# Three columns. Two was tried: it gives the small margins, but 216 px is too narrow for
# the two-line tick labels, which then overlap by 21.9 px and fit_ticks stands them up at
# 42 degrees. Angled two-line labels on a two-bar panel read worse than the margins did.
ax = fig.add_subplot(gs_bot[0, 2]); st.head(fig, ax, "g")
# Frame widened 30 px and shifted 15 px left, so the extra end margins come out of air on
# both sides rather than off the bars. NOTE: unlike figure 2, this figure's offsets arrive
# WHOLE. Figure 2 overflows its canvas and _fit_once scales every axes inward by about
# 0.795 to bring it back, so an offset asked for there lands at four fifths. Figure 1 does
# not overflow, so -15 px of dx is -15 px on the page. Asking for -18.9 and +37.8 here, the
# figure 2 numbers, overshot by exactly the damping.
ax._axes_dx_in = -15.0000 / 300.0
ax._axes_dw_in = 30.0000 / 300.0
# Letter on the block edge, the same rule a, c and f use: flush to this panel's own
# leftmost ink rather than sitting inside it at the axes frame.
ax._panel_head_flush = True
# Narrower than its gridspec cell, so the two bars can come together without the end
# margins opening back up. The arithmetic is in the xlim note further down.
# [cell is now the frame; offset retired with the 3-row grid] ax._axes_dw_in = -94.4 / 300.0
# PROVENANCE, corrected 2026-09-04. An earlier note here said the germline pair had no source
# in the workspace. That was a search failure on my part, not a data gap: germline_lysine.py and
# germline_lysine_nc.py sit in this directory and write their results to .out files, which the
# search that looked for them did not cover.
#
# Both germline values are read from germline_lysine.out at build time rather than pasted, so
# re-running that analysis propagates here instead of silently diverging. Method, from the script:
# dbNSFP 5.4a restricted to rows carrying a MANE annotation (column 29) with an unambiguous
# single-letter amino-acid change, gnomAD v4.1 joint allele frequency (column 428) > 0. The
# per-individual figure is the expected diploid count, sum of 2 x AF, NOT a count of distinct
# variants. The non-cancer control (germline_nc.out, gnomAD v2.1.1 non-cancer exomes) is read too
# and printed beside it, because two independent germline references agreeing is the control this
# panel rests on.
def _read_germline(fn):
    """(acquired lysines per individual, % of missense) from a germline_*.out report."""
    path = find(fn)
    if not path:
        return None, None, "missing"
    txt = open(path).read()
    m1 = re.search(r"ACQUIRED LYSINES\s*:\s*([\d.]+)", txt)
    m2 = re.search(r"lysine fraction\s*:\s*([\d.]+)%", txt)
    if not (m1 and m2):
        return None, None, "unparsed"
    return float(m1.group(1)), float(m2.group(1)), "ok"

germ_n, germ_pct, _gstat = _read_germline("germline_lysine.out")
nc_n, nc_pct, _ncstat = _read_germline("germline_nc.out")
if germ_pct is None:
    notes.append("1g: germline_lysine.out %s, germline bar cannot be drawn" % _gstat)
    germ_pct, germ_n = float("nan"), float("nan")

if "a_union_K" in HEAD:
    som_pct = 100 * HEAD["a_union_K"] / HEAD["a_all_missense"]
else:
    som_pct = float("nan")
    notes.append("1g: union counts unavailable, somatic bar cannot be computed")
HEAD["g_som_pct"] = som_pct
HEAD["g_germ_pct"] = germ_pct
HEAD["g_germ_per_individual"] = germ_n
HEAD["g_germ_nc_pct"] = nc_pct
HEAD["g_germ_nc_per_individual"] = nc_n
HEAD["g_fold"] = som_pct / germ_pct if germ_pct else None
# The somatic bar sits at 0.674, not 1.0, so it and its tick label come in towards the
# germline pair, and both the axes width and xlim are pinned rather than left to the
# automatic margins. Those gave 17.7 px of air inside the y axis and 48.8 px past the
# somatic bar. fit_bars draws each bar at the 16 pt cap, 66.7 px, not at its nominal 0.70
# units, so a margin is measured off the DRAWN edge. 30 px of air plus 33.35 px of half
# bar puts the frame 63.35 px outside each centre, which in the gridspec's own 388.7 px
# axes forced the centres 262 px apart and read as two bars shoved against the walls.
# Closer bars and 30 px margins cannot both be had at a fixed width, so the axes is
# narrowed 94.4 px by _axes_dw_in above instead of surrendering the margins. The centres
# sit 171.6 px apart, which is 8 px above the 163.6 px floor the two tick labels impose.
# That is as close as these bars go while "(gnomAD v4.1)" and "(union)" stay on the
# artwork, and the 8 px is clear air between the labels, not slack.
# 6 Nov: both end margins widened on the PI's instruction, 30 px of air each side where
# there were 30 px... see the xlim note below. The bar centres do not move: the frame is
# widened by the same 30 px the margins take, through _axes_dw_in and _axes_dx_in above,
# so px per unit holds at 254.6 and the centres stay 171.6 px apart, 8 px above the floor
# the two tick labels impose.
_GX = [0, 0.674]
ax.bar(_GX, [germ_pct, som_pct], st.BARW, color=[st.LIGHTGREY, st.BLUE])
for _x, v in zip(_GX, [germ_pct, som_pct]):
    # Body size, matching the in-panel text of panels a to f. sz(6) gave 6.51.
    ax.text(_x, v + 0.12, "%.2f%%" % v, ha="center", fontsize=st.sz(7))
# 0.3078 either side of the pair instead of 0.2488. At 254.6 px per unit that is 78.4 px
# from each centre to the frame, and 33.35 px of that is half a drawn bar, so the air goes
# from 30.0 px to 45.0 px at each end. The span grows from 1.1716 to 1.2896 units and the
# frame grows 30 px with it, which is the only way to open the margins without pulling the
# bars together: at a fixed frame a wider window just shrinks the pixels.
ax.set_xlim(-0.3078, 0.9818)
ax.set_xticks(_GX)
# No explicit size: the house xtick.labelsize is 7.19, as in c to f.
# The source stays on the artwork. "(gnomAD v4.1)" is 223.6 px wide against a 66.7 px
# bar and a tick label is centred on its bar, so these two lines are what set how close
# the bars can sit. The floor is 163.6 px between centres, at which the labels touch.
ax.set_xticklabels(["germline\n(gnomAD v4.1)", "somatic\n(union)"])
ax.set_ylabel("acquired lysines as % of\nall missense")
# 0 to 10 by 2.5, the PI's call, replacing 0 to 8 by 2. Five labels either way.
ax.set_ylim(0, 10); ax.set_yticks([0, 2.5, 5, 7.5, 10])
ax.set_yticklabels([("%g" % v) for v in (0, 2.5, 5, 7.5, 10)])
# 7.25 clears the 6.72% label, whose top sits at about 7.17, and keeps the fold inside the
# 8 the axis now stops at. It was som_pct * 1.20, which lands at 8.07 and would print above
# the frame.
# Top left, over the germline bar's headroom, which is empty above 3.29. Centred at the
# top the long form overlapped the 6.72% label by 5.8%, under the 12% gate and touching on
# the page. Placed in AXES FRACTION, not data units: at x = -0.30 in data it drifted out
# past the frame and onto the "8" tick label at 22% the moment the x limits were retuned
# for the margins. 0.035 of the axes keeps it just inside the frame whatever xlim becomes.
# Bringing the two bars together then carried "6.72%" left with its own bar until it stood
# 1.3 px from this label, so the label is lifted to 0.9635 as well, as high as it goes with
# its top still inside the frame. The top band of this panel is tight in both directions.
ax.text(0.035, 0.9635, "%.2f-fold" % (som_pct / germ_pct), transform=ax.transAxes,
        ha="left", va="center", fontsize=st.sz(7), color=st.BLUE)
# "471 per individual" moved to the caption on 2026-11-05. It is still in the ledger as
# g_germ_per_individual and the caption already carried the fuller form, 471.1 acquired
# lysines among 14,342 missense variants, with the statement that this is an expected
# diploid count and not a count of distinct variants.

st.save(fig, "fig1")
if __import__("os").environ.get("FIG1_FINAL"):
    # The 12.00 in file, the one that goes to CCB-submission/figures-final/.
    # The rule in that folder's README: render the SAME trimmed figure at a finer raster and
    # declare the result 300 dpi. It scales everything at once, type included, and it is not
    # a resample. st.save() has already sized the figure to the trim, so the only thing to
    # choose is the raster: 300 x 12 / native_width_in. Done here rather than by hand so the
    # finalized file can be rebuilt from the script after any panel change.
    import PIL.Image as _PI
    _nat_in = fig.get_size_inches()[0]
    _dpi = 300.0 * 12.0 / _nat_in
    fig.savefig("fig1_12in.png", dpi=_dpi)
    _im = _PI.open("fig1_12in.png")
    _im.save("fig1_12in.png", dpi=(300, 300))      # declare it 300 dpi, so 3600 px is 12.00 in
    print("12 in file: native %.4f in -> raster %.2f dpi -> %d x %d px declared 300 dpi (%.3f x %.3f in)"
          % (_nat_in, _dpi, _im.size[0], _im.size[1], _im.size[0] / 300.0, _im.size[1] / 300.0))
if __import__("os").environ.get("FIG1_PROBE"):
    # geometry probe, env-gated: the frame each panel actually got after the fit damped the
    # declared offsets. Every offset in this file is converged against these, never computed.
    for _a in fig.axes:
        _bb = _a.get_window_extent()
        print("probe %s  x0 %8.1f  x1 %8.1f  w %7.1f  |  xlim %s"
              % (getattr(_a, "_panel_letter", "?"), _bb.x0, _bb.x1, _bb.width, _a.get_xlim()))

fd.head("fig1", HEAD)
for _e in (".png", ".pdf"):          # review copies alongside fig4_v2 / figS_react_v2
    fig.savefig("../results/fig1_v3" + _e, dpi=400, bbox_inches="tight")
def _fmt(v):
    if isinstance(v, float):
        return "%.4g" % v
    if isinstance(v, int):
        return f"{v:,}"
    return str(v)

print("fig1 headline numbers")
for kk in sorted(HEAD):
    print("  %-26s %s" % (kk, _fmt(HEAD[kk])))
print("  germline source   germline_lysine.out (dbNSFP 5.4a, MANE, gnomAD v4.1 joint, ~807k)")
print("  germline control  germline_nc.out (gnomAD v2.1.1 non-cancer exomes, ~118k): %s%% "
      "against %s%%" % (HEAD.get("g_germ_nc_pct"), HEAD.get("g_germ_pct")))
if notes:
    print("  NOTES: " + " | ".join(notes))
print("wrote fig1.pdf and fig1.png")
for n in notes: print("  MISSING:", n)
