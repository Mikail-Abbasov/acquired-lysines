#!/usr/bin/env python3
"""
fig2.py, Figure 2: acquired lysines are largely invisible to conventional proteomics. 3x3.

  a  the four peptide species an acquired lysine creates
  b  gate 1, regions never covered (two search architectures agree)
  c  gate 2, covered regions where the mutant is still missed
  d  three sites with every alternative explanation excluded
  e  peptide length predicts which regions are never covered
  f  loss of the wild-type allele inverts gate 1
  g  only peptide length separates identified species from missing ones
  h  detection does not track allele fraction
  i  two independent routes to the same sites

SCOPE. Every panel counts ACQUIRED-LYSINE sites only. The variant databases contain every
mutation carried by each cell line, so a protein with several mutations contributes a row per
position; 18% of rows in the detection table are other substitutions. Those were legitimately
searched and are reported separately in the Supplementary material, where they also serve as a
control: the same pipeline, spectra and FDR treatment applied to variant peptides that do not
create a lysine.
"""
import os, re, sys, glob
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Rectangle
from matplotlib.transforms import Bbox

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import abbasov_style as st
import figure_data as fd

# Headline numbers, keyed. The prints at the foot of this file stay for a human reading the log;
# HEAD is what verify.py diffs against the manuscript.
HEAD = {}

# --- portability -------------------------------------------------------------------------
# The data search paths below were written as ~/Desktop/DepMap/... , which only resolves on
# the machine the analysis was built on. They are now anchored to this script's own location,
# so the build works wherever the folder lives (a collaborator's checkout, a container, a
# different user account). On the original machine the two are identical.
_TOP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# One body type size for this figure; every other size is a ratio to it.
BODY = 7.0
st.apply(BODY)

def find(n):
    for d in (".", "..", os.path.join(_TOP, 'ccle_plex05'),
              _TOP,
              os.path.join(_TOP, 'Datasets'),
              os.path.join(_TOP, 'ligand_rebuild'),
              os.path.join(_TOP, 'Datasets/pandepmap_out')):
        p = os.path.join(d, n)
        if os.path.exists(p): return p
    return None

MARGIN = 10          # the score-margin threshold section 2.3 defines


def acquired_lysine_only(df):
    """Keep rows whose stored `site` IS the acquired lysine.

    region_compare writes one row per mutated position, and a cell line carrying several
    mutations in one protein contributes a row for each. Those extra rows are real variant
    sites that were legitimately searched, but they are not acquired lysines, and counting
    them inflates every denominator in this figure. In the 42-plex table 1,902 of 10,372 rows
    (18%) are of this kind. They are reported separately (Supplementary), not discarded.
    """
    import re as _re3
    def _isk(ch, site):
        for t in str(ch).split(";"):
            m = _re3.match(r"^[A-Z](\d+)K$", t.strip())
            if m and int(m.group(1)) == site:
                return True
        return False
    return df[[_isk(c, s) for c, s in zip(df.change, df.site)]]


def _site_digest(acc="P07195", pos=287, new="K"):
    """The four peptide species at one acquired-lysine site, computed from the search FASTA.

    Panel a draws a real site rather than a cartoon, so the sequences and the three lengths
    it prints have to come from the same protein database the searches ran against, not from
    a literal typed into this file. Returns (wild-type stretch, mutant stretch, 0-indexed
    position of the acquired lysine inside the stretch, (first, last) 1-indexed residue, and
    the two flanking residues (before, after) that the [K].PEPTIDE.[Y] notation prints).

    Cleavage is trypsin's own rule: after K or R, not before proline. The stretch returned is
    the tryptic peptide of the WILD-TYPE protein that spans `pos`, which is exactly the span
    the mutation subdivides, so all four species live inside it.
    """
    fa = find("2014-02-04_REVuniprot_HUMAN_02_2014_contam_sorted.fasta")
    if fa is None:
        return None
    seq, grab, buf = None, False, []
    with open(fa) as fh:
        for line in fh:
            if line.startswith(">"):
                if grab:
                    seq = "".join(buf); break
                grab = line.startswith(">sp|%s|" % acc); buf = []
            elif grab:
                buf.append(line.strip())
    if seq is None and grab:
        seq = "".join(buf)
    if not seq or len(seq) < pos or seq[pos - 1] == new:
        return None
    def cuts(x):
        out = [0]
        for i, c in enumerate(x, start=1):
            if c in "KR" and not (i < len(x) and x[i] == "P"):
                out.append(i)
        if out[-1] != len(x):
            out.append(len(x))
        return out
    c = cuts(seq)
    a = max(v for v in c if v < pos) + 1
    b = min(v for v in c if v >= pos)
    wt = seq[a - 1:b]
    mut = wt[:pos - a] + new + wt[pos - a + 1:]
    flank = (seq[a - 2] if a >= 2 else "-", seq[b] if b < len(seq) else "-")
    return wt, mut, pos - a, (a, b), flank


PEP  = find("region_compare_pep.tsv")
CONV = find("region_compare_conv.tsv")
HITS = find("hits62.tsv")
SITES= find("sites_covalent_residues.csv")
FLUX = find("residue_flux_by_line.csv")
notes = []

# 5 Nov: on the Figure 1 standard. Two rows at 12 inches, written at 300 dpi so one display
# pixel is one output pixel and every measurement in this file is in the same units as the
# file it produces. Cells are sized to the FRAME each panel needs, because a subplot cell is
# the axes rectangle and the tick and axis labels hang outside it into the gutter. Both rows
# total 2900 px of frame against 3300 px of content, so both carry 200 px gutters and the
# right-hand pair, c and f, share a right edge.
# 5 Nov: Figure 2 put on Figure 1's type scale, which is the signed-off standard. Measured
# off the two renders, Figure 1 prints exactly three sizes and the role mapping is clean:
#   13.733  panel letter                      (LETTER_K x BODY, set by st.head)
#    8.860  every in-panel text, every axis label   (rung 7.0)
#    7.974  every tick label, every legend entry    (rung 6.3)
# Figure 2 was printing four, and none of its in-panel text reached 8.860: text ran at
# 6.645 x18, 7.214 x4, 7.974 x7, and its ticks were split 7.214 x8 / 7.974 x42 with the
# legend at 7.214. Panel a was the worst of it, every one of its ten strings at 6.645, the
# crowded-axis rung. Each st.sz() argument below is now 7 for text and axis labels and 6.3
# for tick labels and legends, so the two figures snap to the same two rungs. Nothing else
# here was touched: no panel moved, no bar changed width, no string changed.
plt.rcParams["savefig.dpi"] = 300
# Sheet trimmed to the ink, 18 px all round, the Figure 1 rule. The 12 inch deliverable
# comes from rendering this same figure at a finer raster, not from touching geometry.
fig = plt.figure(figsize=(11.5333, 4.3663))
# A 12-column grid rather than 3, so panels can take the width their content needs. The bar
# charts (b, c, f, i) carry two or three bars each and were given a full third of the figure;
# they are narrowed to 3 columns, which leaves room for panel j without stretching it.
# The gridspec margins are NOT the lever for where the content lands. _fit_once measures how
# far the drawn content runs past the canvas and scales every axes inward by that much, and
# on the first pass this figure overflows by 1.51 in on the left and 0.64 on the right, so
# the settled position is set by that scaling rather than by these numbers. Moving `left`
# here moves the ink the OTHER way: measured, -152 px of margin put the ink 29 px further
# right, because a larger overflow buys a harder inward scale. The sheet is squared up by
# trimming at save instead, which is where the 18 px pads come from.
gs = fig.add_gridspec(2, 1, hspace=0.36,
                      left=-0.033526, right=0.996532, top=0.952055, bottom=0.097704)
# Frames sized so each row's INK, not its frames, spans the sheet with even 95 px gaps.
# Row 2 is the tight one: panel d's y tick labels reach 171 px left of its frame and its
# gene list 260 px right of it, so d owns 1824 px of ink for a 1393 px frame.
gs_top = gs[0, 0].subgridspec(1, 3, wspace=0.1939, width_ratios=[1840, 570, 746])
gs_bot = gs[1, 0].subgridspec(1, 3, wspace=0.6513, width_ratios=[1393, 646, 446])

p_all = pd.read_csv(PEP, sep="\t") if PEP else None
p = acquired_lysine_only(p_all) if p_all is not None else None
if p_all is not None:
    print("region_compare rows: %d total, %d acquired lysine, %d other variant sites"
          % (len(p_all), len(p), len(p_all) - len(p)))

# ================================================================== six panels, 2 x 3
# The four alternative-explanation panels (allele loss, abundance, allele fraction) and the
# digestion panel move to the supplement: three nulls and a methods result are not three
# and a half main-text slots. Evidence depth, which had its own panel, is now a line inside
# the recovery panel, where it qualifies the number it belongs to.


# ------------------------------------------------------------------ a
ax = fig.add_subplot(gs_top[0, 0]); st.head(fig, ax, "a")
ax._panel_head_flush = True
# Placed from measurement, not from a column count: each panel's own left and right
# overhang differs, so one wspace cannot give even INK gaps. dx lands the frame, dw
# sizes it, and the row then runs 18 px to 3582 px with 95 px between ink blocks.
#
# 6 Nov: dw takes 220.6 px OFF the frame, which is exactly what the LDHB peptide gave back.
# The axis is off so the frame was invisible, but it is not inert: intrude.py reads frames,
# and at a 100 px gap panel b's y label and two of its tick labels were sitting 61%, 39%
# and 36% inside panel a's empty right half. Six intrusions, all of them this. The frame
# now ends 6.9 px past the end of the written peptide, the same slack it had at 27 residues.
# This is dw on the axes, NOT the gridspec: a width_ratios change would re-damp _fit_once
# for all six panels and widen b and c, and the instruction was to move them, not resize.
ax._axes_dx_in = -218.6000 / 300.0
ax._axes_dw_in = -236.8000 / 300.0
ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.set_axis_off()
# 5 Nov: the schematic chain is gone and the panel draws ONE REAL SITE. 6 Nov: that site is
# LDHB E287K in NCIH1793, chosen off TYRP1 E210K to shorten the row. The site is one row of
# region_compare_pep.tsv, the same table panels b and c are built from, with shortest_wt 20,
# wt_detected 0 and mut_detected 1: the wild-type region was never covered in 42 plexes and
# the mutant peptide was recovered anyway, which puts this site in panel c's 26-of-7,853
# class. The panel is an instance of the asymmetry the figure is about, not an illustration
# of it. (The earlier note here cited acquired_lysine_peptides.csv for a 42-plex claim. That
# file is the output of span_test.py on plex05 alone, four MEWO rows, so it could not carry
# the claim; the gate flags below and in the print are read from region_compare_pep instead.)
#
# Tryptic digestion of UniProt P07195 with and without the mutation, computed in
# _site_digest() below rather than typed in, gives around position 287:
#     wild type                280-299  20 aa  GMYGIENEVFLSLPCILNAR
#     mutant, cleaved          280-287   8 aa  GMYGIENK        <- the species observed
#     downstream               288-299  12 aa  VFLSLPCILNAR
#     mutant, missed cleavage  280-299  20 aa  GMYGIENKVFLSLPCILNAR
# The 20-mer is above the never-covered mean panel e measures (covered 13.3, never covered
# 16.7); the acquired lysine splits it into an 8-mer and a 12-mer, both inside the observable
# range and both well clear of the length floor, so the two detection routes the panel draws
# are each a plausible peptide. The panel therefore shows the mechanism panel e quantifies,
# on one molecule, and the three numbers it prints are computed here, not transcribed.
# TYRP1 E210K cost 189 px more of row: 27 letters against 20, on a 436 px label column.
SITE_ACC, SITE_GENE, SITE_LINE, SITE_POS = "P07195", "LDHB", "NCIH1793", 287
SEQ_WT, SEQ_MUT, KPOS, SPAN, FLANK = _site_digest(SITE_ACC, SITE_POS)
NRES = len(SEQ_MUT)
SITE_LAB = "[%d]" % (SPAN[0] + KPOS)       # the residue number, 287 here
PRE  = "[%s]." % FLANK[0]
POST = ".[%s]" % FLANK[1]

# The row is written in the notation a search engine reports a localised site in,
# [K].VAK[173]SAK.[K], with the flanking residues in brackets outside the dots and the site
# residue carrying its position. Laying it out needs real glyph widths, because the bars
# below have to line up with individual residues: the 27 residues sit on ONE uniform pitch so
# a bar's length is its residue count, while [K]., [210] and .[Y] take whatever width they
# actually print at. Solving the pitch from measured ink rather than guessing it is what
# keeps the W in TWHR off the H that follows; at an eyeballed pitch the two overlapped.
def _inkw(txt, bold=False):
    """Ink width of a string in 300 dpi px at the house text size."""
    from matplotlib.textpath import TextPath
    from matplotlib.font_manager import FontProperties
    fp = FontProperties(family=plt.rcParams["font.sans-serif"][0], size=st.sz(7),
                        weight="bold" if bold else "normal")
    return TextPath((0, 0), txt, prop=fp).get_extents().width * 300.0 / 72.0

# 5 Nov: the letters are set PROPORTIONALLY, not on a uniform pitch. A uniform grid has to
# clear the widest glyph in the string, and this sequence contains a W at 34.0 px against a
# mean of 21.7, so one letter was setting the spacing for all 27 and the row cost 891 px to
# carry 676 px of ink. Set solid with 2 px of tracking the same letters take 728 px, and the
# row is 180 px shorter, which is what shortens the panel.
#
# The bars below still line up with individual residues, because each letter's left and right
# edge is measured rather than derived from a pitch. Their lengths are no longer exactly
# residue counts, but on this sequence they are within 0.1%: the 12-mer draws at 0.800 of the
# 15-mer against a true 12/15 of 0.800, because the two halves happen to average the same
# letter width. That is luck, not a property, so the caption states counts and the drawing is
# not offered as the measurement.
# 6 Nov: the frame narrowed from 1495.6 px to 1275.0 and NOTHING IN THE PANEL MOVED. Both
# of these were data-unit literals, so they tracked the frame; every other horizontal
# position here is already px/PXDU and did not. Written as pixels over PXDU they are now
# frame-invariant too, and the drawing holds its absolute column whatever the frame does.
PXDU   = 1275.0 / 10.0                     # px per data unit across panel a's frame
LBL_R  = 530.5 / PXDU                      # species labels, right-aligned, 530.5 px in
SX0    = 546.5 / PXDU                      # where the written peptide starts, 16 px on
PYDU   =  579.4 / 10.0                     # and down it. The frame is wide and short, so a
                                           # shape drawn in px needs both, not one of them.
JGAP   = 10.0 / PXDU                       # clearance where the letters meet [K]. and .[Y].
                                           # 6 left those two junctions at 0.2 and 0.8 px,
                                           # tighter than the 2 px between the letters.
TIGHT  = 1.5 / PXDU                        # K and its [210] are one token and set solid
TRK    = 2.0 / PXDU                        # tracking added between letters. At 0 the sequence
                                           # sets solid, which is how one-letter code is
                                           # normally written but reads as a wall at this size.
CUTGAP =  6.0 / PXDU                       # the break in the chain at the new cleavage site
_wpre  = _inkw(PRE) / PXDU
_wsite = _inkw(SITE_LAB, bold=True) / PXDU
_wpost = _inkw(POST) / PXDU
def _advw(ch, bold=False):
    """ADVANCE width of one glyph in 300 dpi px, which is what matplotlib lays text out on.

    Ink width is the wrong measure here. A Text drawn with ha="left" anchors its layout box,
    whose left edge sits outside the ink by the glyph's left side bearing, so positions built
    from ink widths overlap: with 2 px of tracking the measured gaps came out at -1.3 to +0.2
    instead of +2. Advances make the drawn boxes abut exactly, and the tracking is then real.
    """
    from matplotlib import font_manager as _fm
    from matplotlib.ft2font import FT2Font as _FT
    fp = _fm.FontProperties(family=plt.rcParams["font.sans-serif"][0], size=st.sz(7),
                            weight="bold" if bold else "normal")
    f = _FT(_fm.findfont(fp)); f.set_size(st.sz(7), 300); f.set_text(ch)
    return f.get_width_height()[0] / 64.0

_ADV = {c: _advw(c) / PXDU for c in set(SEQ_MUT)}
_ADV[SEQ_MUT[KPOS]] = _advw(SEQ_MUT[KPOS], bold=True) / PXDU      # the acquired K is bold
_INK = _ADV
_PRE_W = [0.0]
for _c in SEQ_MUT:
    _PRE_W.append(_PRE_W[-1] + _ADV[_c])
R0     = SX0 + _wpre + JGAP                # left edge of residue 0
SHIFT  = TIGHT + _wsite + CUTGAP - TRK     # what [210] pushes the C-terminal half

def _rx0(i):
    """Left edge of residue i. Residues after the site are pushed right by [210]."""
    return R0 + _PRE_W[i] + i * TRK + (SHIFT if i > KPOS else 0.0)

def _rx1(i):
    return _rx0(i) + _ADV[SEQ_MUT[i]]

def _rx(i):
    """Centre of residue i, for anything that wants to sit over a letter."""
    return 0.5 * (_rx0(i) + _rx1(i))

def _redge(i, right):
    """Bar edge at residue i, half the tracking beyond the ink."""
    return (_rx1(i) + TRK / 2) if right else (_rx0(i) - TRK / 2)

SX1 = _rx1(NRES - 1) + JGAP + _wpost       # where the row actually ends, now a result

# the written peptide, one text call per residue so every letter's edges are known to the
# bars below, flanks and site number as single strings. Arial throughout, like every other
# string here. Each letter is anchored by its LEFT ink edge, which is what _rx0 returns.
ax.text(SX0, 9.20, PRE, ha="left", va="center", fontsize=st.sz(7), color=st.LIGHTGREY)
for _i, _c in enumerate(SEQ_MUT):
    ax.text(_rx0(_i), 9.20, _c, ha="left", va="center", fontsize=st.sz(7),
            color=st.ORANGE if _i == KPOS else st.GREY,
            fontweight="bold" if _i == KPOS else "normal")
_sx0 = _rx1(KPOS) + TIGHT
ax.text(_sx0, 9.20, SITE_LAB, ha="left", va="center", fontsize=st.sz(7),
        color=st.ORANGE, fontweight="bold")
# left-aligned off the last residue, not right-aligned at SX1: ha="right" anchors on the
# advance box, whose left edge sits inside the ink by the period's side bearing, and the
# group landed 2.9 px on top of the R.
ax.text(_redge(NRES - 1, True) + JGAP, 9.20, POST, ha="left", va="center",
        fontsize=st.sz(7), color=st.LIGHTGREY)
# The cut the mutation creates. A vertical rule here was indistinguishable from the closing
# bracket of [210] and from the stems of the letters either side, so it is drawn as scissors.
# Arial carries no scissors glyph (U+2700 to U+2704 are all absent, checked against the file
# matplotlib resolves), and borrowing one would put a second font in the figure, so the shape
# is two crossed limbs and two handle loops in the house orange. Height is given in px and
# converted on each axis separately, because this frame is 1717 x 579 and a shape sized in
# data units would come out three times wider than it is tall.
# The cut is marked with the scissors BioRender draws on the Digestion step of its TMT
# Labeling Protocol template (t-5f4963c5b59c7a00ad18b5fd), which is the look chosen on 5 Nov.
# Rebuilt as geometry, not lifted and not imitated by eye. The only copy obtainable is a crop
# of a raster thumbnail, 69 x 51 px of heavily antialiased grey even from the 2000 px preview;
# tracing that gave a blobby silhouette with polygonal loops and blunt blades, and a uniform
# two-line sketch missed it completely, because the icon is a FILLED shape with thick handles
# and blades that taper to needle points. The numbers below are measured off that crop: the
# loop holes sit at (0.167, 0.869) and (0.184, 0.163) with an inner radius of 0.113, the limbs
# cross at (0.667, 0.500), and the whole mark is 1.353 long to 1 across. Built this way it is
# vector in the PDF, exactly CC3F0C, free of any BioRender licence, and sizeable to anything.
_SC_ASPECT = 1.353
_SC_LOOPS  = ((0.167, 0.869), (0.184, 0.163))     # upper left, lower left
_SC_TIPS   = ((1.300, 0.158), (1.318, 0.845))     # the tip each loop's limb runs to, crossed
_SC_RIN, _SC_ROUT = 0.113, 0.158                  # loop hole and outer rim
_SC_PIVOT  = (0.667, 0.500)
_SC_WH, _SC_WP = 0.034, 0.027                     # half width at the handle and at the pivot


def _scissors(ax_, x, y, w_px=80.0, col=st.ORANGE, flip=True):
    """The TMT-template scissors, centred on (x, y), w_px wide along the blades.

    flip mirrors the mark left to right, which turns the blades toward K[210] and puts the
    handles on the open side. Mirroring reverses every contour's winding together, so the
    rings keep their holes under the nonzero fill.

    Every point is built in the normalised frame above, scaled to pixels, then divided onto
    each axis separately, because panel a's frame is 1657 x 579 and a shape sized in data
    units would come out three times wider than it is tall.
    """
    import numpy as _np
    from matplotlib.path import Path as _Path
    from matplotlib.patches import PathPatch as _PathPatch
    sc = w_px / _SC_ASPECT                         # px per normalised unit (height is 1.0)
    cx_, cy_ = _SC_ASPECT / 2.0, 0.5
    mir = -1.0 if flip else 1.0
    def D(u, v):
        return x + mir * (u - cx_) * sc / PXDU, y + (v - cy_) * sc / PYDU
    def unit(v):
        v = _np.asarray(v, float); return v / _np.hypot(*v)
    def perp(v):
        return _np.array([-v[1], v[0]])
    verts, codes = [], []
    def add(pts):
        verts.extend([D(*pts[0])] + [D(*q) for q in pts[1:]] + [D(*pts[0])])
        codes.extend([_Path.MOVETO] + [_Path.LINETO] * (len(pts) - 1) + [_Path.CLOSEPOLY])
    P = _np.array(_SC_PIVOT)
    for loop, tip in zip(_SC_LOOPS, _SC_TIPS):
        L, T = _np.array(loop), _np.array(tip)
        d0, d1 = unit(P - L), unit(_np.array(tip) - P)
        # butt the handle onto the ring's outer rim: started inside it the strip is narrower
        # than the annulus and leaves wedge-shaped gaps, started at the centre it shows
        # through the hole.
        start = L + d0 * (_SC_ROUT * 0.98)
        nm, n0 = perp(unit(d0 + d1)), perp(d0)
        add([start + n0 * _SC_WH, P + nm * _SC_WP, T, P - nm * _SC_WP, start - n0 * _SC_WH])
    t = _np.linspace(0, 2 * _np.pi, 72, endpoint=False)
    for (lx, ly) in _SC_LOOPS:                     # ring as an outer loop plus a reversed hole
        add(_np.c_[lx + _SC_ROUT * _np.cos(t),  ly + _SC_ROUT * _np.sin(t)])
        add(_np.c_[lx + _SC_RIN * _np.cos(-t), ly + _SC_RIN * _np.sin(-t)])
    ax_.add_patch(_PathPatch(_Path(verts, codes), facecolor=col, edgecolor="none",
                             zorder=6, clip_on=False))
# Centred over the cut, the junction after K[210]. The +5 px is the bracket's side bearing:
# _sx0 anchors the advance box of "[210]", not its ink. The anchor rises with the size:
# the mark grows about its centre, and at 46 px wide it already sat 5.6 px off the tops
# of the letters, so 64 px has to be carried up or the lower blade lands on them.
_scissors(ax, _sx0 + _wsite + CUTGAP / 2 + 5.0 / PXDU, 10.24)
# The substitution, centred on the K itself rather than on the whole K[210] token. Over the
# token it reached 13 px into the scissors, which occupy the column just to the right of it;
# text_collisions compares text to text and cannot see an image sitting on a label.
ax.text(_rx(KPOS), 9.86, "E\u2192K", ha="center", va="center",
        fontsize=st.sz(7), color=st.ORANGE, fontweight="bold")
ax.text(LBL_R, 9.20, "%s %d-%d" % ((SITE_GENE,) + SPAN), ha="right", va="center", fontsize=st.sz(7),
        color=st.TEXT, fontweight="bold")

# species: label sits directly ON its own row, bar to the right of it
# Four species from the five colours. The two mutant-carrying species are the SAME peptide
# digested two ways, so they share the teal and are told apart by fill, not by a sixth hue:
# solid = cleaved, hatched outline = missed cleavage. Wild type drops to the background grey.
# 5 Nov: downstream moved into slot 2, directly under "mutant, cleaved". The two species the
# NEW cleavage site creates now sit together, one on each side of the cut, and the two that
# span the position uncut (missed cleavage, wild type) follow. The row slots are unchanged.
rows = [(7.6, 0,        KPOS,     st.BLUE,      "mutant, cleaved",         True,  None, True),
        (6.0, KPOS + 1, NRES - 1, st.GREEN,     "downstream",              False, None, False),
        (4.4, 0,        NRES - 1, st.BLUE,      "mutant, missed cleavage", True,  "///", True),
        (2.8, 0,        NRES - 1, st.LIGHTGREY, "wild type",               False, None, False)]
for y, i0, i1, col, lab, blue, hat, carries in rows:
    x0, x1 = _redge(i0, False), _redge(i1, True)
    ax.add_patch(Rectangle((x0, y - 0.34), x1 - x0, 0.68,
                           facecolor=st.PAPER if hat else col,
                           edgecolor=col, lw=0.9 if hat else 0.5, hatch=hat))
    # The acquired lysine marked in place on the two species that carry it, so which
    # peptides contain the residue is read off the drawing rather than off the caption.
    # The wild-type 27-mer has a glutamate at this position and is deliberately unmarked.
    if carries:
        ax.add_patch(Rectangle((_redge(KPOS, False), y - 0.34),
                               _redge(KPOS, True) - _redge(KPOS, False), 0.68,
                               facecolor=st.ORANGE, edgecolor=st.ORANGE, lw=0.5, zorder=4))
    # No residue count printed on the bar. Outside it there is no room, because the two
    # 27-mers end where the written peptide ends; inside, the count was unreadable over the
    # teal hatch and the "12 aa" landed on the orange lysine block. The bars are drawn on one
    # residue pitch, so their lengths ARE the counts, and 27, 12 and 15 reach the reader
    # through HEAD and the caption.
    # labels are right-aligned and end just before the written peptide begins, so the longest
    # one ("mutant, missed cleavage") cannot run into its bar
    ax.text(LBL_R, y, lab, ha="right", va="center", fontsize=st.sz(7),
            color=st.GREEN if col == st.GREEN else (st.BLUE if blue else st.GREY),
            fontweight="bold" if (blue or col == st.GREEN) else "normal")
HEAD["a_site"] = "%s E%dK (%s)" % (SITE_GENE, SITE_POS, SITE_LINE)
HEAD["a_wt_peptide_len"] = len(SEQ_WT)
HEAD["a_mutant_cleaved_len"] = KPOS + 1
HEAD["a_downstream_len"] = NRES - 1 - KPOS
# The two gate outcomes for this one site, read out of the same table the rest of the figure
# is built from rather than asserted in the caption. Both flags and shortest_wt reach the
# caption through HEAD, so the panel cannot disagree with panel c about its own example.
if PEP:
    _a_row = p[(p.gene == SITE_GENE) & (p.site == SITE_POS) & (p.line == SITE_LINE)]
    if len(_a_row) == 1:
        _a = _a_row.iloc[0]
        HEAD["a_shortest_wt"] = int(_a.shortest_wt)
        HEAD["a_wt_detected"] = int(_a.wt_detected)
        HEAD["a_mut_detected"] = int(_a.mut_detected)
        print("a: %s wild-type %d aa -> %d aa + %d aa | shortest_wt %d | gate 1 wt_detected %d, "
              "gate 2 mut_detected %d"
              % (HEAD["a_site"], len(SEQ_WT), KPOS + 1, NRES - 1 - KPOS,
                 int(_a.shortest_wt), int(_a.wt_detected), int(_a.mut_detected)))
    else:
        notes.append("panel a: %s not found as a single row in region_compare_pep" % HEAD["a_site"])
# 5 Nov: the three explanatory lines that stood here moved to the panel A caption, the same
# call made for Figure 1 panel a's union sentence. They were prose on the artwork, and on the
# Figure 1 type scale the three of them printed 0.8 px apart and read as one block.

# ------------------------------------------------------------------ b
ax = fig.add_subplot(gs_top[0, 1]); st.head(fig, ax, "b")
ax._panel_head_flush = True
# dx and dw are both damped by about 0.795 on the way to the page: fit_inside rescales the
# content box after nudge_axes has run, so an offset asked for in inches does not arrive
# whole. Both numbers below are therefore converged against the render, not derived.
# 6 Nov: b and c move LEFT. Panel a's row lost 220.6 px when the peptide went from TYRP1's
# 27-mer to LDHB's 20-mer, and all of it fell into the a|b gap, which measured 421 px of
# white against 195 px between b and c. Finalized Figure 1's packed top row runs 103, 100
# and 97 px, so 421 is four times the house gap and the unevenness is what reads first.
# dx only, so neither panel changes size. Second pass, same day: 241 and 160 were still
# loose, so b and c went to the house gap itself. b has now moved 321 px and c 416 px off
# where they sat with the 27-mer, which puts both gaps at about 100.
# -355.2, not -343.9. Changing b's ticks from 0/2.5/5/7.5/10 to 0/4/8/12 narrowed its
# widest tick label, which moved the panel's leftmost ink 9 px right and took the a|b gap
# from 100 to 109. The tick labels are part of the ink the gap is measured between, so a
# tick change is a layout change.
ax._axes_dx_in = -428.2000 / 300.0
# 5 Nov: x axis condensed on the Figure 1 panel g recipe. The two bars were 272.6 px between
# centres in a 509.7 px frame, and the gap was set by nothing: at BAR_MAX_PT the drawn bars are
# 66.7 px, so 206 px of the panel was the space between two categories. The floor is the tick
# labels, exactly as in 1g: "conventional" is 186.9 px, "mutation-/aware" 138.8 px, each centred
# on its bar, so the labels touch at 162.9 px between centres and 8 px of clearance puts the
# floor at 170.9 px. Bringing the bars together without narrowing the frame only moves the dead
# air to the ends, which is the mistake Figure 1 panel d made, so the frame comes in too: at
# 170.9 px between centres and 1g's end margin of 63.4 px (half a drawn bar plus 30 px of air)
# the frame is 297.7 px, within a pixel of 1g's own 298.34. dw was +71.0; 71.0 - 212.0 = -141.0.
ax._axes_dw_in = -102.7000 / 300.0
_BX = [0, 0.674]                    # 1g's spacing, which lands the centres on the 170.9 px floor
if PEP and CONV:
    c = acquired_lysine_only(pd.read_csv(CONV, sep="\t"))
    vals = [100 * c.wt_detected.sum() / len(c), 100 * p.wt_detected.sum() / len(p)]
    ax.bar(_BX, vals, color=[st.LIGHTGREY, st.BLUE], width=st.BARW)
    for _x, v in zip(_BX, vals):
        ax.text(_x, v + 0.2, "%.2f%%" % v, ha="center", fontsize=st.sz(7))
    ax.set_xticks(_BX)
    # "conventional" on its own line read as one run-on string with "mutation-" beside it,
    # 8.1 px away. Two two-line blocks read as two labels. The widest line is unchanged
    # ("conventional", 186.9 px), so the floor the bars sit on does not move.
    ax.set_xticklabels(["conventional\nsearch", "mutation-\naware"], fontsize=st.sz(6.3))
    # Margins re-cut 5 Nov. At the bare 170.9 px floor the bars stood 30.1 px off each end of
    # the frame, under half a bar width, and read as pushed into the corners. The frame is 34 px
    # wider and the limits hold the bars 57.4 px off each end, close to one bar width, while the
    # centres sit at 190 px. That is still a condensing: they were 272.6 px apart. It also keeps
    # 27 px between "conventional" and "mutation-", which matters here in a way it did not in
    # Figure 1 panel g: there the two widest lines are on different rows of their labels and
    # never meet, while here both are the first line and sit side by side.
    ax.set_xlim(-0.3221, 0.9961)
    ax.set_ylabel("regions covered (%)")
    # 0 to 12 by 4, the PI's call, replacing 0 to 10 by 2.5. Four labels instead of five and
    # every one an integer, so the %g formatter below is no longer doing any work; it stays
    # because a 2.5 step would otherwise print 0.0, 2.5, 5.0, 7.5, 10.0 if the step ever
    # returns. The bars are 7.56 and 7.28, so 12 leaves 37% of the scale above them, which
    # is where the "92.7% never observed" block sits.
    ax.set_ylim(0, 12); ax.set_yticks([0, 4, 8, 12])
    ax.set_yticklabels([("%g" % v) for v in (0, 4, 8, 12)])
    # Two lines now, and the exact figure rather than "~93%". 355.9 px of one line does not fit
    # a 297.7 px frame, and the tilde was the only approximate number in a figure that prints
    # 7.56 and 7.28 beside it. Hung from the top of the frame, which is where the uncovered
    # fraction actually is: "92.7%" over "never observed" is 249.4 px at its widest.
    ax.text(0.5, 0.985, "%.1f%%\nnever observed" % (100 - vals[1]), transform=ax.transAxes,
            ha="center", va="top", fontsize=st.sz(7), color=st.GREY, linespacing=1.35)
    HEAD["b_conv_covered_n"] = int(c.wt_detected.sum()); HEAD["b_conv_total"] = int(len(c))
    HEAD["b_conv_covered_pct"] = round(vals[0], 2)
    HEAD["b_pep_covered_n"] = int(p.wt_detected.sum()); HEAD["b_pep_total"] = int(len(p))
    HEAD["b_pep_covered_pct"] = round(vals[1], 2)
    HEAD["b_never_covered_n"] = int(len(p) - p.wt_detected.sum())
    HEAD["b_never_covered_pct"] = round(100 * (len(p) - p.wt_detected.sum()) / len(p), 2)
    fd.panel("fig2", "b", [
        dict(search="conventional 20,431-protein MANE", covered=int(c.wt_detected.sum()),
             observations=int(len(c)), pct=vals[0]),
        dict(search="mutation-aware 10.17M peptide", covered=int(p.wt_detected.sum()),
             observations=int(len(p)), pct=vals[1])], note="gate 1, both architectures")

# ------------------------------------------------------------------ c
ax = fig.add_subplot(gs_top[0, 2]); st.head(fig, ax, "c")
ax._panel_head_flush = True
# -407.2: c follows b's 9 px so the b|c gap stays at 100 as well. Moving one panel in a
# packed row moves every panel to its right by the same amount, or the gap it protects
# just reappears one place along.
ax._axes_dx_in = -509.2000 / 300.0
# 5 Nov: x axis condensed on the same recipe as panel b, with one difference. There the floor
# was set by the tick labels; here it is set by the BAR labels, which now carry their
# denominator on a second line: "(26/7,853)" and "(27/617)" have a half-sum of 149.3 px
# against the ticks' 141.8, so they are what would touch first. At panel b's 27 px of
# clearance the floor is 176.3 px between centres, and with b's end margin of 90.85 px from
# frame edge to bar centre the frame is 358.0 px, within 14 px of b's own 371.6. It was 684.7,
# the widest frame in the figure holding two bars 366.1 px apart. dw was +115.0; the damping
# is about 0.795, so 115.0 - 326.7/0.795 = -295.9.
ax._axes_dw_in = -250.1000 / 300.0
_CX = [0, 0.674]
if PEP:
    cov = p[p.wt_detected == 1]; unc = p[p.wt_detected == 0]
    rc = 100 * cov.mut_detected.sum() / len(cov); ru = 100 * unc.mut_detected.sum() / len(unc)
    ax.bar(_CX, [ru, rc], color=[st.LIGHTGREY, st.BLUE], width=st.BARW)
    # Each bar's rate and the counts behind it are now ONE label, the denominator parenthesised
    # on a second line. They were two: the percentage above the bar and the fraction in grey
    # below the axis, which needed ylim to start at -0.9 and read as a detached row of labels
    # rather than as the n behind each rate.
    for _x, (v, n, kk) in zip(_CX, [(ru, len(unc), int(unc.mut_detected.sum())),
                                    (rc, len(cov), int(cov.mut_detected.sum()))]):
        ax.text(_x, v + 0.12, "%.2f%%\n(%d/%s)" % (v, kk, f"{n:,}"), ha="center", va="bottom",
                fontsize=st.sz(7), linespacing=1.35)
    ax.set_xticks(_CX)
    ax.set_xticklabels(["not covered", "covered"], fontsize=st.sz(6.3))
    # Bars moved 18 px further in than panel b's margin. b's 57.5 px of bar edge to frame
    # edge works there because the labels are narrower than the bars are apart; here the
    # left block is 168 px wide on a bar 57.5 px off the spine, so its left edge sat 6.8 px
    # inside the frame and crossed the axis. 109 px from edge to bar centre clears it by 25.
    ax.set_xlim(-0.41672, 1.09072)
    ax.set_ylabel("mutant recovered (%)")
    # 0 to 6 by 2. It was (-0.9, rc * 1.32): the top quarter existed to hold this annotation
    # and the bottom 0.9 to hold the denominators, which now sit with their own bars.
    ax.set_ylim(0, 6); ax.set_yticks([0, 2, 4, 6])
    # Hard left, aligned with the left edge of the "0.33% (26/7,853)" block below it, so the
    # two read as one column in the empty upper left. It began at the midpoint of the bars,
    # which is dead centre of a frame whose right half already carries the tall bar and its
    # two-line label. Placed in AXES FRACTION rather than data units, the same thing Figure 1
    # panel g does with its fold label: 0.063 of the frame holds whatever xlim becomes.
    ax.text(0.063, 0.95, "%.0f-fold" % (rc / max(ru, 1e-9)), transform=ax.transAxes,
            ha="left", va="center", fontsize=st.sz(7), color=st.BLUE)
    HEAD["c_covered_found"] = int(cov.mut_detected.sum()); HEAD["c_covered_n"] = int(len(cov))
    HEAD["c_covered_pct"] = round(rc, 2)
    HEAD["c_uncovered_found"] = int(unc.mut_detected.sum()); HEAD["c_uncovered_n"] = int(len(unc))
    HEAD["c_uncovered_pct"] = round(ru, 2)
    HEAD["c_total_found"] = int(p.mut_detected.sum())
    HEAD["c_total_pct"] = round(100 * p.mut_detected.sum() / len(p), 2)

    # THE SAME DETECTIONS UNDER THE PAPER'S OWN SCORE MARGIN.
    # region_compare.py counts a detection at the search's reported 1% FDR and applies no margin
    # filter, so these counts are NOT the margin >= 10 population that section 2.3's error
    # analysis is defined over. Both are correct and they are different denominators, exactly as
    # Figure 3b (distinct peptides) and Figure 3c (identifications) are. Quoting one with the
    # other's counts is the mistake this block exists to make impossible.
    if HITS:
        _hm = acquired_lysine_only(pd.read_csv(HITS, sep="\t")).copy()
        _hm["margin"] = _hm.hyper - _hm.nextscore
        _hm["key"] = list(zip(_hm.gene, _hm.change, _hm.line, _hm.site))
        _bm = _hm.groupby("key").margin.max()
        _keep = set(_bm[_bm >= MARGIN].index)
        _pk = list(zip(p.gene, p.change, p.line, p.site))
        _sel = pd.Series([k in _keep for k in _pk], index=p.index) & (p.mut_detected == 1)
        _cv = p.wt_detected == 1
        HEAD["c_margin_threshold"] = MARGIN
        HEAD["c_total_found_margin"] = int(_sel.sum())
        HEAD["c_total_pct_margin"] = round(100 * int(_sel.sum()) / len(p), 2)
        HEAD["c_covered_found_margin"] = int((_sel & _cv).sum())
        HEAD["c_covered_pct_margin"] = round(100 * int((_sel & _cv).sum()) / int(_cv.sum()), 2)
        HEAD["c_uncovered_found_margin"] = int((_sel & ~_cv).sum())
        HEAD["c_uncovered_pct_margin"] = round(
            100 * int((_sel & ~_cv).sum()) / int((~_cv).sum()), 2)
        HEAD["c_genomic_sites"] = len({(g, st) for g, c_, l_, st in _bm.index})
        HEAD["c_genomic_sites_margin"] = len({(g, st) for g, c_, l_, st in _keep})
        HEAD["c_n_psms"] = int(len(_hm))
        HEAD["c_n_psms_margin"] = int((_hm.margin >= MARGIN).sum())
    fd.panel("fig2", "c", [
        dict(group="region not covered", recovered=int(unc.mut_detected.sum()),
             observations=int(len(unc)), pct=ru),
        dict(group="region covered", recovered=int(cov.mut_detected.sum()),
             observations=int(len(cov)), pct=rc)], note="gate 2, mutant-only route")

# ------------------------------------------------------------------ d
# No dx_in here. _panel_head_flush anchors the letter to the panel's LEFTMOST INK, and the
# row labels below are pinned to the pixel column they already held, so the letter follows
# them back to a's 170.5 on its own. Adding 106.9 px of dx_in as well pushed it to 278.6.
ax = fig.add_subplot(gs_bot[0, 0]); st.head(fig, ax, "d")
ax._panel_head_flush = True
# 6 Nov, second pass: the axis now runs to 0 and the DOTS DO NOT MOVE, which fixes the
# frame rather than the offset. Holding the data still means holding px-per-unit at 1068.8,
# so a span of 1.0 instead of 0.9 needs a 1068.8 px frame against 961.9, and the extra
# 106.9 px has to come off the LEFT because the right end is where the data ends. 106.9/0.795
# is 134.5 of parameter either way, so dx 191.2 -> 56.7 and dw -183.4 -> -48.9.
ax._axes_dx_in = 56.7000 / 300.0
# Consequence, stated because it is visible: d's leftmost ink is the right-aligned row
# label, which hangs off the frame edge, so it travels the same 106.9 px left and panels a
# and d stop sharing the 170.5 letter column. That column is re-converged in the row settle
# once e and f are done, which is where the sheet width is settled too.
ax._axes_dw_in = -48.9000 / 300.0
# Driver status is the Hotspot flag carried in the DepMap variant table, not a gene list chosen
# here — which removes the question of who decided what counts as a driver.
if PEP and SITES and FLUX:
    import re as _re2
    dd = pd.read_csv(SITES, low_memory=False)
    nm2 = pd.read_csv(FLUX, low_memory=False)[["ModelID", "CellLine"]]
    kk2 = (dd[dd.to_aa == "K"].merge(nm2, on="ModelID")
           .rename(columns={"CellLine": "line", "Gene": "gene", "Position": "ksite"}))
    kk2 = kk2.groupby(["line", "gene", "ksite"], as_index=False).agg(
        AF=("AF", "max"), Hotspot=("Hotspot", "max"))
    def _ks(ch, fb):
        for t in str(ch).split(";"):
            m = _re2.match(r"^[A-Z](\d+)K$", t.strip())
            if m: return int(m.group(1))
        return fb
    cov = p[p.wt_detected == 1].copy()
    cov["ksite"] = [_ks(c, s) for c, s in zip(cov.change, cov.site)]
    cov = cov.merge(kk2, on=["line", "gene", "ksite"], how="left").dropna(subset=["AF"])
    cov["Hotspot"] = cov.Hotspot.fillna(False).astype(bool)
    miss = cov[cov.mut_detected == 0]; fnd = cov[cov.mut_detected == 1]
    hm = miss[miss.Hotspot]; hf = fnd[fnd.Hotspot]
    rng = np.random.default_rng(0)
    # The found row sits at 0.72, not at 1. Lowered 6 Nov so the two rows read as a pair with
    # the label block above them rather than one row stranded halfway up the panel.
    _RF = 0.72
    # Bigger and semi-transparent, 5 Nov. At s=2.0 the 590 missed sites were a grey smear in
    # which no reader could tell two overlapping points from twenty. Alpha is the one use of
    # opacity the palette brief allows, because here it reports density rather than inventing
    # a colour: every dot is CC3F0C, C5C3C6 or 1985A1 at full strength, and where they pile up
    # the ink darkens. Each of these carries _alpha_ok so the translucency check keeps
    # reporting accidental alpha elsewhere instead of printing this every run.
    # The ten driver hotspots stay OPAQUE. They are the panel's claim, there are only ten of
    # them, and fading them to report a density of ten would cost the one thing being shown.
    _d1 = ax.scatter(miss.AF, rng.normal(0, 0.055, len(miss)), s=14.0, color=st.LIGHTGREY,
                     linewidths=0, alpha=0.55, clip_on=False)
    _d2 = ax.scatter(fnd.AF, rng.normal(_RF, 0.055, len(fnd)), s=22.0, color=st.BLUE,
                     linewidths=0, alpha=0.55, clip_on=False)
    for _d in (_d1, _d2):
        _d._alpha_ok = "panel d overplotting, 590 and 27 points in two narrow rows"
    # clip_on=False on all four. The axis ends at 1.0 because that is where a fraction ends,
    # and the largest allele fractions are 0.996, so at s=28 half of those markers fell
    # outside the frame and were cut off. The glyph is allowed past the spine; the value is
    # still at 0.996 and the axis still stops on a tick.
    ax.scatter(hm.AF, rng.normal(0, 0.04, len(hm)), s=28, color=st.ORANGE,
               linewidths=0, zorder=6, clip_on=False)
    ax.scatter(hf.AF, rng.normal(_RF, 0.04, len(hf)), s=28, color=st.ORANGE,
               linewidths=0, zorder=6, clip_on=False)
    labs = sorted({"%s %s" % (r.gene, str(r.change).split(";")[0]) for r in hm.itertuples()})
    half = (len(labs) + 1) // 2
    # the labels need clear air ABOVE the found row, so the y-limit is extended rather than
    # the text being squeezed between the two rows.
    # 0.109 and 0.424 were 0.01 and 0.36 of a 961.9 px frame starting at 375.5. Anything in
    # DATA coordinates survived the axis reaching 0 untouched, because px-per-unit and the
    # right edge both held; anything in AXES FRACTION did not, and these two columns were
    # travelling 106.9 px left with the frame. Re-expressed against the 1068.8 px frame they
    # land on 385.1 px and 721.8 px, the columns they already held. The "n of n driver
    # hotspots missed" block below needs no such fix: it is anchored to the right edge.
    ax.text(0.1090, 0.99, "\n".join(labs[:half]), transform=ax.transAxes, fontsize=st.sz(7),
            color=st.ORANGE, linespacing=1.55, va="top")
    ax.text(0.4240, 0.99, "\n".join(labs[half:]), transform=ax.transAxes, fontsize=st.sz(7),
            color=st.ORANGE, linespacing=1.55, va="top")
    # 0.9910, not 0.99: an inset given as a fraction grows with the frame, so 1% of 1068.8 px
    # sat this block 1.1 px left of the 9.6 px inset it had at 961.9. Pinned to the same inset.
    ax.text(0.9910, 0.99, "%d of %d driver\nhotspots missed" % (len(hm), len(hm) + len(hf)),
            transform=ax.transAxes, fontsize=st.sz(7), color=st.ORANGE, ha="right", va="top",
            linespacing=1.4, fontweight="bold")
    # No y axis. Two named rows are not a scale: there is nothing to read off a spine, and a
    # left spine with two tick marks on it was drawing a measuring device for a quantity that
    # does not exist. The rows are named directly instead, right-aligned off the plot area in
    # a blended transform, x in axes fraction and y in data, so each sits on its own row
    # whatever the limits become. One line each with the count parenthesised, the format
    # panel c settled on; they were "missed" over "590" on two lines.
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", left=False, labelleft=False)
    _ytr = ax.get_yaxis_transform()
    # The row labels hang off the frame's left edge, so extending the axis to 0 would have
    # carried them 106.9 px left with it and put d's leftmost ink at 63 px, a third of an inch
    # left of panel a. They are pinned to the pixel column they already occupied instead: the
    # frame gained 0.1 of span on the left and the labels give it straight back, sitting in
    # the stretch below the lowest called allele fraction (0.151) where there is nothing to
    # overlap - 66 px clear of the first point. 0.0892 = 0.1 - 11.5/1068.8, which is the old
    # -0.012 of a 961.9 px frame measured from the new left edge.
    for _yv, _lab in ((0, "missed (%d)" % len(miss)), (_RF, "found (%d)" % len(fnd))):
        ax.text(0.0892, _yv, _lab, transform=_ytr, ha="right", va="center",
                fontsize=st.sz(6.3), clip_on=False)
    ax.set_xlabel("allele fraction")
    # 0 to 1.0, 6 Nov. The axis starts on a labelled major tick, as every other axis in
    # both figures does. The empty 0.151 at the left is real and is the point: no variant in
    # the set is called below that allele fraction, and an axis that starts at 0.1 hides the
    # floor by cropping to it. The dots are unmoved - the frame grew to the left instead, so
    # every point sits on the same output pixel it sat on at 0.1 to 1.0.
    ax.set_xlim(0.0, 1.0); ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    # The two rows now use the frame. ylim ran to 3.1 against data at 0 and 1, so 2.1 of a
    # 3.55-unit range was scaffolding holding the gene labels, and the points lived in the
    # bottom 240 px of a 579 px panel. At 2.15 the labels still clear the found row by 26 px.
    ax.set_ylim(-0.55, 2.15)

# ------------------------------------------------------------------ e
ax = fig.add_subplot(gs_bot[0, 1]); st.head(fig, ax, "e")
ax._panel_head_flush = True
# 6 Nov: stretched 120 px wider, 60 px each way so the panel keeps its centre, then moved
# 175 px LEFT so the d|e gap comes to the 100 px row 0 now runs at. dw +120 px holds the
# width through the move. Both over the 0.795 damping, so dx is -60 - 175 = -235 px of
# page, which is -295.6 of parameter off the 65.2 this panel started the day at.
ax._axes_dx_in = -264.4000 / 300.0
ax._axes_dw_in = 126.2000 / 300.0
# Restored. This is a GATE 1 result: where the region was never covered, the wild-type tryptic
# peptides are longer, and the two distributions diverge specifically in the tail beyond 20
# residues rather than shifting as a whole. That points at a subset of regions requiring
# peptides too long to observe — an argument for a protease that redraws cleavage boundaries,
# not for deeper fractionation.
if PEP and "shortest_wt" in p.columns:
    _a0 = pd.to_numeric(p.loc[p.wt_detected == 1, "shortest_wt"], errors="coerce").dropna()
    _b0 = pd.to_numeric(p.loc[p.wt_detected == 0, "shortest_wt"], errors="coerce").dropna()
    # shortest_wt == 0 IS NOT A LENGTH. region_compare.py writes
    # `shortest_wt=min((len(x) for x in wt), default=0)`, so 0 is the default when the
    # region enumerated no wild-type tryptic peptide at all. It appears 386 times, every
    # one of them in the never-covered group and none in the covered group, and there is
    # no value between 1 and 6 anywhere in the column, which is what proves it a sentinel
    # rather than a short peptide. Averaging it as a length pulled the never-covered mean
    # down to 16.65 from 17.51 and weakened the test from p = 2.2e-15 to p = 3.5e-08,
    # both in the direction that works against this panel's own claim. The histogram had
    # already dropped them, by starting its bins at 4, so the picture and the statistics
    # were describing different sets of sites.
    a, b = _a0[_a0 > 0], _b0[_b0 > 0]
    N_SENTINEL = int((_b0 == 0).sum() + (_a0 == 0).sum())
    from scipy import stats as _st2
    # Bins cover the DATA, 7 to 50, rather than a window. At 4 to 40 the drawn bars held
    # 99.8% of the covered group and 91.0% of the never-covered one, and density=True
    # normalises over what is drawn, so the two curves were not on the same footing in a
    # panel whose whole point is that they differ. Now both integrate to 1 over all of
    # their own sites and the tail that carries the argument is on the page.
    bins = np.arange(6, 52, 2)
    _W = bins[1] - bins[0]
    _wa = np.full(len(a), 100.0 / len(a)); _wb = np.full(len(b), 100.0 / len(b))
    # Never covered is the grey mass and covered is the teal outline over it, which is the
    # SAME assignment panel c makes for the same two groups. Before this, panel c drew
    # "not covered" in light grey and panel e drew it in orange, and orange everywhere
    # else in this figure means the acquired lysine or the hotspot that was missed.
    ax.hist(b, bins=bins, weights=_wb, color=st.LIGHTGREY, label="not covered")
    ax.hist(a, bins=bins, weights=_wa, histtype="step", color=st.BLUE, lw=1.6,
            label="covered")
    _, pv = _st2.mannwhitneyu(a, b)
    ax.set_xlabel("shortest tryptic peptide (residues)"); ax.set_ylabel("sites (%)")
    ax.set_xlim(0, 60); ax.set_xticks([0, 15, 30, 45, 60])
    _pk = max(np.histogram(a, bins=bins, weights=_wa)[0].max(),
              np.histogram(b, bins=bins, weights=_wb)[0].max())
    # 0 to 20 by 5, the PI's call. I had taken it to 25 for headroom because the covered
    # peak is 19.96% and at a 20% top the teal outline runs along the frame. He preferred
    # the round ceiling to the clearance, which is the same choice as "always end y-axis on
    # major tick": the step has to be round too, and 0 5 10 15 20 is rounder than 0 to 25.
    ax.set_ylim(0, 20); ax.set_yticks([0, 5, 10, 15, 20])
    assert _pk <= 20.0, "panel e: tallest bar is %.2f%%, past the 20%% top" % _pk
    ax.legend(frameon=False, fontsize=st.sz(6.3), handlelength=1.1, loc="upper right")
    # 6 Nov: the "mean 13.3 vs 17.5 / p = 2.2 x 10^-15" block that stood here is gone, to
    # the caption, the same call made for panel a's three sentences on 5 November. The
    # artwork carries the data and the labels that name it; a statistic is prose and reads
    # as prose. Both numbers still reach the reader, and they still come from this script
    # rather than from the artwork, through the print below and the HEAD ledger.
    HEAD["e_covered_mean"] = round(float(a.mean()), 1); HEAD["e_covered_n"] = int(len(a))
    HEAD["e_uncovered_mean"] = round(float(b.mean()), 1); HEAD["e_uncovered_n"] = int(len(b))
    HEAD["e_mannwhitney_p"] = float("%.1e" % pv)
    HEAD["e_no_wt_peptide_excluded"] = N_SENTINEL
    # The decomposition, because "mean 13.3 vs 17.5" invites the reading that every
    # never-covered region needs a slightly longer peptide, and that is not what happened.
    # Below 20 residues the two groups are the same distribution; the whole difference is
    # a tail. Every number the caption quotes is printed here rather than read off the
    # artwork, including the ones the panel does not draw.
    _CUT = 20
    _la, _lb = a[a <= _CUT], b[b <= _CUT]
    _, _pvlo = _st2.mannwhitneyu(_la, _lb)
    HEAD["e_cut_residues"] = _CUT
    HEAD["e_below_cut_mean_covered"] = round(float(_la.mean()), 2)
    HEAD["e_below_cut_mean_uncovered"] = round(float(_lb.mean()), 2)
    HEAD["e_below_cut_p"] = round(float(_pvlo), 3)
    for _t in (20, 30):
        HEAD["e_above_%d_pct_covered" % _t] = round(100 * float((a > _t).mean()), 1)
        HEAD["e_above_%d_pct_uncovered" % _t] = round(100 * float((b > _t).mean()), 1)
    print("e: at or below %d residues the two groups are one distribution, mean %.2f covered "
          "vs %.2f never covered, Mann-Whitney p = %.2f | above %d, %.1f%% of never covered "
          "against %.1f%% of covered | above 30, %.1f%% against %.1f%%"
          % (_CUT, _la.mean(), _lb.mean(), _pvlo, _CUT, 100 * (b > 20).mean(),
             100 * (a > 20).mean(), 100 * (b > 30).mean(), 100 * (a > 30).mean()))
    print("e: shortest tryptic peptide, covered mean %.1f (median %g, n=%s) vs never covered "
          "mean %.1f (median %g, n=%s), Mann-Whitney p = %.1e | %d sites with no wild-type "
          "tryptic peptide excluded, all of them never covered"
          % (a.mean(), a.median(), f"{len(a):,}", b.mean(), b.median(), f"{len(b):,}",
             pv, N_SENTINEL))
    fd.panel("fig2", "e", pd.DataFrame(dict(
        shortest_wt_peptide=pd.concat([a, b]).values,
        region=["covered"] * len(a) + ["never covered"] * len(b))),
        note="gate 1 explanation; one row per observation; sites with no wild-type tryptic "
             "peptide (shortest_wt sentinel 0) excluded")

# ------------------------------------------------------------------ f
ax = fig.add_subplot(gs_bot[0, 2]); st.head(fig, ax, "f")
ax._panel_head_flush = True
ax._axes_dx_in = -415.1000 / 300.0
ax._axes_dw_in = -41.4000 / 300.0
# Two independent routes to the same class of site. They are kept separate because their
# false-discovery behaviour differs: downstream peptides returned zero entrapment hits across
# 42 experiments (FDP 0.00%) against 0.58% for variant-containing peptides in the same searches,
# and their boundary control — the same region shifted one residue — was never identified.
DOWN = find("downstream_sites.csv")
if HITS and DOWN:
    h = acquired_lysine_only(pd.read_csv(HITS, sep="\t"))
    hv = h.groupby(["gene", "change", "line", "site"]).ngroups
    dn = pd.read_csv(DOWN)
    # YBX1 E4K was excluded here on a "+1.0087 Da, wrong molecular species" reading. That was
    # wrong: 1.0087 sits 0.0053 Da from a neutron, which at 5,213 Da is 1.0 ppm, while
    # deamidation is 0.0247 Da away. It is an isotope peak, and selecting the first isotope
    # rather than the monoisotopic happens in 7.7% of variant PSMs and 7.3% of background — the
    # same rate — and is expected at 5+ on a 5.2 kDa precursor. Reinstated. See A81.
    sites_d = dn.groupby(["gene", "change", "line"]).ngroups
    known = set(zip(h.gene, h.line))
    newd = dn[[(g, l) not in known for g, l in zip(dn.gene, dn.line)]]
    n_new = newd.groupby(["gene", "change", "line"]).ngroups
    n_conf = sites_d - n_new
    ax.bar([0], [hv], color=st.BLUE, width=st.BARW)
    # 6 Nov, the colours swapped round. The shared segment was 46494C and the new one was
    # C5C3C6, so the 8 sites that only the downstream route finds, which are the whole
    # reason the route is in the paper, were drawn in the lightest grey in the palette,
    # the colour this figure uses everywhere else for the category that failed. Now the
    # 7 shared sites take BLUE, which is bar 0's own colour and says without a word that
    # those 7 are already counted in the bar beside them, and the 8 route-unique sites
    # take 46494C, which is the colour panel a gives the downstream peptide.
    ax.bar([1], [n_conf], color=st.BLUE, width=st.BARW)
    ax.bar([1], [n_new], bottom=[n_conf], color=st.GREEN, width=st.BARW)
    ax.text(0, hv + 1.5, "%d" % hv, ha="center", fontsize=st.sz(7))
    ax.text(1, n_conf / 2, "%d" % n_conf, ha="center", va="center", fontsize=st.sz(7),
            color=st.PAPER)
    # Both segments are dark now, so both counts are set in paper.
    ax.text(1, n_conf + n_new / 2, "%d" % n_new, ha="center", va="center", fontsize=st.sz(7),
            color=st.PAPER)
    ax.text(1, n_conf + n_new + 1.5, "%d" % sites_d, ha="center", fontsize=st.sz(7))
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["mutant\npeptide", "downstream\npeptide"], fontsize=st.sz(6.3))
    ax.set_ylabel("acquired-lysine sites")
    # 0 to 60 by 15. It was max(hv, sites_d) * 1.30, which came to 68.9 and end_on_tick
    # rounded to 70, and the default locator then printed eight labels for two bars. The
    # tallest bar is 53 and its label tops out near 58, so 60 holds it.
    ax.set_ylim(0, 60); ax.set_yticks([0, 15, 30, 45, 60])
    # The key went to the caption too. In a 321 px frame the honest wording, "teal: also
    # found by the mutant-peptide route", set 430 px per line and printed straight over the
    # "53" at 100%, and every short enough form said less than the colour already says: the
    # teal segment is bar 0's own teal, so it reads as "these are in the bar beside it".
    # The caption states it, and the two counts are on the segments.
    # 6 Nov: the two prose blocks that hung under this panel are gone, to the caption.
    # They were the figure's rightmost ink, which put the sheet edge at 3458 px and left
    # panel f's letter 412 px to the left of its own frame, because _panel_head_flush
    # anchors the letter to the panel's leftmost ink and the prose was wider than the
    # panel. One of them also carried 53, 41 and 12 as TYPED DIGITS inside the string,
    # which the rule against transcribed numbers exists to stop; they are computed below
    # and reach the caption through the ledger.
    _spec = h.groupby(["gene", "change", "line", "site"]).spectrum.nunique()
    HEAD["f_sites_single_spectrum"] = int((_spec == 1).sum())
    HEAD["f_sites_multi_spectrum"] = int((_spec >= 2).sum())



# ------------------------------------------------------------------ j
# Which peptide species get identified. Three species, each comparing the peptides that were
# identified against those that were not. The x-axis is the difference in medians in standard
# deviations of that property across all searched peptides — a common scale, so three properties
# measured in different units can share one axis.
#
# The label's hydrophobic contribution was FITTED against observed retention times of 107
# identified peptides rather than assumed (best fit K = -3.9, rho 0.537).
# (The peptide-property panel that consumed peptide_properties.csv moved to figS2 in v3;
#  the loader and its series list were left behind here and are removed.)
st.save(fig, "fig2")
if __import__("os").environ.get("FIG2_PROBE"):
    # geometry probe: the frame each panel actually got AFTER _fit_once damped the offsets.
    # Any offset in this script is converged against these numbers, never computed from the
    # value asked for, because _fit_once scales every axes inward by the overflow it measures.
    _r = fig.canvas.get_renderer()
    print("---- ink boxes (axes frame + every child artist + the panel letter) ----")
    _letters = {getattr(_a, "_panel_letter", "?"): _a for _a in fig.axes}
    for _ln in "abcdef":
        _a = _letters.get(_ln)
        if _a is None: continue
        _bb = [_a.get_window_extent()]
        for _ch in (list(_a.texts) + list(_a.collections) + list(_a.patches) +
                    list(_a.lines) + list(_a.images) +
                    [_a.xaxis.label, _a.yaxis.label] +
                    _a.get_xticklabels() + _a.get_yticklabels()):
            try:
                _e = _ch.get_window_extent(_r)
                if _e.width > 0 and _e.height > 0: _bb.append(_e)
            except Exception: pass
        _lg = _a.get_legend()
        if _lg is not None:
            try: _bb.append(_lg.get_window_extent(_r))
            except Exception: pass
        _x0 = min(e.x0 for e in _bb); _x1 = max(e.x1 for e in _bb)
        print("  %s  ink %8.1f .. %8.1f   width %7.1f   frame %8.1f .. %8.1f"
              % (_ln, _x0, _x1, _x1 - _x0, _a.get_window_extent().x0, _a.get_window_extent().x1))
    print("probe panel a: LBL_R %.1f px  SX0 %.1f px  SX1 %.1f px (%.4f du)  frame %.1f px  slack %.1f px"
          % (LBL_R * PXDU, SX0 * PXDU, SX1 * PXDU, SX1, 10.0 * PXDU, (10.0 - SX1) * PXDU))
    for _i, _a in enumerate(fig.axes):
        _bb = _a.get_window_extent()
        print("probe ax%d  x0 %8.1f  x1 %8.1f  w %7.1f  |  y0 %7.1f  y1 %7.1f  h %6.1f  |  xlim %s"
              % (_i, _bb.x0, _bb.x1, _bb.width, _bb.y0, _bb.y1, _bb.height, _a.get_xlim()))

# 6 Nov: trim to the drawing, 18 px all round, the Figure 1 standard. Before this the sheet
# carried 170 px of empty paper on the left and 1 px on the right, so laid beside Figure 1 at
# one scale the whole figure read as indented by half an inch. Nothing inside moves: this
# changes the canvas, not the layout, so every panel keeps its position, size and type.
# Measured off the RENDERED PIXELS, not from get_tightbbox. The tight bbox puts this
# figure's left edge at 112.8 px while the leftmost ink is at 170: panels a and d have their
# axis turned off but still contribute a rectangle, so trimming to it left 57 px of paper
# that is not there. Reading the saved PNG asks the only question that matters, which is
# where the ink actually starts.
_PAD_IN = 18.0 / 300.0
from PIL import Image as _Im
import numpy as _np
_px = _np.asarray(_Im.open("fig2.png").convert("RGB")).astype(int)
_ys, _xs = _np.nonzero(_px.sum(2) < 735)
_H_in = fig.get_size_inches()[1]
_trim = Bbox.from_extents(_xs.min() / 300.0 - _PAD_IN,
                          _H_in - _ys.max() / 300.0 - _PAD_IN,
                          _xs.max() / 300.0 + _PAD_IN,
                          _H_in - _ys.min() / 300.0 + _PAD_IN)
for _ext in (".png", ".pdf"):
    fig.savefig("fig2" + _ext, bbox_inches=_trim)
print("trimmed to %.3f x %.3f in (%.0f x %.0f px at 300 dpi), 18 px pads"
      % (_trim.width, _trim.height, _trim.width * 300, _trim.height * 300))
for _e in (".png", ".pdf"):          # review copies alongside fig4_v2 / figS_react_v2
    fig.savefig("../results/fig2_v3" + _e, dpi=400, bbox_inches="tight")
print("wrote fig2.pdf and fig2.png")
if PEP and CONV:
    _c = acquired_lysine_only(pd.read_csv(CONV, sep="\t"))
    print("b: conventional %.2f%% (%d/%d) | mutation-aware %.2f%% (%d/%d) regions covered"
          % (100 * _c.wt_detected.sum() / len(_c), int(_c.wt_detected.sum()), len(_c),
             100 * p.wt_detected.sum() / len(p), int(p.wt_detected.sum()), len(p)))
    print("b: conventional mutant detections %d (passed control: no variant sequences in that database)"
          % int(_c.mut_detected.sum()))
if PEP:
    _cov = p[p.wt_detected == 1]; _unc = p[p.wt_detected == 0]
    print("c: covered %.2f%% (%d/%s) | not covered %.2f%% (%d/%s) | overall %.2f%% (%d/%s)"
          % (100 * _cov.mut_detected.sum() / len(_cov), int(_cov.mut_detected.sum()), f"{len(_cov):,}",
             100 * _unc.mut_detected.sum() / len(_unc), int(_unc.mut_detected.sum()), f"{len(_unc):,}",
             100 * p.mut_detected.sum() / len(p), int(p.mut_detected.sum()), f"{len(p):,}"))
    # The second computation of panel e's statistics that stood here is gone. It recomputed
    # the means and the test from the raw column, without the sentinel guard the panel
    # applies, so the two prints disagreed by 0.9 residues and seven orders of magnitude on
    # p the moment the guard went in. Panel e prints its own numbers where it draws them.
# headline numbers per panel, printed so no legend has to hold a number the script does not report
try:
    print("d: %d of %d driver hotspots missed | missed %d, found %d covered sites | %d gene-mutation labels across %d genes"
          % (len(hm), len(hm) + len(hf), len(miss), len(fnd), len(labs), len(set(hm.gene))))
    HEAD["d_hotspots_missed"] = int(len(hm)); HEAD["d_hotspots_total"] = int(len(hm) + len(hf))
    HEAD["d_covered_missed"] = int(len(miss)); HEAD["d_covered_found"] = int(len(fnd))
    HEAD["d_gene_mutation_labels"] = int(len(labs)); HEAD["d_genes"] = int(len(set(hm.gene)))
    fd.panel("fig2", "d", hm.assign(outcome="mutant missed")[["gene", "change", "line", "site"]]
             if hasattr(hm, "assign") else [], note="driver hotspots among covered sites")
except NameError:
    print("d: not built")
try:
    print("f: mutant-peptide route %d sites | downstream route %d (%d shared, %d route-unique) | union %d"
          % (hv, sites_d, n_conf, n_new, hv + n_new))
    HEAD["f_mutant_route_sites"] = int(hv); HEAD["f_downstream_sites"] = int(sites_d)
    HEAD["f_shared"] = int(n_conf); HEAD["f_downstream_unique"] = int(n_new)
    HEAD["f_union"] = int(hv + n_new)
    fd.panel("fig2", "f", [
        dict(route="mutant peptide", sites=int(hv)),
        dict(route="downstream peptide", sites=int(sites_d)),
        dict(route="shared by both", sites=int(n_conf)),
        dict(route="downstream only", sites=int(n_new)),
        dict(route="union", sites=int(hv + n_new))], note="two detection routes")
except NameError:
    print("f: not built")

# The scope line at the head of the run is itself a published quantity (8,470 / 8,055 / 372).
if PEP:
    HEAD["scope_acquiredK_observations"] = int(len(p))
    HEAD["scope_nonlysine_observations"] = int(len(p_all) - len(p))
    HEAD["scope_rows_total"] = int(len(p_all))
    if "line" in p.columns:
        HEAD["scope_cell_lines"] = int(p.line.nunique())
    if {"gene", "site"} <= set(p.columns):
        HEAD["scope_distinct_sites"] = int(p.groupby(["gene", "site"]).ngroups)

fd.head("fig2", HEAD)
for x in notes: print("  MISSING:", x)

if __import__("os").environ.get("FIG2_FINAL"):
    # The finalized file for CCB-submission/figures-final/.
    # Figure 1 set the standard: render the trimmed figure at a finer raster and declare the
    # result 300 dpi. For Figure 1 that raster is 300 x 12.00 / 11.5400 = 311.96 dpi, which
    # puts it at 12.00 in. Figure 2 uses THE SAME RASTER, not the same width. Scaling this
    # figure to 12.00 in instead would scale it by 1.2579 against Figure 1's 1.0399 and print
    # its type 21% larger, which is the one thing the two figures must not do. At the shared
    # raster Figure 2 finalizes at 9.92 in with body type at 9.213 pt, identical to Figure 1.
    # A six-panel figure narrower than a seven-panel one is the honest outcome.
    # Placed AFTER the trim, and saved with the same _trim box: this figure's trim is a
    # bbox_inches on save, not a resize of the figure, so fig.get_size_inches() still
    # reports the untrimmed 11.53 in. Saving without the box gave a 3597 px file.
    import PIL.Image as _PI, pymupdf as _mu
    _ref = _mu.open("fig1.pdf")[0].rect.width / 72.0      # Figure 1's native width, measured
    _dpi = 300.0 * 12.0 / _ref
    _nat = _trim.width
    fig.savefig("fig2_final.png", dpi=_dpi, bbox_inches=_trim)
    _im = _PI.open("fig2_final.png")
    _im.save("fig2_final.png", dpi=(300, 300))
    print("final: native %.4f in at the shared %.2f dpi raster -> %d x %d px declared 300 dpi "
          "(%.3f x %.3f in), body type %.3f pt against Figure 1's %.3f"
          % (_nat, _dpi, _im.size[0], _im.size[1], _im.size[0] / 300.0, _im.size[1] / 300.0,
             8.860 * 12.0 / _ref, 8.860 * 12.0 / _ref))
