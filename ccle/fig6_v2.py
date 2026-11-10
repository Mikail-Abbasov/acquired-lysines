#!/usr/bin/env python3
"""
Figure 6, Engagement is a property of the lysine, not of the molecule.

DESIGN RULE FOR THIS FIGURE (set 2026-09-02, after review)
No sentences inside panels. A panel carries axis labels, tick labels, series labels and numbers
attached to marks, nothing a reader has to stop and read. Everything explanatory belongs in the
legend. The previous version had a three- to five-line grey caption in five of the six panels;
they are gone, and FIGURE_LEGENDS_fig6.md carries what they said.

PANEL a, the assay, then the design, as one flow
  Two different things were being labelled with two headings and drawn the same way, so the panel
  read as one idea stated twice. It is now left to right:
      competition  ->  occupancy distribution with the 60% cut  ->  within one compound  ->  pooled
  WHAT IS MEASURED is occupancy: one number per (lysine, compound), 100 x (1 - compound/DMSO),
  the share of probe labelling the compound competes away. THE DESIGN is what is done with those
  numbers: a compound's engaged sites are compared with ITS OWN below-cut sites, standardised
  inside that compound, and only then pooled across compounds.
  The distribution is every measurement in all three screens. An earlier draft drew one compound
  (27d), picked by sorting on engaged-site count, selection on the displayed quantity. 27d is the
  99th percentile for engagement rate, so it implied a hit rate ~20x the truth.

PANEL d, the null as a distribution, not a count
  Three bars, one of them zero, asked the reader to read an absence off blank space. Every
  correlation is now plotted: the cloud contracts from |rho| 0.45 pooled to 0.24 within NatChem to
  0.18 across the 173 slopes, and the significant points stop existing. The within-NatChem values
  are reproduced by step37b_save_within.py, which step 37 computes and never saved.

PANEL e, a real axis
  The plotted quantity is the Euclidean length of a K-component vector of SD-unit mean
  differences. It is not itself an SD; with K = 15 a uniform 0.1 SD shift reads 0.39. Divided by
  sqrt(K) it is the average shift per feature, which is judgeable. K is read from step 40's own
  feature list so a changed feature set rescales the axis instead of silently misreporting it.

PANEL f, the confound AND its consequence
  Left: molecular weight all but names the screen. Right: the descriptors that separate the
  screens are the descriptors that produced pooled positives (rho 0.69, p 0.027, step 41).
  This replaces "every one of the 21 pooled hits involves a SIZE descriptor", which was wrong
  twice, it called rotatable-bond count and fsp3 "size", and it described step 31 (150 pairs)
  while Fig. 6d plots step 37 (84 pairs). Both return 21 by coincidence; step 37's 21 include six
  from logP.
"""
import os, re, sys, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
from matplotlib.patches import Rectangle, FancyArrowPatch
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import abbasov_style as st
import figure_data as fd
# One body type size for this figure; every other size is a ratio to it.
BODY = 7.0
st.apply(BODY)
# 300, not the style module's 600. THIS SCRIPT NEVER SET IT, so its PNG came out at twice the
# canvas and every pixel measurement taken off it was out by a factor of two. Figures 1 to 5
# all pin it; Figure 4 and Figure 5 each carried this same defect before they came through.
plt.rcParams["savefig.dpi"] = 300

def find(n):
    for d in (".", "..", "ligand_rebuild", "../ligand_rebuild"):
        p = os.path.join(d, n)
        if os.path.exists(p):
            return p
    raise FileNotFoundError(n)

R37 = open(find("step37_report.txt")).read()
R38 = open(find("step38_report.txt")).read()
R40 = open(find("step40_report.txt")).read()
POOL = pd.read_csv(find("step37_within_compound.csv"))
NAT  = pd.read_csv(find("step37_within_NatChem.csv"))
SLOP = pd.read_csv(find("step38_tiers.csv"))
CONF = pd.read_csv(find("step41_confound_strength.csv"))

NICE = lambda f: (f.replace("_8A", " 8Å").replace("_10A", " 10Å").replace("_", " ")
                   .replace("rsa tien", "exposure").replace("propka pKa", "predicted pKa")
                   .replace("nz sasa", "NZ exposure").replace("sasa sidechain", "side-chain exposure")
                   .replace("dist centroid", "distance to centroid").replace("n lys", "nearby lysines"))

# ---------------------------------------------------------------- parse
Q1 = [(m.group(1), float(m.group(2)), float(m.group(3)), int(m.group(4)), float(m.group(5)),
       bool(m.group(6)))
      for m in re.finditer(r"^\s{2,}(\w+)\s+(-?[\d.]+)\s+([\d.]+)\s+(\d+)\s+([\d.e+-]+)\s*(\*)?$",
                           R37[:R37.index("Q2")], re.M)]
N37 = Q1[0][3]
NCONS = sum(1 for q in Q1 if q[5])

TIERS, hdr = [], None
for line in R38.splitlines():
    if line.strip().startswith("tier "):
        hdr = line.split()[2:]
    m = re.match(r"^\s{2,}([\d>=-]+)\s+(\d+)\s+((?:-?[\d.]+\s*)+)$", line)
    if m and hdr:
        TIERS.append((m.group(1), int(m.group(2)), [float(v) for v in m.group(3).split()]))
TIERS = [t for t in TIERS if len(t[2]) == len(hdr)]

NPOOL = int(re.search(r"(\d+) compounds matched to a structure", R37).group(1))
NNAT  = int(re.search(r"NatChem\s+(\d+) compounds \|", R37).group(1))
NSLOP = int(re.search(r"(\d+) compounds with slopes and structures", R38).group(1))

CUTS = []
for m in re.finditer(r"CUT AT (\d+)% — (\d+) sites with zero liganding compounds\n(.*?)-> (NOT )?monotonic",
                     R40, re.S):
    rows = [(r.group(1), int(r.group(2)), float(r.group(3)), float(r.group(4)))
            for r in re.finditer(r"^\s{2,}([\d>=-]+)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+[\d.]+\s*\*?$",
                                 m.group(3), re.M)]
    CUTS.append((int(m.group(1)), rows, m.group(4) is None))

SFLIST = [x.strip().strip("'\"") for x in
          re.search(r"SF=\[(.*?)\]", open(find("step40_threshold_free.py")).read(),
                    re.S).group(1).split(",")]
KFEAT = len(SFLIST)
RT = np.sqrt(KFEAT)

occ = pd.read_csv(find("occ_long.csv"), low_memory=False)
occ["base"] = [re.match(r"^(\d+[a-z])", str(c)).group(1)
               if (d == "NatChem" and re.match(r"^(\d+[a-z])", str(c))) else str(c)
               for c, d in zip(occ.compound, occ.dataset)]
OCCD = {ds: g.occupancy_pct.dropna().values for ds, g in occ.groupby("dataset")}
NOBS = int(sum(len(v) for v in OCCD.values()))
NHIT = int(sum((v >= 60).sum() for v in OCCD.values()))
CPD = occ.groupby(["dataset", "base"]).occupancy_pct.agg(n="size", lig=lambda v: (v >= 60).sum())

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
from scipy import stats as sstats
RDLogger.DisableLog("rdApp.*")
_x = pd.ExcelFile(find("SMILES.xlsx")); _rows = []
for _sh in _x.sheet_names:
    for _r in _x.parse(_sh).dropna(subset=["SMILES"]).itertuples():
        _m = Chem.MolFromSmiles(str(_r.SMILES))
        if _m:
            _rows.append(dict(dataset=_sh, MW=Descriptors.MolWt(_m)))
MW = pd.DataFrame(_rows)
KW = sstats.kruskal(*[v.values for _, v in MW.groupby("dataset").MW])
RHO_CONF, P_CONF = sstats.spearmanr(CONF.eta2, CONF.hits)

assert len(Q1) == 14 and len(TIERS) == 10 and len(CUTS) == 5, (len(Q1), len(TIERS), len(CUTS))

HEAD = {}

# ---- a: the assay
fd.panel("fig6", "a", [dict(dataset=_ds, occupancy_bin_low=float(_e),
                            occupancy_bin_high=float(_e) + 4, n_measurements=int(_n))
                       for _ds, _v in OCCD.items()
                       for _e, _n in zip(np.arange(0, 100, 4),
                                         np.histogram(_v, bins=np.arange(0, 102, 4))[0])],
         note="every site x compound occupancy measurement, binned; the 60% cut is drawn on it")
HEAD["a_n_measurements"] = NOBS
HEAD["a_n_engaged"] = NHIT
HEAD["a_pct_engaged"] = round(100.0 * NHIT / NOBS, 2)
HEAD["a_median_quantified_per_compound"] = int(CPD.n.median())
HEAD["a_median_engaged_per_compound"] = int(CPD.lig.median())

# ---- b: the fourteen site features, within compound
fd.panel("fig6", "b", [dict(feature=_q[0], shift_sd=_q[1], sd=_q[2], n_compounds=_q[3],
                            q=_q[4], consistent=bool(_q[5])) for _q in Q1],
         note="engaged vs the same compound's below-cut sites, standardised within compound")
HEAD["b_n_consistent"] = NCONS
HEAD["b_n_features"] = len(Q1)
HEAD["b_n_compounds"] = N37
for _q in Q1:
    if not _q[5]:
        HEAD["b_exception_%s_shift" % _q[0]] = round(_q[1], 3)
        HEAD["b_exception_%s_q" % _q[0]] = round(_q[4], 3)

# ---- c: the occupancy gradient
fd.panel("fig6", "c", [dict(tier=_t[0], n_observations=_t[1],
                            **{_h: _v for _h, _v in zip(hdr, _t[2])}) for _t in TIERS],
         note="shift from the 0-10%% tier, SD, per site feature")
HEAD["c_n_observations"] = int(sum(_t[1] for _t in TIERS))
HEAD["c_n_tiers"] = len(TIERS)
HEAD["c_n_lowest_tier"] = int(TIERS[0][1])
HEAD["c_n_highest_tier"] = int(TIERS[-1][1])

# ---- d: does any compound descriptor predict site preference?
fd.panel("fig6", "d", [dict(column=_lab.split("\n")[0], descriptor_pair=_i,
                            rho=float(_r.rho), q=float(_r.q))
                       for _lab, _df in (("pooled", POOL), ("NatChem only", NAT),
                                         ("NatChem slopes", SLOP))
                       for _i, _r in enumerate(_df.itertuples())],
         note="every compound-descriptor x site-feature correlation, all three columns")
for _lab, _df, _n in (("pooled", POOL, NPOOL), ("natchem", NAT, NNAT), ("slopes", SLOP, NSLOP)):
    HEAD["d_%s_significant" % _lab] = int((_df.q < 0.05).sum())
    HEAD["d_%s_pairs" % _lab] = int(len(_df))
    HEAD["d_%s_compounds" % _lab] = int(_n)
    HEAD["d_%s_expected_at_0.05" % _lab] = round(0.05 * len(_df), 1)

# ---- e: breadth, at five occupancy cuts
fd.panel("fig6", "e", [dict(cut_pct=_c, compounds_engaging=_r[0], n_sites=_r[1],
                            observed_shift_per_feature=_r[2] / RT,
                            shuffled_shift_per_feature=_r[3] / RT, monotonic=_m)
                       for _c, _rows, _m in CUTS for _r in _rows],
         note="Euclidean length of the K-feature shift vector divided by sqrt(K = %d)" % KFEAT)
_c60 = [c for c in CUTS if c[0] == 60][0]
HEAD["e_broadest_shift_per_feature"] = round(_c60[1][-1][2] / RT, 2)
HEAD["e_broadest_shuffled"] = round(_c60[1][-1][3] / RT, 2)
HEAD["e_n_cuts"] = len(CUTS)
HEAD["e_k_features"] = KFEAT
HEAD["e_cuts_not_monotonic"] = [c[0] for c in CUTS if not c[2]]

# ---- f: the confound
fd.panel("fig6", "f", CONF, note="how strongly each descriptor separates the three screens "
                                 "(eta2) against how many pooled false positives it produced")
for _ds, _v in MW.groupby("dataset").MW:
    HEAD["f_mw_median_%s" % _ds] = round(float(_v.median()), 0)
HEAD["f_kruskal_p"] = "%.1e" % KW[1]
HEAD["f_confound_rho"] = round(float(RHO_CONF), 2)
HEAD["f_confound_p"] = round(float(P_CONF), 3)
HEAD["f_n_descriptors"] = int(len(CONF))

# ---------------------------------------------------------------- figure
# ALL FOUR MARGINS IN INCHES, AND THE CANVAS SOLVED FROM THE PANEL BLOCK, as Figures 1 to 5
# are built. Fractional margins shrink with the canvas and walk the figure into _fit_once,
# which scales every axes at once and takes the y axes off the family's 483.6 px. _PS is the
# sum of the three panel widths and is HELD, so changing wspace closes the gaps instead of
# widening all six panels. ONE gridspec, not two: the old pair pinned each row to its own
# fractional band, which is what left the two rows on different panel heights.
_ML, _MR = 1.60, 1.20          # left and right margins, inches
_MB, _MT = 1.60, 0.80          # bottom and top margins, inches
_PS = 8.245                    # sum of the three panel widths, held. Chosen so the trimmed
                               # sheet is 3461 px, which is fig1.png, fig4.png and fig5.png, so this
                               # figure finalizes at Figure 1's 12.00 in on the shared raster.
_WS, _HS = 0.42, 0.72          # column and row spacing, fractions of the average panel
_YAX = 1.55020                 # the family y axis, 483.6 px at the 311.96 raster
_PW = _PS * (1 + 2 * _WS / 3)
_PH = (2 + _HS) * _YAX
_FW, _FH = _PW + _ML + _MR, _PH + _MB + _MT
fig = plt.figure(figsize=(_FW, _FH))
# (text artist, bar patch) pairs whose x is settled only after st.save has fitted the bar
# thickness and re-centred the groups.
_SNAP = []
_GS = fig.add_gridspec(2, 3, hspace=_HS, wspace=_WS,
                       left=_ML / _FW, right=1 - _MR / _FW,
                       top=1 - _MT / _FH, bottom=_MB / _FH,
                       width_ratios=[1.62, 1.0, 1.0])


# ================================================================ a
# THE VISUAL SYSTEM. Three stages on one band: every stage sits between y = 2.4 and y = 8.4, is
# introduced by a stage label in the same style at y = 9.4, and is joined to the next by the same
# arrow at the same height. Nothing is placed by eye. Panel-a units are ~5.4 x per y, so widths
# below are in x units and heights in y units and they are NOT visually proportional.
ax = fig.add_subplot(_GS[0, 0]); st.head(fig, ax, "a")
ax._panel_head_flush = True
# A slides LEFT to put its letter on D's. Flush letters sit on their panel's leftmost ink, and
# d carries a two-line y-axis title and its tick numbers while a, axis off, starts at its own
# frame, so A printed 192 px right of D.
ax._axes_dx_in = -192.0 / 300.0
ax.set_xlim(0, 100); ax.set_ylim(2.15, 9.75); ax.set_axis_off()

# TWO RAILS, ONE SET OF LYSINES. Side by side, the two strips were two pictures a reader had to
# believe were the same eighteen objects. Stacked, with a thread from each square on the top
# rail to the same square on the bottom one, the panel SHOWS it: sort by occupancy, which makes
# the two groups, then sort the same squares by how buried they are, and watch the teal ones
# sweep right while the grey ones do not. The thread carries its square's colour, so the five
# engaged lysines are the five lines that travel.
# Where occupancy comes from, 100 x (1 - compound / DMSO), is in the caption.
SQW, SQH = 1.90, 0.85
X0, XW = 6.0, 90.0
# The top rail dropped 0.35 to open a second line above it. Six numeric ticks fill the
# first one, and the rail still needs to say what the numbers are.
# The frame is 483.6 px whatever ylim says, so height comes from drawing BELOW it, into the
# hspace where b and c put their tick labels and x titles. The bottom rail dropped a unit,
# which lengthens every thread by a third, and the stack under it follows; the key sits
# 0.95 units under the frame, about 60 px, which is less overhang than b's angled labels.
TOP, BOT = 8.05, 4.00
_px = lambda v: X0 + v * XW / 100.0

# occupancy, then burial, for the same eighteen lysines. Fixed values, not a random generator:
# nothing in a cartoon should look sampled. The two groups OVERLAP on burial, because the shift
# panel b measures is a few tenths of an SD, not a separation.
LYS = [(3, 8), (7, 16), (11, 23), (16, 31), (21, 13), (25, 40), (30, 35), (35, 47), (39, 20),
       (44, 54), (48, 43), (53, 60), (57, 27),
       (63, 46), (70, 58), (77, 68), (85, 79), (93, 88)]
for _occ, _bur in LYS:
    _c = st.BLUE if _occ >= 60 else st.LIGHTGREY
    ax.plot([_px(_occ) + SQW / 2, _px(_bur) + SQW / 2], [TOP - SQH / 2, BOT + SQH / 2],
            color=_c, lw=0.5, zorder=1)
for _occ, _bur in LYS:
    _c = st.BLUE if _occ >= 60 else st.LIGHTGREY
    for _x, _y in ((_px(_occ), TOP), (_px(_bur), BOT)):
        ax.add_patch(Rectangle((_x, _y - SQH / 2), SQW, SQH, facecolor=_c, lw=0, zorder=2))

# ---- top rail: sorted by occupancy, which is what makes the two groups
ax.plot([X0, X0 + XW], [TOP, TOP], color=st.GREY, lw=0.7, zorder=0)
ax.plot([_px(60), _px(60)], [TOP - 0.45, TOP + 0.45], color=st.LIGHTGREY, lw=1.0,
        ls=(0, (1.2, 2.0)), zorder=0)
# A ROUND STEP, with the cut on a tick and both ends on one. Not 6c's 0, 30, 60, >=90: that
# axis is ten occupancy bins and its last one collects everything at 90 and above, while this
# rail is continuous, so a lysine at 93% sits at 93 and a >=90 tick would claim a bin that is
# not there. Six numeric ticks leave no room for the unit on this line, and only one line fits
# above the rail, so "occupancy" moved into the key at the foot.
for _v in (0, 20, 40, 60, 80, 100):
    ax.text(_px(_v) + SQW / 2, TOP + 0.55, "%d" % _v, fontsize=st.sz(6.3), color=st.GREY,
            ha="center", va="bottom")
# THE WORDS PANEL C USES, and they also name the object. Panel c's x axis already reads
# "occupancy of the lysine (%)", and since the two-rail rebuild dropped the "one
# compound's lysines" title this panel had stopped saying what a square is at all.
ax.text(_px(50) + SQW / 2, 9.35, "occupancy of the lysine (%)", fontsize=st.sz(7),
        color=st.GREY, ha="center", va="bottom")

# ---- bottom rail: the same squares, sorted by one site feature
ax.plot([X0, X0 + XW], [BOT, BOT], color=st.GREY, lw=0.7, zorder=0)
ax.text(X0, BOT - 0.72, "exposed", fontsize=st.sz(6.3), color=st.GREY, ha="left", va="top")
ax.text(X0 + XW, BOT - 0.72, "buried", fontsize=st.sz(6.3), color=st.GREY, ha="right", va="top")

# ---- the shift between the two groups on that feature, which is one point in b
_g = [b for o, b in LYS if o < 60]; _e = [b for o, b in LYS if o >= 60]
_gm, _em = sum(_g) / len(_g), sum(_e) / len(_e)
ax.annotate("", xy=(_px(_em) + SQW / 2, 2.60), xytext=(_px(_gm) + SQW / 2, 2.60),
            arrowprops=dict(arrowstyle="<->", color=st.ORANGE, lw=1.1, shrinkA=0, shrinkB=0))
ax.text((_px(_gm) + _px(_em)) / 2 + SQW / 2, 2.35, "engaged vs below-cut (SD)",
        fontsize=st.sz(6.3), color=st.ORANGE, ha="center", va="top")

# ---- which square is which, on one line at the foot
_kx = X0
# The unit went back to the axis title, so the key does not carry it twice.
for _c, _t, _b in ((st.BLUE, "engaged, ≥ 60%", True),
                   (st.LIGHTGREY, "below the cut", False)):
    # Dropped again: at 1.20 the grey swatch sat inside the orange label's box. The gate is
    # blind to that, a swatch being a patch and not text, so it is measured by eye.
    ax.add_patch(Rectangle((_kx, 0.50), SQW, SQH, facecolor=_c, lw=0, clip_on=False))
    ax.text(_kx + SQW + 1.0, 0.92, _t, fontsize=st.sz(6.3), va="center",
            color=st.BLUE if _b else st.GREY, fontweight="bold" if _b else "normal")
    _kx += 34.0      # measured off the render, not estimated. The collision gate cannot see
                     # a swatch landing inside a word, a swatch being a patch, not text.

# ================================================================ b
ax = fig.add_subplot(_GS[0, 1]); st.head(fig, ax, "b")
ax._panel_head_flush = True
ax._axes_dx_in = -137.0 / 300.0    # row 0 re-evened after a was rebuilt     # row 0 evened at 247 px
# WIDER BY 15% OF ITS COLUMN, taken out of row 0's gaps. Fourteen names angled at 45
# degrees sat on a 48.8 px tick pitch against Figure 4's 73 for fifteen of them, so the
# ends of the longer ones ran up to the start of the next. _axes_dw_in holds x0, so the
# panel grows to the right.
ax._axes_dw_in = 0.15 * _PS / 3.62
# TRANSPOSED, the move Figure 4 panels b and c make for the same reason. Fourteen feature
# names down a 483.6 px y axis is 34.5 px a row against a 33 px label, so they touched; angled
# along the x axis at 45 degrees they have the full width. Features run most positive on the
# left, as Figure 4's do, which leaves the bottom left empty for the count.
Q = sorted(Q1, key=lambda q: -q[1])
for i, (f_, sh, sd, n_, q_, sig) in enumerate(Q):
    se = sd / np.sqrt(n_)
    c = st.BLUE if sig else st.ORANGE
    # CAPPED, as Figure 5 panels a and b are. A bare line through a marker reads as a
    # connector rather than as an interval, which is what got asked about in 5a.
    ax.errorbar([i], [sh], yerr=[[1.96 * se], [1.96 * se]], fmt="none", ecolor=c,
                elinewidth=1.3, capsize=2.8, capthick=1.3, zorder=2)
    ax.plot([i], [sh], "o", color=c, ms=4.5, mew=0, zorder=3)
# The family's zero rule, copied from Figure 4 panels b, c, d and f rather than invented:
# light grey, 0.8 pt, dotted on a (1.2, 2.0) pattern. It was solid slate here, which reads
# as a drawn series rather than as the reference the values are measured against.
ax.axhline(0, color=st.LIGHTGREY, lw=0.8, ls=(0, (1.2, 2.0)))
ax.set_xticks(np.arange(len(Q)))
ax.set_xticklabels([NICE(q[0]) for q in Q], fontsize=st.sz(6.3), rotation=45, ha="right",
                   rotation_mode="anchor")
# The exception carries its colour into its own name. One feature of the fourteen points
# the other way, and a rust mark among thirteen teal ones is only findable if the reader
# can see which feature it is without counting across to the axis.
for _tl, _q in zip(ax.get_xticklabels(), Q):
    if not _q[5]:
        _tl.set_color(st.ORANGE)
ax.set_ylabel("engaged vs below-cut\nsites (SD)", fontsize=st.sz(7), linespacing=1.4)
ax.set_xlim(-0.7, len(Q) - 0.3)
# ON A TICK AT BOTH ENDS. The shifts run -0.45 to +0.33 and their intervals a little past
# that, so -0.6 to 0.6 on a step of 0.3 holds everything and reads round.
ax.set_ylim(-0.6, 0.6); ax.set_yticks([-0.6, -0.3, 0.0, 0.3, 0.6])
ax.text(0.02, 0.02, "%d of %d\nn = %d compounds" % (NCONS, len(Q1), N37), transform=ax.transAxes,
        fontsize=st.sz(6.3), color=st.BLUE, va="bottom", linespacing=1.4)
for s in ("top", "right"): ax.spines[s].set_visible(False)

# ================================================================ c
ax = fig.add_subplot(_GS[0, 2]); st.head(fig, ax, "c")
ax._panel_head_flush = True
ax._axes_dx_in = 20.0 / 300.0      # row 0 evened after the frame came in
LAB = {"rsa_tien": "exposure", "nz_sasa": "NZ exposure", "crowding_": "crowding 8Å",
       "atoms_10A": "atoms 10Å", "hydrophob": "hydrophobic 8Å", "aromatic_": "aromatic 8Å"}
COLS = [st.ORANGE, st.ORANGE, st.BLUE, st.BLUE, st.GREEN, st.GREEN]
MK = ["o", "s", "o", "s", "o", "s"]
xs = np.arange(len(TIERS))
for j, h in enumerate(hdr):
    ax.plot(xs, [t[2][j] for t in TIERS], "-", color=COLS[j], lw=1.0,
            ls="-" if MK[j] == "o" else "--")
    # clip_on=False: the frame now ends ON the first and last tier, so a marker centred
    # there would lose half of itself to the axes box.
    ax.plot(xs, [t[2][j] for t in TIERS], MK[j], color=COLS[j], ms=4.5, mew=0,
            clip_on=False)

# COMPUTED, NOT SET. 0.10 was measured against the old type on the old axis; the ladder took
# the label to 7.974 pt and the tick rule took the axis to -1 to 1, and six labels at a fixed
# 0.10 overlapped again. One line box is fontsize/72 x 1.25 inches over a _YAX-inch panel,
# times the span in data units. The set is then put back on its own mean, so the spreading
# does not walk the whole group up off the top of the frame.
_GAP_C = (st.sz(6.3) / 72.0 * 1.25) / _YAX * 2.0

def _spread(vals, gap):
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    out = list(vals)
    for k in range(1, len(order)):
        a, b = order[k - 1], order[k]
        if out[b] - out[a] < gap:
            out[b] = out[a] + gap
    _sh = sum(out) / len(out) - sum(vals) / len(vals)
    return [o - _sh for o in out]

_ends = _spread([TIERS[-1][2][j] for j in range(len(hdr))], _GAP_C)
for j, h in enumerate(hdr):
    ax.plot([len(TIERS) - 1, len(TIERS) - 0.62], [TIERS[-1][2][j], _ends[j]],
            color=COLS[j], lw=0.5)
    ax.text(len(TIERS) - 0.55, _ends[j], LAB.get(h, h), fontsize=st.sz(6.3), color=COLS[j], va="center")
# The family's zero rule, as panel b now carries: light grey, 0.8 pt, dotted on (1.2, 2.0).
ax.axhline(0, color=st.LIGHTGREY, lw=0.8, ls=(0, (1.2, 2.0)))
# AND THE 60% RULE IS LIGHT GREY DOTTED TOO, not rust dashed. Rust marks the two falling
# features in this panel, so a rust rule made the colour mean two things inside one panel,
# the same fault Figure 5 panel b's model line had before it went grey. Its label is gone
# with it: caption (C) says what the line is.
ax.axvline(6.0, color=st.LIGHTGREY, lw=1.0, ls=(0, (1.2, 2.0)), zorder=1)
# The tier names ("0-10", "10-20", ...) do not fit side by side at this width. The lower edge
# of each tier says the same thing and reads horizontally.
# EVERY THIRD TIER, 0 to >= 90 by 30. Ten labels on a 590 px frame put "80" against ">=90";
# four of them read at a glance and the ticks stay on every tier.
ax.set_xticks(xs)
ax.set_xticklabels([t[0].split("-")[0].replace(">=", "≥") if i % 3 == 0 else ""
                    for i, t in enumerate(TIERS)], fontsize=st.sz(6.3))
ax.set_xlabel("occupancy of the lysine (%)", fontsize=st.sz(7))
ax.set_ylabel("shift from the\n0–10% tier (SD)", fontsize=st.sz(7), linespacing=1.4)
# ON A TICK AT BOTH ENDS, and the rail past the last tier is gone. The axis ran to 12.6
# against a last tier at 9, and the 2.6 units it held were for the six direct labels at 9.45;
# they are outside the frame now. The frame comes in by the same 9.4/13.0, so the ten tiers
# keep their pitch and every label lands on the page where it already was.
# FROM THE FIRST TIER TO THE LAST, frame to frame, as panel e now runs. It had 0.4 of a unit
# of lead-in on the left as well as the rail on the right. The frame comes in by the same
# 9.0/13.0, so the ten tiers keep their pitch and every direct label lands where it was.
ax.set_xlim(0, len(TIERS) - 1.0)
ax._axes_dw_in = (9.0 / 13.0 - 1.0) * _PS / 3.62
ax.set_ylim(-1.0, 1.0); ax.set_yticks([-1.0, -0.5, 0.0, 0.5, 1.0])

# ================================================================ d
ax = fig.add_subplot(_GS[1, 0]); st.head(fig, ax, "d")
ax._panel_head_flush = True
ax._axes_dw_in = -520.0 / 300.0    # the two-line tick labels no longer need the width
# THE SECOND AND THIRD COLUMNS ARE THE SAME SCREEN. Citing one of them and not the other read
# as though they came from different screens, so the contrast is named on the axis the three
# columns actually differ on, pooled across screens against within one screen against slopes
# inside that same screen. The caption says which screen, once, for both.
DCOL = [("pooled\n(n = %d)" % NPOOL, POOL),
        ("one screen\n(n = %d)" % NNAT, NAT),
        ("slopes\n(n = %d)" % NSLOP, SLOP)]
rng = np.random.default_rng(7)
for i, (lab, df) in enumerate(DCOL):
    sig = (df.q < 0.05).values
    jx = i + rng.normal(0, 0.075, len(df))
    # ONE ENCODING CHANNEL. Size carried significance alongside colour, 2.1 against 2.7, and
    # both sat under the house 4.5. Colour carries it alone now, at the house marker, and the
    # transparency is what lets 84 points share a column at that size.
    ax.plot(jx[~sig], df.rho.values[~sig], "o", ms=4.5, mew=0, alpha=0.55,
            color=st.LIGHTGREY)
    ax.plot(jx[sig], df.rho.values[sig], "o", ms=4.5, mew=0, alpha=0.55, color=st.ORANGE)
    # Above the frame, as Figure 4 panel e's p row is: the counts are a header for the column
    # below them, not a value on this scale.
    # THE HEADER HAS TO NAME ITS UNIT. A bare "21 of 84" over a column whose tick label reads
    # "82 compounds" puts two denominators in one column with nothing to tell them apart. A
    # point is one descriptor against one site preference, and 82 is the compounds inside each
    # of those tests, which is the vocabulary the caption already uses.
    # TWO LINES, because the header is what bounds the panel width once the tick labels are
    # short: "21 of 84 pairs" is 234 px on one line and 134 on two, and the column pitch has
    # to clear the wider of the two.
    ax.text(i, 0.63, "%d of %d\npairs" % (int(sig.sum()), len(df)), ha="center",
            va="bottom", linespacing=1.35,
            fontsize=st.sz(6.3), color=st.ORANGE if sig.sum() else st.GREY)
# A REFERENCE RULE IS FURNITURE. Solid slate at 0.8 made the zero line the darkest mark in
# the panel, competing with the points it exists to measure against. Light grey dotted and
# behind them, as in b and c.
ax.axhline(0, color=st.LIGHTGREY, lw=0.8, ls=(0, (1.2, 2.0)), zorder=1)
ax.set_xticks(range(3)); ax.set_xticklabels([c[0] for c in DCOL], fontsize=st.sz(6.3), linespacing=1.35)
ax.set_ylabel("property vs site\npreference (ρ)", fontsize=st.sz(7), linespacing=1.4)
# ON A TICK AT BOTH ENDS. It ran -0.60 to 0.72 to hold the counts that now sit above the
# frame. The correlations run -0.45 to +0.45, so -0.6 to 0.6 on a step of 0.3 holds them.
ax.set_ylim(-0.6, 0.6); ax.set_yticks([-0.6, -0.3, 0.0, 0.3, 0.6])
ax.set_xlim(-0.5, 2.5)
for s in ("top", "right"): ax.spines[s].set_visible(False)

# ================================================================ e
ax = fig.add_subplot(_GS[1, 1]); st.head(fig, ax, "e")
ax._panel_head_flush = True
ax._axes_dx_in = -208.0 / 300.0    # row 1 re-evened after d came in      # row 1 evened
# RUST IS NOT FREE IN THIS FIGURE. It already carries the 60% rule in a and c, the falling
# features in c, q < 0.05 in d, and the descriptors that produced false positives in f; a
# fifth meaning here, "not strictly monotonic", is one more than one colour can hold. Colour
# now says one thing only, which cut it is, and the dash says whether that cut is monotonic.
for cut, rows, mono in CUTS:
    c = st.BLUE if cut == 60 else st.GREY
    ls_ = "-" if mono else (0, (3.0, 1.6))
    lw_ = 1.3 if cut == 60 else 0.9
    # clip_on=False: the frame now ends ON the last tier, so a marker centred there would
    # lose its right half to the axes box.
    ax.plot(np.arange(len(rows)), [r[2] / RT for r in rows], marker="o", ls=ls_, color=c,
            lw=lw_, ms=5.0 if cut == 60 else 4.0, mew=0, clip_on=False,
            zorder=3 if cut == 60 else 2)
    ax.text(len(rows) - 0.92, rows[-1][2] / RT, "  %d%%" % cut, fontsize=st.sz(6.3), color=c,
            va="center", fontweight="bold" if cut == 60 else "normal")
base = CUTS[2][1]
# BEHIND THE DATA. Drawn after the five cuts it sat on top of them, so a 1.4 pt light grey
# line crossed every series at the left where they all converge on it. It is the null the
# others are read against, so it belongs underneath them.
ax.plot(np.arange(len(base)), [r[3] / RT for r in base], "-", color=st.LIGHTGREY, lw=1.4,
        zorder=1)
ax.text(len(base) - 0.92, base[-1][3] / RT, "  shuffled", fontsize=st.sz(6.3), color=st.GREY,
        va="center")
ax.set_xticks(np.arange(len(base)))
ax.set_xticklabels([r[0].replace(">=", "≥") for r in base], fontsize=st.sz(6.3))
ax.set_xlabel("compounds that engage the site", fontsize=st.sz(7))
ax.set_ylabel("shift from a\nnever-engaged site (SD)", fontsize=st.sz(7), linespacing=1.4)
# FROM THE FIRST TIER TO THE LAST, frame to frame. It ran -0.3 to 6.2 against tiers at 0 to
# 4: 0.3 of a unit of lead-in on the left and 2.2 units of rail on the right, the rail holding
# the six direct labels at 4.08, which are outside the frame now. The frame comes in by the
# same 4.0/6.5, so the five tiers keep their pitch and every label lands on the page where it
# already was.
ax.set_xlim(0, len(base) - 1.0)
ax._axes_dw_in = (4.0 / 6.5 - 1.0) * _PS / 3.62
ax.set_ylim(0, 1.0); ax.set_yticks([0, 0.25, 0.50, 0.75, 1.0])
for s in ("top", "right"): ax.spines[s].set_visible(False)

# ================================================================ f  (two tiers in one cell)
fsub = _GS[1, 2].subgridspec(2, 1, hspace=0.72, height_ratios=[1.0, 1.28])
ax = fig.add_subplot(fsub[0, 0]); st.head(fig, ax, "f")
ax._panel_head_flush = True
ax._axes_dx_in = 37.0 / 300.0      # row 1 evened
ORD = ["CyMPSA", "Meroterpenoids", "NatChem"]
COL3 = {"CyMPSA": st.GREEN, "Meroterpenoids": st.BLUE, "NatChem": st.ORANGE}   # as 5b
CITEF = {"CyMPSA": "Ryan 2024", "Meroterpenoids": "Bracken 2024", "NatChem": "Abbasov 2021"}
rg = np.random.default_rng(1)
for i, ds in enumerate(ORD):
    v = MW[MW.dataset == ds].MW.values
    # NO MEDIAN RULE. At the house marker a row of 200 compounds is a solid band, and a rule
    # drawn through it was a mark competing with the data it summarised. Transparent dots show
    # the pile-up instead, and the median is printed at the end of the row.
    ax.plot(v, i + rg.normal(0, 0.08, len(v)), "o", color=COL3[ds], ms=4.5, mew=0,
            alpha=0.55, zorder=2)
    # RIGHT OF THE ROW, not above its own tick. Three rows share about 156 px here, so a
    # 33 px number set 0.42 of a row above its tick lands in the row above: 386 printed on
    # the points of the screen above it. Right-aligned at the axis end it sits past every
    # point in every row, 610 Da being the furthest.
    ax.text(695, i, "%.0f" % np.median(v), ha="right", va="center",
            fontsize=st.sz(6.3), color=COL3[ds])
ax.set_yticks(range(3)); ax.set_yticklabels([CITEF[d_] for d_ in ORD], fontsize=st.sz(6.3))
ax.set_xlabel("molecular weight (Da)", fontsize=st.sz(7), labelpad=1.5)
# ROUND LIMITS ON A ROUND STEP. Weights run 150 to 610, so 100 to 700 every 200 holds them
# on four ticks; 100 to 640 answered with seven.
ax.set_xlim(100, 700); ax.set_xticks([100, 300, 500, 700]); ax.set_ylim(-0.75, 2.95)
ax.tick_params(labelsize=st.sz(6.3))
for s in ("top", "right"): ax.spines[s].set_visible(False)

ax = fig.add_subplot(fsub[1, 0])
ax._axes_dx_in = 37.0 / 300.0      # the other tier of f moves with it
C = CONF.sort_values("eta2", ascending=False)
for _, r in C.iterrows():
    ax.plot([r.eta2], [r.hits], "o", ms=4.5, mew=0, clip_on=False, alpha=0.55,
            color=st.ORANGE if r.hits else st.GREY)
# EACH LABEL ON ITS OWN DOT. fsp3 and rotatable sit centred above theirs. MW at 0.2545 and
# heavy at 0.2524 are 0.002 apart, about 6 px here against an 18 px marker, so one label
# cannot be assigned to one dot from where it sits; they are split either side of the pair,
# heavy atoms off the left dot and MW off the right one. TPSA is drawn with rho below, on the
# same baseline.
# Offsets are POINTS, not pixels: at this raster one point is 4.2 px, and the marker is
# 18.8 px across. 4 points clears the dot by 7 px, which reads as sitting on it; the 7 and 9
# the labels first carried were 29 and 38 px and read as floating.
LBL = {"fsp3":      ("fsp3",        0,  4, "center", "bottom"),
       "rotatable": ("rotatable",   0,  4, "center", "bottom"),
       "heavy":     ("heavy atoms", -10.5, -8.0, "right", "top"),
       "MW":        ("MW",            10.5, -8.0, "left",  "top")}
# A LEADER STARTS AT THE LABEL'S CENTRE UNLESS TOLD OTHERWISE, so the same offset gave a
# 172 px arrow on "heavy atoms" and a 110 px one on "MW": the longer the word, the longer its
# arrow. It is drawn on its own now, as an empty annotation anchored at the offset itself, so
# the leader is 55 px at 37 degrees whatever the word is and wherever the word is put.
LEADER = ("heavy", "MW")
# The word then closes on its own leader along that same 37 degrees. 1.5 points, 5 px, put
# the shaft start inside the word; 0.75 points, 2.5 px, leaves the two just clear. The step is
# applied to both labels identically, so they stay level.
# The horizontal step differs because the ink does, not because the labels differ in height:
# the final "s" of "heavy atoms" has a right side bearing, so its box edge stands off its ink,
# while the "M" of "MW" starts on its box edge. An equal step put one word clear of its shaft
# and the other on top of it. The VERTICAL step is identical, which is what holds the two
# labels level.
_CLOSE = {"heavy": (0.60, 0.46), "MW": (-0.40, 0.46)}
# The pair gets leaders, set below and angled. Horizontal ones 9 points out were a head and no
# shaft, so they read as > and < punctuation on the words rather than as arrows. From below,
# the left leader rises to heavy at 0.2524 and the right one to MW at 0.2545 without crossing,
# which is the only way 0.002 of separation can be assigned by eye. Below rather than above
# because the dots sit 2 units under the frame top and a label there would overrun it.
for _, r in C.iterrows():
    if r.descriptor not in LBL:
        continue
    _t, dx, dy, _ha, _va = LBL[r.descriptor]
    if r.descriptor in LEADER:
        # shrinkB CLEARS THE MARKER, it does not clear the point, and it is taken out of the
        # arrow rather than added to it. The 4.5 pt dot has a 2.25 pt radius, so 2 stopped the
        # head inside the ink and 4 stood it off by 11 px on a 55 px leader. 3.0 leaves 3.1 px.
        ax.annotate("", (r.eta2, r.hits), textcoords="offset points", xytext=(dx, dy),
                    ha="center", va="center",
                    arrowprops=dict(arrowstyle="->", color=st.ORANGE, lw=0.8, shrinkA=0,
                                    shrinkB=3.0, relpos=(0.5, 0.5)))
        _cx, _cy = _CLOSE[r.descriptor]
        dx = dx - _cx if dx > 0 else dx + _cx
        dy = dy + _cy
    ax.annotate(_t, (r.eta2, r.hits), textcoords="offset points", xytext=(dx, dy),
                fontsize=st.sz(6.3), color=st.ORANGE, fontweight="bold", ha=_ha, va=_va)
ax.set_xlabel("separates the three screens (η²)", fontsize=st.sz(7), labelpad=1.5)
ax.set_ylabel("false\npositives", fontsize=st.sz(7), linespacing=1.4)
# ON A TICK AT BOTH ENDS. It ran -1.6 to 9.6, so the bottom sat 1.6 below the lowest tick
# for clearance under the descriptors that produced no false positive at all. Those points sit
# ON zero now, drawn whole by clip_on=False, which is what the axis is for.
_XL = (0.05, 0.30)
ax.set_xlim(*_XL); ax.set_xticks([0.05, 0.10, 0.15, 0.20, 0.25, 0.30])
ax.set_ylim(0, 9); ax.set_yticks([0, 3, 6, 9])
ax.tick_params(labelsize=st.sz(6.3))
ax.text(0.97, 0.05, "ρ = %.2f" % RHO_CONF, transform=ax.transAxes, ha="right", va="bottom",
        fontsize=st.sz(6.3), color=st.GREY)
# TPSA KEEPS ITS LABEL. It is the highest ranked descriptor that produced no false positive,
# which is the counterexample the caption rests on, not an arbitrary naming. Anchored at rho's
# own baseline in axes fractions so the two sit on one line, and centred over its own dot.
_tpsa = float(C.loc[C.descriptor == "TPSA", "eta2"].iloc[0])
ax.text((_tpsa - _XL[0]) / (_XL[1] - _XL[0]), 0.05, "TPSA", transform=ax.transAxes,
        ha="center", va="bottom", fontsize=st.sz(6.3), color=st.GREY)
for s in ("top", "right"): ax.spines[s].set_visible(False)

st.save(fig, "fig6")
# The bar fit has run, so every value label can now be put on its own bar.
for _t, _q in _SNAP:
    _t.set_position((_q.get_x() + _q.get_width() / 2.0, _t.get_position()[1]))
# Trim to the drawing, 18 px all round, the figures-final standard Figure 1 set. Measured off
# the RENDERED PIXELS rather than get_tightbbox, because an axis-off panel still contributes
# its rectangle to the tight bbox and panel a has its axis off.
_PAD_IN = 18.0 / 300.0
from PIL import Image as _Im
_px = np.asarray(_Im.open("fig6.png").convert("RGB")).astype(int)
_ys, _xs = np.nonzero(_px.sum(2) < 735)
_H_in = fig.get_size_inches()[1]
_trim = Bbox.from_extents(_xs.min() / 300.0 - _PAD_IN,
                          _H_in - _ys.max() / 300.0 - _PAD_IN,
                          _xs.max() / 300.0 + _PAD_IN,
                          _H_in - _ys.min() / 300.0 + _PAD_IN)
for _ext in (".png", ".pdf"):
    fig.savefig("fig6" + _ext, bbox_inches=_trim)
print("trimmed to %.3f x %.3f in (%.0f x %.0f px at 300 dpi), 18 px pads"
      % (_trim.width, _trim.height, _trim.width * 300, _trim.height * 300))

if __import__("os").environ.get("FIG6_FINAL"):
    # THE FAMILY SHARES A RASTER, NOT A WIDTH. Figure 1 at 12.00 in fixes it at
    # 300 x 12.00 / 11.5400 = 311.96 dpi; every figure renders at that raster and is declared
    # 300 dpi, so all of them print 9.213 pt body type whatever their width.
    import PIL.Image as _PI, pymupdf as _mu
    _ref = _mu.open("fig1.pdf")[0].rect.width / 72.0
    _dpi = 300.0 * 12.0 / _ref
    fig.savefig("fig6_final.png", dpi=_dpi, bbox_inches=_trim)
    _im = _PI.open("fig6_final.png")
    _im.save("fig6_final.png", dpi=(300, 300))
    print("final: native %.4f in at the shared %.2f dpi raster -> %d x %d px declared 300 dpi "
          "(%.3f x %.3f in), body type %.3f pt against Figure 1's %.3f"
          % (_trim.width, _dpi, _im.size[0], _im.size[1], _im.size[0] / 300.0,
             _im.size[1] / 300.0, 8.860 * 12.0 / _ref, 8.860 * 12.0 / _ref))
fd.head("fig6", HEAD)
print("a: %s measurements, %s (%.2f%%) at >=60%%; median compound %s quantified, %d engaged"
      % (f"{NOBS:,}", f"{NHIT:,}", 100.0 * NHIT / NOBS, f"{int(CPD.n.median()):,}",
         int(CPD.lig.median())))
print("b: %d of %d consistent, n = %d compounds" % (NCONS, len(Q1), N37))
print("c: %s site-compound observations across %d tiers" % (f"{sum(t[1] for t in TIERS):,}", len(TIERS)))
print("d: pooled %d/%d (%d cpds) | NatChem %d/%d (%d) | slopes %d/%d (%d)"
      % ((POOL.q < 0.05).sum(), len(POOL), NPOOL, (NAT.q < 0.05).sum(), len(NAT), NNAT,
         (SLOP.q < 0.05).sum(), len(SLOP), NSLOP))
print("e: K = %d, axis / sqrt(K) = %.3f" % (KFEAT, RT))
print("d2: expected @0.05  pooled %.1f  NatChem %.1f  slopes %.1f" %
      (0.05 * len(POOL), 0.05 * len(NAT), 0.05 * len(SLOP)))
print("d3: max|rho|  pooled %.3f  NatChem %.3f  slopes %.3f" %
      (POOL.rho.abs().max(), NAT.rho.abs().max(), SLOP.rho.abs().max()))
print("b2: feature list (n=%d): %s" % (len(Q1), ", ".join(q[0] for q in Q1)))
print("e2: KFEAT features (n=%d): %s" % (KFEAT, ", ".join(SFLIST)))
print("f: MW medians %s | KW p %.2e | confound rho %.3f p %.4f"
      % ({k: round(v, 1) for k, v in MW.groupby("dataset").MW.median().items()}, KW[1],
         RHO_CONF, P_CONF))
