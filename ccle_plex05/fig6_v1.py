#!/usr/bin/env python3
"""
Figure 6 — Engagement is a property of the lysine, not of the molecule.

DESIGN RULE FOR THIS FIGURE (set 2026-09-02, after review)
No sentences inside panels. A panel carries axis labels, tick labels, series labels and numbers
attached to marks — nothing a reader has to stop and read. Everything explanatory belongs in the
legend. The previous version had a three- to five-line grey caption in five of the six panels;
they are gone, and FIGURE_LEGENDS_fig6.md carries what they said.

PANEL a — the assay, then the design, as one flow
  Two different things were being labelled with two headings and drawn the same way, so the panel
  read as one idea stated twice. It is now left to right:
      competition  ->  occupancy distribution with the 60% cut  ->  within one compound  ->  pooled
  WHAT IS MEASURED is occupancy: one number per (lysine, compound), 100 x (1 - compound/DMSO),
  the share of probe labelling the compound competes away. THE DESIGN is what is done with those
  numbers: a compound's engaged sites are compared with ITS OWN below-cut sites, standardised
  inside that compound, and only then pooled across compounds.
  The distribution is every measurement in all three screens. An earlier draft drew one compound
  (27d), picked by sorting on engaged-site count — selection on the displayed quantity. 27d is the
  99th percentile for engagement rate, so it implied a hit rate ~20x the truth.

PANEL d — the null as a distribution, not a count
  Three bars, one of them zero, asked the reader to read an absence off blank space. Every
  correlation is now plotted: the cloud contracts from |rho| 0.45 pooled to 0.24 within NatChem to
  0.18 across the 173 slopes, and the significant points stop existing. The within-NatChem values
  are reproduced by step37b_save_within.py, which step 37 computes and never saved.

PANEL e — a real axis
  The plotted quantity is the Euclidean length of a K-component vector of SD-unit mean
  differences. It is not itself an SD; with K = 15 a uniform 0.1 SD shift reads 0.39. Divided by
  sqrt(K) it is the average shift per feature, which is judgeable. K is read from step 40's own
  feature list so a changed feature set rescales the axis instead of silently misreporting it.

PANEL f — the confound AND its consequence
  Left: molecular weight all but names the screen. Right: the descriptors that separate the
  screens are the descriptors that produced pooled positives (rho 0.69, p 0.027, step 41).
  This replaces "every one of the 21 pooled hits involves a SIZE descriptor", which was wrong
  twice — it called rotatable-bond count and fsp3 "size", and it described step 31 (150 pairs)
  while Fig. 6d plots step 37 (84 pairs). Both return 21 by coincidence; step 37's 21 include six
  from logP.
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
st.apply()

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
                   .replace("nz sasa", "NZ exposure").replace("sasa sidechain", "side-chain exp.")
                   .replace("dist centroid", "dist. to centroid").replace("n lys", "nearby lysines"))

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
fig = plt.figure(figsize=(st.DOUBLE, st.DOUBLE * 0.62))
gtop = fig.add_gridspec(1, 3, left=0.070, right=0.980, top=0.945, bottom=0.585,
                        wspace=0.46, width_ratios=[1.62, 1.0, 1.0])
gbot = fig.add_gridspec(1, 3, left=0.070, right=0.980, top=0.455, bottom=0.085, wspace=0.46)

# ================================================================ a
# THE VISUAL SYSTEM. Three stages on one band: every stage sits between y = 2.4 and y = 8.4, is
# introduced by a stage label in the same style at y = 9.4, and is joined to the next by the same
# arrow at the same height. Nothing is placed by eye. Panel-a units are ~5.4 x per y, so widths
# below are in x units and heights in y units and they are NOT visually proportional.
ax = fig.add_subplot(gtop[0, 0]); st.panel_label(ax, "a", dx=-0.10, dy=1.09)
ax.set_xlim(0, 100); ax.set_ylim(0, 10); ax.set_axis_off()
BAND_LO, BAND_HI, LABEL_Y, ARROW_Y = 2.4, 8.4, 9.35, 5.4
STAGE = [(0, 24, "compete"), (30, 62, "split at 60%"), (68, 100, "compare within")]
for x0, x1, lab in STAGE:
    ax.text((x0 + x1) / 2, LABEL_Y, lab, fontsize=5.0, color=st.GREY, fontweight="bold",
            ha="center", va="center")
for (x0, x1, _), (nx0, _, _) in zip(STAGE[:-1], STAGE[1:]):
    ax.annotate("", xy=(nx0 - 1.5, ARROW_Y), xytext=(x1 + 1.5, ARROW_Y),
                arrowprops=dict(arrowstyle="-|>", color=st.LIGHTGREY, lw=1.1))

# stage 1 — probe signal with and without compound, on a shared baseline
BASE_Y = 4.2
for k, (lab, h, c) in enumerate([("DMSO", 3.6, st.GREY), ("+ cmpd", 1.25, st.BLUE)]):
    ax.add_patch(Rectangle((1.5 + k * 10.0, BASE_Y), 8.0, h, facecolor=c, lw=0))
    ax.text(5.5 + k * 10.0, BASE_Y - 0.30, lab, fontsize=4.6, color=c, ha="center", va="top")
ax.plot([0.5, 20.5], [BASE_Y, BASE_Y], color=st.GREY, lw=0.7)
# The bracket is the whole point of the stage; the word for it is on the formula line below,
# so labelling the arrow as well would be saying it twice in a panel with no room to spare.
ax.annotate("", xy=(22.0, BASE_Y + 1.25), xytext=(22.0, BASE_Y + 3.6),
            arrowprops=dict(arrowstyle="<->", color=st.ORANGE, lw=1.0, shrinkA=0, shrinkB=0))

# stage 2 — every measurement, and where the cut falls on it
axi = ax.inset_axes([0.315, 0.255, 0.290, 0.545], transform=ax.transAxes)
bins = np.arange(0, 102, 4)
DSC = {"NatChem": (st.BLUE, "-"), "CyMPSA": (st.ORANGE, (0, (3.0, 1.4))),
       "Meroterpenoid": (st.GREEN, (0, (1.2, 1.2)))}
for _ds, _v in OCCD.items():
    _c, _ls = DSC[_ds]
    _h, _e = np.histogram(_v, bins=bins)
    axi.step(_e[:-1], _h / _h.sum(), where="post", color=_c, lw=1.0, ls=_ls, label=_ds)
axi.axvline(60, color=st.ORANGE, lw=1.0, ls=(0, (3, 1.6)))
axi.set_yscale("log")
axi.set_xlabel("occupancy (%)", fontsize=5.2, labelpad=1.2)
axi.set_ylabel("share", fontsize=5.2, labelpad=1.2)
axi.tick_params(labelsize=4.3, pad=1.0)
axi.set_xticks([0, 60, 100]); axi.set_yticks([1e0, 1e-2, 1e-4])
axi.legend(frameon=False, fontsize=4.2, loc="upper right", handlelength=1.4, borderaxespad=0.1,
           labelspacing=0.22)
for _s in ("top", "right"): axi.spines[_s].set_visible(False)
ax.text(46.0, 0.35, "occupancy = 100 × (1 − cmpd/DMSO)", fontsize=4.6, color=st.GREY,
        ha="center", va="bottom")

# stage 3 — the two groups a compound supplies, and the shift between them
SW, SH, GAP = 1.55, 0.42, 0.62
for i_, (n_, col, y_) in enumerate([(5, st.BLUE, 7.5), (13, st.LIGHTGREY, 6.0)]):
    for k in range(n_):
        ax.add_patch(Rectangle((68.5 + k * (SW + GAP), y_), SW, SH, facecolor=col, lw=0))
ax.text(68.5, 7.25, "engaged, ≥ 60%", fontsize=4.6, color=st.BLUE, va="top", fontweight="bold")
ax.text(68.5, 5.75, "quantified, below 60%", fontsize=4.6, color=st.GREY, va="top")
ax.annotate("", xy=(72.0, 4.35), xytext=(72.0, 5.15),
            arrowprops=dict(arrowstyle="-|>", color=st.GREY, lw=0.9))
ax.text(73.6, 4.72, "standardise inside\nthis compound", fontsize=4.6, color=st.GREY,
        va="center", linespacing=1.45)
rng = np.random.default_rng(4)
for k in range(20):
    ax.plot([69.0 + k * 1.55], [2.95 + rng.normal(0, 0.36)], "o", ms=2.2, color=st.BLUE, mew=0)
ax.plot([68.5, 100.0], [2.95, 2.95], color=st.GREY, lw=0.7)
ax.text(68.5, 2.05, "one point per compound  →  b", fontsize=4.7, color=st.GREY, va="top")

# ================================================================ b
ax = fig.add_subplot(gtop[0, 1]); st.panel_label(ax, "b", dx=-0.40, dy=1.09)
Q = sorted(Q1, key=lambda q: q[1])
for i, (f_, sh, sd, n_, q_, sig) in enumerate(Q):
    se = sd / np.sqrt(n_)
    c = st.BLUE if sig else st.ORANGE
    ax.plot([sh - 1.96 * se, sh + 1.96 * se], [i, i], color=c, lw=1.3, solid_capstyle="round")
    ax.plot([sh], [i], "o", color=c, ms=3.4, mew=0)
ax.axvline(0, color=st.GREY, lw=0.8)
ax.set_yticks(np.arange(len(Q))); ax.set_yticklabels([NICE(q[0]) for q in Q], fontsize=5.0)
ax.set_xlabel("engaged vs same compound's\nbelow-cut sites (SD)", fontsize=6.0, linespacing=1.4)
ax.set_xlim(-0.55, 0.48); ax.set_ylim(-0.7, len(Q) + 0.3)
# Top left is empty: the highest rows carry the largest positive shifts, which sit right.
ax.text(0.02, 0.98, "%d of %d\nn = %d compounds" % (NCONS, len(Q1), N37), transform=ax.transAxes,
        fontsize=5.4, color=st.BLUE, fontweight="bold", va="top", linespacing=1.4)
for s in ("top", "right"): ax.spines[s].set_visible(False)

# ================================================================ c
ax = fig.add_subplot(gtop[0, 2]); st.panel_label(ax, "c", dx=-0.28, dy=1.09)
LAB = {"rsa_tien": "exposure", "nz_sasa": "NZ exposure", "crowding_": "crowding 8Å",
       "atoms_10A": "atoms 10Å", "hydrophob": "hydrophobic 8Å", "aromatic_": "aromatic 8Å"}
COLS = [st.ORANGE, st.ORANGE, st.BLUE, st.BLUE, st.GREEN, st.GREEN]
MK = ["o", "s", "o", "s", "o", "s"]
xs = np.arange(len(TIERS))
for j, h in enumerate(hdr):
    ax.plot(xs, [t[2][j] for t in TIERS], "-", color=COLS[j], lw=1.0,
            ls="-" if MK[j] == "o" else "--")
    ax.plot(xs, [t[2][j] for t in TIERS], MK[j], color=COLS[j], ms=2.4, mew=0)

def _spread(vals, gap=0.10):
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    out = list(vals)
    for k in range(1, len(order)):
        a, b = order[k - 1], order[k]
        if out[b] - out[a] < gap:
            out[b] = out[a] + gap
    return out

_ends = _spread([TIERS[-1][2][j] for j in range(len(hdr))])
for j, h in enumerate(hdr):
    ax.plot([len(TIERS) - 1, len(TIERS) - 0.62], [TIERS[-1][2][j], _ends[j]],
            color=COLS[j], lw=0.5)
    ax.text(len(TIERS) - 0.55, _ends[j], LAB.get(h, h), fontsize=4.6, color=COLS[j], va="center")
ax.axhline(0, color=st.LIGHTGREY, lw=0.8)
ax.axvline(5.5, color=st.ORANGE, lw=1.0, ls=(0, (3, 1.6)))
ax.text(5.35, -0.94, "60%", ha="right", fontsize=5.0, color=st.ORANGE, fontweight="bold")
# The tier names ("0-10", "10-20", ...) do not fit side by side at this width. The lower edge
# of each tier says the same thing and reads horizontally.
ax.set_xticks(xs)
ax.set_xticklabels([t[0].split("-")[0].replace(">=", "≥") for t in TIERS], fontsize=5.0)
ax.set_xlabel("occupancy of the lysine (%)", fontsize=6.0)
ax.set_ylabel("shift from the 0–10% tier (SD)", fontsize=6.0)
ax.set_xlim(-0.4, len(TIERS) + 2.6); ax.set_ylim(-1.05, 0.95)

# ================================================================ d
ax = fig.add_subplot(gbot[0, 0]); st.panel_label(ax, "d", dx=-0.24, dy=1.09)
DCOL = [("pooled\n%d compounds" % NPOOL, POOL), ("NatChem only\n%d compounds" % NNAT, NAT),
        ("NatChem slopes\n%d compounds" % NSLOP, SLOP)]
rng = np.random.default_rng(7)
for i, (lab, df) in enumerate(DCOL):
    sig = (df.q < 0.05).values
    jx = i + rng.normal(0, 0.075, len(df))
    ax.plot(jx[~sig], df.rho.values[~sig], "o", ms=2.1, mew=0, color=st.LIGHTGREY)
    ax.plot(jx[sig], df.rho.values[sig], "o", ms=2.7, mew=0, color=st.ORANGE)
    ax.text(i, 0.615, "%d of %d" % (int(sig.sum()), len(df)), ha="center", fontsize=5.6,
            fontweight="bold", color=st.ORANGE if sig.sum() else st.BLUE)
ax.axhline(0, color=st.GREY, lw=0.8)
ax.set_xticks(range(3)); ax.set_xticklabels([c[0] for c in DCOL], fontsize=5.0, linespacing=1.35)
ax.set_ylabel("compound property vs site\npreference (Spearman ρ)", fontsize=6.0, linespacing=1.4)
ax.set_ylim(-0.60, 0.72); ax.set_xlim(-0.5, 2.5)
for s in ("top", "right"): ax.spines[s].set_visible(False)

# ================================================================ e
ax = fig.add_subplot(gbot[0, 1]); st.panel_label(ax, "e", dx=-0.26, dy=1.09)
for cut, rows, mono in CUTS:
    c = st.ORANGE if not mono else (st.BLUE if cut == 60 else st.GREY)
    lw_ = 1.3 if cut == 60 else (1.1 if not mono else 0.8)
    ax.plot(np.arange(len(rows)), [r[2] / RT for r in rows], "-o", color=c, lw=lw_,
            ms=2.8 if cut == 60 else 2.2, mew=0)
    ax.text(len(rows) - 0.92, rows[-1][2] / RT, "  %d%%" % cut, fontsize=4.9, color=c,
            va="center", fontweight="bold" if (cut == 60 or not mono) else "normal")
base = CUTS[2][1]
ax.plot(np.arange(len(base)), [r[3] / RT for r in base], "-", color=st.LIGHTGREY, lw=1.4)
ax.text(len(base) - 0.92, base[-1][3] / RT, "  shuffled", fontsize=4.9, color=st.GREY,
        va="center")
ax.set_xticks(np.arange(len(base)))
ax.set_xticklabels([r[0].replace(">=", "≥") for r in base], fontsize=5.2)
ax.set_xlabel("compounds that engage the site", fontsize=6.0)
ax.set_ylabel("average shift per site feature,\nvs a never-engaged site (SD)", fontsize=6.0,
              linespacing=1.4)
ax.set_xlim(-0.3, len(base) + 1.2)
for s in ("top", "right"): ax.spines[s].set_visible(False)

# ================================================================ f  (two tiers in one cell)
fsub = gbot[0, 2].subgridspec(2, 1, hspace=0.72, height_ratios=[1.0, 1.28])
ax = fig.add_subplot(fsub[0, 0]); st.panel_label(ax, "f", dx=-0.26, dy=1.20)
ORD = ["CyMPSA", "Meroterpenoids", "NatChem"]
COL3 = {"CyMPSA": st.ORANGE, "Meroterpenoids": st.GREEN, "NatChem": st.BLUE}
rg = np.random.default_rng(1)
for i, ds in enumerate(ORD):
    v = MW[MW.dataset == ds].MW.values
    ax.plot(v, i + rg.normal(0, 0.10, len(v)), "o", color=COL3[ds], ms=1.7, mew=0)
    ax.plot([np.median(v)] * 2, [i - 0.30, i + 0.30], color=COL3[ds], lw=1.7,
            solid_capstyle="round")
    ax.text(np.median(v), i + 0.42, "%.0f" % np.median(v), ha="center", fontsize=5.2,
            fontweight="bold", color=COL3[ds])
ax.set_yticks(range(3)); ax.set_yticklabels(ORD, fontsize=5.0)
ax.set_xlabel("molecular weight (Da)", fontsize=5.8, labelpad=1.5)
ax.set_xlim(100, 640); ax.set_ylim(-0.75, 2.95)
ax.tick_params(labelsize=4.8)
for s in ("top", "right"): ax.spines[s].set_visible(False)

ax = fig.add_subplot(fsub[1, 0])
C = CONF.sort_values("eta2", ascending=False)
for _, r in C.iterrows():
    ax.plot([r.eta2], [r.hits], "o", ms=4.0 if r.hits else 2.8, mew=0,
            color=st.ORANGE if r.hits else st.GREY)
for _, r in C.iterrows():
    if r.hits > 0 or r.descriptor == "TPSA":
        dx, dy = {"MW": (5, 3), "heavy": (5, -9), "rotatable": (-5, 4), "fsp3": (5, 2),
                  "TPSA": (5, 1)}[r.descriptor]
        ax.annotate(r.descriptor, (r.eta2, r.hits), textcoords="offset points", xytext=(dx, dy),
                    fontsize=4.8, color=st.ORANGE if r.hits else st.GREY,
                    fontweight="bold" if r.hits else "normal",
                    ha="right" if dx < 0 else "left")
ax.set_xlabel("separates the three screens (η²)", fontsize=5.8, labelpad=1.5)
ax.set_ylabel("pooled false\npositives (of 15)", fontsize=5.8, linespacing=1.4)
ax.set_xlim(0.035, 0.30); ax.set_ylim(-1.6, 9.6)
ax.tick_params(labelsize=4.8)
ax.text(0.97, 0.05, "ρ = %.2f" % RHO_CONF, transform=ax.transAxes, ha="right", va="bottom",
        fontsize=5.4, color=st.GREY, fontweight="bold")
for s in ("top", "right"): ax.spines[s].set_visible(False)

st.save(fig, "fig6")
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
