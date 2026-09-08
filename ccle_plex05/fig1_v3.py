#!/usr/bin/env python3
"""
fig1.py — Figure 1: cancer mutation creates a large, driver-enriched class of new lysines.

  a  union assembly: three variant sources -> 295,411 acquired-lysine sites
  b  substitution ranking, E->K first of 363 observed pairs
  c  gnomAD allele-frequency bands, with "not assessed" shown separately
  d  driver enrichment among Hess-assessed positions
  e  mutational-signature composition of driver hotspots
  f  acquired-lysine burden per tumour by cancer type, with the uveal control

SCOPE, which the legend must state: panels a-c are union-wide (DepMap + TCGA + COSMIC) and
panel f is TCGA; panels d and e are DepMap only, because the Hess driver and signature
annotations exist only there.

Panels d and e reproduce the logic of driver_figure.py exactly rather than re-deriving it — a
driver enrichment depends entirely on how the denominator is defined.
"""
import os, re, sys
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch
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

st.apply()

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

fig = plt.figure(figsize=(st.DOUBLE, st.DOUBLE * 1.05))
gs = fig.add_gridspec(3, 12, hspace=0.44, wspace=0.55,
                      left=0.055, right=0.985, top=0.955, bottom=0.055,
                      height_ratios=[1.38, 1.0, 1.0])

# ------------------------------------------------------------------ a
ax = fig.add_subplot(gs[0, 0:4]); st.panel_label(ax, "a", dx=-0.10, dy=1.05)
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
    v = venn3(subsets=reg, set_labels=("DepMap", "TCGA", "COSMIC"), ax=ax,
              set_colors=(st.LIGHTGREY, st.LIGHTGREY, st.LIGHTGREY), alpha=1.0)
    for key in ("100", "010", "001", "110", "101", "011", "111"):
        pt = v.get_patch_by_id(key)
        if pt is not None:
            pt.set_facecolor(st.LIGHTGREY); pt.set_alpha(1.0); pt.set_edgecolor("none")
    for circ, col in zip(venn3_circles(subsets=reg, ax=ax, lw=1.0),
                         (st.BLUE, st.GREEN, st.ORANGE)):
        circ.set_edgecolor(col)
    for t, col in zip(v.set_labels or [], (st.BLUE, st.GREEN, st.ORANGE)):
        if t: t.set_fontsize(6.0); t.set_color(col); t.set_fontweight("bold")
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
    for key in ("100", "010", "001", "110", "101", "011", "111"):
        t = v.get_label_by_id(key)
        if t is None: continue
        if key in KEEP:
            t.set_fontsize(5.4); t.set_color(st.GREY)
            t.set_text(f"{int(t.get_text()):,}")
            if key in OUT:
                x0, y0 = t.get_position()
                dx, dy, ha, ox, oy = OUT[key]
                t.set_position((dx, dy)); t.set_ha(ha)
                ax.annotate("", xy=(x0, y0), xytext=(dx + ox, dy + oy),
                            arrowprops=dict(arrowstyle="-", color=st.LIGHTGREY, lw=0.7,
                                            shrinkA=0, shrinkB=0), zorder=0)
        else:
            t.set_visible(False)
    # venn3 fits its layout to the axes, so the COSMIC set label lands hard against the
    # bottom edge and collided with the caption. Open space below the diagram instead of
    # moving the caption further away from the panel it belongs to.
    _y0, _y1 = ax.get_ylim()
    ax.set_ylim(_y0 - 0.20 * (_y1 - _y0), _y1 + 0.13 * (_y1 - _y0))
    ax.text(0.5, -0.02, "union  %s acquired lysines" % f"{len(k):,}",
            transform=ax.transAxes, ha="center", fontsize=6.8, fontweight="bold")
    ax.text(0.5, -0.15, "%s (%.0f%%) found in more than one source; small overlap\n"
            "regions are unlabelled \u2014 see Table S1" % (f"{multi:,}", 100*multi/len(k)),
            transform=ax.transAxes, ha="center", fontsize=5.4, color=st.GREY, linespacing=1.4)
ax.set_axis_off()

# ------------------------------------------------------------------ b
# The full codon route map, promoted from figS2 panel b. Spans 8 columns (~4.3 in), which is
# what its route labels need. The compact variant built for the old 1.73 in cell is retired.
ax = fig.add_subplot(gs[0, 4:12]); st.panel_label(ax, "b", dx=-0.048, dy=1.05)
ax.set_title("Every single-base route into and out of lysine", loc="left", pad=5, fontsize=7.0)
ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")

_IN = [("glutamate  E",  74.3, "C>T",          True),
       ("arginine  R",   11.2, "C>T",          True),
       ("glutamine  Q",   5.5, "C>A",          False),
       ("asparagine  N",  5.4, "C>A/G, T>A/G", False),
       ("threonine  T",   2.4, "C>A",          False),
       ("Met, Ile",       1.2, "T>A",          False)]
_OUT = [("N  asparagine", 45.5, "C>A/G, T>A/G"), ("T  threonine", 16.0, "T>G"),
        ("E  glutamate",  11.2, "T>C"),          ("R  arginine",  10.8, "T>C"),
        ("Q  glutamine",   6.3, "T>G"),          ("stop codon",    4.8, "T>A"),
        ("Met, Ile",       5.4, "T>A")]
KX0, KX1, KY0, KY1 = 4.30, 5.85, 4.55, 6.25
ax.add_patch(FancyBboxPatch((KX0, KY0), KX1-KX0, KY1-KY0,
                            boxstyle="round,pad=0.05,rounding_size=0.13",
                            facecolor="white", edgecolor=st.GREEN, lw=1.3))
ax.text((KX0+KX1)/2, KY0+(KY1-KY0)*0.65, "LYSINE", ha="center", va="center", fontsize=6.6,
        fontweight="bold", color=st.GREEN)
ax.text((KX0+KX1)/2, KY0+(KY1-KY0)*0.24, "AAA   AAG", ha="center", va="center", fontsize=5.4,
        color=st.INK)
_lw = lambda sh: 0.7 + 0.070 * sh
for (lab, sh, cls, ct), yy in zip(_IN, np.linspace(9.15, 2.85, len(_IN))):
    col = st.ORANGE if ct else st.LIGHTGREY
    ax.annotate("", xy=(KX0-0.05, (KY0+KY1)/2), xytext=(2.98, yy),
                arrowprops=dict(arrowstyle="-|>", color=col, lw=_lw(sh), shrinkA=1, shrinkB=1,
                                mutation_scale=7,
                                connectionstyle="arc3,rad=%.2f" % (-0.11 if yy > 5.4 else 0.11)))
    ax.text(2.86, yy, lab, ha="right", va="center", fontsize=5.4,
            fontweight="bold" if ct else "normal", color=st.INK if ct else st.GREY)
    ax.text(2.86, yy-0.44, "%.1f%%  %s" % (sh, cls), ha="right", va="center", fontsize=4.4,
            color=st.ORANGE if ct else st.GREY)
for (lab, sh, cls), yy in zip(_OUT, np.linspace(9.15, 2.55, len(_OUT))):
    ax.annotate("", xy=(7.10, yy), xytext=(KX1+0.05, (KY0+KY1)/2),
                arrowprops=dict(arrowstyle="-|>", color=st.BLUE, lw=_lw(sh), shrinkA=1,
                                shrinkB=1, mutation_scale=7,
                                connectionstyle="arc3,rad=%.2f" % (0.11 if yy > 5.4 else -0.11)))
    ax.text(7.22, yy, lab, ha="left", va="center", fontsize=5.4, color=st.GREY)
    ax.text(7.22, yy-0.44, "%.1f%%  %s" % (sh, cls), ha="left", va="center", fontsize=4.4,
            color=st.BLUE)
ax.text(2.86, 9.92, "what becomes a lysine", ha="right", va="center", fontsize=5.6,
        fontweight="bold", color=st.INK)
ax.text(7.22, 9.92, "what a lysine becomes", ha="left", va="center", fontsize=5.6,
        fontweight="bold", color=st.INK)
ax.text(2.86, 1.70, "14 codon routes\n4 of them C>T", ha="right", va="top", fontsize=5.2,
        color=st.ORANGE, fontweight="bold", linespacing=1.5)
ax.text(7.22, 1.45, "16 codon routes\nnone of them C>T", ha="left", va="top", fontsize=5.2,
        color=st.BLUE, fontweight="bold", linespacing=1.5)
ax.text(5.05, 3.85, "AAA carries no G or C.\nAAG's only G\u2192A gives AAA,\nwhich is still lysine.",
        ha="center", va="top", fontsize=5.0, color=st.INK, linespacing=1.55)
ax.text(5.05, 1.45, "orange = G:C\u2192A:T, half of all\ncancer missense mutation", ha="center",
        va="top", fontsize=4.8, color=st.GREY, linespacing=1.5)

# ------------------------------------------------------------------ c
ax = fig.add_subplot(gs[1, 0:4]); st.panel_label(ax, "c", dx=-0.16, dy=1.07)
if UNIONG:
    g = pd.read_csv(UNIONG, sep="\t", usecols=["band"])
    order = ["absent", "lt_1e-5", "1e-5_1e-4", "1e-4_1e-3", "gt_1e-3", "not_assessed"]
    show = ["absent", "<10$^{-5}$", "10$^{-5}$-\n10$^{-4}$", "10$^{-4}$-\n10$^{-3}$",
            ">10$^{-3}$", "not\nassessed"]
    v = np.array([(g.band == b).sum() for b in order])
    # A five-step lightness ramp was tried here and rejected on sight: the one chromatic step
    # in the middle of it pulled the eye to the 10^-5-10^-4 bin, which is not the point. The
    # point is the annotation, so colour carries the annotation: teal = the two rare bins that
    # make up the 91%, background grey = the three common bins, rust = never assessed.
    cols = [st.BLUE, st.BLUE, st.LIGHTGREY, st.LIGHTGREY, st.LIGHTGREY, st.ORANGE]
    ax.bar(range(len(v)), v / 1000, color=cols, width=0.76,
           edgecolor=st.GREY, linewidth=0.5)
    ax.set_xticks(range(len(v))); ax.set_xticklabels(show, fontsize=5.4)
    ax.set_ylabel("sites (thousands)")
    assessed = v.sum() - v[-1]; som = v[0] + v[1]
    HEAD["c_bands"] = dict(zip(order, [int(x) for x in v]))
    HEAD["c_pct_somatic_of_assessed"] = 100 * som / assessed
    HEAD["c_pct_not_assessed"] = 100 * v[5] / v.sum()
    ax.text(0.52, 0.97, "%.1f%% of assessed sites\nabsent or <10$^{-5}$"
            % (100 * som / assessed), transform=ax.transAxes, ha="center", va="top",
            fontsize=6, color=st.GREY, linespacing=1.4)
    ax.text(5, v[5]/1000 + v.max()/1000*0.05, "%.1f%%" % (100*v[5]/v.sum()),
            ha="center", fontsize=5.6, color=st.ORANGE, fontweight="bold")
else:
    notes.append("1c: union_sites_gnomad_v2.tsv not found")

# ------------------------------------------------------------------ d  (DepMap scope)
ax = fig.add_subplot(gs[1, 4:8]); st.panel_label(ax, "d", dx=-0.16, dy=1.07)
if PARQ:
    v = pd.read_parquet(PARQ)
    m = v[v.vclass == "missense"].dropna(subset=["alt_aa"]).copy()
    m["hess"] = m.hess_driver.astype(str).str.lower().isin(["true", "yes", "y", "1"])
    ass = m[m.hess_driver.notna()]
    K = ass[ass.alt_aa == "K"]; O = ass[ass.alt_aa != "K"]
    kr = 100 * K.hess.sum() / len(K); orr = 100 * O.hess.sum() / len(O)
    OR, P = stats.fisher_exact([[K.hess.sum(), len(K) - K.hess.sum()],
                                [O.hess.sum(), len(O) - O.hess.sum()]])
    ax.bar([0, 1], [orr, kr], color=[st.LIGHTGREY, st.BLUE], width=0.58)
    for i, val in enumerate([orr, kr]):
        ax.text(i, val * 1.03, "%.3f%%" % val, ha="center", fontsize=6.2, fontweight="bold")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["other\nsubstitutions", "acquired\nlysines"],
                                              fontsize=6.2)
    ax.set_ylabel("driver hotspots (%)"); ax.set_ylim(0, kr * 1.34)
    HEAD["d_pct_driver_K"] = kr; HEAD["d_pct_driver_other"] = orr
    HEAD["d_OR"] = OR; HEAD["d_p"] = P
    HEAD["d_n_K"] = len(K); HEAD["d_n_other"] = len(O)
    mant = P / (10 ** np.floor(np.log10(P))); expo = int(np.floor(np.log10(P)))
    ax.text(0.5, kr * 1.22, "OR %.2f, p = %.1f \u00d7 10$^{%d}$" % (OR, mant, expo),
            ha="center", fontsize=6, color="k")
else:
    notes.append("1d: depmap_variants_annotated.parquet not found")

# ------------------------------------------------------------------ e  (DepMap scope)
ax = fig.add_subplot(gs[1, 8:12]); st.panel_label(ax, "e", dx=-0.16, dy=1.07)
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
    x = np.arange(len(order)); w = 0.38
    ax.bar(x - w/2, pk, w, color=st.BLUE, label="acquired lysine")
    ax.bar(x + w/2, po, w, color=st.LIGHTGREY, label="other")
    # 30 degrees let POLE and APOBEC overlap by 31% of the smaller box (found by
    # figtools/overlap_check.py, not by eye). 42 degrees and a hair smaller clears it.
    ax.set_xticks(x); ax.set_xticklabels(order, fontsize=5.7, rotation=42, ha="right",
                                         rotation_mode="anchor")
    ax.set_ylabel("driver hotspots (%)")
    ax.set_ylim(0, max(max(pk), max(po)) * 1.30)
    ax.legend(frameon=False, fontsize=5.8, loc="upper left", handlelength=1.1,
              bbox_to_anchor=(0.0, 1.02))
    for i in range(len(order)):
        r = pk[i] / po[i] if po[i] else 0
        # label every ratio that departs from parity, in both directions. Depletion in
        # smoking and CpG is part of the argument — those processes do not generate the
        # transitions that create lysine — so labelling only enrichments would be selective.
        if r and (r >= 1.25 or r <= 0.8):
            ax.text(i, max(pk[i], po[i]) * 1.04, "%.1f\u00d7" % r, ha="center",
                    fontsize=5.6, fontweight="bold",
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
ax = fig.add_subplot(gs[2, 0:9]); st.panel_label(ax, "f", dx=-0.045, dy=1.06)
if BURD:
    T = pd.read_csv(BURD).sort_values("median", ascending=False)
    HEAD["f_n_lineages"] = len(T)
    HEAD["f_SKCM"] = float(T.loc[T.cancer == "SKCM", "median"].iloc[0]) if (T.cancer == "SKCM").any() else None
    HEAD["f_UVM"] = float(T.loc[T.cancer == "UVM", "median"].iloc[0]) if (T.cancer == "UVM").any() else None
    HEAD["f_SKCM_fold"] = float(T.loc[T.cancer == "SKCM", "fold_vs_pan"].iloc[0]) if "fold_vs_pan" in T.columns and (T.cancer == "SKCM").any() else None
    x = np.arange(len(T))
    cols = [st.ORANGE if c == "SKCM" else (st.GREEN if c == "UVM" else st.BLUE)
            for c in T.cancer]
    ax.bar(x, T["median"], color=cols, width=0.72)
    ax.set_xticks(x); ax.set_xticklabels(T.cancer, fontsize=4.4, rotation=90)
    ax.set_ylabel("median acquired lysines\nper tumour", fontsize=6.0)
    ax.tick_params(axis="y", labelsize=5.6)
    ax.set_ylim(0, T["median"].max() * 1.22)
    for c, lab, xo, yf, col in (("SKCM", "melanoma, 13\u00d7 pan-TCGA", 2.4, 0.80, st.ORANGE),
                                ("UVM", "uveal melanoma \u2014 not UV-driven,\nmedian 0",
                                 -7.5, 0.36, st.GREEN)):
        if c in set(T.cancer):
            i = list(T.cancer).index(c)
            ax.annotate(lab, xy=(i, T["median"].iloc[i]),
                        xytext=(i + xo, T["median"].max() * yf), fontsize=5.2, color=col,
                        arrowprops=dict(arrowstyle="-", color=col, lw=0.6), linespacing=1.45)
    ax.text(0.99, 0.97, "%d lineages, %s tumours" % (len(T), f"{int(T.n.sum()):,}"),
            transform=ax.transAxes, ha="right", va="top", fontsize=5.2, color=st.GREY)
else:
    notes.append("1f: tcga_burden_by_lineage.csv not found")

# ------------------------------------------------------------------ g
# Promoted from figS1 panel d. This is the control that excludes the cheapest alternative to
# the whole premise — that acquired lysines are simply what mutation does everywhere, and the
# same excess would appear in healthy people. It does not: somatic mutation creates them at
# twice the germline rate. Nothing else in this figure addresses that.
ax = fig.add_subplot(gs[2, 9:12]); st.panel_label(ax, "g", dx=-0.30, dy=1.06)
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
    notes.append("1g: germline_lysine.out %s — germline bar cannot be drawn" % _gstat)
    germ_pct, germ_n = float("nan"), float("nan")

if "a_union_K" in HEAD:
    som_pct = 100 * HEAD["a_union_K"] / HEAD["a_all_missense"]
else:
    som_pct = float("nan")
    notes.append("1g: union counts unavailable — somatic bar cannot be computed")
HEAD["g_som_pct"] = som_pct
HEAD["g_germ_pct"] = germ_pct
HEAD["g_germ_per_individual"] = germ_n
HEAD["g_germ_nc_pct"] = nc_pct
HEAD["g_germ_nc_per_individual"] = nc_n
HEAD["g_fold"] = som_pct / germ_pct if germ_pct else None
ax.bar([0, 1], [germ_pct, som_pct], 0.5, color=[st.LIGHTGREY, st.BLUE])
for i, v in enumerate([germ_pct, som_pct]):
    ax.text(i, v + 0.12, "%.2f%%" % v, ha="center", fontsize=6, fontweight="bold")
ax.set_xticks([0, 1])
ax.set_xticklabels(["germline\n(gnomAD v4.1)", "somatic\n(union)"], fontsize=5.6)
ax.set_ylabel("acquired lysines as % of\nall missense", fontsize=6.0)
ax.tick_params(axis="y", labelsize=5.6)
ax.set_ylim(0, som_pct * 1.35)
ax.text(0.5, som_pct * 1.20, "%.2f-fold" % (som_pct / germ_pct), ha="center",
        fontsize=6.4, color=st.BLUE, fontweight="bold")
ax.text(0, -0.30, "%.0f per individual" % germ_n, ha="center",
        transform=ax.get_xaxis_transform(), va="top", fontsize=5.0, color=st.GREY)

st.save(fig, "fig1")
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
