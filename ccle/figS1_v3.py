#!/usr/bin/env python3
"""
Supplementary figure for Figure 1 — the acquired-lysine class, supporting analyses.

a  the full substitution matrix, all 20x20 pairs, of which 363 are observed and 17 never occur
b  observed against trinucleotide-weighted expectation: E->K leads because of mutational
   context, not in spite of it
c  gnomAD population frequency, with 'not assessed' kept as its own category
d  germline against somatic, the control for whether acquired lysines are simply what mutation
   does everywhere
e  burden per tumour across every TCGA lineage, with UVM as the causal control

NOTE ON A CORRECTED NUMBER
The germline figure is 471.1 acquired lysines per individual (3.29% of missense), read from
germline_lysine.out. An earlier summary carried 572 and 3.23%; the file is the authority.
"""
import os, sys, re
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt
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

def find(n):
    # ligand_rebuild added 2026-09-04: the two panels consolidated in from figS2 read their
    # step48 outputs from there.
    for d in (".", "..", os.path.join(_TOP, 'ccle'),
              os.path.join(_TOP, 'ligand_rebuild')):
        p = os.path.join(d, n)
        if os.path.exists(p): return p
    return None

AA = list("ACDEFGHIKLMNPQRSTVWY")
# CONSOLIDATED 2026-09-04. The CODON SUPPLEMENT built on 2026-09-03 no longer exists; its
# script is in ligand_rebuild/superseded/figS2_codon_routes.py.
#
# NAME WARNING: that retired script was called figS2_codon_routes.py and this comment used to
# say "figS2 is gone", which is false and misleading. The manuscript's figS2 is figS2_v3.py,
# the supplement to FIGURE 2 (non-lysine control, peptide length, missed cleavage, score margin,
# spectrum quality, annotated spectrum). It is untouched and still built. Nothing about Figure 2
# was consolidated. Only the codon supplement was.
#
# What happened to its three panels: the chord was DROPPED — panel e here already carries the
# net-flux-by-residue message on a more accurate encoding (bar length beats ribbon width) and
# across every substitution class rather than C>T alone, and the routing the chord uniquely
# showed is now delivered for lysine, better, by Figure 1b. The other two moved here as f and g.
# Every support panel for Figure 1 is therefore in this one figure, the SI count returns to
# figS1-figS8, and no renumbering is needed.
CLASSES = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
SP = pd.read_csv(find("step48_spectra.csv")).set_index("set")
CT = pd.read_csv(find("step48c_by_cancer_type.csv")).sort_values("bg_CT")
fig = plt.figure(figsize=(st.DOUBLE, st.DOUBLE * 1.30))
gs = fig.add_gridspec(4, 12, hspace=0.66, wspace=2.4,
                      left=0.06, right=0.985, top=0.955, bottom=0.055)

# HEADLINE NUMBERS, added 2026-09-04, for the same reason as fig1: this printed only its panel
# count, so its legend would have been written from prose rather than from the build.
HEAD = {}

# ---------------------------------------------------------------- a: substitution matrix
ax = fig.add_subplot(gs[0, 0:4]); st.panel_label(ax, "a", dx=-0.20, dy=1.07)
# Observed over expected rather than raw counts. Counts show only that some substitutions are
# commoner than others, which nobody doubts; this shows which are commoner THAN THE MUTATIONAL
# PROCESS PREDICTS. Expected is computed over the 19,212 MANE genes whose coding sequence could
# be reconstructed; observed spans the full union, so the two sides are not quite the same gene
# set — 99.8% overlap by gene count.
# Normalised over the WHOLE matrix, not per row: each cell is that substitution's share of all
# 4.4 M missense sites, so cells are comparable across rows and E->K's rank is visible. Row
# normalisation showed each residue's own profile but made every row's darkest cell look alike.
#
# Direct counts from the union. An earlier version showed observed/expected, which required
# deriving a mutational spectrum from the same substitutions being tested; that circularity made
# individual cells unreliable (arginine's enrichment survived independent checks, tyrosine's did
# not), and the enrichment framing answered a reviewer objection rather than the question the
# panel is for.
M = pd.read_csv(find("union_subst_matrix.csv"))
grid = np.full((20, 20), np.nan)
for r in M.itertuples():
    if r.ref_aa in AA and r.alt_aa in AA:
        grid[AA.index(r.ref_aa), AA.index(r.alt_aa)] = r.n_sites
pct = 100 * grid / np.nansum(grid)
# Light to dark, anchored only on the supplied five plus the page: paper -> light grey ->
# teal -> charcoal. magma_r belonged to no other figure and its purple/black high end sat
# outside the palette entirely; the tint and shade of the teal that replaced it were also
# outside it, just less obviously.
_CMAP = st.seq_cmap()
im = ax.imshow(pct, cmap=_CMAP, aspect="equal", vmin=0, vmax=np.nanmax(pct))
ax.set_xticks(range(20)); ax.set_xticklabels(AA, fontsize=4.0)
ax.set_yticks(range(20)); ax.set_yticklabels(AA, fontsize=4.0)
ax.set_xlabel("becomes", fontsize=6); ax.set_ylabel("from", fontsize=6)
ei, ki = AA.index("E"), AA.index("K")
HEAD["a_EK_sites"] = int(grid[ei, ki]); HEAD["a_total_subs"] = int(np.nansum(grid))
HEAD["a_EK_share_pct"] = round(100.0 * grid[ei, ki] / np.nansum(grid), 2)
HEAD["a_pairs_observed"] = int(np.nansum(grid > 0))
ax.add_patch(plt.Rectangle((ki - 0.5, ei - 0.5), 1, 1, fill=False, edgecolor=st.ORANGE, lw=1.4))
flat = pct[~np.isnan(pct)]
rank = int((flat > pct[ei, ki]).sum()) + 1
ax.text(ki + 0.9, ei - 0.1, "E\u2192K\n%.2f%%, rank %d" % (pct[ei, ki], rank),
        fontsize=4.8, color=st.ORANGE, va="center", linespacing=1.3)
cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
cb.set_label("% of all missense sites", fontsize=5.2)
cb.ax.tick_params(labelsize=4.6)
ax.text(0.5, -0.34, "whole matrix sums to 100%; white = never observed",
        transform=ax.transAxes, ha="center", fontsize=4.8, color=st.GREY)

# ---------------------------------------------------------------- b: codon context
ax = fig.add_subplot(gs[0, 5:8]); st.panel_label(ax, "b", dx=-0.28, dy=1.07)
C = find("codon_context96.csv")
if C:
    D = pd.read_csv(C).nlargest(12, "observed")
    y = np.arange(len(D))[::-1]
    cols = [st.ORANGE if (r.ref == "E" and r.alt == "K") else
            (st.BLUE if r.ratio > 2 else st.LIGHTGREY) for r in D.itertuples()]
    ax.barh(y, D.ratio, height=0.62, color=cols)
    ax.axvline(1, color=st.GREY, lw=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels(["%s\u2192%s" % (r.ref, r.alt) for r in D.itertuples()], fontsize=5.0)
    ax.set_xlabel("observed / expected", fontsize=5.8)
    ax.set_xlim(0, 4.0)
    HEAD["b_top12"] = {"%s->%s" % (r.ref, r.alt): round(float(r.ratio), 3) for r in D.itertuples()}
    ek = D[(D.ref == "E") & (D.alt == "K")]
    if len(ek):
        HEAD["b_EK_obs_over_exp"] = round(float(ek.ratio.iloc[0]), 3)
    if len(ek):
        ax.text(1.25, y[list(D.index).index(ek.index[0])], "E\u2192K = %.2f" % ek.ratio.iloc[0],
                fontsize=5.2, color=st.ORANGE, va="center", fontweight="bold")
    # below the axes entirely: every interior region is occupied by a bar or its label
    ax.text(0.5, -0.30, "expectation from trinucleotide context and codon usage; "
            "R\u2192Q/H/C are the CpG hotspots",
            transform=ax.transAxes, ha="center", va="top", fontsize=4.8, color=st.GREY)

# ---------------------------------------------------------------- c: gnomAD
ax = fig.add_subplot(gs[0, 9:12]); st.panel_label(ax, "c", dx=-0.28, dy=1.07)
G = find("union_sites_gnomad_v2.tsv")
if G:
    g = pd.read_csv(G, sep="\t", low_memory=False)
    b = g.band.fillna("not assessed").value_counts()
    order = [k for k in ("absent", "<1e-5", "1e-5 to 1e-4", "1e-4 to 1e-3", ">1e-3")
             if k in b.index] + [k for k in b.index if k not in
                                 ("absent", "<1e-5", "1e-5 to 1e-4", "1e-4 to 1e-3", ">1e-3")]
    vals = [100 * b[k] / b.sum() for k in order]
    cols = [st.BLUE if k in ("absent", "<1e-5") else
            (st.LIGHTGREY if "assess" not in str(k) else st.GREY) for k in order]
    PRETTY = {"absent": "absent", "lt_1e-5": "< 10$^{-5}$", "1e-5_1e-4": "10$^{-5}$\u201310$^{-4}$",
              "1e-4_1e-3": "10$^{-4}$\u201310$^{-3}$", "gt_1e-3": "> 10$^{-3}$",
              "not_assessed": "not assessed"}
    ax.barh(range(len(order))[::-1], vals, height=0.6, color=cols)
    ax.set_yticks(range(len(order))[::-1])
    ax.set_yticklabels([PRETTY.get(str(k), str(k)) for k in order], fontsize=5.0)
    for i, v in enumerate(vals):
        ax.text(v + 1, len(order) - 1 - i, "%.1f%%" % v, va="center", fontsize=5.0)
    ax.set_xlabel("% of acquired-lysine sites", fontsize=5.8)
    ax.set_xlim(0, max(vals) * 1.34)
    # the two low bands as a fraction of ASSESSED sites, which is how the figure is quoted
    HEAD["c_band_pct"] = {str(k): round(float(v), 2) for k, v in zip(order, vals)}
    low = sum(v for k, v in zip(order, vals) if str(k) in ("absent", "lt_1e-5"))
    na = sum(v for k, v in zip(order, vals) if "assess" in str(k))
    if na < 100:
        ax.text(0.97, 0.06, "%.1f%% of ASSESSED sites are\nabsent or below 10$^{-5}$"
                % (100 * low / (100 - na)), transform=ax.transAxes, ha="right",
                fontsize=5.0, color=st.BLUE, linespacing=1.5)

# ---------------------------------------------------------------- d: germline vs somatic
# PANELS d AND e MOVED TO FIGURE 1 (2026-09-03).
#   d  germline vs somatic  -> Figure 1g. It is the control that excludes the cheapest
#      alternative to the premise, and nothing in the main figure addressed it.
#   e  burden across every lineage -> Figure 1f, replacing the eight-lineage version. The UVM
#      control is the point of the panel and a truncated list cannot carry it.
# The two panels below are relettered d and e and moved up into the vacated row.

# ---------------------------------------------------------------- d: which residues are gained
ax = fig.add_subplot(gs[1, 0:6]); st.panel_label(ax, "d", dx=-0.12, dy=1.07)
# Lysine is 4th, not 1st. Stating that plainly is better than inviting the question: the
# argument was never that lysine is uniquely acquired, but that it is the most-acquired residue
# covalent chemistry can address. Leucine, serine and valine gains are as numerous and inert.
gain = M.groupby("alt_aa").n_sites.sum().reindex(AA).fillna(0)
lost = M.groupby("ref_aa").n_sites.sum().reindex(AA).fillna(0)
o = gain.sort_values(ascending=False)
HEAD["d_gain_rank"] = [str(a) for a in o.index[:12]]
HEAD["d_K_gained"] = int(gain.get("K", 0)); HEAD["d_C_gained"] = int(gain.get("C", 0))
HEAD["d_K_rank"] = int(list(o.index).index("K") + 1)
HEAD["d_C_rank"] = int(list(o.index).index("C") + 1)
cols = [st.ORANGE if a == "K" else (st.GREEN if a == "C" else st.LIGHTGREY) for a in o.index]
ax.bar(range(20), o.values / 1000, color=cols, width=0.72)
ax.set_xticks(range(20)); ax.set_xticklabels(o.index, fontsize=5.0)
ax.set_ylabel("sites gained (thousands)", fontsize=5.8)
ax.set_ylim(0, o.max() / 1000 * 1.22)
for a, c in (("K", st.ORANGE), ("C", st.GREEN)):
    i = list(o.index).index(a)
    ax.text(i, o[a] / 1000 + 8, a, ha="center", fontsize=6, fontweight="bold", color=c)
ax.text(0.98, 0.95, "lysine 4th at %sk; cysteine 11th at %sk\n"
        "\u2014 %.2f-fold more acquired lysines than cysteines"
        % (int(gain['K'] / 1000), int(gain['C'] / 1000), gain['K'] / gain['C']),
        transform=ax.transAxes, ha="right", va="top", fontsize=5.2, color=st.GREY,
        linespacing=1.5)

# ---------------------------------------------------------------- g: net flux
ax = fig.add_subplot(gs[1, 6:12]); st.panel_label(ax, "e", dx=-0.12, dy=1.07)
net = (gain - lost).sort_values(ascending=False)
HEAD["e_K_lost"] = int(lost.get("K", 0)); HEAD["e_C_lost"] = int(lost.get("C", 0))
HEAD["e_K_net"] = int(net.get("K", 0)); HEAD["e_C_net"] = int(net.get("C", 0))
HEAD["e_net_rank"] = [str(a) for a in net.index[:6]]
cols = [st.ORANGE if a == "K" else (st.GREEN if a == "C" else
        (st.BLUE if v > 0 else st.LIGHTGREY)) for a, v in net.items()]
ax.bar(range(20), net.values / 1000, color=cols, width=0.72)
ax.axhline(0, color=st.GREY, lw=0.8)
ax.set_xticks(range(20)); ax.set_xticklabels(net.index, fontsize=5.0)
ax.set_ylabel("gained \u2212 lost (thousands)", fontsize=5.8)
ax.text(0.98, 0.95, "cysteine is rarely LOST (%sk vs %sk), so it nets\n"
        "higher on far smaller numbers. Lysine turns over more."
        % (int(lost['C'] / 1000), int(lost['K'] / 1000)),
        transform=ax.transAxes, ha="right", va="top", fontsize=5.2, color=st.GREY,
        linespacing=1.5)

# ---------------------------------------------------------------- f: the six classes
# From the retired figS2. The evidence behind Figure 1b's claim: gains are one substitution
# class, losses are every other one.
ax = fig.add_subplot(gs[2, 0:4]); st.panel_label(ax, "f", dx=-0.24, dy=1.07)
ax.set_title("Gains are one class.\nLosses are every other one.", loc="left", pad=5, fontsize=6.6)
HEAD["f_gain_by_class"] = {c: round(float(SP.loc["gain", c]), 2) for c in CLASSES}
HEAD["f_loss_by_class"] = {c: round(float(SP.loc["loss", c]), 2) for c in CLASSES}
HEAD["f_allmis_by_class"] = {c: round(float(SP.loc["all missense", c]), 2) for c in CLASSES}
x = np.arange(len(CLASSES)); w = 0.27
ax.bar(x - w, [SP.loc["all missense", c] for c in CLASSES], w, color=st.LIGHTGREY,
       label="all missense")
ax.bar(x, [SP.loc["gain", c] for c in CLASSES], w, color=st.ORANGE, label="lysine gained")
ax.bar(x + w, [SP.loc["loss", c] for c in CLASSES], w, color=st.BLUE, label="lysine lost")
ax.set_xticks(x); ax.set_xticklabels(CLASSES, fontsize=5.0)
ax.set_ylabel("share of substitutions (%)", fontsize=5.8)
ax.set_ylim(0, 118); ax.tick_params(labelsize=5.4)
ax.legend(fontsize=4.8, loc="upper left", frameon=False, handlelength=0.8,
          handletextpad=0.3, borderpad=0.12, bbox_to_anchor=(-0.03, 1.04))
ax.annotate("", xy=(2.0, 89.0), xytext=(2.0, 106.0),
            arrowprops=dict(arrowstyle="-|>", color=st.ORANGE, lw=0.9, mutation_scale=5))
ax.text(2.20, 107.0, "85.5%", fontsize=5.4, color=st.ORANGE, fontweight="bold", va="top")
ax.annotate("", xy=(2.30, 5.0), xytext=(3.05, 30.0),
            arrowprops=dict(arrowstyle="-|>", color=st.BLUE, lw=0.9, mutation_scale=5))
ax.text(3.12, 32.0, "0 of\n60,240", fontsize=5.2, color=st.BLUE, fontweight="bold",
        va="bottom", linespacing=1.4)
ax.text(0.5, -0.26, "missense only, so the three sets are like for like;\nnonsense adds 3,020 "
        "further losses, all T>A", transform=ax.transAxes, ha="center", va="top", fontsize=4.8,
        color=st.GREY, linespacing=1.5)

# ---------------------------------------------------------------- g: the signature control
# From the retired figS2. The control that separates the codon explanation from the cheaper
# one — that C>T simply dominates cancer mutation, so the gain figure is inherited from the
# background. That predicts the asymmetry should weaken where C>T is not dominant. It does not.
ax = fig.add_subplot(gs[2, 4:12]); st.panel_label(ax, "g", dx=-0.085, dy=1.07)
ax.set_title("The asymmetry is not inherited from the background spectrum", loc="left", pad=5,
             fontsize=6.6)
xx = np.arange(len(CT))
ax.plot(xx, CT.bg_CT, "-", color=st.LIGHTGREY, lw=1.1, zorder=1)
ax.plot(xx, CT.bg_CT, "s", ms=2.2, color=st.GREY, zorder=3, label="all mutations")
ax.plot(xx, CT.gain_CT, "-", color=st.ORANGE, lw=0.9, zorder=2)
ax.plot(xx, CT.gain_CT, "o", ms=2.6, color=st.ORANGE, zorder=4, label="K gained")
ax.plot(xx, CT.loss_CT, "o", ms=2.6, color=st.BLUE, zorder=4, label="K lost")
ax.set_xticks(xx)
ax.set_xticklabels(CT.cancer, fontsize=4.6, rotation=45, ha="right", rotation_mode="anchor")
ax.set_ylabel("C>T share (%)", fontsize=5.8)
ax.set_ylim(-6, 120); ax.tick_params(axis="y", labelsize=5.4)
ax.legend(fontsize=4.8, loc="upper left", frameon=False, handlelength=0.9, handletextpad=0.3,
          borderpad=0.12, ncol=3, bbox_to_anchor=(-0.012, 1.05), columnspacing=0.6)
ax.text(len(CT)/2 - 0.5, 13.0, "lysine lost: 0.0% in every type", fontsize=5.0, color=st.BLUE,
        ha="center", va="center", fontweight="bold")
ax.text(0.5, -0.40, "Ordered by each cancer's own C>T content, tobacco-driven lung (LUAD, 22%) "
        "to ultraviolet melanoma (SKCM, 86%).\nGains follow it; losses do not.",
        transform=ax.transAxes, ha="center", va="top", fontsize=4.8, color=st.GREY,
        linespacing=1.5)

# ============================================ h  added 2026-09-07
# Figure 1g reports one residue's somatic:germline ratio. This panel reports all twenty, computed
# by step54_germline_by_residue.py on the identical pass over dbNSFP 5.4a that produced Figure 1g's
# lysine figure -- same file, same MANE filter, same gnomAD v4.1 joint frequency field, same
# expected-diploid-count definition. The lysine value it returns is 471.1 per individual, which is
# Figure 1g's number to one decimal place; that agreement is the check that generalising the
# calculation did not change it.
#
# WHY THE WHOLE DISTRIBUTION AND NOT JUST LYSINE. Panel f shows all 32 lineages because a truncated
# ranking cannot carry its own control. The same applies here: 2.05 means nothing without the
# nineteen values it sits among. It also supplies the comparator the field actually uses -- cysteine,
# on identical terms rather than from a differently computed literature figure.
#
# WHAT THE PANEL IS NOT. The ratio is not a measure of selection in cancer, and the caption says so.
# Germline shares are depressed by purifying selection, so residues that are disruptive to introduce
# are purged from inherited variation and score high here for a reason unrelated to cancer -- which
# is the most likely reading of F, Y and W at the top. Selection is Figure 1d's odds ratio of 1.68,
# and the context-corrected enrichment is panel b's 1.12-fold.
ax = fig.add_subplot(gs[3, 0:12]); st.panel_label(ax, "h", dx=-0.055, dy=1.07)
_R = find("germline_by_residue.csv")
if _R:
    RB = pd.read_csv(_R).sort_values("ratio", ascending=False).reset_index(drop=True)
    COL = {"K": st.ORANGE, "C": st.BLUE}
    cols = [COL.get(r.res, st.LIGHTGREY) for r in RB.itertuples()]
    x = np.arange(len(RB))
    ax.bar(x, RB.ratio, color=cols, width=0.72)
    ax.axhline(1.0, color=st.GREY, lw=0.9, ls="--")
    ax.text(len(RB) - 0.4, 1.05, "created just as often by inheritance", fontsize=5.2,
            color=st.GREY, ha="right", va="bottom")
    for i, r in enumerate(RB.itertuples()):
        if r.res in COL:
            ax.text(i, r.ratio + 0.05, "%.2f" % r.ratio, ha="center", fontsize=6,
                    fontweight="bold", color=COL[r.res])
            ax.text(i, -0.16, {"K": "lysine", "C": "cysteine"}[r.res], ha="center", fontsize=5.4,
                    color=COL[r.res], fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(RB.res, fontsize=6)
    ax.set_xlabel("residue the mutation creates", fontsize=6)
    ax.set_ylabel("times more often created by cancer\nthan by inherited variation", fontsize=6)
    ax.set_ylim(-0.28, max(RB.ratio) * 1.16)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    _k = RB[RB.res == "K"].iloc[0]; _c = RB[RB.res == "C"].iloc[0]
    HEAD["h_K_ratio"] = round(float(_k.ratio), 2)
    HEAD["h_C_ratio"] = round(float(_c.ratio), 2)
    HEAD["h_K_rank"] = int(RB[RB.res == "K"].index[0]) + 1
    HEAD["h_C_rank"] = int(RB[RB.res == "C"].index[0]) + 1
    # Split rather than stored as one formatted sentence: a display string is not a checkable
    # quantity, and verify.py should be diffing the numbers the caption quotes.
    HEAD["h_ratio_min"] = round(float(RB.ratio.min()), 2)
    HEAD["h_ratio_min_res"] = str(RB.iloc[-1].res)
    HEAD["h_ratio_max"] = round(float(RB.ratio.max()), 2)
    HEAD["h_ratio_max_res"] = str(RB.iloc[0].res)
    HEAD["h_ratio_median"] = round(float(RB.ratio.median()), 2)
    HEAD["h_below_one"] = int((RB.ratio < 1).sum())
    fd.panel("figS1", "h", RB, note="somatic:germline share, all 20 residues, step54 output")
else:
    ax.set_axis_off()
    ax.text(0.5, 0.5, "germline_by_residue.csv not found \u2014 run step54_germline_by_residue.py",
            ha="center", va="center", fontsize=6, color=st.GREY)

st.save(fig, "figS1")
fd.head("figS1", HEAD)
for _e in (".png", ".pdf"):          # review copies alongside fig4_v2 / figS_react_v2
    fig.savefig("../results/figS1_v3" + _e, dpi=400, bbox_inches="tight")
# the germline/somatic and per-lineage summaries moved to Figure 1 with their panels
def _fmt(v):
    if isinstance(v, float):
        return "%.4g" % v
    if isinstance(v, int):
        return f"{v:,}"
    return str(v)


print("figS1 headline numbers")
for kk in sorted(HEAD):
    print("  %-22s %s" % (kk, _fmt(HEAD[kk])))
print("figS1: 8 panels (a-h)")
