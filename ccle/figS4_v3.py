#!/usr/bin/env python3
"""
Supplementary figure for Figure 3, how the entrapment set was built, and what the controls did.

a  THE FOUR-WAY DESIGN. FDRBench builds four matched sets per plex: target (the acquired-lysine
   peptides), paired target (entrapment, sequences that cannot be present), decoy, and paired
   decoy. All four are matched on length and composition, so a hit in the entrapment set is a
   direct observation of a false positive at the same difficulty as a real one.

b  MATCHING, demonstrated rather than asserted: length distributions of the four sets.

c  WHAT THE SHUFFLE PRESERVES. Each entrapment peptide is its target's own residues permuted,
   with the C-terminal K or R kept so it remains tryptic. Identical residues means identical
   mass, so the search engine cannot separate target from entrapment at the precursor level at
   all, every bit of discrimination must come from fragment ions. That is what makes this a
   fair control rather than an easy one.

d  THE DOWNSTREAM BOUNDARY CONTROL across 42 experiments. A downstream peptide begins
   immediately after the acquired lysine and carries no mutated residue; its only claim to being
   mutant-derived is that its N-terminal boundary exists solely because the mutation created a
   tryptic site there. The boundary control is the same region shifted one residue, tryptic in
   neither protein. If the search were matching sequences indiscriminately, both would be
   identified at similar rates.

e  PER-PLEX FDP, to show whether the pooled figure is stable or driven by a few experiments.

f  NOMINAL AGAINST OBSERVED. The reason the design carries decoys as well as entrapment: with
   both, you can ask whether decoy-based FDR is accurate. It is not, nominal q reaches 0.60%
   in the plex-2 search where the entrapment set puts the observed proportion at 0.99%, an
   underestimate of about 1.7-fold. This is the empirical case for not trusting nominal FDR
   alone in variant space, which panel a of Figure 3 otherwise only asserts.

g  ALL FOUR SETS on one score axis, showing entrapment behaves as a population that cannot be
   present should.

WHAT IS NOT HERE. The channel test, which asks whether a variant peptide's reporter signal
sits in the cell line carrying the mutation, and which rejected RAB2B, PFKM, ACBD3 and MYH10 ,
covers only 6 acquired-lysine sites across 4 genes in the data currently on disk. The main
searches ran with quantification disabled, and TMT-Integrator cannot read the mzBIN files they
produced. Presenting four genes as evidence that the filter works across a list of 61 would
overstate it.
"""
import os, sys, re, glob
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
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
BODY = 7.0
st.apply(BODY)

HEAD = {}       # headline numbers per panel, captured where they are computed
MARGIN = 10     # score-margin threshold -- one definition for the axis line and the labels

def find(n):
    for d in (".", "..", os.path.join(_TOP, 'ccle_plex05')):
        p = os.path.join(d, n)
        if os.path.exists(p): return p
    return None

VP = sorted(glob.glob("plex*_varpair.txt"))
_parts = []
for _k, _f in enumerate(VP[:8]):
    _d = pd.read_csv(_f, sep="\t")
    _d["src"] = _k                      # pair_index restarts in every file
    _parts.append(_d)
V = pd.concat(_parts, ignore_index=True)
V["L"] = V.sequence.str.len()
S = pd.read_csv(find("downstream_summary.csv"))
print("varpair rows read: %s from %d plexes" % (f"{len(V):,}", min(8, len(VP))))

fig = plt.figure(figsize=(st.DOUBLE, st.SINGLE * 2.76))
gs = fig.add_gridspec(4, 12, hspace=0.50, wspace=2.4,
                      left=0.065, right=0.985, top=0.955, bottom=0.065)

# ---------------------------------------------------------------- a: the design
ax = fig.add_subplot(gs[0, 0:4]); st.head(fig, ax, "a")
ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.set_axis_off()
n_each = int(V.peptide_type.value_counts().min())
cells = [(0.6, 6.0, "target", st.BLUE, "acquired-lysine peptides\nthat may be present"),
         (5.3, 6.0, "paired target", st.ORANGE, "matched sequences that\nCANNOT be present"),
         (0.6, 2.2, "decoy", st.GREY, "reversed target"),
         (5.3, 2.2, "paired decoy", st.GREY, "reversed entrapment")]
for x, y, name, col, sub in cells:
    ax.add_patch(FancyBboxPatch((x, y), 4.1, 2.9, boxstyle="round,pad=0.08",
                                facecolor=st.PAPER, edgecolor=col, lw=1.1))
    ax.text(x + 0.3, y + 2.35, name, fontsize=st.sz(5.8), fontweight="bold", color=col)
    ax.text(x + 0.3, y + 1.35, sub, fontsize=st.sz(5), color=st.GREY, va="center", linespacing=1.5)
ax.text(5.0, 0.9, "%s peptides in each set. The entrapment peptide is its target's OWN\n"
        "residues permuted, C-terminal K or R kept, identical mass, still tryptic.\n"
        "Every entrapment identification is an OBSERVED false positive." % f"{n_each:,}",
        fontsize=st.sz(5), color=st.GREY, ha="center", linespacing=1.5)

HEAD["a_n_varpair_rows"] = int(len(V))
HEAD["a_n_plexes"] = int(min(8, len(VP)))
HEAD["a_n_peptides_per_set"] = int(n_each)
fd.panel("figS4", "a", [dict(peptide_type=str(_t), n_peptides=int(_n))
                        for _t, _n in V.peptide_type.value_counts().items()],
         note="the four matched sets FDRBench builds per plex; one row per set")

# ---------------------------------------------------------------- b: length matching
ax = fig.add_subplot(gs[0, 5:8]); st.head(fig, ax, "b")
bins = np.arange(6, 52, 2)
_brows = []
for t, col, ls in (("target", st.BLUE, "-"), ("p_target", st.ORANGE, "--")):
    v = V.loc[V.peptide_type == t, "L"]
    ax.hist(v, bins=bins, density=True, histtype="step", lw=1.3, color=col, ls=ls,
            label="%s (median %.0f)" % ("entrapment" if t == "p_target" else t, v.median()))
    _lab = "entrapment" if t == "p_target" else t
    HEAD["b_median_length_%s_aa" % _lab] = int(round(float(v.median())))
    HEAD["b_n_%s_peptides" % _lab] = int(len(v))
    _dn, _ = np.histogram(v, bins=bins, density=True)
    _ct, _ = np.histogram(v, bins=bins)
    _brows += [dict(series=_lab, bin_left=float(bins[_i]), bin_right=float(bins[_i + 1]),
                    density=float(_dn[_i]), n=int(_ct[_i])) for _i in range(len(bins) - 1)]
fd.panel("figS4", "b", _brows, note="one row per histogram bin per set")
ax.set_xlabel("peptide length (residues)", fontsize=st.sz(6))
ax.set_ylabel("density", fontsize=st.sz(6))
ax.legend(frameon=False, fontsize=st.sz(5.0), loc="upper right", handlelength=1.4)
ax.text(0.97, 0.55, "the two distributions are\nidentical by construction",
        transform=ax.transAxes, ha="right", fontsize=st.sz(5), color=st.GREY, linespacing=1.4)

# ---------------------------------------------------------------- c: what the shuffle preserves
ax = fig.add_subplot(gs[0, 9:12]); st.head(fig, ax, "c")
# The entrapment peptide is its target's OWN residues in a different order, with the C-terminal
# K or R kept so it stays tryptic. An earlier version of this panel plotted target against
# entrapment amino-acid composition and got r = 1.0000 — true by construction, and therefore
# not evidence of anything.
#
# What the permutation does establish is that the two are identical at the PRECURSOR level:
# same residues means same mass. The search engine cannot separate them by mass at all, so
# every bit of discrimination has to come from fragment ions.
# Pair on plex AND pair index: the index restarts in every varpair file, so joining on it
# alone collides across plexes. Indexing a Series by a repeated label also made en[i][-1] a
# label lookup rather than "last character", which is what raised KeyError: -1.
V["plex"] = V.get("plex", pd.Series([""] * len(V)))
if "src" not in V.columns:
    V["src"] = 0
P = V[V.peptide_type.isin(["target", "p_target"])].copy()
P["key"] = P.src.astype(str) + "_" + P.peptide_pair_index.astype(str)
piv = P.pivot_table(index="key", columns="peptide_type", values="sequence",
                    aggfunc="first").dropna()
same_res = np.mean([sorted(a) == sorted(b)
                    for a, b in zip(piv["target"].head(2000), piv["p_target"].head(2000))])
ct_tg = piv["target"].str[-1].value_counts(normalize=True)
ct_res = piv["p_target"].str[-1].value_counts(normalize=True)
res = sorted(set(ct_tg.index) | set(ct_res.index), key=lambda a: -ct_tg.get(a, 0))[:6]
xx = np.arange(len(res)); w = st.grp(2)
ax.bar(xx - w/2, [100 * ct_tg.get(a, 0) for a in res], w, color=st.BLUE, label="target")
ax.bar(xx + w/2, [100 * ct_res.get(a, 0) for a in res], w, color=st.ORANGE,
       label="entrapment")
ax.set_xticks(xx); ax.set_xticklabels(res, fontsize=st.sz(6))
ax.set_xlabel("C-terminal residue", fontsize=st.sz(6))
ax.set_ylabel("% of peptides", fontsize=st.sz(6))
ax.legend(frameon=False, fontsize=st.sz(5.0), loc="upper right", handlelength=1.1)
# This four-line block hung below the panel and printed into the next panel's "19 ids"
# at 81%. It is a definition of what the shuffle preserves, which legend (C) now carries.

HEAD["c_pct_pairs_same_residues"] = round(float(100 * same_res), 1)
HEAD["c_n_pairs_checked"] = int(min(2000, len(piv)))
HEAD["c_n_pairs"] = int(len(piv))
for _a in res:
    HEAD["c_cterm_pct_target_%s" % _a] = round(float(100 * ct_tg.get(_a, 0)), 2)
    HEAD["c_cterm_pct_entrapment_%s" % _a] = round(float(100 * ct_res.get(_a, 0)), 2)
fd.panel("figS4", "c", [dict(cterm_residue=_a, series=_lab, pct=float(100 * _src.get(_a, 0)))
                        for _a in res
                        for _lab, _src in (("target", ct_tg), ("entrapment", ct_res))],
         note="C-terminal residue composition; one row per bar")

# ---------------------------------------------------------------- d: boundary control
ax = fig.add_subplot(gs[1, 0:5]); st.head(fig, ax, "d")
cats = [("downstream\ntarget", int(S.down.sum()), int(S.down_all.sum()), st.GREEN),
        ("boundary\ncontrol", int(S.bound.sum()), int(S.bound_all.sum()), st.ORANGE),
        ("entrapment", int(S.down_entrap.sum()), int(S.down_entrap_all.sum()), st.LIGHTGREY)]
x = np.arange(3); w = st.grp(2)
# "all" was the same bar at 45% opacity; it is now the same bar as an outline.
ax.bar(x - w/2, [c[2] for c in cats], w, facecolor=st.PAPER,
       edgecolor=[c[3] for c in cats], linewidth=1.2, label="all identifications")
ax.bar(x + w/2, [c[1] for c in cats], w, color=[c[3] for c in cats],
       label="margin \u2265 %d" % MARGIN)
for i, c in enumerate(cats):
    ax.text(i - w/2, c[2] + 1.2, "%d" % c[2], ha="center", fontsize=st.sz(5.6), color=st.GREY)
    ax.text(i + w/2, c[1] + 1.2, "%d" % c[1], ha="center", fontsize=st.sz(6.2))
ax.set_xticks(x); ax.set_xticklabels([c[0] for c in cats], fontsize=st.sz(5.6))
ax.set_ylabel("identifications across %d experiments" % len(S), fontsize=st.sz(6))
ax.set_ylim(0, max(c[2] for c in cats) * 1.35)
ax.legend(frameon=False, fontsize=st.sz(5.0), loc="upper right", handlelength=1.1)
HEAD["d_n_experiments"] = int(len(S))
HEAD["d_downstream_target_all"] = int(S.down_all.sum())
HEAD["d_downstream_target_above_margin"] = int(S.down.sum())
HEAD["d_boundary_control_all"] = int(S.bound_all.sum())
HEAD["d_boundary_control_above_margin"] = int(S.bound.sum())
HEAD["d_downstream_entrapment_all"] = int(S.down_entrap_all.sum())
HEAD["d_downstream_entrapment_above_margin"] = int(S.down_entrap.sum())
fd.panel("figS4", "d", [dict(category=_c[0].replace(chr(10), " "), n_all=int(_c[2]),
                             n_above_margin=int(_c[1])) for _c in cats],
         note="identifications across all experiments, and those surviving the score margin")
ax.text(0.0, -0.24, "the boundary control is the same region shifted one residue, tryptic in "
        "neither the wild-type nor the mutant protein",
        transform=ax.transAxes, ha="left", va="top", fontsize=st.sz(5), color=st.GREY)

# ---------------------------------------------------------------- e: per-plex FDP
ax = fig.add_subplot(gs[1, 6:12]); st.head(fig, ax, "e")
Sp = S[S["var"] > 0].copy()
Sp["fdp"] = 100 * Sp.var_entrap / (Sp["var"] + Sp.var_entrap)
Sp = Sp.sort_values("var", ascending=False)
xx = np.arange(len(Sp))
# The two tallest bars are plexes with FEW identifications, where a single entrapment hit is a
# large proportion. Bar width is scaled by identification count so that is visible rather than
# leaving a reader to read them as problem experiments.
wdt = st.BARW * (0.43 + 0.57 * (Sp["var"] / Sp["var"].max()))
ax.bar(xx, Sp.fdp, color=[st.ORANGE if f > 0 else st.LIGHTGREY for f in Sp.fdp], width=wdt)
pooled = 100 * S.var_entrap.sum() / (S["var"].sum() + S.var_entrap.sum())
ax.axhline(pooled, color=st.BLUE, lw=1.0, ls="--")
ax.text(len(Sp) * 0.62, pooled + 0.4, "pooled %.2f%%" % pooled, fontsize=st.sz(5.2), color=st.BLUE)
ax.set_xticks(xx); ax.set_xticklabels(["%02d" % p for p in Sp.plex], fontsize=st.sz(5), rotation=90)
ax.set_xlabel("plex, ordered by identifications", fontsize=st.sz(6))
ax.set_ylabel("false-discovery proportion (%)", fontsize=st.sz(6))
_hi = Sp.nlargest(2, "fdp")
ax.text(0.02, 0.97, "%d of %d experiments contributed no\nentrapment identification."
        % (int((Sp.fdp == 0).sum()), len(Sp)),
        transform=ax.transAxes, ha="left", va="top", fontsize=st.sz(5), color=st.GREY,
        linespacing=1.5)
# Bar width encodes identification count, but the tall bars are the FEW-identification plexes,
# so width alone was a hairline and carried nothing to the reader. The count is labelled.
for _i in _hi.index:
    _k = int(np.where(Sp.index == _i)[0][0])
    ax.text(_k, Sp.fdp.loc[_i] + 0.18, "%d ids" % int(Sp["var"].loc[_i]), ha="center",
            fontsize=st.sz(5), color=st.ORANGE)
HEAD["e_n_plexes"] = int(len(Sp))
HEAD["e_n_plexes_no_entrapment"] = int((Sp.fdp == 0).sum())
HEAD["e_fdp_min_pct"] = round(float(Sp.fdp.min()), 2)
HEAD["e_fdp_max_pct"] = round(float(Sp.fdp.max()), 2)
HEAD["e_pooled_fdp_pct"] = round(float(pooled), 2)
HEAD["e_n_variant_identifications"] = int(S["var"].sum())
HEAD["e_n_entrapment_identifications"] = int(S.var_entrap.sum())
fd.panel("figS4", "e", Sp[["plex", "var", "var_entrap", "fdp"]].rename(
             columns={"var": "n_variant", "var_entrap": "n_entrapment", "fdp": "fdp_pct"}),
         note="one row per plex, ordered by identifications; the bars of the panel")

# ---------------------------------------------------------------- f: nominal vs observed FDR
ax = fig.add_subplot(gs[2, 0:5]); st.head(fig, ax, "f")
# The reason the four-way design carries decoys AS WELL AS entrapment. With both, you can ask
# whether decoy-based FDR is ACCURATE: nominal q comes from the decoys, observed FDP comes from
# counting entrapment identifications, and the two should agree. Panels a-e show the entrapment
# set is fair; this shows what it reveals.
F = find("fdp_p2.csv")
if F:
    W = pd.read_csv(F).sort_values("q_value")
    q = 100 * pd.to_numeric(W.q_value, errors="coerce")
    o = 100 * pd.to_numeric(W.paired_fdp, errors="coerce")
    ax.plot(q, o, color=st.BLUE, lw=1.4)
    lim = max(q.max(), o.max()) * 1.08
    ax.plot([0, lim], [0, lim], color=st.LIGHTGREY, lw=0.9, ls="--")
    ax.set_xlim(0, lim); ax.set_ylim(0, lim)
    ax.set_xlabel("nominal q-value from decoys (%)", fontsize=st.sz(6))
    ax.set_ylabel("observed FDP from entrapment (%)", fontsize=st.sz(6))
    ratio = o.max() / q.max() if q.max() else np.nan
    ax.text(0.05, 0.92, "decoy FDR UNDERESTIMATES\nthe observed proportion\nby %.1f-fold"
            % ratio, transform=ax.transAxes, fontsize=st.sz(5.4), color=st.BLUE, va="top",
            linespacing=1.5, fontweight="bold")
    ax.text(0.97, 0.10, "plex 2, %s identifications\nnominal reaches %.2f%%, observed %.2f%%"
            % (f"{len(W):,}", q.max(), o.max()), transform=ax.transAxes, ha="right",
            fontsize=st.sz(5), color=st.GREY, linespacing=1.45)
    ax.text(0.5, -0.26, "the diagonal is where the two would agree",
            transform=ax.transAxes, ha="center", va="top", fontsize=st.sz(5), color=st.GREY)
    HEAD["f_n_identifications_plex2"] = int(len(W))
    HEAD["f_nominal_q_max_pct"] = round(float(q.max()), 2)
    HEAD["f_observed_fdp_max_pct"] = round(float(o.max()), 2)
    HEAD["f_underestimate_fold"] = round(float(ratio), 1)
    fd.panel("figS4", "f", pd.DataFrame({"nominal_q_pct": q.values,
                                         "observed_fdp_pct": o.values}).drop_duplicates(),
             note="the plotted curve, one row per distinct vertex on plex 2")

# ---------------------------------------------------------------- g: all four sets
ax = fig.add_subplot(gs[2, 6:12]); st.head(fig, ax, "g")
# Entrapment should behave like decoy — both are populations that cannot be present. If the
# entrapment set scored like the targets it would be too easy; if it scored below the decoys it
# would be too hard.
CACHE = "fig3_margins.csv"
if os.path.exists(CACHE):
    M4 = pd.read_csv(CACHE)
    bins = np.linspace(0, 45, 34)
    for lab, key, col in (("target", "target", st.BLUE),
                          ("entrapment", "p_target", st.ORANGE)):
        v = M4.loc[M4.label == key, "margin"]
        if len(v):
            ax.hist(v, bins=bins, density=True, histtype="step", lw=1.4, color=col,
                    label="%s (%s)" % (lab, f"{len(v):,}"))
    ax.axvline(MARGIN, color=st.GREY, lw=1.0, ls="--")
    ax.set_xlabel("score margin", fontsize=st.sz(6))
    ax.set_ylabel("density", fontsize=st.sz(6))
    ax.legend(frameon=False, fontsize=st.sz(5.2), loc="upper right", handlelength=1.3)
    ax.text(0.5, -0.26, "entrapment peptides are identified, and score like a population that "
            "cannot be present should",
            transform=ax.transAxes, ha="center", va="top", fontsize=st.sz(5), color=st.GREY)
    HEAD["g_margin_threshold"] = int(MARGIN)
    _grows = []
    for _lab, _key in (("target", "target"), ("entrapment", "p_target")):
        _v = M4.loc[M4.label == _key, "margin"]
        HEAD["g_n_%s_peptides" % _lab] = int(len(_v))
        if len(_v):
            HEAD["g_median_margin_%s" % _lab] = round(float(_v.median()), 2)
        _dn, _ = np.histogram(_v, bins=bins, density=True)
        _ct, _ = np.histogram(_v, bins=bins)
        _grows += [dict(series=_lab, bin_left=float(bins[_i]), bin_right=float(bins[_i + 1]),
                        density=float(_dn[_i]), n=int(_ct[_i])) for _i in range(len(bins) - 1)]
    fd.panel("figS4", "g", _grows, note="one row per histogram bin per set; score margin")


# ---------------------------------------------------------------- h: moved from Figure 3
# The field's own estimator on our data. It runs over the WHOLE plex-2 search, not the
# variant subset, so it is NOT comparable with the 0.58% in Fig. 3c — which is exactly why
# it sits here rather than beside that number.
# "not concentrated in one experiment" overstates what six events can show — with six hits you
# would expect them scattered by chance. The count is the honest headline.

# ---------------------------------------------------------------- e: published estimator
ax = fig.add_subplot(gs[3, 0:6]); st.head(fig, ax, "h")
F = find("fdp_p2.csv")
if F:
    W = pd.read_csv(F)
    for col, lab, col_c in (("combined_fdp", "combined", st.BLUE),
                            ("paired_fdp", "paired", st.GREEN),
                            ("lower_bound_fdp", "lower bound", st.LIGHTGREY)):
        if col in W.columns:
            v = pd.to_numeric(W[col], errors="coerce").dropna().sort_values()
            ax.plot(np.arange(1, len(v) + 1), 100 * v.values, lw=1.2, color=col_c,
                    ls="--" if col == "paired_fdp" else "-", label=lab)
    ax.set_xlabel("identifications, ranked by score", fontsize=st.sz(6))
    ax.set_ylabel("estimated FDP (%)", fontsize=st.sz(6))
    ax.legend(frameon=False, fontsize=st.sz(5.2), loc="lower right", handlelength=1.6)
    # NOT comparable to our variant FDP: this estimator runs over the WHOLE plex-2 search,
    # ~80,000 identifications, while our 0.58% is the variant subset across all plexes. An
    # earlier version placed the two side by side as though they measured the same thing.
    ax.text(0.02, 0.98, "Wen et al. 2025 estimator over the whole\nplex-2 search "
            "(%s identifications), not the\nvariant subset, so this is NOT comparable\n"
            "with the %.2f%% in c." % (f"{len(W):,}", pooled),
            transform=ax.transAxes, ha="left", va="top", fontsize=st.sz(5), color=st.GREY,
            linespacing=1.45)
    HEAD["h_n_identifications_plex2"] = int(len(W))
    HEAD["h_pooled_variant_fdp_pct_quoted"] = round(float(pooled), 2)
    _hrows = []
    for _col, _lab in (("combined_fdp", "combined"), ("paired_fdp", "paired"),
                       ("lower_bound_fdp", "lower bound")):
        if _col in W.columns:
            _v = pd.to_numeric(W[_col], errors="coerce").dropna().sort_values().values
            HEAD["h_max_%s_pct" % _lab.replace(" ", "_")] = round(float(100 * _v[-1]), 3)
            _u, _ix = np.unique(_v, return_index=True)
            _hrows += [dict(series=_lab, rank=int(_ix[_j] + 1),
                            estimated_fdp_pct=float(100 * _u[_j])) for _j in range(len(_u))]
    fd.panel("figS4", "h", _hrows,
             note="the three plotted curves as step vertices: rank at which each distinct value is first reached")

st.save(fig, "figS4")
fd.head("figS4", HEAD)
for _e in (".png", ".pdf"):          # review copies alongside fig4_v2 / figS_react_v2
    fig.savefig("../results/figS4_v3" + _e, dpi=400, bbox_inches="tight")
print("entrapment sets %s each | %.0f%% share residues with their target | pooled FDP %.2f%%"
      % (f"{n_each:,}", 100 * same_res, pooled))
try:    print("  d boundary control (all identifications -> surviving margin >= 10):",
              ", ".join("%s %d -> %d" % (c[0].replace(chr(10), " "), c[2], c[1]) for c in cats))
except Exception as _e: print("  d: %s" % _e)
try:    print("  e per-plex FDP: %d plexes, range %.2f%%-%.2f%%, pooled %.2f%%" % (len(Sp), Sp.fdp.min(), Sp.fdp.max(), pooled))
except Exception as _e: print("  e: %s" % _e)
try:    print("  f nominal vs observed: nominal reaches %.2f%%, observed %.2f%% on plex 2 (%s ids), %.1f-fold" % (q.max(), o.max(), f"{len(W):,}", o.max()/q.max()))
except Exception as _e: print("  f: %s" % _e)
try:    print("  h Wen et al. 2025 estimator over %s plex-2 identifications" % f"{len(W):,}")
except Exception as _e: print("  h: %s" % _e)
