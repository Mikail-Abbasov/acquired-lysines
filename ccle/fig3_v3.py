#!/usr/bin/env python3
"""
Figure 3 — how the false-discovery rate was controlled, and what it cost.

THE PROBLEM THIS FIGURE ANSWERS
Standard target-decoy FDR assumes the candidate population has a reasonable prior of being
present. A variant database does not: most of the 8,470 acquired-lysine peptides are absent from
any given cell line, so a 1% FDR computed against the whole search is meaningless for the
variant subset. Two things were done about it — a two-stage search that establishes the
reference population first, and a paired entrapment set that measures the false-discovery
proportion directly rather than estimating it.

a  the design: two-stage search, and what each stage is for
b  the score margin dose-response — how FDP falls as the margin requirement rises. This is what
   justifies the >= 10 threshold; without it the threshold is arbitrary
c  the three populations at margin >= 10, side by side
d  per-plex consistency: is the entrapment rate stable across 42 experiments, or driven by a
   few outliers?
e  the published estimator (Wen et al. 2025) applied to plex 2, for comparison with our own
   count-based figure

ENTRAPMENT DESIGN. Each variant peptide is paired with an entrapment peptide of matched length
and composition that cannot be present. An entrapment identification is therefore a direct
observation of a false positive, not an estimate.
"""
import os, sys, re, csv, glob
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

st.apply()

def find(n):
    for d in (".", "..", os.path.join(_TOP, 'ccle')):
        p = os.path.join(d, n)
        if os.path.exists(p): return p
    return None

# ---------------------------------------------------------------- gather margins and labels
CACHE = "fig3_margins.csv"
if os.path.exists(CACHE):
    D = pd.read_csv(CACHE)
    print("margins from cache: %s rows" % f"{len(D):,}")
else:
    rows = []
    for pf in sorted(glob.glob("plex*_varpair.txt")):
        plex = pf.split("_")[0]
        g = glob.glob("out_%s_pep/*/psm.tsv" % plex)
        if not g:
            continue
        lab = {}
        with open(pf) as fh:
            fh.readline()
            for ln in fh:
                f = ln.split("\t", 5)
                if len(f) > 3 and f[3] in ("target", "p_target"):
                    lab[f[0]] = f[3]
        best = {}
        for r in csv.DictReader(open(g[0]), delimiter="\t"):
            p = r.get("Peptide", "").strip()
            if p not in lab:
                continue
            try:
                d = float(r.get("Hyperscore") or 0) - float(r.get("Nextscore") or 0)
            except ValueError:
                continue
            if p not in best or d > best[p][0]:
                best[p] = (d, lab[p])
        for p, (d, l) in best.items():
            rows.append(dict(plex=plex, peptide=p, margin=d, label=l))
    D = pd.DataFrame(rows)
    # Never write an empty cache. In a deposit without the raw search output this branch finds no
    # psm.tsv, and an empty fig3_margins.csv here is worse than none: it is written once, then
    # every later run reads it and dies on "No columns to parse from file" with no clue why.
    if D.empty:
        raise SystemExit(
            "fig3: no PSM tables found, so no margins could be computed, and an empty cache was "
            "NOT written.\n       This figure needs either fig3_margins.csv (deposited) or the "
            "out_plex*_pep/*/psm.tsv search output.")
    D.to_csv(CACHE, index=False)
    print("margins computed: %s rows across %d plexes" % (f"{len(D):,}", D.plex.nunique()))

S = pd.read_csv(find("downstream_summary.csv"))

fig = plt.figure(figsize=(st.DOUBLE, st.SINGLE * 2.30))
gs = fig.add_gridspec(3, 12, hspace=0.58, wspace=2.4,
                      left=0.065, right=0.985, top=0.955, bottom=0.075)

# ============================================== six panels
# Panel d of the previous version — the Wen et al. estimator over the whole plex-2 search —
# stated in its own annotation that the populations differ and the figure is not comparable
# with the 0.58% next to it. A main-text panel showing a number that cannot be compared to
# our own belongs in the supplement, and is now figS_fdr h.


# ---------------------------------------------------------------- a: the design
ax = fig.add_subplot(gs[0, 0:5]); st.panel_label(ax, "a", dx=-0.22, dy=1.09)
ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.set_axis_off()
boxes = [(7.9, "stage 1 \u2014 reference space", st.LIGHTGREY,
          "MSFragger 4.4.1 \u00b7 Percolator \u00b7 Philosopher\n"
          "MANE v1.2 + contaminants, reversed decoys\ntarget-decoy FDR behaves here"),
         (4.9, "stage 2 \u2014 variant space", st.BLUE,
          "same engine, database restricted to\nacquired-lysine peptides not explained\n"
          "by stage 1 \u00b7 8,470 sites, DepMap 26Q1"),
         (1.9, "paired entrapment", st.ORANGE,
          "FDRBench 1.1.1 \u00b7 shuffled partners of\nmatched length and composition\n"
          "every hit is an observed false positive")]
for y, title, col, sub in boxes:
    ax.add_patch(FancyBboxPatch((0.4, y), 9.2, 2.05, boxstyle="round,pad=0.08",
                                facecolor=st.PAPER, edgecolor=col, lw=1.1))
    ax.text(0.8, y + 1.68, title, fontsize=5.8, fontweight="bold", color=col)
    ax.text(0.8, y + 0.78, sub, fontsize=4.5, color=st.GREY, va="center", linespacing=1.5)
for y in (7.7, 4.7):
    ax.annotate("", xy=(5.0, y - 0.35), xytext=(5.0, y),
                arrowprops=dict(arrowstyle="-|>", color=st.GREY, lw=0.9))
ax.text(5.0, 0.75, "a variable-modification search pools all of this into one\n"
        "FDR calculation, which the low variant prior breaks",
        fontsize=4.5, color=st.GREY, ha="center", linespacing=1.45)

# ---------------------------------------------------------------- b: margin dose-response
ax = fig.add_subplot(gs[0, 6:12]); st.panel_label(ax, "b", dx=-0.24, dy=1.09)
th = np.arange(0, 31, 1)
fdp, n_id = [], []
for t in th:
    sub = D[D.margin >= t]
    nt = int((sub.label == "target").sum()); ne = int((sub.label == "p_target").sum())
    fdp.append(100 * ne / max(nt + ne, 1)); n_id.append(nt)
HEAD = {}
_ne_lad = []
for _t in th:
    _s = D[D.margin >= _t]
    _ne_lad.append(int((_s.label == "p_target").sum()))
fd.panel("fig3", "b", [dict(margin_threshold=int(_t), fdp_pct=_f,
                            n_variant_peptides=int(_n), n_entrapment_peptides=int(_e))
                       for _t, _f, _n, _e in zip(th, fdp, n_id, _ne_lad)],
         note="distinct peptides, best PSM each")
HEAD["b_n_rows"] = len(D)
for _t in (0, 2, 5, 10, 15, 20):
    _i = list(th).index(_t)
    HEAD["b_fdp_margin%d_pct" % _t] = round(fdp[_i], 2)
    HEAD["b_n_variant_margin%d" % _t] = int(n_id[_i])
    HEAD["b_n_entrapment_margin%d" % _t] = int(_ne_lad[_i])
ax.plot(th, fdp, color=st.ORANGE, lw=1.4)
ax.axvline(10, color=st.BLUE, lw=1.0, ls="--")
ax.set_xlabel("score margin threshold", fontsize=6)
ax.set_ylabel("false-discovery proportion (%)", fontsize=6, color=st.ORANGE)
ax.tick_params(axis="y", colors=st.ORANGE)
ax2 = ax.twinx()
ax2.plot(th, n_id, color=st.GREY, lw=1.0, ls=":")
ax2.set_ylabel("identifications retained", fontsize=6, color=st.GREY)
ax2.tick_params(axis="y", colors=st.GREY); ax2.spines["right"].set_visible(True)
i10 = list(th).index(10)
# This curve counts DISTINCT PEPTIDES, best PSM each, which is what the entrapment pairing is
# defined over. downstream_summary counts identifications per plex and gives a larger number
# (1,033) for the same threshold. Both are correct; they are not the same denominator, and the
# panel says so rather than letting the two numbers look like a contradiction.
ax.text(11, max(fdp) * 0.72, "margin \u2265 10\nFDP %.2f%%\n%s distinct peptides"
        % (fdp[i10], f"{n_id[i10]:,}"), fontsize=5.0, color=st.BLUE, linespacing=1.45)
ax.text(0.5, -0.30, "distinct peptides, best PSM each", transform=ax.transAxes,
        ha="center", va="top", fontsize=4.6, color=st.GREY)

# ---------------------------------------------------------------- c: per-plex consistency
ax = fig.add_subplot(gs[1, 0:12]); st.panel_label(ax, "c", dx=-0.14, dy=1.09)
Sp = S.sort_values("var", ascending=False)
xx = np.arange(len(Sp))
# grouped, not stacked, and log-scaled: one entrapment hit against 79 identifications is
# invisible as a sliver at the base of a linear bar
# DOTS, NOT BARS. A bar encodes its value as a length, and on a log axis length is not
# proportional to value: drawn as bars, one entrapment hit rose to about 40% of the height of a
# 25-identification bar, which reads as a large false-discovery share against a true 0.58%.
# A point encodes position only, so the log axis stays honest.
mask = (Sp.var_entrap > 0).values
fd.panel("fig3", "c", Sp[["plex", "var", "var_entrap"]].rename(
    columns={"var": "n_variant_identifications", "var_entrap": "n_entrapment_identifications"}),
    note="identifications at margin >= 10, per plex")
HEAD["c_n_entrapment_identifications"] = int(S.var_entrap.sum())
HEAD["c_n_variant_identifications"] = int(S["var"].sum())
HEAD["c_fdp_pooled_pct"] = round(100 * S.var_entrap.sum() / max(S["var"].sum(), 1), 2)
HEAD["c_n_plexes_with_entrapment"] = int((S.var_entrap > 0).sum())
HEAD["c_n_plexes"] = int(len(S))
ax.scatter(xx, Sp["var"], s=7, color=st.GREY, zorder=3, label="variant identifications")
ax.scatter(xx[mask], Sp.var_entrap[mask], s=16, color=st.ORANGE, zorder=4,
           label="entrapment hits")
ax.set_yscale("log")
ax.set_ylim(0.75, Sp["var"].max() * 2.6)
for i in np.where(mask)[0]:
    ax.text(xx[i], Sp.var_entrap.iloc[i] * 1.45, "%d" % Sp.var_entrap.iloc[i],
            ha="center", fontsize=4.4, color=st.ORANGE, fontweight="bold")
ax.set_xticks(xx); ax.set_xticklabels(["%02d" % p for p in Sp.plex], fontsize=3.6, rotation=90)
ax.set_xlabel("plex", fontsize=6)
ax.set_ylabel("identifications at margin \u2265 10", fontsize=6)
ax.legend(frameon=False, fontsize=5.0, loc="upper right", bbox_to_anchor=(1.0, 0.84),
          handlelength=1.1, scatterpoints=1)
# identifications per plex, from downstream_summary. Panel g counts DISTINCT PEPTIDES and gives
# 7 — the same events under a different denominator, stated in both panels so the two numbers
# do not read as a contradiction.
ax.text(0.985, 0.99, "%d in %s identifications \u2014 FDP %.2f%%\n"
        "falling in %d of %d plexes"
        % (int(S.var_entrap.sum()), f"{int(S['var'].sum()):,}",
           100 * S.var_entrap.sum() / max(S['var'].sum(), 1),
           int((S.var_entrap > 0).sum()), len(S)),
        transform=ax.transAxes, ha="right", va="top", fontsize=5.0, color=st.GREY,
        linespacing=1.45)


# ---------------------------------------------------------------- d: SpecQuality
ax = fig.add_subplot(gs[2, 0:4]); st.panel_label(ax, "d", dx=-0.22, dy=1.09)
# SpecQuality scores a spectrum on its own properties, with no reference to the search engine,
# the database or the FDR.
#
# SQS_GM CANNOT BE USED HERE. It is a geometric mean across features, and the hit spectra were
# scored from an MGF export in which the Intensity feature is zero (median 0.000 against 6.818
# for the mzML background). A single zero collapses a geometric mean, which is the entire
# 1.31-vs-6.96 difference an earlier version of this panel showed — a file-format artifact, not
# a quality difference. PeakCount and SNR came through the export intact and are compared here.
import glob as _g
SQ = None
for c in ("../SpecQuality", os.path.join(_TOP, 'SpecQuality')):
    if os.path.isdir(c):
        SQ = c; break
if SQ:
    H = pd.read_csv(os.path.join(SQ, "hit_spectra_mgf_features.tsv"), sep="\t")
    bg = []
    for f in sorted(_g.glob(os.path.join(SQ, "sq_out", "*_features.tsv")))[:20]:
        try:
            bg.append(pd.read_csv(f, sep="\t", usecols=["PeakCount", "SNR", "Intensity"]))
        except Exception:
            continue
    B = pd.concat(bg, ignore_index=True) if bg else pd.DataFrame()
    feats = [("PeakCount", "peak count"), ("SNR", "signal-to-noise")]
    xs = np.arange(len(feats)); w = 0.34
    hv = [pd.to_numeric(H[f], errors="coerce").dropna() for f, _ in feats]
    bv = [pd.to_numeric(B[f], errors="coerce").dropna() for f, _ in feats]
    # each feature on its own scale: peak count is in the hundreds, SNR near 0.2
    for k, ((f, lab), a, b) in enumerate(zip(feats, hv, bv)):
        norm = b.median() if b.median() else 1
        ax.bar(k - w/2, a.median() / norm, w, color=st.BLUE)
        ax.bar(k + w/2, 1.0, w, color=st.LIGHTGREY)
        # the bar is a ratio, so the label on it is the ratio; the raw medians, which are on a
        # different scale entirely (hundreds vs 0.2), go underneath as text
        ax.text(k - w/2, a.median() / norm + 0.03, "%.2f\u00d7" % (a.median() / norm),
                ha="center", fontsize=5.0, color=st.BLUE)
        ax.text(k + w/2, 1.03, "1.00\u00d7", ha="center", fontsize=5.0, color=st.GREY)
        ax.text(k, -0.14, "medians %.3g vs %.3g" % (a.median(), b.median()), ha="center",
                va="top", transform=ax.get_xaxis_transform(), fontsize=4.4, color=st.GREY)
    ax.axhline(1.0, color=st.GREY, lw=0.8, ls="--")
    ax.set_xticks(xs); ax.set_xticklabels([l for _, l in feats], fontsize=5.8)
    ax.set_ylabel("median, relative to background", fontsize=6)
    ax.set_ylim(0, 1.45)
    ax.text(0.02, 0.96, "acquired-lysine hits (%d)" % len(H), transform=ax.transAxes,
            fontsize=5.0, color=st.BLUE, va="top")
    ax.text(0.02, 0.89, "all spectra (%s)" % f"{len(B):,}", transform=ax.transAxes,
            fontsize=5.0, color=st.GREY, va="top")
    ax.text(0.5, -0.26, "SQS$_{GM}$ omitted: the MGF export zeroed its Intensity term",
            transform=ax.transAxes, ha="center", va="top", fontsize=4.6, color=st.GREY)
    _drows = []
    for (f, lab), a, b in zip(feats, hv, bv):
        _drows.append(dict(feature=lab, hit_median=float(a.median()),
                           background_median=float(b.median()),
                           ratio_hit_over_background=float(a.median() / (b.median() or 1))))
        HEAD["d_%s_ratio" % lab.replace("-", "_").replace(" ", "_")] = \
            round(float(a.median() / (b.median() or 1)), 2)
    fd.panel("fig3", "d", _drows, note="SpecQuality medians, hits vs background")
    HEAD["d_n_hit_spectra"] = int(len(H))
    HEAD["d_n_background_spectra"] = int(len(B))

# ---------------------------------------------------------------- e: precursor mass error
ax = fig.add_subplot(gs[2, 5:8]); st.panel_label(ax, "e", dx=-0.28, dy=1.09)
# Observed minus calculated, from the PSM tables. An earlier version plotted the fractional part
# of the observed mass alone, which is a property of the peptide's own composition and says
# nothing about error. A cluster at +1 Da is the signature that withdrew YBX1 E4K: one neutron,
# a deamidation or isotope error rather than the peptide claimed.
CACHE_DM = "fig3_deltamass.csv"
if os.path.exists(CACHE_DM):
    DM = pd.read_csv(CACHE_DM)
else:
    rows = []
    for f in sorted(_g.glob("out_plex*_pep/*/psm.tsv")):
        try:
            q = pd.read_csv(f, sep="\t", low_memory=False,
                            usecols=lambda c: c in ("Protein", "Delta Mass"))
        except Exception:
            continue
        if "Delta Mass" not in q.columns:
            continue
        pr = q.get("Protein", pd.Series([""] * len(q))).astype(str)
        rows.append(pd.DataFrame(dict(
            dm=pd.to_numeric(q["Delta Mass"], errors="coerce"),
            grp=np.where(pr.str.startswith("VAR|"), "variant", "background"))))
    DM = pd.concat(rows, ignore_index=True).dropna()
    DM.to_csv(CACHE_DM, index=False)
v = DM.loc[DM.grp == "variant", "dm"]
b = DM.loc[DM.grp == "background", "dm"].sample(min(200000, int((DM.grp == "background").sum())),
                                                random_state=0)
# Zoomed to the monoisotopic population. Over the full +/-1.5 Da range 91% of the mass sits in
# one bin at zero and the plot shows nothing; here the actual accuracy is visible, and variant
# errors sit on top of background rather than merely falling in the same bins.
_v1 = 100 * ((v - 1.00336).abs() < 0.05).mean()
_b1 = 100 * ((b - 1.00336).abs() < 0.05).mean()
_vz, _bz = v[v.abs() < 0.06], b[b.abs() < 0.06]
HEAD["e_pct_variant_within_0.05_da"] = round(100 * float((v.abs() < 0.05).mean()), 1)
HEAD["e_median_variant_da"] = round(float(_vz.median()), 4)
HEAD["e_median_background_da"] = round(float(_bz.median()), 4)
HEAD["e_plus1_variant_pct"] = round(float(_v1), 1)
HEAD["e_plus1_background_pct"] = round(float(_b1), 1)
HEAD["e_n_variant_psms"] = int(len(v))
bins = np.linspace(-0.06, 0.06, 49)
ax.hist(_bz, bins=bins, density=True, color=st.LIGHTGREY,
        label="background (%s)" % f"{len(b):,}")
ax.hist(_vz, bins=bins, density=True, histtype="step", lw=1.4, color=st.BLUE,
        label="variant (%s)" % f"{len(v):,}")
ax.axvline(0, color=st.GREY, lw=0.8, ls="--")
_eb = np.linspace(-0.06, 0.06, 49)
_hv_, _ = np.histogram(_vz, bins=_eb); _hb_, _ = np.histogram(_bz, bins=_eb)
fd.panel("fig3", "e", [dict(bin_left_da=float(_eb[_k]), bin_right_da=float(_eb[_k + 1]),
                            n_variant=int(_hv_[_k]), n_background=int(_hb_[_k]))
                       for _k in range(len(_eb) - 1)],
         note="precursor mass error, monoisotopic window; background is a 200k sample")
ax.set_xlabel("observed \u2212 calculated mass (Da)", fontsize=6)
ax.set_ylabel("density", fontsize=6)
ax.legend(frameon=False, fontsize=5.0, loc="upper left", handlelength=1.1)
ax.text(0.98, 0.70, "%.0f%% within 0.05 Da\nmedian %.4f vs %.4f Da"
        % (100 * (v.abs() < 0.05).mean(), _vz.median(), _bz.median()),
        transform=ax.transAxes, ha="right", va="top", fontsize=4.8, color=st.GREY,
        linespacing=1.45)
ax.text(0.0, -0.30, "the remaining %.1f%% of variant PSMs and %.1f%% of background sit at "
        "+1 Da \u2014 the first isotope\nselected instead of the monoisotopic peak, at the "
        "same rate in both" % (_v1, _b1),
        transform=ax.transAxes, ha="left", va="top", fontsize=4.4, color=st.GREY,
        linespacing=1.45)

# ---------------------------------------------------------------- f: paired score separation
ax = fig.add_subplot(gs[2, 9:12]); st.panel_label(ax, "f", dx=-0.28, dy=1.09)
# The empirical version of the argument. An earlier panel here plotted false positives computed
# FROM the nominal FDR and then expressed them as a fraction of that same total, which is
# algebraically just the FDR — two flat lines proving nothing.
#
# This instead uses the pairing directly: every variant peptide has exactly one entrapment
# partner (peptide_pair_index in the varpair files). Where the two score distributions overlap
# is where an identification cannot be distinguished from a peptide known to be absent.
tg = D.loc[D.label == "target", "margin"]
en = D.loc[D.label == "p_target", "margin"]
bins = np.linspace(0, 45, 34)
# Two translucent fills overlap into a third colour. Outlines overlap into nothing.
ax.hist(tg, bins=bins, density=True, histtype="step", color=st.BLUE, lw=1.3,
        label="variant targets (%s)" % f"{len(tg):,}")
ax.hist(en, bins=bins, density=True, histtype="step", color=st.ORANGE, lw=1.3,
        ls=(0, (2.6, 1.3)), label="entrapment (%d)" % len(en))
ax.axvline(10, color=st.GREY, lw=1.0, ls="--")
ax.set_xlabel("score margin", fontsize=6)
ax.set_ylabel("density", fontsize=6)
ax.legend(frameon=False, fontsize=5.0, loc="upper right", handlelength=1.1)
above = (en >= 10).sum()
fd.panel("fig3", "f", D[["plex", "peptide", "margin", "label"]],
         note="every distinct peptide with its best score margin; label target = variant, "
              "p_target = entrapment")
HEAD["f_n_entrapment_above_margin10"] = int(above)
HEAD["f_n_entrapment_peptides"] = int(len(en))
HEAD["f_n_variant_peptides"] = int(len(tg))
HEAD["f_pct_entrapment_below_margin10"] = round(100 * float((en < 10).mean()), 1)
ax.text(1.0, -0.46, "%d of %d entrapment peptides clear the threshold; %.0f%% of the "
        "entrapment mass sits below it.\nPanel c counts identifications rather than peptides, "
        "giving %d for the same events."
        % (above, len(en), 100 * (en < 10).mean(), int(S.var_entrap.sum())),
        transform=ax.transAxes, ha="right", va="top", fontsize=4.4, color=st.GREY,
        linespacing=1.45)
st.save(fig, "fig3")
fd.head("fig3", HEAD)
for _e in (".png", ".pdf"):          # review copies alongside fig4_v2 / figS_react_v2
    fig.savefig("../results/fig3_v3" + _e, dpi=400, bbox_inches="tight")
print("FDP at margin>=10: %.2f%% | variant %d, entrapment %d"
      % (fdp[i10], int(S["var"].sum()), int(S.var_entrap.sum())))
try:    print("  e precursor mass: %.0f%% within 0.05 Da | medians %.4f vs %.4f Da | +1 Da tail %.1f%% variant vs %.1f%% background"
              % (100 * (v.abs() < 0.05).mean(), _vz.median(), _bz.median(), _v1, _b1))
except Exception as _e: print("  e: %s" % _e)
print("  b margin ladder (distinct peptides, best PSM each):")
for _t in (0, 2, 5, 10, 15, 20):
    _sub = D[D.margin >= _t]
    _nt = int((_sub.label == "target").sum()); _ne = int((_sub.label == "p_target").sum())
    print("     margin >= %-2d  FDP %5.2f%%   variant %5d   entrapment %3d"
          % (_t, 100 * _ne / max(_nt + _ne, 1), _nt, _ne))
print("  c per-plex (identifications, downstream_summary): %d entrapment in %s identifications "
      "= %.2f%%; entrapment appears in %d of %d plexes"
      % (int(S.var_entrap.sum()), f"{int(S['var'].sum()):,}",
         100 * S.var_entrap.sum() / max(S["var"].sum(), 1),
         int((S.var_entrap > 0).sum()), len(S)))
print("  NOTE two denominators: panel b counts DISTINCT PEPTIDES, panel c counts "
      "IDENTIFICATIONS. Not interchangeable.")
