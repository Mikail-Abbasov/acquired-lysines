#!/usr/bin/env python3
"""
Figure S6 — support for Figure 5: everything the central comparison rests on.

Figure 5a claims a 0.106 AUC gap between what a structural model achieves (0.617) and how
reproducible the labels are (0.723). Four things have to hold for that claim, and each gets a
panel here so a reviewer can check rather than trust:

  a  WHY NatChem IS THE COMPARATOR.  Three datasets, three ceilings, not interchangeable. The
     ceiling rises with how deeply each half is screened, and the two shallow screens return
     0.88-0.94 at engaged rates of 1-3%, where a few dozen positives make AUC unstable. NatChem
     is the only screen deep enough to estimate a ceiling at a realistic engaged rate.
  b  THE COMPOSITE BARELY BEATS ITS BEST COMPONENT: +0.019 AUC over distance-to-centroid. Shown,
     not buried — it is part of the same argument, not an embarrassment to it.
  c  THE MODEL IS NOT OVERFIT.  Train 0.619, test 0.617, optimism +0.002, held out by protein.
     The gap in Figure 5a is not the model failing to fit; it is the model fitting all there is.
  d  MORE FEATURES WOULD NOT HELP.  14 features carry about 5 independent properties.
  e  THE BREADTH GRADIENT SURVIVES ITS OWN CONFOUND.  A site tested against few chemotypes
     cannot show high breadth. Holding testing depth constant strengthens most of the gradient
     and flips three features — reported either way.
  f  A GENERAL COMPOSITE IS DEFENSIBLE.  Two unrelated chemotypes prefer the same features
     (0 of 13 differ significantly) and engage the same sites (OR 6.05, p 2e-11).
"""
import os, re, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import abbasov_style as st
import figure_data as fd
st.apply()
from adjustText import adjust_text

def find(n):
    for d in (".", "..", "ligand_rebuild", "../ligand_rebuild"):
        p = os.path.join(d, n)
        if os.path.exists(p):
            return p
    raise FileNotFoundError(n)

R10 = open(find("step10_report.txt")).read()

# The number of protein-held-out splits and of half-splits are stated in the analysis scripts,
# not in their reports, and were typed into this figure and its caption. Parse them from source
# so the caption cannot drift from the run. N_REP is env-overridable and its default is what is
# reported; step29's report does not record the value actually used.
def _const(path, name, default):
    try:
        t = open(find(path)).read()
    except Exception:
        return default
    m = re.search(name + r"\s*=\s*(?:int\(os\.environ\.get\(\s*\"[^\"]+\"\s*,\s*)?(\d+)", t)
    return int(m.group(1)) if m else default

N_SPLITS = _const("step10_composite.py", "N_SPLITS", 20)
N_HALFSPLITS = _const("step29_ceiling_fast.py", "N_REP", 50)
R29 = open(find("step29_report.txt")).read()
R16 = open(find("step16_report.txt")).read()
R13 = open(find("step13_report.txt")).read()
# read for panel g: the space's dimensionality and the size of the permutation null were
# typed into the panel, and they are stated here by the step that produced step15_pca.csv
R15 = open(find("step15_report.txt")).read()
# regenerated 2026-09-02: step10_composite.py reports an AUC but saves no per-site scores,
# so the ROC had to be reproduced from the same protocol (see step10_roc_report.txt)
ROC_S = pd.read_csv(find("step10_roc_splits.csv"))
ROC_P = pd.read_csv(find("step10_roc_pooled.csv"))
COEF  = pd.read_csv(find("step10_composite.csv"))
S15   = pd.read_csv(find("step15_pca.csv"))

NICE = lambda f: (f.replace("_8A", " 8Å").replace("_10A", " 10Å").replace("_", " ")
                   .replace("rsa tien", "exposure").replace("propka pKa", "predicted pKa")
                   .replace("nz sasa", "NZ exposure").replace("sasa sidechain", "side-chain exp.")
                   .replace("dist centroid", "dist. to centroid").replace("n lys", "nearby lysines")
                   .replace("dist nearest pocket", "dist. to pocket")
                   .replace("hbond partners NZ", "H-bond partners")
                   .replace("plddt site", "pLDDT").replace("plddt 8A", "pLDDT 8Å"))

# ---------------------------------------------------------------- parse
CEIL, cur = {}, None
for line in R29.splitlines():
    m = re.match(r"^(\w+) — (\d+) sites \| compounds per site: median (\d+)", line)
    if m:
        cur = m.group(1); CEIL[cur] = dict(n=int(m.group(2)), med=int(m.group(3)), rows=[])
    m = re.match(r"^\s*(\d+) compounds/half \|\s*(\d+) sites \|\s*([\d.]+)% engaged \| "
                 r"ceiling AUC ([\d.]+) \[([\d.]+), ([\d.]+)\]", line)
    if m and cur:
        CEIL[cur]["rows"].append(tuple(float(x) for x in m.groups()))

def auc_row(name):
    m = re.search(re.escape(name) + r"\s+([\d.]+)\s+\[([\d.]+), ([\d.]+)\]", R10)
    return tuple(float(x) for x in m.groups())
SINGLES = [(f,) + auc_row(f) for f in
           ("dist_centroid", "plddt_site", "propka_pKa", "hydrophobic_8A", "rsa_tien")]
COMP, TRAIN = auc_row("COMPOSITE (test)"), auc_row("composite (train)")

PCS = [(int(m.group(1)), float(m.group(2)), float(m.group(3)))
       for m in re.finditer(r"PC(\d)\s+([\d.]+)%\s+cumulative\s+([\d.]+)%", R10)]

# The bare pattern also matches a row of the MULTIVARIATE table further down ("1  1529  0.491
# 0.120  0.0000"), so the feature name must actually be a feature: a letter-led token with an
# underscore or a known bare name. 21 rows, not 22.
RHO = [(m.group(1), float(m.group(2)), float(m.group(3)))
       for m in re.finditer(r"^\s*([A-Za-z]\w*)\s+(-?[\d.]+)\s+[\de.+-]+\s+(-?[\d.]+)\s+[\de.+-]+\s*\*?$",
                            R16, re.M)]

T1 = [(m.group(1), float(m.group(2)), float(m.group(3)), m.group(4), float(m.group(5)))
      for m in re.finditer(r"^\s*(\w+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(yes|NO)\s+([\d.]+)\s*$",
                           R13, re.M)]
T2 = dict(both=int(re.search(r"engaged by BOTH\s+: (\d+)", R13).group(1)),
          exp=float(re.search(r"expected by chance ([\d.]+)", R13).group(1)),
          cy=int(re.search(r"engaged by CyMPSA\s+: (\d+)", R13).group(1)),
          me=int(re.search(r"engaged by Meroterpenoid : (\d+)", R13).group(1)),
          n=int(re.search(r"measured by BOTH chemotypes: (\d+)", R13).group(1)),
          orv=float(re.search(r"odds ratio ([\d.]+)", R13).group(1)))
assert len(PCS) >= 6 and len(RHO) == 21 and len(T1) == 13, (len(PCS), len(RHO), len(T1))

HEAD = {}   # headline numbers per panel, captured where they are computed

# ---------------------------------------------------------------- figure
fig = plt.figure(figsize=(st.DOUBLE, st.DOUBLE * 1.02))
gs = fig.add_gridspec(4, 3, hspace=0.62, wspace=0.52,
                      left=0.070, right=0.980, top=0.955, bottom=0.052,
                      width_ratios=[1.0, 0.92, 0.92],
                      height_ratios=[1.0, 0.60, 0.60, 0.52])

# ---------------------------------------------------------------- a: ceiling per dataset
ax = fig.add_subplot(gs[0, 0:2]); st.panel_label(ax, "a", dx=-0.10, dy=1.07)
STYLE = {"NatChem": (st.ORANGE, "o"), "CyMPSA": (st.GREY, "s"), "Meroterpenoid": (st.GREEN, "^")}
for ds, d in CEIL.items():
    c, mk = STYLE[ds]
    xs = [r[0] for r in d["rows"]]; ys = [r[3] for r in d["rows"]]
    for r in d["rows"]:
        ax.plot([r[0], r[0]], [r[4], r[5]], color=c, lw=0.9, solid_capstyle="round")
    ax.plot(xs, ys, "-", color=c, lw=1.1)
    ax.plot(xs, ys, mk, color=c, ms=3.6, mew=0,
            label="%s — median %d compounds/site, %.1f%% of sites engaged\nat the deepest split"
                  % (ds, d["med"], d["rows"][-1][2]))
    # % engaged moved into the legend: CyMPSA and Meroterpenoid both end at 5 compounds/half
    # with ceilings 0.012 apart, so point labels landed on top of each other
ax.axhline(COMP[0], color=st.BLUE, lw=1.0, ls=(0, (3, 1.6)))
ax.text(1.78, COMP[0] + 0.012, "the model, %.3f" % COMP[0], va="bottom", ha="left",
        fontsize=5.4, color=st.BLUE)
ax.set_xscale("log"); ax.set_xticks([2, 3, 5, 10, 20])
ax.set_xticklabels([2, 3, 5, 10, 20], fontsize=6)
ax.set_xlabel("compounds in each half (the deeper the split, the less noisy the label)",
              fontsize=6.2)
ax.set_ylabel("reproducibility ceiling (AUC)", fontsize=6.2)
ax.set_xlim(1.7, 30); ax.set_ylim(0.55, 1.05)
ax.legend(frameon=False, fontsize=5.0, loc="upper right", handlelength=1.3, ncol=1,
          handletextpad=0.4, labelspacing=0.8, borderaxespad=0.2)
# All five numbers in this note were typed, and one of them was wrong: the engaged-rate range
# read 0.6–2.6%, but 2.6% is NatChem's rate at five compounds per half, not a shallow screen's.
# The two shallow screens run 0.6–2.3%. Every number now comes from CEIL.
_SH  = [r for _d in CEIL if _d != "NatChem" for r in CEIL[_d]["rows"]]
_SHC = [r[3] for r in _SH]; _SHE = [r[2] for r in _SH]
_NCR = CEIL["NatChem"]["rows"][-1]
ax.text(0.30, 0.55, "The two shallow screens reach %.2f–%.2f, but at engaged rates of\n"
        "%.1f–%.1f%%, where a few dozen positives make an AUC unstable.\nNatChem is the only screen "
        "deep enough to measure a ceiling\nat a realistic engaged rate, and its ceiling is still "
        "rising at\n%g compounds per half — so %.3f is a lower bound."
        % (min(_SHC), max(_SHC), min(_SHE), max(_SHE), _NCR[0], _NCR[3]),
        transform=ax.transAxes, fontsize=5.0, color=st.GREY, va="top", linespacing=1.55)
HEAD["a_model_auc"] = float(COMP[0])
HEAD["a_shallow_ceiling_min"] = round(float(min(_SHC)), 3)
HEAD["a_shallow_ceiling_max"] = round(float(max(_SHC)), 3)
HEAD["a_shallow_engaged_pct_min"] = round(float(min(_SHE)), 2)
HEAD["a_shallow_engaged_pct_max"] = round(float(max(_SHE)), 2)
HEAD["a_natchem_deepest_split"] = int(_NCR[0])
HEAD["a_natchem_ceiling_at_deepest"] = float(_NCR[3])
fd.panel("figS6", "a",
         [dict(dataset=_d, compounds_per_half=_r[0], sites=int(_r[1]), engaged_pct=_r[2],
               ceiling_auc=_r[3], ci_lo=_r[4], ci_hi=_r[5],
               median_compounds_per_site=CEIL[_d]["med"])
          for _d in CEIL for _r in CEIL[_d]["rows"]],
         note="one point per dataset per split depth; intervals over 50 random half-splits")

# ---------------------------------------------------------------- e: breadth, depth-adjusted
ax = fig.add_subplot(gs[1:3, 2]); st.panel_label(ax, "e", dx=-0.46, dy=1.03)
R = sorted(RHO, key=lambda r: r[2])
y = np.arange(len(R))
for i, (f_, raw, adj) in enumerate(R):
    ax.plot([raw, adj], [i, i], color=st.LIGHTGREY, lw=1.6, solid_capstyle="round", zorder=1)
    ax.plot([raw], [i], "o", color=st.LIGHTGREY, ms=3.2, mec=st.GREY, mew=0.5, zorder=2)
    ax.plot([adj], [i], "o", color=st.BLUE if adj > 0 else st.ORANGE, ms=3.6, mew=0, zorder=3)
ax.axvline(0, color=st.GREY, lw=0.8)
ax.set_yticks(y); ax.set_yticklabels([NICE(r[0]) for r in R], fontsize=5.0)
ax.set_xlabel("correlation with how many chemotypes\nengage the site (Spearman ρ)",
              fontsize=6.0, linespacing=1.4)
ax.set_ylim(-3.2, len(R) + 1.5)
ax.plot([], [], "o", color=st.LIGHTGREY, mec=st.GREY, mew=0.5, ms=3.2, label="raw")
ax.plot([], [], "o", color=st.BLUE, ms=3.6, label="holding testing depth constant")
ax.legend(frameon=False, fontsize=5.0, loc="lower right", handletextpad=0.3, borderaxespad=0.3)
# "three features" was typed; it is the number of features whose sign changes under adjustment.
_FLIP = sum(1 for _f, _raw, _adj in RHO if np.sign(_raw) != np.sign(_adj))
ax.text(0.02, 0.995, "a site tested against few chemotypes cannot show\nhigh breadth; holding "
        "depth constant strengthens most\nof the gradient and flips %d features outright" % _FLIP,
        transform=ax.transAxes, fontsize=4.9, color=st.GREY, va="top", linespacing=1.5)
for s in ("top", "right"): ax.spines[s].set_visible(False)
HEAD["e_n_breadth_features"] = int(len(RHO))
HEAD["e_n_features_sign_flipped"] = int(_FLIP)
fd.panel("figS6", "e",
         [dict(feature=_f, rho_raw=_raw, rho_depth_adjusted=_adj,
               sign_flipped=bool(np.sign(_raw) != np.sign(_adj))) for _f, _raw, _adj in R],
         note="one dumbbell per breadth feature, raw and holding testing depth constant")

# ---------------------------------------------------------------- b: single features
ax = fig.add_subplot(gs[0, 2]); st.panel_label(ax, "b", dx=-0.26, dy=1.10)
rows = [(NICE(f), v, lo, hi) for f, v, lo, hi in SINGLES][::-1]
y = np.arange(len(rows))
for yy, (lab, v, lo, hi) in zip(y, rows):
    ax.plot([lo, hi], [yy, yy], color=st.LIGHTGREY, lw=1.5, solid_capstyle="round")
    ax.plot([v], [yy], "o", color=st.GREY, ms=3.6, mew=0)
ax.plot([COMP[1], COMP[2]], [len(rows)] * 2, color=st.BLUE, lw=1.8, solid_capstyle="round")
ax.plot([COMP[0]], [len(rows)], "o", color=st.BLUE, ms=4.6, mew=0)
ax.set_yticks(list(y) + [len(rows)])
ax.set_yticklabels([r[0] for r in rows] + ["all %d together" % len(COEF)], fontsize=5.4)
ax.get_yticklabels()[-1].set_color(st.BLUE); ax.get_yticklabels()[-1].set_fontweight("bold")
ax.axvline(0.5, color=st.GREY, lw=0.8, ls=(0, (2, 1.6)))
ax.set_xlabel("AUC, held out by protein", fontsize=6.2)
ax.set_xlim(0.36, 0.70); ax.set_ylim(-2.4, len(rows) + 0.6)
ax.text(0.5, 0.035, "the composite gains only +%.3f over the best single\nfeature — the %d "
        "features are not independent" % (COMP[0] - SINGLES[0][1], len(COEF)),
        transform=ax.transAxes, ha="center", va="bottom", fontsize=5.0, color=st.GREY,
        linespacing=1.5)
for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
ax.tick_params(axis="y", length=0)
HEAD["b_composite_auc"] = float(COMP[0])
HEAD["b_composite_ci_lo"] = float(COMP[1])
HEAD["b_composite_ci_hi"] = float(COMP[2])
HEAD["b_best_single_feature"] = SINGLES[0][0]
HEAD["b_best_single_auc"] = float(SINGLES[0][1])
HEAD["b_composite_gain_over_best_single"] = round(float(COMP[0] - SINGLES[0][1]), 3)
HEAD["b_n_features_in_composite"] = int(len(COEF))
fd.panel("figS6", "b",
         [dict(feature=_f, auc=_v, ci_lo=_lo, ci_hi=_hi) for _f, _v, _lo, _hi in SINGLES]
         + [dict(feature="COMPOSITE (all %d)" % len(COEF), auc=COMP[0],
                 ci_lo=COMP[1], ci_hi=COMP[2])],
         note="AUC held out by protein; one row per single feature plus the composite")

# ---------------------------------------------------------------- d: train vs held-out
ax = fig.add_subplot(gs[1, 1]); st.panel_label(ax, "d", dx=-0.26, dy=1.10)
for i, (lab, v, c) in enumerate([("training fit", TRAIN, st.LIGHTGREY),
                                 ("held-out proteins", COMP, st.BLUE)]):
    ax.plot([v[1], v[2]], [i, i], color=c, lw=2.0, solid_capstyle="round")
    ax.plot([v[0]], [i], "o", color=c, ms=4.6, mew=0)
    ax.text(v[0], i + 0.22, "%.3f" % v[0], ha="center", fontsize=6.0, fontweight="bold", color=c)
ax.set_yticks([0, 1]); ax.set_yticklabels(["training fit", "held-out\nproteins"], fontsize=5.6,
                                          linespacing=1.35)
ax.set_xlabel("AUC", fontsize=6.2)
ax.set_xlim(0.575, 0.66); ax.set_ylim(-0.7, 1.75)
ax.text(0.5, 0.94, "optimism = %+.3f AUC" % (TRAIN[0] - COMP[0]), transform=ax.transAxes,
        ha="center", fontsize=6.0, color=st.GREEN, fontweight="bold")
ax.text(0.5, -0.42, "the model is not overfit — it is fitting all\nthere is to fit, and that is "
        "the point", transform=ax.transAxes, ha="center", va="top", fontsize=5.2,
        color=st.GREY, linespacing=1.5)
for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
ax.tick_params(axis="y", length=0)
HEAD["d_train_auc"] = float(TRAIN[0])
HEAD["d_heldout_auc"] = float(COMP[0])
HEAD["d_optimism"] = round(float(TRAIN[0] - COMP[0]), 3)
fd.panel("figS6", "d",
         [dict(fit="training", auc=TRAIN[0], ci_lo=TRAIN[1], ci_hi=TRAIN[2]),
          dict(fit="held-out proteins", auc=COMP[0], ci_lo=COMP[1], ci_hi=COMP[2])],
         note="intervals are the 2.5th to 97.5th percentile over the protein-held-out splits")

# ---------------------------------------------------------------- f: two chemotypes agree
ax = fig.add_subplot(gs[2, 1]); st.panel_label(ax, "f", dx=-0.30, dy=1.07)
xs = [t[1] for t in T1]; ys = [t[2] for t in T1]
same = [t[3] == "yes" for t in T1]
# "7 of 13", "6 of 13" and "0 of 13" were all typed. They are sums over T1.
_NSAME = int(sum(same)); _NDIFF = int(sum(1 for t in T1 if t[4] < 0.05))
ax.axhline(0, color=st.LIGHTGREY, lw=0.8); ax.axvline(0, color=st.LIGHTGREY, lw=0.8)
lim = max(max(map(abs, xs)), max(map(abs, ys))) * 1.30
ax.plot([-lim, lim], [-lim, lim], color=st.LIGHTGREY, lw=0.8, ls=(0, (2, 2)), zorder=0)
ax.scatter([x for x, s_ in zip(xs, same) if s_], [y_ for y_, s_ in zip(ys, same) if s_],
           s=16, color=st.BLUE, linewidths=0, zorder=3,
           label="same direction (%d of %d)" % (_NSAME, len(T1)))
ax.scatter([x for x, s_ in zip(xs, same) if not s_], [y_ for y_, s_ in zip(ys, same) if not s_],
           s=16, facecolor=st.PAPER, edgecolor=st.ORANGE, linewidths=0.9, zorder=3,
           label="opposite direction (%d of %d)" % (len(T1) - _NSAME, len(T1)))
# Thirteen labels do not fit in this cell and no amount of nudging makes them: adjustText was
# still returning overlaps at 3.9 pt. Name only the three features that disagree MOST between
# the two chemotypes — the ones a reader would actually look for — and leave the rest as points.
# Every per-feature value is in panel c and in step13_report.txt.
_gap = sorted(range(len(T1)), key=lambda k: -abs(xs[k] - ys[k]))[:3]
txt = [ax.text(xs[k], ys[k], NICE(T1[k][0]), fontsize=4.4, color=st.GREY) for k in _gap]
adjust_text(txt, ax=ax, arrowprops=dict(arrowstyle="-", color=st.GREY, lw=0.35),
            expand=(1.5, 1.8))
ax.set_xlim(-lim * 1.15, lim * 1.15); ax.set_ylim(-lim * 1.35, lim * 1.35)
ax.set_xlabel("effect in CyMPSA (d)", fontsize=6.2)
ax.set_ylabel("effect in meroterpenoids (d)", fontsize=5.8, labelpad=1.5)
ax.legend(frameon=False, fontsize=4.7, loc="upper left", markerscale=0.9,
          handletextpad=0.3, labelspacing=0.35, borderaxespad=0.3)
ax.text(0.5, -0.30, "%d of %d features differ significantly between the two chemotypes, and the "
        "same sites\nare engaged by both (%d of %s shared, %.1f expected, OR %.2f) — so one "
        "general\ncomposite is defensible. The three most discordant features are named."
        % (_NDIFF, len(T1), T2["both"], f"{T2['n']:,}", T2["exp"], T2["orv"]),
        transform=ax.transAxes, ha="center", va="top", fontsize=5.0, color=st.GREY,
        linespacing=1.5)
HEAD["f_n_chemotype_features"] = int(len(T1))
HEAD["f_n_same_direction"] = _NSAME
HEAD["f_n_opposite_direction"] = int(len(T1) - _NSAME)
HEAD["f_n_differing_significantly"] = _NDIFF
HEAD["f_n_sites_engaged_by_both"] = int(T2["both"])
HEAD["f_n_sites_measured_by_both"] = int(T2["n"])
HEAD["f_sites_both_expected_by_chance"] = float(T2["exp"])
HEAD["f_odds_ratio"] = float(T2["orv"])
fd.panel("figS6", "f",
         [dict(feature=_t[0], effect_cympsa=_t[1], effect_meroterpenoid=_t[2],
               same_direction=(_t[3] == "yes"), q=_t[4]) for _t in T1],
         note="one point per feature; the two chemotypes' effects against each other")


# ---------------------------------------------------------------- c: every coefficient
# An AUC below 0.5 is an INVERTED feature, not a weak one: distance to centroid reads 0.402,
# which step 10 reports as its flipped value, 0.598. Showing both makes the sign structure
# visible — the two largest weights belong to features barely better than chance on their own.
ax = fig.add_subplot(gs[1:3, 0]); st.panel_label(ax, "c", dx=-0.44, dy=1.03)
CF = COEF.copy()
CF["flip"] = np.maximum(CF.single_auc, 1 - CF.single_auc)
CF = CF.sort_values("mean_coef")
_y = np.arange(len(CF))
ax.barh(_y, CF.mean_coef, height=0.66,
        color=[st.BLUE if v > 0 else st.ORANGE for v in CF.mean_coef])
ax.axvline(0, color=st.GREY, lw=0.8)
ax.set_yticks(_y); ax.set_yticklabels([NICE(f_) for f_ in CF.feature], fontsize=5.0)
ax.set_xlabel("weight in the composite\n(mean over %d splits)" % N_SPLITS,
              fontsize=6.0, linespacing=1.4)
ax.set_xlim(-0.40, 0.30); ax.set_ylim(-3.6, len(CF) - 0.35)
_axt = ax.twiny()
_axt.plot(CF.flip, _y, "o", color=st.GREEN, ms=3.2, mew=0, zorder=5)
_axt.axvline(0.5, color=st.LIGHTGREY, lw=0.8, ls=(0, (2, 2)))
_axt.set_xlim(0.36, 0.72); _axt.set_ylim(ax.get_ylim())
_axt.set_xlabel("AUC on its own", fontsize=6.0, color=st.GREEN)
_axt.tick_params(axis="x", colors=st.GREEN, labelsize=5.2)
ax.spines["right"].set_visible(False)
# The feature and both AUCs were typed. They are the row of CF with the lowest single AUC.
_INV = CF.loc[CF.single_auc.idxmin()]
ax.text(0.985, 0.02, "an AUC below 0.5 is an inverted feature,\nnot a weak one: %s reads\n"
        "%.3f, which is %.3f with the sign flipped"
        % (NICE(_INV.feature), _INV.single_auc, _INV.flip),
        transform=ax.transAxes, ha="right", va="bottom", fontsize=4.8, color=st.GREY,
        linespacing=1.5)
HEAD["c_n_coefficients"] = int(len(CF))
HEAD["c_n_splits"] = N_SPLITS
HEAD["a_n_half_splits"] = N_HALFSPLITS
HEAD["c_n_principal_components"] = int(len(PCS))
HEAD["c_inverted_feature"] = str(_INV.feature)
HEAD["c_inverted_feature_auc"] = round(float(_INV.single_auc), 3)
HEAD["c_inverted_feature_auc_flipped"] = round(float(_INV.flip), 3)
fd.panel("figS6", "c",
         [dict(feature=_r.feature, mean_coef=float(_r.mean_coef),
               single_auc=float(_r.single_auc), single_auc_flipped=float(_r.flip))
          for _r in CF.itertuples()],
         note="one bar per coefficient, with its own AUC on the twin axis")

# ---------------------------------------------------------------- g: across datasets — a LIMIT
# This is the only cross-dataset evidence in the paper and it is the weakest thing in it. Panel g
# is within ONE dataset; step 15 is what happens when the question is asked across two, and the
# honest answer is that the group carrying the claim is 25 sites in 21 dimensions and two of the
# six pairs do not clear their null. It sits here as a stated limitation, not as reinforcement.
ax = fig.add_subplot(gs[3, 0:3]); st.panel_label(ax, "g", dx=-0.30, dy=1.05)
S = S15.iloc[::-1].reset_index(drop=True)
_y = np.arange(len(S))
LBL = {"SHARED": "engaged by both", "CYMPSA_ONLY": "CyMPSA only",
       "MERO_ONLY": "meroterpenoid only", "NEITHER": "neither"}
for i_, r_ in S.iterrows():
    sig = r_.p < 0.05
    ax.plot([r_.null_mean, r_.distance], [i_, i_], color=st.LIGHTGREY, lw=1.6,
            solid_capstyle="round", zorder=1)
    ax.plot([r_.null_mean], [i_], "o", color=st.LIGHTGREY, mec=st.GREY, mew=0.6, ms=3.0, zorder=2)
    ax.plot([r_.distance], [i_], "o", color=st.BLUE if sig else st.ORANGE, ms=4.0, mew=0, zorder=3)
    ax.text(r_.distance + 0.09, i_, "p %s   n %d vs %s"
            % (("= %.3f" % r_.p).rstrip("0").rstrip(".") if r_.p > 0.001 else "\u2264 0.001",
               r_.n_a, f"{int(r_.n_b):,}"),
            va="center", fontsize=4.5, color=st.BLUE if sig else st.ORANGE)
ax.set_yticks(_y)
ax.set_yticklabels(["%s vs %s" % (LBL[p.split(" vs ")[0]], LBL[p.split(" vs ")[1]])
                    for p in S.pair], fontsize=5.0)
# 21 dimensions, 2,000 shuffles, 25 sites and "two of six" were all typed into this panel.
# The first two are stated by step 15; the last two are columns of the table being plotted.
_NDIM   = int(re.search(r"full (\d+)-dimensional", R15).group(1))
_NSHUF  = int(re.search(r"\((\d+) shuffles", R15).group(1))
_NSHARE = int(S.loc[S.pair.str.startswith("SHARED"), "n_a"].min())
_NFAIL  = int((S.p >= 0.05).sum())
ax.set_xlabel("distance between group centres in %d-dimensional structural space" % _NDIM,
              fontsize=6.0)
ax.set_xlim(0, 4.6); ax.set_ylim(-0.85, len(S) - 0.25)
ax.text(0.985, 0.03, "grey = null over %s shuffles, teal = clears it, rust = does not.\n"
        "THE WEAKEST EVIDENCE IN THE PAPER, AND FIGURE 5e DEPENDS ON IT: the group\n"
        "carrying the cross-dataset claim is %d sites in %d dimensions, and %d of %d pairs fall\n"
        "short of their null. Figure 5e repeats this test along the breadth gradient WITHIN one\n"
        "dataset; this is what happens when it is asked across two."
        % (f"{_NSHUF:,}", _NSHARE, _NDIM, _NFAIL, len(S)),
        transform=ax.transAxes, ha="right", fontsize=4.8, color=st.GREY, va="bottom",
        linespacing=1.5)
HEAD["g_n_pairs"] = int(len(S))
HEAD["g_n_pairs_short_of_null"] = _NFAIL
HEAD["g_n_shared_sites"] = _NSHARE
HEAD["g_n_dimensions"] = _NDIM
HEAD["g_n_shuffles"] = _NSHUF
fd.panel("figS6", "g",
         [dict(pair=_r.pair, distance=float(_r.distance), null_mean=float(_r.null_mean),
               p=float(_r.p), n_a=int(_r.n_a), n_b=int(_r.n_b),
               clears_null=bool(_r.p < 0.05)) for _r in S.itertuples()],
         note="one row per group pair; distance against its permutation null")
for _sp in ("top", "right", "left"): ax.spines[_sp].set_visible(False)
ax.tick_params(axis="y", length=0)

st.save(fig, "figS6")
fd.head("figS6", HEAD)
for _e in (".png", ".pdf"):
    fig.savefig("../results/figS6_v1" + _e, dpi=400, bbox_inches="tight")
print("ceiling datasets:", {k: len(v["rows"]) for k, v in CEIL.items()})
print("PCs:", len(PCS), "| breadth features:", len(RHO), "| chemotype features:", len(T1))
print("composite %.3f vs best single %.3f (+%.3f) | train %.3f, optimism %+.3f"
      % (COMP[0], SINGLES[0][1], COMP[0]-SINGLES[0][1], TRAIN[0], TRAIN[0]-COMP[0]))
