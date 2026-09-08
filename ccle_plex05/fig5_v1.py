#!/usr/bin/env python3
"""
Figure 5 — Ligandability is reproducible, and structurally unpredictable.

THE CENTRAL COMPARISON (panel a)
A structural model of ligandability reaches AUC 0.617 (95% CI 0.594-0.640), cross-validated by
PROTEIN so no test protein is ever seen in training. On its own that number says nothing: an AUC
is only interpretable against how reproducible the labels themselves are. Step 29 measures that
directly — split the compounds tested at each site into halves, and ask how well engagement in
one half predicts engagement in the other. That is the best any site-level classifier could do.

WHICH CEILING IS THE FAIR ONE
Three datasets give three ceilings, and they are not interchangeable:
  NatChem        median 41 compounds/site, 17.5% of sites engaged -> deepest split 20/half
  CyMPSA         median 10 compounds/site,  2.9% engaged          -> deepest split  5/half
  Meroterpenoid  median 10 compounds/site,  1.5% engaged          -> deepest split  5/half
Only NatChem is screened deeply enough to estimate a ceiling at a realistic engaged rate, and its
ceiling RISES with split depth (0.628 at 2/half to 0.723 at 20/half) because a deeper half is a
less noisy label. The shallow screens return higher numbers (0.88-0.94) at engaged rates of 1-3%,
where a handful of positives makes AUC both unstable and easy to inflate. NatChem at 20/half,
0.723 [0.706, 0.744], is therefore the binding comparator, and step 29 is explicit that all of
these are LOWER bounds. figS6 shows all three so the choice can be checked rather than trusted.

A CONTRADICTION IN THE REPO, RESOLVED HERE
step10_report.txt still carries an interpretation guide from an earlier pass: "the Part 1 ceiling
analysis put the maximum achievable AUC at 0.55-0.65 ... an AUC in that band means the model is at
the ceiling, not that it is weak." Read that way, 0.617 sits AT the ceiling and the conclusion
inverts — there would be nothing left to measure. Step 29 supersedes it with an empirical,
per-dataset estimate, and the gap is real. The stale guide should be struck from step10_report.
"""
import os, re, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import abbasov_style as st
import figure_data as fd
st.apply()

def find(n):
    for d in (".", "..", "ligand_rebuild", "../ligand_rebuild"):
        p = os.path.join(d, n)
        if os.path.exists(p):
            return p
    raise FileNotFoundError(n)          # a hard failure, never a silent blank panel

R10 = open(find("step10_report.txt")).read()
R29 = open(find("step29_report.txt")).read()
R29B = open(find("step29b_report.txt")).read()
R16 = open(find("step16_report.txt")).read()
R13 = open(find("step13_report.txt")).read()
ROC_S = pd.read_csv(find("step10_roc_splits.csv"))
ROC_P = pd.read_csv(find("step10_roc_pooled.csv"))
PCS = [(int(m.group(1)), float(m.group(2)), float(m.group(3)))
       for m in re.finditer(r"PC(\d)\s+([\d.]+)%\s+cumulative\s+([\d.]+)%", R10)]

# ---------------------------------------------------------------- parse, then assert
def auc_row(txt, name):
    m = re.search(re.escape(name) + r"\s+([\d.]+)\s+\[([\d.]+), ([\d.]+)\]", txt)
    return tuple(float(x) for x in m.groups())

COMP = auc_row(R10, "COMPOSITE (test)")
TRAIN = auc_row(R10, "composite (train)")
BEST = auc_row(R10, "dist_centroid")
SINGLES = [(f, auc_row(R10, f)[0]) for f in
           ("dist_centroid", "plddt_site", "propka_pKa", "hydrophobic_8A", "rsa_tien")]

# ceiling curves: dataset -> [(compounds per half, n sites, % engaged, auc, lo, hi)]
CEIL, cur = {}, None
for line in R29.splitlines():
    m = re.match(r"^(\w+) — (\d+) sites \| compounds per site: median (\d+)", line)
    if m:
        cur = m.group(1); CEIL[cur] = dict(n=int(m.group(2)), med=int(m.group(3)), rows=[])
    m = re.match(r"^\s*(\d+) compounds/half \|\s*(\d+) sites \|\s*([\d.]+)% engaged \| "
                 r"ceiling AUC ([\d.]+) \[([\d.]+), ([\d.]+)\]", line)
    if m and cur:
        CEIL[cur]["rows"].append(tuple(float(x) for x in m.groups()))

# ceiling by structural-confidence band
BANDS, cur = {}, None
for line in R29B.splitlines():
    if re.match(r"^(NatChem|CyMPSA|Meroterpenoid)\s*$", line.strip()):
        cur = line.strip(); BANDS[cur] = []
    m = re.match(r"^\s*(< 50\s+disordered|50-70 low|70-90 confident|>= 90 very high|no structure)"
                 r"\s+n=\s*(\d+) \| deepest split\s+(\d+)/half \((\d+) sites\) \| "
                 r"ceiling ([\d.]+) \[([\d.]+), ([\d.]+)\] \|\s*([\d.]+)% engaged", line)
    if m and cur:
        BANDS[cur].append((m.group(1).strip(), int(m.group(2)), int(m.group(3)), int(m.group(4)),
                           float(m.group(5)), float(m.group(6)), float(m.group(7)),
                           float(m.group(8))))

# breadth: multivariate distance from the 0-chemotype group
BRE = []
for line in R16.splitlines():
    m = re.match(r"^\s*(\d|\d-\d|>=\d)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)", line)
    if m and "chemotypes TESTED" not in line:
        BRE.append((m.group(1), int(m.group(2)), float(m.group(3)), float(m.group(4))))
BRE = [b for b in BRE if b[0] in ("1", "2-3", "4-5", ">=6")]

# enantiomer pairs
ENA = []
for m in re.finditer(r"\((.)\)-([\w ]+?)\s+vs\s+\((.)\)-([\w ]+?)\s+n=\s*(\d+)\n"
                     r"\s+both\s+(\d+) \(expected ([\d.]+)\) \| only\(\+\)\s+(\d+) \| "
                     r"only\(-\)\s+(\d+) \| r = ([\d.]+)", R13):
    ENA.append(dict(name=m.group(2).strip(), n=int(m.group(5)), both=int(m.group(6)),
                    exp=float(m.group(7)), plus=int(m.group(8)), minus=int(m.group(9)),
                    r=float(m.group(10))))

NAT = CEIL["NatChem"]["rows"]
NAT_DEEP = NAT[-1]                      # 20 compounds/half
assert abs(COMP[0] - 0.617) < 1e-9 and abs(NAT_DEEP[3] - 0.723) < 1e-9, (COMP, NAT_DEEP)
assert len(BRE) == 4 and len(ENA) == 2, (BRE, ENA)
GAP = NAT_DEEP[3] - COMP[0]

HEAD = {
    "a_composite_auc": COMP[0], "a_composite_lo": COMP[1], "a_composite_hi": COMP[2],
    "a_train_auc": TRAIN[0], "a_optimism": round(TRAIN[0] - COMP[0], 3),
    "a_ceiling_auc": NAT_DEEP[3], "a_ceiling_lo": NAT_DEEP[4], "a_ceiling_hi": NAT_DEEP[5],
    "a_gap": round(GAP, 3),
    "a_ceiling_shallowest": NAT[0][3], "a_compounds_per_half_deepest": int(NAT_DEEP[0]),
    "a_best_single_auc": BEST[0],
    "a_composite_over_best_single": round(COMP[0] - BEST[0], 3),
}
fd.panel("fig5", "a", [dict(row="composite (test)", auc=COMP[0], lo=COMP[1], hi=COMP[2]),
                       dict(row="composite (train)", auc=TRAIN[0], lo=TRAIN[1], hi=TRAIN[2]),
                       dict(row="best single feature (dist_centroid)",
                            auc=BEST[0], lo=BEST[1], hi=BEST[2])] +
         [dict(row="NatChem ceiling, %d compounds/half" % int(r[0]), auc=r[3], lo=r[4], hi=r[5])
          for r in NAT],
         note="AUC and 2.5-97.5 percentile interval; ceiling series is the split-half ladder")

fd.panel("fig5", "b", [dict(dataset=ds, band=r[0], n_sites=r[1], deepest_split_per_half=r[2],
                            n_sites_at_split=r[3], ceiling_auc=r[4], lo=r[5], hi=r[6],
                            pct_engaged=r[7])
                       for ds, v in BANDS.items() for r in v],
         note="reproducibility ceiling within AlphaFold-confidence bands")
_nb = {r[0]: r for r in BANDS["NatChem"]}
HEAD["b_natchem_ceiling_min"] = round(min(r[4] for r in BANDS["NatChem"]), 3)
HEAD["b_natchem_ceiling_max"] = round(max(r[4] for r in BANDS["NatChem"]), 3)
for _k, _lab in ((">= 90 very high", "plddt90"),):
    if _k in _nb:
        HEAD["b_natchem_ceiling_%s" % _lab] = round(_nb[_k][4], 3)
        HEAD["b_natchem_n_sites_%s" % _lab] = int(_nb[_k][1])
for _ds in BANDS:
    HEAD["b_n_bands_%s" % _ds] = len(BANDS[_ds])

fd.panel("fig5", "d", [dict(component=p_[0], variance_pct=p_[1], cumulative_pct=p_[2])
                       for p_ in PCS], note="principal components of the 14 features")
HEAD["d_pc1_variance_pct"] = PCS[0][1]
HEAD["d_n_components_to_80pct"] = int(next(p_[0] for p_ in PCS if p_[2] >= 80))
HEAD["d_n_features"] = 14

fd.panel("fig5", "e", [dict(chemotypes=b[0], n_sites=b[1], observed_distance=b[2],
                            null_distance=b[3]) for b in BRE],
         note="multivariate distance from a never-engaged site; null over 2,000 label shuffles")
for _b in BRE:
    _t = _b[0].replace(">=", "ge").replace("-", "_")
    HEAD["e_observed_%s" % _t] = _b[2]
    HEAD["e_null_%s" % _t] = _b[3]
    HEAD["e_n_sites_%s" % _t] = _b[1]

fd.panel("fig5", "f", [dict(pair=e["name"], n_sites_measured=e["n"], engaged_by_both=e["both"],
                            expected_by_chance=e["exp"], only_plus=e["plus"],
                            only_minus=e["minus"], correlation=e["r"]) for e in ENA],
         note="enantiomer pairs with enough shared coverage to compare")
for _e in ENA:
    _t = _e["name"].replace(" ", "_")
    HEAD["f_%s_both" % _t] = _e["both"]
    HEAD["f_%s_one_hand" % _t] = _e["plus"] + _e["minus"]
    HEAD["f_%s_expected" % _t] = _e["exp"]
    HEAD["f_%s_n_sites" % _t] = _e["n"]
HEAD["f_n_pairs"] = len(ENA)

# ---------------------------------------------------------------- figure
fig = plt.figure(figsize=(st.DOUBLE, st.DOUBLE * 0.63))
gs = fig.add_gridspec(2, 3, hspace=0.58, wspace=0.42,
                      left=0.078, right=0.982, top=0.930, bottom=0.090,
                      width_ratios=[1.22, 1.0, 1.0])

# ---------------------------------------------------------------- a: the central comparison
ax = fig.add_subplot(gs[0, 0]); st.panel_label(ax, "a", dx=-0.20, dy=1.08)
# The two comparator rows are set in the darker neutral, not the light one: their value labels
# sit partly over the light-grey gap band, where light grey on light grey is unreadable.
rows = [("chance", 0.50, None, None, st.GREY),
        ("best single feature\n(distance to centroid)", BEST[0], BEST[1], BEST[2], st.GREY),
        ("all 14 features together\ncross-validated by protein", COMP[0], COMP[1], COMP[2], st.BLUE),
        ("how reproducible the labels are\n(NatChem, by screen depth)",
         NAT_DEEP[3], NAT_DEEP[4], NAT_DEEP[5], st.ORANGE)]
y = np.arange(len(rows))[::-1]
ax.axvspan(COMP[0], NAT_DEEP[3], color=st.LIGHTGREY, zorder=0)
for yy, (lab, v, lo, hi, c) in zip(y, rows):
    if lo is not None:
        ax.plot([lo, hi], [yy, yy], color=c, lw=1.4, solid_capstyle="round", zorder=3)
    ax.plot([v], [yy], "o", color=c, ms=5.0, mew=0, zorder=4)
    ax.text(v, yy + 0.26, "%.3f" % v, ha="center", fontsize=5.8, fontweight="bold", color=c)
# The ceiling is not one number: it rises with how deeply each half is screened, from 0.628 at
# 2 compounds per half to 0.723 at 20, and it is still climbing at the deepest split available.
# Putting the whole series on this axis lets the comparator justify itself here rather than only
# in figS6 a, and shows 0.723 for what it is — a lower bound.
_cy = y[-1]
# The rail sits inside the light-grey gap band, so it has to be the darker neutral or it
# disappears into the band.
ax.plot([r[3] for r in NAT], [_cy] * len(NAT), "-", color=st.GREY, lw=1.0, zorder=2)
for r in NAT[:-1]:
    ax.plot([r[3]], [_cy], "o", color=st.ORANGE, ms=2.2, mew=0, zorder=3)
ax.text(NAT[0][3] - 0.008, _cy - 0.30, "2", ha="center", fontsize=4.8, color=st.ORANGE)
ax.text(NAT_DEEP[3], _cy - 0.30, "20 compounds per half", ha="center", fontsize=4.8,
        color=st.ORANGE)
ax.text((NAT[0][3] + NAT_DEEP[3]) / 2, _cy - 0.55, "still rising at the deepest split —\n"
        "0.723 is a lower bound", ha="center", va="top", fontsize=4.9, color=st.ORANGE,
        linespacing=1.4)
ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=5.6, linespacing=1.35)
ax.set_xlabel("how well a site can be sorted into engaged or not (AUC)", fontsize=6.2)
# dist_centroid's CI runs from 0.375 — below chance — because a single feature is unstable
# across protein-held-out splits. Truncating it would flatter the comparison, so the axis
# opens far enough to show it.
# The two-line gap statement below the arrow needs its own room; at -2.00 its second line
# ran along the x-axis spine.
ax.set_xlim(0.355, 0.79); ax.set_ylim(-2.62, len(rows) - 0.30)
ax.annotate("", xy=(NAT_DEEP[3], -1.30), xytext=(COMP[0], -1.30),
            arrowprops=dict(arrowstyle="<->", color=st.GREY, lw=0.9))
ax.text((COMP[0] + NAT_DEEP[3]) / 2, -1.44,
        "%.3f of reproducible signal the structure does not reach\n"
        "(and the gap is not overfitting: optimism %+.3f AUC)"
        % (GAP, TRAIN[0] - COMP[0]),
        ha="center", va="top", fontsize=5.6, color=st.GREY, fontweight="bold",
        linespacing=1.4)
for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
ax.tick_params(axis="y", length=0)

# ---------------------------------------------------------------- b: ceiling vs structural order
ax = fig.add_subplot(gs[0, 1]); st.panel_label(ax, "b", dx=-0.16, dy=1.08)
ORDER = ["< 50  disordered", "50-70 low", "70-90 confident", ">= 90 very high", "no structure"]
NICE = ["< 50\ndisordered", "50–70\nlow", "70–90\nconfident", "≥ 90\nvery high", "no\nstructure"]
# (colour, marker, line style). The third slot used to be an opacity.
STYLE = {"NatChem": (st.ORANGE, "o", "-"), "CyMPSA": (st.GREY, "s", (0, (3.0, 1.4))),
         "Meroterpenoid": (st.GREEN, "^", (0, (1.2, 1.2)))}
for ds, rowset in BANDS.items():
    c, mk, a = STYLE[ds]
    xs = [ORDER.index(r[0]) for r in rowset if r[0] in ORDER]
    ys = [r[4] for r in rowset if r[0] in ORDER]
    los = [r[5] for r in rowset if r[0] in ORDER]
    his = [r[6] for r in rowset if r[0] in ORDER]
    # a was an opacity per dataset; it is now a dash pattern, so the three series separate
    # without any of them being drawn in a colour that is not in the palette.
    for x_, lo_, hi_ in zip(xs, los, his):
        ax.plot([x_, x_], [lo_, hi_], color=c, lw=0.9, solid_capstyle="round")
    ax.plot(xs, ys, color=c, lw=1.0, ls=a)
    ax.plot(xs, ys, mk, color=c, ms=3.4, mew=0, label=ds)
ax.axhline(COMP[0], color=st.BLUE, lw=1.0, ls=(0, (3, 1.6)))
ax.text(-0.42, COMP[0] + 0.012, "the model, 0.617", va="bottom", ha="left", fontsize=5.4,
        color=st.BLUE)
ax.set_xticks(range(5)); ax.set_xticklabels(NICE, fontsize=5.2, linespacing=1.3)
ax.set_xlabel("AlphaFold confidence at the lysine", fontsize=6.2)
ax.set_ylabel("reproducibility ceiling (AUC)", fontsize=6.2)
ax.set_xlim(-0.5, 4.62); ax.set_ylim(0.36, 1.05)
ax.legend(frameon=False, fontsize=5.2, loc="center left", bbox_to_anchor=(0.005, 0.60),
          handlelength=1.4, handletextpad=0.4, borderaxespad=0.0)
ax.text(0.5, 0.02, "the ceiling does not depend on structural order, so it cannot be dismissed\n"
        "as an artefact of disordered sites", transform=ax.transAxes, ha="center",
        fontsize=5.2, color=st.GREY, va="bottom", linespacing=1.45)

# ---------------------------------------------------------------- e: breadth across chemotypes
ax = fig.add_subplot(gs[1, 1]); st.panel_label(ax, "e", dx=-0.22, dy=1.08)
x = np.arange(len(BRE)); w = 0.36
obs = [b[2] for b in BRE]; nul = [b[3] for b in BRE]
ax.bar(x - w/2, obs, w, color=st.BLUE, label="observed")
ax.bar(x + w/2, nul, w, color=st.LIGHTGREY, label="2,000 label shuffles")
for i, b in enumerate(BRE):
    ax.text(i - w/2, b[2] + 0.06, "%.2f" % b[2], ha="center", fontsize=5.4,
            fontweight="bold", color=st.BLUE)
    ax.text(i, -0.20, "n=%d" % b[1], ha="center", fontsize=5.0, color=st.GREY)
ax.set_xticks(x); ax.set_xticklabels([b[0].replace(">=", "≥") for b in BRE], fontsize=6)
ax.set_xlabel("chemotypes that engage the site (of 34 tested)", fontsize=6.2)
ax.set_ylabel("how far the site's environment sits\nfrom a never-engaged site", fontsize=6.2,
              linespacing=1.4)
ax.set_ylim(-0.32, max(obs) * 1.52)
ax.legend(frameon=False, fontsize=5.4, loc="upper left", handlelength=1.1,
          bbox_to_anchor=(0.0, 0.80))
ax.text(0.98, 0.99, "what marks a site is not WHICH chemistry\nengages it but HOW MANY can — "
        "every bin\np < 0.0005 against 2,000 label shuffles",
        transform=ax.transAxes, ha="right", va="top", fontsize=5.4, color=st.GREY,
        linespacing=1.5)

# ---------------------------------------------------------------- f: enantiomers
ax = fig.add_subplot(gs[1, 2]); st.panel_label(ax, "f", dx=-0.18, dy=1.08)
lab, both, one, exp = [], [], [], []
for e in sorted(ENA, key=lambda z: -z["both"]):
    lab.append("%s\n(%s sites measured)" % (e["name"], f"{e['n']:,}"))
    both.append(e["both"]); one.append(e["plus"] + e["minus"]); exp.append(e["exp"])
x = np.arange(len(lab)); w = 0.34
ax.bar(x - w/2, both, w, color=st.BLUE, label="engaged by both hands")
ax.bar(x + w/2, one, w, color=st.ORANGE, label="engaged by one hand only")
for i in range(len(lab)):
    ax.text(i - w/2, both[i] + 0.45, "%d" % both[i], ha="center", fontsize=6.0,
            fontweight="bold", color=st.BLUE)
    ax.text(i + w/2, one[i] + 0.45, "%d" % one[i], ha="center", fontsize=6.0,
            fontweight="bold", color=st.ORANGE)
    # A rule at 0.1 sites is invisible next to a bar of 19, so the expectation is stated
    # rather than drawn — the whole point is that it is two orders of magnitude smaller.
    # inside the bar only where the bar is tall enough to hold three lines
    # Set inside the bar this ran wider than the bar itself, and the overflow was white on
    # white paper. Both expectations now sit above their bar, above the count.
    ax.text(i - w/2 + 0.05, both[i] + 3.1, "vs %.1f expected\nby chance" % exp[i], ha="center",
            va="bottom", fontsize=5.0, color=st.BLUE, linespacing=1.4)
ax.set_xticks(x); ax.set_xticklabels(lab, fontsize=5.6, linespacing=1.35)
ax.set_ylabel("sites engaged at ≥60% occupancy", fontsize=6.2)
ax.set_ylim(0, max(both) * 1.95)
ax.legend(frameon=False, fontsize=5.4, loc="upper right", handlelength=1.1,
          bbox_to_anchor=(1.0, 0.66), borderaxespad=0.0)
ax.text(0.02, 0.98, "same chemistry, opposite three-dimensional shape.\n"
        "The pair overwhelmingly engages the SAME sites, so\n"
        "engagement is a property of the lysine — but the sites\n"
        "engaged by one hand only cannot be intrinsic reactivity.",
        transform=ax.transAxes, fontsize=5.2, color=st.GREY, va="top", linespacing=1.5)


# ---------------------------------------------------------------- c: is it weak everywhere?
# The AUC says the model is weak on average. A reviewer's next hope is that it is sharp at high
# specificity — useful for picking a shortlist even if useless overall. It is not: at 5% false
# positives it recovers 11% of engaged sites, at 10% it recovers 19%. Predictions regenerated
# 2026-09-02 (step10_composite.py saves none); the re-run reproduces 0.617 exactly.
ax = fig.add_subplot(gs[0, 2]); st.panel_label(ax, "c", dx=-0.22, dy=1.08)
ax.plot([0, 1], [0, 1], color=st.LIGHTGREY, lw=0.9, ls=(0, (3, 2)), zorder=1)
ax.fill_between(ROC_S.fpr, ROC_S.lo_tpr, ROC_S.hi_tpr, color=st.LIGHTGREY, zorder=2,
                label="95% of 20 held-out splits")
ax.plot(ROC_P.fpr, ROC_P.tpr, color=st.BLUE, lw=1.6, zorder=4, label="every site predicted once")
for _fp in (0.05, 0.10):
    _tp = float(np.interp(_fp, ROC_P.fpr, ROC_P.tpr))
    ax.plot([_fp], [_tp], "o", color=st.ORANGE, ms=4.0, mew=0, zorder=6)
    # The two operating points are only 0.08 apart in recovery; put one label under its
    # point and one over so they read as two numbers rather than one stack.
    _lo = _fp < 0.075
    ax.text(_fp - 0.030, _tp + (-0.055 if _lo else 0.055), "%.0f%%" % (100*_tp), fontsize=5.4,
            color=st.ORANGE, ha="right", va="top" if _lo else "bottom", fontweight="bold")
fd.panel("fig5", "c", [dict(false_positive_rate=float(_f),
                            engaged_recovered=float(np.interp(_f, ROC_P.fpr, ROC_P.tpr)))
                       for _f in np.linspace(0, 1, 201)],
         note="pooled ROC of the composite, every site predicted once")
for _fpv in (0.05, 0.10):
    HEAD["c_recovered_at_%dpct_fpr" % int(100 * _fpv)] = round(
        100 * float(np.interp(_fpv, ROC_P.fpr, ROC_P.tpr)), 0)
ax.set_xlabel("false positives (of sites never engaged)", fontsize=6.2)
ax.set_ylabel("engaged sites recovered", fontsize=6.2)
ax.set_xlim(0, 1); ax.set_ylim(0, 1.02)
ax.legend(frameon=False, fontsize=5.2, loc="lower right", handlelength=1.2,
          handletextpad=0.4, borderaxespad=0.4)
ax.text(0.035, 0.965, "there is no operating point where this model\nis useful: 11% of engaged "
        "sites recovered at 5%\nfalse positives, 19% at 10%",
        transform=ax.transAxes, va="top", fontsize=5.2, color=st.GREY, linespacing=1.5)

# ---------------------------------------------------------------- d: would more features help?
# The reflex objection to a modest AUC is that the feature set is too small. It is not a
# feature-count problem: 14 features carry about five independent properties.
ax = fig.add_subplot(gs[1, 0]); st.panel_label(ax, "d", dx=-0.20, dy=1.08)
_x = [p[0] for p in PCS]; _ind = [p[1] for p in PCS]; _cum = [p[2] for p in PCS]
ax.bar(_x, _ind, color=st.LIGHTGREY, width=0.66, label="each property on its own")
ax.plot(_x, _cum, "-o", color=st.BLUE, lw=1.3, ms=3.4, mew=0, label="running total")
ax.axhline(80, color=st.ORANGE, lw=1.0, ls=(0, (3, 1.6)))
ax.text(len(_x) + 0.12, 81.5, " 80%", va="bottom", fontsize=5.4, color=st.ORANGE)
ax.set_xticks(_x); ax.set_xticklabels(_x, fontsize=6)
ax.set_xlabel("independent properties kept", fontsize=6.2)
ax.set_ylabel("of everything the 14 features\nmeasure, how much is captured (%)",
              fontsize=6.2, linespacing=1.4)
ax.set_xlim(0.4, len(_x) + 1.05); ax.set_ylim(0, 112)
ax.legend(frameon=False, fontsize=5.4, loc="center right", handlelength=1.2)
ax.text(0.03, 0.985, "14 features, about five independent properties.\nThe gap in a is not a "
        "feature-count problem.", transform=ax.transAxes, va="top", fontsize=5.2,
        color=st.GREY, linespacing=1.5)

st.save(fig, "fig5")
fd.head("fig5", HEAD)
for _e in (".png", ".pdf"):
    fig.savefig("../results/fig5_v1" + _e, dpi=400, bbox_inches="tight")
print("composite AUC %.3f [%.3f, %.3f] | NatChem ceiling %.3f [%.3f, %.3f] | gap %.3f"
      % (COMP + NAT_DEEP[3:6] + (GAP,)))
print("ceiling curves parsed:", {k: len(v["rows"]) for k, v in CEIL.items()})
print("bands parsed:", {k: len(v) for k, v in BANDS.items()})
print("breadth bins:", BRE)
print("enantiomer pairs:", [(e["name"], e["both"], e["plus"] + e["minus"], e["exp"]) for e in ENA])
print("c: operating points  " + " | ".join(
      "%.0f%% FPR -> %.1f%% recovered" % (100*f, 100*float(np.interp(f, ROC_P.fpr, ROC_P.tpr)))
      for f in (0.05, 0.10)))
print("d: PCs %d | individual %% %s | cumulative %% %s"
      % (len(PCS), [round(p[1], 1) for p in PCS], [round(p[2], 1) for p in PCS]))
print("d: PCs to clear 80%%: %d" % next(p[0] for p in PCS if p[2] >= 80))
_bb = [(ds, r[0], r[4], r[5], r[6], r[1]) for ds, v in BANDS.items() for r in v]
print("b: ceiling per band %.3f-%.3f over %d band-dataset points" %
      (min(b[2] for b in _bb), max(b[2] for b in _bb), len(_bb)))
for b in _bb:
    print("     %-14s %-18s ceiling %.3f [%.3f, %.3f]  n=%d" % b)
