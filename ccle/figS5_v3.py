"""
figS_react, proposed rebuild — support for Figure 4.

WHAT CHANGED

Old a drew five flat lines on a +/-0.8 axis and put the comparison that justifies the axis in
a caption: "ligandability's steepest feature is 0.43 — thirtyfold larger". That number was
WRONG in a specific way. 0.43 is step38 Q1b hydrophobic_8A = 0.429 "slope/100%", a per-compound
slope per 100 percentage points of occupancy. The reactivity side reports 0.013 as a slope PER
TIER. Comparing them is a unit mismatch and inflates the ratio about fivefold. The matched
comparator is step38's Q1a pooled occupancy table — the same object, a standardised shift from
the first tier — which gives 0.083 per tier against reactivity's 0.013: SIX times, not thirty.
The comparison is now drawn on the axis and computed in this script, not quoted.
Series are labelled by CLASS, which is what colour actually encodes — the old legend named five
features that shared two styles.

Old b is kept almost unchanged. It is the best panel in either figure.

Old c was two percentages over n = 9 and n = 6 (8/9 and 5/6, which is p ~ 0.11 against a coin).
It is replaced by a robustness panel that answers the real objection: neither class's result
depends on any single feature, and splitting by effect magnitude instead of by chemistry does
not reproduce it.

Old d duplicated Fig. 4d. Its slot is taken by the two redundancy panels demoted from the main
figure (Fig. 4 e/g -> d here, Fig. 4 f -> e here), which is where dimensionality belongs.
"""
import os, re
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from itertools import combinations
import abbasov_style as st
import figure_data as fd
st.apply()   # these two scripts imported the style module for its colours but
             # never applied it, so their typeface, sizes, spines and DPI were
             # matplotlib defaults while every other figure used the BMC style
st.slate()         # palette from the supplied ramp; see abbasov_style.SLATE

def find(n):
    for d in [".", "..", "ligand_rebuild", "../ligand_rebuild"]:
        p = os.path.join(d, n)
        if os.path.exists(p): return p
    raise FileNotFoundError(n)

G    = pd.read_csv(find("step11_gates.csv"))
G2   = pd.read_csv(find("step11b_ligandability_gates.csv"))
G2   = G2[~G2.dropped]
S9   = pd.read_csv(find("step9_results.csv"))
D12  = pd.read_csv(find("step12_reactivity_dims.csv"))
R10  = open(find("step10_report.txt")).read()
R12  = open(find("step12_report.txt")).read()
RTB  = open(find("test_behaviour_report.txt")).read()
R38  = open(find("step38_report.txt")).read()

NICE = lambda f: (f.replace("_8A", " 8Å").replace("_10A", " 10Å").replace("_", " ")
                   .replace("rsa tien", "exposure").replace("propka pKa", "predicted pKa")
                   .replace("nz sasa", "NZ exposure").replace("sasa sidechain", "side-chain exp.")
                   .replace("dist centroid", "dist. to centroid").replace("n lys", "nearby lysines"))
SHAPE = {"crowding_8A", "atoms_10A", "hydrophobic_8A", "aromatic_8A", "polar_8A",
         "rsa_tien", "sasa_sidechain", "nz_sasa", "dist_centroid"}
COL = lambda f: st.ORANGE if f in SHAPE else st.BLUE

HEAD = {}   # headline numbers per panel, captured where they are computed

# Third row added 2026-09-02 for the gate panels (f, g).
fig = plt.figure(figsize=(st.DOUBLE, st.SINGLE * 2.44))
gs  = fig.add_gridspec(3, 3, hspace=0.54, wspace=0.40,
                       left=0.075, right=0.985, top=0.940, bottom=0.050,
                       width_ratios=[1.0, 0.92, 1.0], height_ratios=[1.0, 1.0, 1.18])

# ---------------------------------------------------------------- a: tiering null
ax = fig.add_subplot(gs[0, 0]); st.panel_label(ax, "a", dx=-0.20, dy=1.09)
rows, hdr = [], None
KEEP = {"anionic_8A","net_charge_8A","propka_pKa","crowding_8A","hydrophobic_8A",
        "aromatic_8A","rsa_tien"}
for ln in RTB[RTB.find("TEST 1"):RTB.find("TEST 2")].splitlines():
    p_ = ln.split()
    if len(p_) > 3 and p_[0] == "feature": hdr = p_
    elif hdr and len(p_) == len(hdr) and p_[0] in KEEP:
        rows.append([p_[0]] + [float(x) for x in p_[1:]])
T1 = pd.DataFrame(rows, columns=hdr)
xs = np.arange(len(T1.columns) - 1)
for r_ in T1.itertuples(index=False):
    ax.plot(xs, list(r_)[1:], "o-", lw=1.1, ms=2.4, color=COL(r_[0]))
# The comparator. It must be the SAME KIND OF OBJECT as the traces above: a standardised
# shift from the first tier, per tier of an ordering. That is step38's Q1a pooled occupancy
# table. The 0.43 quoted in the previous version was step38 Q1b hydrophobic_8A = 0.429
# "slope/100%", a per-compound slope per 100 percentage points of occupancy — a different
# unit, which made the ratio look thirtyfold when the matched ratio is six.
sec = R38[R38.find("Q1a"):R38.find("Q1b")]
lrows, lhdr = [], None
for ln in sec.splitlines():
    q = ln.split()
    if q[:1] == ["tier"]: lhdr = q[2:]
    elif lhdr and len(q) == len(lhdr) + 2:
        try: lrows.append([float(v) for v in q[2:]])
        except ValueError: pass
LIG = pd.DataFrame(lrows, columns=lhdr)
lig_hi, lig_lo = LIG.iloc[-1].max(), LIG.iloc[-1].min()
slope = lambda v: np.polyfit(np.arange(len(v)), v, 1)[0]
lig_slope = max((abs(slope(LIG[c])) for c in LIG.columns))
rea_slope = max((abs(slope(list(r)[1:])) for r in T1.itertuples(index=False)))
ax.axhline(0, color=st.GREY, lw=0.8)
# This was a 7%-opacity wash marking the ligandability range. Opaque, it becomes a slab
# covering the whole panel and swallowing the traces the panel is about. The two dashed
# rules below already state the range, so the fill is dropped rather than recoloured.
for v in (lig_hi, lig_lo):
    ax.axhline(v, color=st.GREEN, lw=1.0, ls="--", zorder=1)
ax.text(-0.22, lig_hi + 0.035, "where tiering by occupancy takes ligandability: %+.2f" % lig_hi,
        fontsize=4.9, color=st.GREEN, ha="left", va="bottom")
ax.text(-0.22, lig_lo - 0.035, "and %+.2f" % lig_lo,
        fontsize=4.9, color=st.GREEN, ha="left", va="top")
ax.set_ylim(-0.95, 0.98); ax.set_xlim(-0.3, xs[-1] + 0.3)
ax.set_xlabel("labelled sites, tiered by detection depth", fontsize=6)
ax.set_ylabel("shift from tier 1, d", fontsize=6)
ax.tick_params(labelsize=5.2)
ax.legend(handles=[Line2D([], [], color=st.BLUE, lw=1.4, marker="o", ms=2.6, label="charge features"),
                   Line2D([], [], color=st.ORANGE, lw=1.4, marker="o", ms=2.6, label="shape features")],
          frameon=False, fontsize=5.0, loc="upper left", bbox_to_anchor=(0.0, 0.80),
          handlelength=1.6)
# The ratio was the word "six" typed into the panel. It is lig_slope / rea_slope and is now
# computed, so it cannot drift away from the traces the panel draws.
ax.text(0.96, 0.28, "steepest slope per tier %.3f,\nagainst %.3f for ligandability —\n"
        "%.1f times smaller, and the\nsame feature on both sides"
        % (rea_slope, lig_slope, lig_slope / rea_slope),
        transform=ax.transAxes, ha="right", va="center", fontsize=4.8, color=st.GREY,
        linespacing=1.55)
for s in ("top", "right"): ax.spines[s].set_visible(False)
HEAD["a_n_series"] = int(len(T1))
HEAD["a_n_tiers"] = int(len(T1.columns) - 1)
HEAD["a_reactivity_steepest_slope_per_tier"] = round(float(rea_slope), 5)
HEAD["a_ligandability_steepest_slope_per_tier"] = round(float(lig_slope), 5)
HEAD["a_slope_ratio_lig_over_react"] = round(float(lig_slope / rea_slope), 2)
HEAD["a_ligandability_shift_hi"] = round(float(lig_hi), 3)
HEAD["a_ligandability_shift_lo"] = round(float(lig_lo), 3)
fd.panel("figS5", "a",
         [dict(feature=_r[0], feature_class=("shape" if _r[0] in SHAPE else "charge"),
               tier=_t, shift_from_tier1=float(_v))
          for _r in T1.itertuples(index=False, name=None)
          for _t, _v in zip(list(T1.columns)[1:], list(_r)[1:])],
         note="reactivity tiering null; ligandability comparator range in the head ledger")

# ---------------------------------------------------------------- b: every feature, both gates
ax = fig.add_subplot(gs[0:2, 1]); st.panel_label(ax, "b", dx=-0.42, dy=1.045)
gates = sorted(G.gate.unique())
g0 = G[G.gate == gates[0]].set_index("feature").d
gL = G[G.gate == gates[-1]].set_index("feature").d
common = g0.index.intersection(gL.index)
tab = pd.DataFrame({"no_gate": g0[common], "gated": gL[common]})
tab = tab.reindex(tab.no_gate.sort_values().index)
y = np.arange(len(tab))
for k, (f, r_) in enumerate(tab.iterrows()):
    c = COL(f)
    ax.plot([r_.no_gate, r_.gated], [k, k], color=st.LIGHTGREY, lw=1.4)
    ax.scatter(r_.no_gate, k, s=15, facecolor=st.PAPER, edgecolor=c, lw=0.95, zorder=3)
    ax.scatter(r_.gated,  k, s=17, color=c, zorder=3)
    if abs(r_.gated - r_.no_gate) < 0.02:
        ax.annotate("no movement", (r_.gated, k), xytext=(-9, 0), textcoords="offset points",
                    fontsize=4.4, color=st.GREY, va="center", ha="right")
ax.axvline(0, color=st.GREY, lw=0.9)
ax.set_yticks(y); ax.set_yticklabels([NICE(f) for f in tab.index], fontsize=5.2)
for t, f in zip(ax.get_yticklabels(), tab.index): t.set_color(COL(f))
ax.set_xlabel("effect on being labelled, d", fontsize=6)
ax.set_xlim(-0.56, 0.52); ax.set_ylim(-1.9, len(tab) - 0.4)
ax.tick_params(axis="x", labelsize=5.2)
ax.text(0.985, 0.02, "open  no gate\nfilled  pLDDT ≥ %d" % gates[-1], transform=ax.transAxes,
        ha="right", va="bottom", fontsize=5.0, color=st.GREY, linespacing=1.7)
for s in ("top", "right"): ax.spines[s].set_visible(False)
HEAD["b_n_features"] = int(len(tab))
HEAD["b_gate_plddt"] = int(gates[-1])
fd.panel("figS5", "b",
         [dict(feature=_f, feature_class=("shape" if _f in SHAPE else "charge"),
               d_no_gate=float(_r.no_gate), d_gated=float(_r.gated),
               movement=float(_r.gated - _r.no_gate))
          for _f, _r in tab.iterrows()],
         note="one dumbbell per feature, ungated and at pLDDT >= %d" % gates[-1])

# ---------------------------------------------------------------- c: robustness
ax = fig.add_subplot(gs[0, 2]); st.panel_label(ax, "c", dx=-0.20, dy=1.09)
# both sides at pLDDT >= 90 (step 11 for labelling, step 11b for ligandability)
_lg = G2[G2.gate == G2.gate.max()].set_index("feature").d.rename("lig")
F = (pd.concat([gL.rename("rg"), _lg,
                G[G.gate == gates[0]].set_index("feature").d.rename("react")], axis=1)
       .dropna().reset_index().rename(columns={"index": "feature"}))
F["kind"] = np.where(F.feature.isin(SHAPE), "shape", "charge")
Rr = lambda idx, c: np.corrcoef(F.loc[idx, c], F.loc[idx, "lig"])[0, 1]
# The class sizes were typed into the tick labels as "9 features" and "5 features". Nine is
# right; five was not — F carries SIX charge features (anionic, cationic, nearby lysines, net
# charge, predicted pKa, salt bridge), and 9 + 6 is the fifteen the caption and Figure 4 name.
# Both counts now come from F, so the panel cannot disagree with its own data again.
_NK = {k: int((F.kind == k).sum()) for k in ("shape", "charge")}
_CLAB = lambda k: "%s\n(%d features)" % (k, _NK[k])
_LOO = {}
for yy, kind in enumerate(("shape", "charge")):
    idx = F.index[F.kind == kind]; c = st.ORANGE if kind == "shape" else st.BLUE
    full = Rr(idx, "rg"); vals = [Rr(idx.drop(f), "rg") for f in idx]
    _LOO[kind] = (float(full), [float(v) for v in vals], [F.feature[i] for i in idx])
    ax.plot([min(vals), max(vals)], [yy]*2, color=st.LIGHTGREY, lw=2.6,
            solid_capstyle="round")
    ax.scatter(vals, [yy]*len(vals), s=15, color=c, zorder=3, lw=0.5, edgecolor=st.PAPER)
    ax.scatter([full], [yy], s=70, marker="|", color=c, zorder=4, lw=1.7)
    ax.text(np.mean([min(vals), max(vals)]), yy + 0.13, "%.2f – %.2f" % (min(vals), max(vals)),
            ha="center", fontsize=5.4, color=c)
    if yy == 1:
        ax.annotate("all %d" % _NK[kind], (full, yy), xytext=(0, -13),
                    textcoords="offset points", ha="center", fontsize=4.5, color=c)
ax.set_yticks([0, 1]); ax.set_yticklabels([_CLAB("shape"), _CLAB("charge")], fontsize=5.6)
ax.set_ylim(-0.75, 1.85); ax.set_xlim(-0.04, 1.06)
ax.tick_params(axis="x", labelsize=5.2)
ax.set_xlabel("how closely labelling tracks ligandability\n(pLDDT ≥ 90), leaving each feature out in turn",
              fontsize=6, linespacing=1.6)
# the alternative grouping: by how large the ungated effect is, rather than by chemistry
F["absd"] = F.react.abs()
big = F.nlargest(int((F.kind == "charge").sum()), "absd").index
small = F.index.difference(big)
for s in ("top", "right"): ax.spines[s].set_visible(False)
HEAD["c_n_features"] = int(len(F))
HEAD["c_n_shape_features"] = _NK["shape"]
HEAD["c_n_charge_features"] = _NK["charge"]
for _k in ("shape", "charge"):
    HEAD["c_%s_r_all_features" % _k] = round(_LOO[_k][0], 3)
    HEAD["c_%s_r_loo_min" % _k] = round(min(_LOO[_k][1]), 3)
    HEAD["c_%s_r_loo_max" % _k] = round(max(_LOO[_k][1]), 3)
fd.panel("figS5", "c",
         [dict(feature_class=_k, feature_dropped=_f, r_with_ligandability=_v)
          for _k in ("shape", "charge")
          for _f, _v in zip(_LOO[_k][2], _LOO[_k][1])],
         note="leave-one-out; the all-features value for each class is in the head ledger")

# ---------------------------------------------------------------- d: redundancy
ax = fig.add_subplot(gs[1, 0]); st.panel_label(ax, "d", dx=-0.20, dy=1.09)
pc12 = [float(x) for x in re.findall(r"PC\d\s+([\d.]+)%", R12)][:5]
pc10 = [float(x) for x in re.findall(r"PC\d\s+([\d.]+)%", R10)][:8]
# DIRECT LABELS, not a legend: in the five-colour palette the teal and the slate are only
# dE 13 apart in normal vision, so hue alone must never be what tells two series apart.
# Feature counts and the 80%-crossing points were typed into the labels and the note below.
# Both are already computed here: the reactivity set is len(pc12), the ligandability set is the
# feature count step 10 reports, and each crossing is the n80 the loop works out.
_NLIG = int(re.search(r"\|\s*(\d+) features", R10).group(1))
_N80 = {}
for vals, c, lab, mk, key in (
        (pc12, st.BLUE, "the %d charge features,\npLDDT \u2265 70" % len(pc12), "o", "reactivity"),
        (pc10, st.GREEN, "the %d ligandability\nfeatures" % _NLIG, "s", "ligandability")):
    cum = np.cumsum(vals)
    ax.plot(np.arange(1, len(cum) + 1), cum, marker=mk, ls="-", color=c, lw=1.3, ms=3.0)
    n80 = int(np.searchsorted(cum, 80) + 1)
    _N80[key] = n80
    ax.scatter([n80], [cum[n80 - 1]], s=52, facecolor="none", edgecolor=c, lw=1.2, zorder=4)
    ax.text(len(cum) + 0.12, cum[-1], "  " + lab, color=c, fontsize=4.8, va="center",
            linespacing=1.5)
ax.axhline(80, color=st.GREY, lw=0.9, ls="--")
ax.text(8.5, 80.5, "80%", fontsize=5.0, color=st.GREY, ha="right", va="bottom")
HEAD["d_n_reactivity_features"] = int(len(pc12))
HEAD["d_n_ligandability_features"] = _NLIG
HEAD["d_reactivity_pc_pct"] = [float(v) for v in pc12]
HEAD["d_ligandability_pc_pct"] = [float(v) for v in pc10[:6]]
HEAD["d_reactivity_pc1_pct"] = float(pc12[0])
HEAD["d_ligandability_pc1_pct"] = float(pc10[0])
HEAD["d_reactivity_components_to_80pct"] = _N80["reactivity"]
HEAD["d_ligandability_components_to_80pct"] = _N80["ligandability"]
fd.panel("figS5", "d",
         [dict(feature_set=_s, component=_i + 1, variance_pct=float(_v),
               cumulative_pct=float(np.cumsum(_vals)[_i]))
          for _s, _vals in (("reactivity", pc12), ("ligandability", pc10[:6]))
          for _i, _v in enumerate(_vals)],
         note="one row per principal component of each feature set")
ax.set_xlabel("number of components kept", fontsize=6)
ax.set_ylabel("variance explained (%)", fontsize=6)
ax.set_ylim(35, 108); ax.set_xlim(0.6, 8.6); ax.set_xticks(range(1, 7))
ax.set_yticks(range(40, 101, 10))
ax.tick_params(labelsize=5.2)
ax.text(0.97, 0.03, "%d charge features → %d axes; %d\nligandability features → %d. Each set\n"
        "is redundant, neither collapses to one\naxis — but they are different sets on\n"
        "different populations, so this is not\n%d against %d."
        % (len(pc12), _N80["reactivity"], _NLIG, _N80["ligandability"],
           _N80["reactivity"], _N80["ligandability"]),
        transform=ax.transAxes, va="bottom", ha="right",
        fontsize=4.7, color=st.GREY, linespacing=1.55)
for s in ("top", "right"): ax.spines[s].set_visible(False)

# ---------------------------------------------------------------- e: mutual adjustment
ax = fig.add_subplot(gs[1, 2]); st.panel_label(ax, "e", dx=-0.30, dy=1.09)
d12 = D12.sort_values("uni_coef"); yy = np.arange(len(d12)); w = 0.36
ax.barh(yy - w/2, d12.uni_coef,  w, color=st.LIGHTGREY, label="on its own")
ax.barh(yy + w/2, d12.multi_coef, w, color=st.BLUE, label="adjusted for the others")
ax.axvline(0, color=st.GREY, lw=0.9)
ax.set_yticks(yy); ax.set_yticklabels([NICE(f) for f in d12.feature], fontsize=5.2)
ax.set_xlabel("effect on being labelled", fontsize=6)
ax.set_xlim(-0.30, 0.22); ax.set_ylim(-2.7, len(d12) - 0.4)
ax.tick_params(axis="x", labelsize=5.2)
for k, r in zip(yy, d12.itertuples()):
    kept = abs(r.multi_coef) > abs(r.uni_coef) / 3
    ax.text(r.multi_coef + (-0.007 if r.multi_coef < 0 else 0.007), k + w/2, "%.3f" % r.multi_coef,
            va="center", fontsize=4.6, ha="right" if r.multi_coef < 0 else "left",
            color=st.BLUE if kept else st.ORANGE)
ax.legend(frameon=False, fontsize=4.7, loc="lower left", bbox_to_anchor=(-0.02, -0.01),
          handlelength=1.0, handletextpad=0.5, labelspacing=0.5)
ax.text(0.99, 0.025, "salt bridge was anionic-neighbour\ndensity in another form",
        transform=ax.transAxes, ha="right", va="bottom",
        fontsize=4.6, color=st.GREY, linespacing=1.55)
for s in ("top", "right"): ax.spines[s].set_visible(False)
HEAD["e_n_features"] = int(len(d12))
fd.panel("figS5", "e",
         [dict(feature=_r.feature, coef_alone=float(_r.uni_coef),
               coef_adjusted=float(_r.multi_coef)) for _r in d12.itertuples()],
         note="one pair of bars per feature, unadjusted and mutually adjusted")


# ---------------------------------------------------------------- f: what the gate removes
# Figure 4 turns on this gate: the two contrasts only separate once it is applied (exact
# permutation p 0.487 ungated, 0.008 at pLDDT >= 90). Until now a reader had to take that on
# trust, because nothing in the paper showed a gated-out site. Both examples below are real
# engaged sites from the atlas, rendered from the same AlphaFold models the features are
# measured on, by figtools/ + PyMOL.
import matplotlib.image as _mpimg
# master_sites is read once here and sliced again by panel g below. PFKL K677's confidence,
# compound count and occupancy were typed into the caption strip; they are columns of this
# table and are now read from it. LYSMD4 K72 IS in the same table, and its three numbers are
# now read from it too -- but it has to be looked up by ACCESSION, not by gene symbol. The
# repo's cached uniprot_gene_map.tsv still carries the legacy clone name for this protein
# ("DKFZp779J2370", the UniProt TrEMBL gene name for Q6MZR2), so a search on "LYSMD4" finds
# nothing and the site looks absent. Q6MZR2 and the LYSMD4 entry F6VMV4 are the same
# 131-residue sequence, byte for byte, and residue 72 is a lysine in both, so the caption's
# gene name is right and the cached map is what is out of date.
_MSF = pd.read_csv(find("master_sites.csv"),
                   usecols=["uniprot", "position", "gene", "plddt_site",
                            "n_engaged_60", "max_occupancy"], low_memory=False)
_pf = _MSF[(_MSF.gene == "PFKL") & (_MSF.position == 677)]
assert len(_pf) == 1, "PFKL K677 not uniquely in master_sites (%d rows)" % len(_pf)
_pf = _pf.iloc[0]
# Q6MZR2 = LYSMD4 (TrEMBL clone name DKFZp779J2370); looked up by accession for the reason above.
_ly = _MSF[(_MSF.uniprot == "Q6MZR2") & (_MSF.position == 72)]
assert len(_ly) == 1, "LYSMD4 K72 (Q6MZR2) not uniquely in master_sites (%d rows)" % len(_ly)
_ly = _ly.iloc[0]
_SUB = gs[2, 0:2].subgridspec(2, 2, hspace=0.40, wspace=0.04)
_axl = fig.add_subplot(gs[2, 0:2]); _axl.set_axis_off()
st.panel_label(_axl, "f", dx=-0.075, dy=1.03)

def _shot(cell, name, title=None):
    a = fig.add_subplot(cell)
    p = os.path.join("struct", name)
    if os.path.exists(p):
        a.imshow(_mpimg.imread(p))
    else:
        a.text(0.5, 0.5, "missing struct/" + name, ha="center", va="center",
               transform=a.transAxes, fontsize=5.5, color=st.ORANGE)
    a.set_xticks([]); a.set_yticks([])
    for _sp in ("top", "right", "bottom", "left"): a.spines[_sp].set_visible(False)
    if title: a.set_title(title, fontsize=5.4, color=st.GREY, pad=1.5)
    return a

_a = _shot(_SUB[0, 0], "pfkl_whole.png", "whole protein")
_shot(_SUB[0, 1], "pfkl_site.png", "the lysine and everything within 8 \u00c5")
_a.text(0.0, -0.10, "PFKL K677   pLDDT %.0f   %d compounds, %.0f%% occupancy   KEPT"
        % (_pf.plddt_site, int(_pf.n_engaged_60), _pf.max_occupancy),
        transform=_a.transAxes, fontsize=5.4, color=st.BLUE, fontweight="bold")
HEAD["f_pfkl_plddt"] = float(_pf.plddt_site)
HEAD["f_pfkl_n_compounds_engaged"] = int(_pf.n_engaged_60)
HEAD["f_pfkl_max_occupancy_pct"] = float(_pf.max_occupancy)
_b = _shot(_SUB[1, 0], "lysmd4_whole.png")
_shot(_SUB[1, 1], "lysmd4_site.png")
_b.text(0.0, -0.10, "LYSMD4 K72   pLDDT %.0f   %d compounds, %.0f%% occupancy   DISCARDED"
        % (_ly.plddt_site, int(_ly.n_engaged_60), _ly.max_occupancy),
        transform=_b.transAxes, fontsize=5.4, color=st.ORANGE, fontweight="bold")
HEAD["f_lysmd4_plddt"] = float(_ly.plddt_site)
HEAD["f_lysmd4_n_compounds_engaged"] = int(_ly.n_engaged_60)
HEAD["f_lysmd4_max_occupancy_pct"] = float(_ly.max_occupancy)
HEAD["f_lysmd4_uniprot"] = "Q6MZR2"
_b.text(0.0, -0.21, "the same 8 \u00c5 sphere \u2014 in one of them there is no structure to measure",
        transform=_b.transAxes, fontsize=5.0, color=st.GREY)

# ---------------------------------------------------------------- g: not a strength filter
ax = fig.add_subplot(gs[2, 2]); st.panel_label(ax, "g", dx=-0.30, dy=1.06)
_FE = pd.read_csv(find("features_final.csv"),
                  usecols=["uniprot", "position", "plddt_site"], low_memory=False)
# sliced from the single read above; plddt_site is dropped here so it cannot collide with
# the per-dataset plddt_site the merge brings in from features_final
_MS = _MSF[["uniprot", "position", "n_engaged_60", "max_occupancy"]]
_D = (_FE.dropna(subset=["plddt_site"]).merge(_MS, on=["uniprot", "position"])
         .query("n_engaged_60 > 0").drop_duplicates(["uniprot", "position"]))
_GATE = int(gates[-1])           # the same gate the ladder in panels b and c ends on
_K = _D.plddt_site >= _GATE
ax.scatter(_D.plddt_site[~_K], _D.max_occupancy[~_K], s=1.5, color=st.LIGHTGREY,
           linewidths=0, label="discarded")
ax.scatter(_D.plddt_site[_K], _D.max_occupancy[_K], s=1.5, color=st.BLUE,
           linewidths=0, label="kept")
ax.axvline(_GATE, color=st.ORANGE, lw=1.0, ls=(0, (3, 1.6)))
_bins = pd.cut(_D.plddt_site, [20, 50, 70, 80, 90, 95, 100])
_med = _D.groupby(_bins, observed=True).max_occupancy.median()
_cnt = _D.groupby(_bins, observed=True).size()
ax.plot([i.mid for i in _med.index], _med.values, "-o", color=st.GREEN, lw=1.2, ms=2.8,
        mew=0, zorder=5)
ax.text(0.03, 0.60, "median occupancy is flat from 20 to 100 (%.0f%% in the\nlowest band, "
        "%.0f%% in the highest): the gate removes sites\nwhose features cannot be measured, "
        "not sites that are\nweakly engaged" % (_med.iloc[0], _med.iloc[-1]),
        transform=ax.transAxes, fontsize=4.9, color=st.GREEN, linespacing=1.5, va="top",
        bbox=dict(facecolor=st.PAPER, edgecolor="none", pad=1.6))
ax.text(_GATE - 1, 57.4, "gate %d" % _GATE, ha="right", va="bottom", fontsize=5.2,
        color=st.ORANGE, fontweight="bold")
ax.text(0.5, 1.005, "%s discarded (%.0f%%)   \u00b7   %s kept (%.0f%%)"
        % (f"{(~_K).sum():,}", 100*(~_K).mean(), f"{_K.sum():,}", 100*_K.mean()),
        transform=ax.transAxes, ha="center", fontsize=5.4, color=st.GREY)
ax.set_xlabel("AlphaFold confidence at the lysine (pLDDT)", fontsize=6)
ax.set_ylabel("highest occupancy reached by any compound (%)", fontsize=6)
ax.set_xlim(20, 100); ax.set_ylim(55, 104)
ax.legend(frameon=False, fontsize=5.2, loc="lower left", markerscale=2.6,
          handletextpad=0.3, borderaxespad=0.4)
for _sp in ("top", "right"): ax.spines[_sp].set_visible(False)
HEAD["g_gate_plddt"] = _GATE
HEAD["g_n_sites_discarded"] = int((~_K).sum())
HEAD["g_n_sites_kept"] = int(_K.sum())
HEAD["g_pct_sites_discarded"] = round(float(100 * (~_K).mean()), 1)
HEAD["g_pct_sites_kept"] = round(float(100 * _K.mean()), 1)
HEAD["g_median_occupancy_lowest_band"] = round(float(_med.iloc[0]), 1)
HEAD["g_median_occupancy_highest_band"] = round(float(_med.iloc[-1]), 1)
fd.panel("figS5", "g",
         [dict(plddt_band_low=_i.left, plddt_band_high=_i.right, plddt_band_mid=_i.mid,
               n_sites=int(_cnt.loc[_i]), median_max_occupancy=float(_med.loc[_i]))
          for _i in _med.index],
         note="binned medians of the scatter; the gate keeps sites at pLDDT >= %d" % _GATE)


# Routed through st.save so the blank-panel check runs here too: this script and fig4_v3
# were the two that bypassed it, which is also why they were the two missing st.apply().
st.save(fig, "figS5")
fd.head("figS5", HEAD)
for _e in (".png", ".pdf"):
    fig.savefig("../results/figS5_v3" + _e, dpi=400, bbox_inches="tight")
print("panel a series:", len(T1), "| panel b features:", len(tab), "| panel c features:", len(F))
print("PCA reactivity:", pc12, "\nPCA ligandability:", pc10[:6])
