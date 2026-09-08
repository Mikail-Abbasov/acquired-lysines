#!/usr/bin/env python3
"""
Figure S7 — support for Figure 6.

Figure 6 claims that engagement is a property of the lysine and not of the molecule, resting on a
null: no compound descriptor predicts which lysines a compound prefers. A null needs three things
to be believed — that it is not an artefact of the design, that it is not an artefact of the 60%
threshold, and that it is well powered — and it needs its boundary stated. Those are panels a-f.

  a  THE CONFOUND, NAMED. Pooled across datasets, 21 of 150 descriptor x site-feature pairs clear
     q < 0.05. Every one of them involves a SIZE descriptor — molecular weight, heavy-atom count,
     rotatable bonds, fsp3 — and molecular weight is almost a label for which screen a compound
     came from. Within dataset: 0 of 150.
  b  The per-compound version of Fig. 6b, across occupancy rather than at one cut: 10 of 12
     features strengthen with occupancy in 173 compounds.
  c  No threshold at all. Site features against three continuous occupancy measures; 14, 12 and
     15 of 15 features remain significant.
  d  The decile profile, again threshold-free.
  e  The compound-count gradient at five cuts, WITH the group sizes — 40% and 80% are the two
     that break monotonicity, and the panel shows why: n collapses at the extremes.
  f  WHAT THIS NULL CANNOT SEE. Not one of the ten descriptors or six shape terms encodes charge,
     and every shape term is achiral by construction — NPR1 and NPR2 are identical for
     enantiomers (step 32). So Fig. 6d is silent on exactly the two axes Figs. 5f and 7 test.
"""
import os, re, sys, warnings
warnings.filterwarnings("ignore")
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
    raise FileNotFoundError(n)

R38 = open(find("step38_report.txt")).read()
R40 = open(find("step40_report.txt")).read()
# Panel a is built from step 37 / step 38 -- the same three tables Figure 6d plots -- so the
# SI panel and the main-text claim cannot disagree. See the note above panel a.
POOL = pd.read_csv(find("step37_within_compound.csv"))
NATC = pd.read_csv(find("step37_within_NatChem.csv"))
SLOP = pd.read_csv(find("step38_tiers.csv"))
C40 = pd.read_csv(find("step40_continuous.csv"))
S37 = pd.read_csv(find("step37_within_compound.csv"))

NICE = lambda f: (f.replace("_8A", " 8Å").replace("_10A", " 10Å").replace("_", " ")
                   .replace("rsa tien", "exposure").replace("propka pKa", "predicted pKa")
                   .replace("nz sasa", "NZ exposure").replace("sasa sidechain", "side-chain exp.")
                   .replace("dist centroid", "dist. to centroid").replace("n lys", "nearby lysines")
                   .replace("dist nearest pocket", "dist. to pocket")
                   .replace("hbond partners NZ", "H-bond partners")
                   .replace("plddt site", "pLDDT").replace("plddt 8A", "pLDDT 8Å"))

# 173 was typed in three places -- this regex, the standard-error denominator in panel b and
# panel b's own caption. Step 38 states it once, so it is read from there and used in all three.
_NCMP = int(re.search(r"compounds with slopes across >=\d+ tiers: (\d+)", R38).group(1))
SLOPE = [(m.group(1), float(m.group(2)), float(m.group(3)), float(m.group(4)), bool(m.group(5)))
         for m in re.finditer(r"^\s{2,}(\w+)\s+(-?[\d.]+)\s+([\d.]+)\s+%d\s+([\d.e+-]+)\s*(\*)?$"
                              % _NCMP, R38, re.M)]
DEC, dhdr = [], None
_blk = R40[R40.index("DECILE PROFILE"):R40.index("THRESHOLD-VARYING")]
for line in _blk.splitlines():
    if "rsa_tien" in line and "crowding" in line:
        dhdr = line.split()
    m = re.match(r"^\s{2,}(\d+)\s+(\d+)\s+([\d.]+)\s+((?:-?[\d.]+\s*)+)$", line)
    if m and dhdr:
        DEC.append((int(m.group(1)), int(m.group(2)), float(m.group(3)),
                    [float(v) for v in m.group(4).split()]))
CUTS = []
for m in re.finditer(r"CUT AT (\d+)% — (\d+) sites with zero liganding compounds\n(.*?)-> (NOT )?monotonic",
                     R40, re.S):
    rows = [(r.group(1), int(r.group(2)), float(r.group(3)), float(r.group(4)))
            for r in re.finditer(r"^\s{2,}([\d>=-]+)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+[\d.]+\s*\*?$",
                                 m.group(3), re.M)]
    CUTS.append((int(m.group(1)), rows, m.group(4) is None))
assert len(SLOPE) == 12 and len(DEC) == 10 and len(CUTS) == 5, (len(SLOPE), len(DEC), len(CUTS))

HEAD = {}   # headline numbers per panel, captured where they are computed

fig = plt.figure(figsize=(st.DOUBLE, st.DOUBLE * 0.90))
gs = fig.add_gridspec(3, 3, hspace=0.66, wspace=0.46,
                      left=0.072, right=0.982, top=0.945, bottom=0.058,
                      width_ratios=[1.0, 1.0, 1.0], height_ratios=[1.0, 1.0, 0.82])

# ---------------------------------------------------------------- a: name the confound
ax = fig.add_subplot(gs[0, 0]); st.panel_label(ax, "a", dx=-0.30, dy=1.07)
# REBUILT ON STEP 37, 2026-09-07.
# This panel used to plot step 31: 150 pairs over 10 descriptors -- a DIFFERENT analysis from the
# one Figure 6d plots (step 37, 84 pairs over 6 descriptors). Both returned 21 significant pairs
# by coincidence, which invited a reader to take them for one number. Worse, step 31's count
# depended on a MIN_SITES threshold the panel never disclosed: 21 of 150 at >= 20 engaged sites,
# 16 of 150 at >= 15 (step 37's own rule), and 4 of 150 -- below the 7.5 expected by chance, with
# no size descriptor among them -- at the script default of 5. The panel now names the confound
# using the same analysis the main text relies on, at step 37's single stated inclusion rule.
SIG = POOL[POOL.q < 0.05]
cnt = SIG["shape"].value_counts()
allc = POOL["shape"].value_counts()
order = list(cnt.index) + [d for d in allc.index if d not in cnt.index]
y = np.arange(len(order))[::-1]
# Rust means only what it can mean here -- this descriptor produced pooled positives. WHY those
# descriptors is answered quantitatively in Figure 6f (step 41).
for yy, d in zip(y, order):
    n = int(cnt.get(d, 0))
    ax.barh([yy], [n], height=0.66, color=st.ORANGE if n > 0 else st.LIGHTGREY)
    ax.text(-0.12, yy, str(d).replace("_", " "), ha="right", va="center", fontsize=5.2,
            color=st.GREY, fontweight="bold" if n > 0 else "normal")
ax.set_yticks([])
ax.set_xlabel("pairs clearing q < 0.05 when the three\ndatasets are pooled",
              fontsize=6.0, linespacing=1.4)
ax.set_xlim(0, max(cnt.values) * 1.18); ax.set_ylim(-1.2, len(order) - 0.4)
ax.text(0.97, 0.05, "%d of %d pooled (%d compounds)\n%d of %d within NatChem (%d)\n"
        "%d of %d across slopes (%d)"
        % (len(SIG), len(POOL), int(POOL.n.iloc[0]),
           int((NATC.q < 0.05).sum()), len(NATC), int(NATC.n_compounds.iloc[0]),
           int((SLOP.q < 0.05).sum()), len(SLOP), int(SLOP.n.iloc[0])),
        transform=ax.transAxes, ha="right", va="bottom", fontsize=5.0, color=st.ORANGE,
        fontweight="bold", linespacing=1.4)
for s_ in ("top", "right", "left"): ax.spines[s_].set_visible(False)
HEAD["a_n_pairs_tested"] = int(len(POOL))
HEAD["a_n_pairs_significant_pooled"] = int(len(SIG))
HEAD["a_n_compounds_pooled"] = int(POOL.n.iloc[0])
HEAD["a_n_descriptors"] = int(len(allc))
HEAD["a_n_descriptors_with_hits"] = int(len(cnt))
HEAD["a_pairs_by_descriptor"] = {str(k): int(v) for k, v in cnt.items()}
HEAD["a_n_significant_within_natchem"] = int((NATC.q < 0.05).sum())
HEAD["a_n_pairs_within_natchem"] = int(len(NATC))
HEAD["a_n_compounds_within_natchem"] = int(NATC.n_compounds.iloc[0])
HEAD["a_n_significant_slopes"] = int((SLOP.q < 0.05).sum())
HEAD["a_n_pairs_slopes"] = int(len(SLOP))
HEAD["a_n_compounds_slopes"] = int(SLOP.n.iloc[0])
fd.panel("figS7", "a",
         [dict(descriptor=_d, pairs_significant=int(cnt.get(_d, 0)),
               pairs_tested=int(allc.get(_d, 0))) for _d in order],
         note="step 37 pooled: one bar per compound descriptor, pairs clearing q < 0.05")

# ---------------------------------------------------------------- b: slopes across occupancy
ax = fig.add_subplot(gs[0, 1]); st.panel_label(ax, "b", dx=-0.34, dy=1.07)
S = sorted(SLOPE, key=lambda r: r[1])
y = np.arange(len(S))
for i, (f_, sl, sd, q_, sig) in enumerate(S):
    se = sd / np.sqrt(_NCMP)
    c = st.BLUE if sig else st.ORANGE
    ax.plot([sl - 1.96 * se, sl + 1.96 * se], [i, i], color=c, lw=1.3, solid_capstyle="round")
    ax.plot([sl], [i], "o", color=c, ms=3.6, mew=0)
ax.axvline(0, color=st.GREY, lw=0.8)
ax.set_yticks(y); ax.set_yticklabels([NICE(r[0]) for r in S], fontsize=5.2)
ax.set_xlabel("how much the preference strengthens\nper 100 points of occupancy",
              fontsize=6.0, linespacing=1.4)
ax.set_xlim(-0.55, 0.55); ax.set_ylim(-0.6, len(S) + 0.3)
# Top left is empty: the strongest slopes are positive and sit right.
ax.text(0.02, 0.98, "%d of %d\nn = %d compounds"
        % (sum(1 for r in SLOPE if r[4]), len(SLOPE), _NCMP),
        transform=ax.transAxes, ha="left", va="top", fontsize=5.4, color=st.BLUE,
        fontweight="bold", linespacing=1.4)
for s in ("top", "right"): ax.spines[s].set_visible(False)
HEAD["b_n_features"] = int(len(SLOPE))
HEAD["b_n_features_strengthening"] = int(sum(1 for r in SLOPE if r[4]))
HEAD["b_n_features_not_strengthening"] = int(len(SLOPE) - sum(1 for r in SLOPE if r[4]))
HEAD["b_n_compounds"] = _NCMP
fd.panel("figS7", "b",
         [dict(feature=_f, slope_per_100pct_occupancy=_sl, sd=_sd,
               se=_sd / np.sqrt(_NCMP), q=_q, strengthens=bool(_sig))
          for _f, _sl, _sd, _q, _sig in S],
         note="one point per site feature, mean per-compound slope over %d compounds" % _NCMP)

# ---------------------------------------------------------------- c: no threshold at all
ax = fig.add_subplot(gs[0, 2]); st.panel_label(ax, "c", dx=-0.34, dy=1.07)
MEAS = ["mean_occ", "max_occ", "sum_occ_norm"]
LBLM = {"mean_occ": "mean occupancy", "max_occ": "highest occupancy",
        "sum_occ_norm": "total occupancy"}
feats = (C40[C40.measure == "mean_occ"].sort_values("rho").feature.tolist())
y = np.arange(len(feats))
MKS = {"mean_occ": ("o", st.BLUE), "max_occ": ("s", st.ORANGE), "sum_occ_norm": ("^", st.GREEN)}
for m_ in MEAS:
    d = C40[C40.measure == m_].set_index("feature")
    mk, c = MKS[m_]
    ax.plot([d.rho.get(f, np.nan) for f in feats], y, mk, color=c, ms=3.0, mew=0,
            label="%s (%d of %d)" % (LBLM[m_], int((d.q < 0.05).sum()), len(d)))
ax.axvline(0, color=st.GREY, lw=0.8)
ax.set_yticks(y); ax.set_yticklabels([NICE(f) for f in feats], fontsize=4.6)
ax.set_xlabel("correlation with occupancy (Spearman ρ)", fontsize=6.0)
ax.set_ylim(-6.3, len(feats) + 0.6)   # blank strip below for the legend; no caption above
ax.legend(frameon=False, fontsize=4.9, loc="lower left", handletextpad=0.3, borderaxespad=0.3)
# The caption used to sit bottom-right, directly on the legend. The top-left corner is empty
# because features are sorted by rho, so the highest rows are the ones furthest right.

for s in ("top", "right"): ax.spines[s].set_visible(False)
HEAD["c_n_features"] = int(len(feats))
for _m in MEAS:
    _d = C40[C40.measure == _m]
    HEAD["c_n_significant_" + _m] = int((_d.q < 0.05).sum())
fd.panel("figS7", "c",
         [dict(measure=_r.measure, feature=_r.feature, rho=float(_r.rho), q=float(_r.q),
               significant=bool(_r.q < 0.05))
          for _r in C40[C40.measure.isin(MEAS)].itertuples()],
         note="Spearman rho against three threshold-free occupancy measures")

# ---------------------------------------------------------------- d: decile profile
ax = fig.add_subplot(gs[1, 0]); st.panel_label(ax, "d", dx=-0.24, dy=1.07)
DLB = {"rsa_tien": "exposure", "crowding_": "crowding 8Å", "atoms_10A": "atoms 10Å",
       "propka_pK": "predicted pKa", "in_pocket": "in pocket"}
DC = [st.ORANGE, st.BLUE, st.BLUE, st.ORANGE, st.GREEN]
xs = [d[2] for d in DEC]
for j, h in enumerate(dhdr):
    ax.plot(xs, [d[3][j] for d in DEC], "-o", color=DC[j], lw=1.0, ms=2.4, mew=0,
            ls="-" if j % 2 == 0 else "--", label=DLB.get(h, h))
ax.axhline(0, color=st.LIGHTGREY, lw=0.8)
ax.set_xlabel("median highest occupancy in the decile (%)", fontsize=6.0)
ax.set_ylabel("shift from the lowest decile", fontsize=6.0)
ax.legend(frameon=False, fontsize=4.9, loc="lower left", handlelength=1.4, ncol=2,
          columnspacing=0.9, borderaxespad=0.3)
ax.set_ylim(-0.34, 0.32)
HEAD["d_n_deciles"] = int(len(DEC))
HEAD["d_n_features"] = int(len(dhdr))
fd.panel("figS7", "d",
         [dict(decile=_d[0], n=_d[1], median_max_occupancy_pct=_d[2],
               feature=_h, shift_from_lowest_decile=_d[3][_j])
          for _d in DEC for _j, _h in enumerate(dhdr)],
         note="one row per feature per occupancy decile")


# ---------------------------------------------------------------- e: where n collapses
ax = fig.add_subplot(gs[1, 1]); st.panel_label(ax, "e", dx=-0.26, dy=1.07)
xs = np.arange(len(CUTS[2][1]))
for cut, rows, mono in CUTS:
    # rust = breaks monotonicity in Fig. 6e; blue = the 60% cut the paper uses;
    # slate = the other monotonic cuts, thinner. No transparency, so every pixel is palette.
    c = st.ORANGE if not mono else (st.BLUE if cut == 60 else st.GREY)
    lw_ = 1.3 if cut == 60 else (1.1 if not mono else 0.8)
    ax.plot(xs, [r[1] for r in rows], "-o", color=c, lw=lw_, ms=2.9 if cut == 60 else 2.3, mew=0)
    ax.text(len(rows) - 0.9, rows[-1][1], "  %d%%" % cut, fontsize=4.9, color=c, va="center",
            fontweight="bold" if (not mono or cut == 60) else "normal")
_NREF = 50                       # the reference line this panel draws, named once
ax.axhline(_NREF, color=st.GREY, lw=0.9, ls=(0, (3, 1.6)))
ax.text(0.02, _NREF + 2, "%d sites" % _NREF, fontsize=4.9, color=st.GREY, va="bottom")
ax.set_yscale("log")
ax.set_xticks(xs); ax.set_xticklabels([r[0].replace(">=", "≥") for r in CUTS[2][1]], fontsize=5.4)
ax.set_xlabel("compounds that engage the site", fontsize=6.0)
ax.set_ylabel("sites in the group", fontsize=6.0)
ax.set_xlim(-0.3, len(xs) + 0.7)
# Open a clear decade below the data so the caption has its own space instead of lying across
# the 70% and 80% end labels, which fall to the bottom right of the log axis.
_lo, _hi = ax.get_ylim()
ax.set_ylim(_lo / 9.0, _hi * 1.10)
HEAD["e_n_cuts"] = int(len(CUTS))
HEAD["e_n_cuts_not_monotonic"] = int(sum(1 for _c in CUTS if not _c[2]))
HEAD["e_cuts_not_monotonic"] = [int(_c[0]) for _c in CUTS if not _c[2]]
HEAD["e_smallest_group_n"] = int(min(_r[1] for _c in CUTS for _r in _c[1]))
fd.panel("figS7", "e",
         [dict(occupancy_cut_pct=_cut, compounds_engaging=_r[0], n_sites=_r[1],
               distance=_r[2], null=_r[3], cut_is_monotonic=_mono)
          for _cut, _rows, _mono in CUTS for _r in _rows],
         note="group sizes behind each of the five occupancy cuts in Figure 6e")


# ---------------------------------------------------------------- f: what the null cannot see
# This panel used to be a table of chips with two rows reading "nothing tested" — an assertion in
# the shape of a figure. Both claims are now measured (step 42) and plotted.
CHG = pd.read_csv(find("step42_charge_information.csv"))
CHI = pd.read_csv(find("step42_chirality.csv"))
NPR = pd.read_csv(find("step42_npr_spread.csv"))

fsub = gs[1, 2].subgridspec(2, 1, hspace=0.85, height_ratios=[1.25, 1.0])
ax = fig.add_subplot(fsub[0, 0]); st.panel_label(ax, "f", dx=-0.26, dy=1.18)
C = CHG.sort_values("rho_mw")
yy = np.arange(len(C))
ax.barh(yy - 0.19, C.rho_mw, height=0.36, color=st.LIGHTGREY)
ax.barh(yy + 0.19, C.rho_charge, height=0.36, color=st.ORANGE)
ax.set_yticks(yy); ax.set_yticklabels([d.replace("_", " ") for d in C.descriptor], fontsize=4.6)
ax.set_xlabel("|ρ| with molecular weight (grey)\nand with net charge (rust)", fontsize=5.6,
              linespacing=1.4)
ax.set_xlim(0, 1.04); ax.set_ylim(-0.8, len(C) - 0.2)
ax.tick_params(labelsize=4.6)
for s_ in ("top", "right"): ax.spines[s_].set_visible(False)

ax = fig.add_subplot(fsub[1, 0])
GRP = [("reflected", CHI[CHI.kind == "reflected"].dNPR1.values, st.BLUE),
       ("re-embedded", CHI[CHI.kind == "re-embedded"].dNPR1.values, st.GREY),
       ("between compounds", np.abs(NPR.NPR1.values - NPR.NPR1.mean()), st.ORANGE)]
rgc = np.random.default_rng(3)
for i_, (lab, v, c) in enumerate(GRP):
    ax.plot(v, i_ + rgc.normal(0, 0.10, len(v)), "o", ms=2.2, mew=0, color=c)
    ax.plot([np.median(v)] * 2, [i_ - 0.28, i_ + 0.28], color=c, lw=1.6, solid_capstyle="round")
_REFMAX = float(CHI[CHI.kind == "reflected"].dNPR1.max())
ax.text(0.013, 0.0, "exactly %g" % _REFMAX, fontsize=5.0, color=st.BLUE, fontweight="bold",
        va="center", ha="left")
ax.set_yticks(range(3)); ax.set_yticklabels([g[0] for g in GRP], fontsize=4.8)
ax.set_xlabel("change in NPR1", fontsize=5.6)
ax.set_xlim(-0.012, 0.33); ax.set_ylim(-0.7, 2.7)
ax.tick_params(labelsize=4.6)
for s_ in ("top", "right"): ax.spines[s_].set_visible(False)
HEAD["f_n_descriptors"] = int(len(CHG))
HEAD["f_charge_max_abs_rho"] = round(float(CHG.rho_charge.max()), 3)
HEAD["f_mw_median_abs_rho"] = round(float(CHG.rho_mw.median()), 3)
HEAD["f_chirality_reflected_max_dnpr1"] = round(_REFMAX, 3)
HEAD["f_chirality_reembedded_max_dnpr1"] = round(
    float(CHI[CHI.kind == "re-embedded"].dNPR1.max()), 3)
HEAD["f_between_compound_npr1_sd"] = round(float(NPR.NPR1.std()), 3)
fd.panel("figS7", "f",
         [dict(part="charge", label=_r.descriptor, measure="abs_rho_with_molecular_weight",
               value=float(_r.rho_mw)) for _r in CHG.itertuples()]
         + [dict(part="charge", label=_r.descriptor, measure="abs_rho_with_net_charge",
                 value=float(_r.rho_charge)) for _r in CHG.itertuples()]
         + [dict(part="chirality", label=_g, measure="change_in_NPR1", value=float(_v))
            for _g, _vals, _c in GRP for _v in _vals],
         note="upper axis: charge information in the tested descriptors; "
              "lower axis: NPR1 change under reflection, re-embedding and between compounds")

st.save(fig, "figS7")
fd.head("figS7", HEAD)
for _e in (".png", ".pdf"):
    fig.savefig("../results/figS7_v1" + _e, dpi=400, bbox_inches="tight")
print("a  step37 pooled: %d of %d pairs over %d compounds | by descriptor: %s"
      % (len(SIG), len(POOL), int(POOL.n.iloc[0]),
         SIG["shape"].value_counts().to_dict()))
print("   nulls: %d of %d within NatChem (%d compounds) | %d of %d across slopes (%d)"
      % (int((NATC.q < 0.05).sum()), len(NATC), int(NATC.n_compounds.iloc[0]),
         int((SLOP.q < 0.05).sum()), len(SLOP), int(SLOP.n.iloc[0])))
print("slopes %d (%d significant) | deciles %d | cuts %s"
      % (len(SLOPE), sum(1 for r in SLOPE if r[4]), len(DEC), [c[0] for c in CUTS]))
print("continuous measures:", {m: int((C40[C40.measure == m].q < 0.05).sum()) for m in MEAS})
print("f charge: max |rho| %.3f vs median |rho| MW %.3f"
      % (CHG.rho_charge.max(), CHG.rho_mw.median()))
print("f chirality: reflected max %.2e | re-embedded max %.3f | between-compound sd %.3f"
      % (CHI[CHI.kind == "reflected"].dNPR1.max(),
         CHI[CHI.kind == "re-embedded"].dNPR1.max(), NPR.NPR1.std()))
