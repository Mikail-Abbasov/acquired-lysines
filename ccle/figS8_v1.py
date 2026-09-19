#!/usr/bin/env python3
"""
Figure S8 — support for Figure 7.

  a  Each probe against unlabelled lysines, all 22 features. Figure 7b shows only the
     difference between probes; this shows what each one prefers, so a reader can see that the
     two probes agree almost everywhere and part company on one axis.
  b  The enclosure measure is a distinct axis: mostly solvent accessibility, but a fifth of it
     is pocket geometry that accessibility does not carry.
  c  POOLING REVERSES THE SIGN. Site-level pooling makes labelled lysines look LESS enclosed
     (d = -0.03); within protein they are MORE enclosed (d = +0.073). Same data. This is the
     hazard the within-protein and within-compound designs exist to avoid, and it is plotted
     rather than described because it nearly reached a figure.
  d  What bounds the claim: if NASA is substantially ionised at pH 7.4 then both probes are
     anionic, and the measured contrast is a lower bound on the true reagent effect.
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
        if os.path.exists(p): return p
    raise FileNotFoundError(n)

T22 = pd.read_csv(find("step22_two_probe.csv"))
R43 = open(find("step43b_report.txt")).read()
S43 = pd.read_csv(find("step43b_site_vs_unlabelled.csv"))
E = pd.read_csv(find("step43_enclosure.csv"))
P = pd.read_csv(find("step19_pockets.csv"),
                usecols=["uniprot","position","rsa_tien","in_pocket","plddt_site","labelled"])
M = E.merge(P, on=["uniprot","position"]); M = M[M.plddt_site >= 70].dropna(subset=["enclosure","rsa_tien"])

CHARGE = ["anionic_8A","cationic_8A","net_charge_8A","propka_pKa","salt_bridge",
          "hbond_partners_NZ","n_lys_8A"]
POCKET = ["in_pocket","dist_nearest_pocket","pocket_volume","pocket_score","pocket_drugscore"]
FAM = lambda f: "charge" if f in CHARGE else ("pocket" if f in POCKET else "shape")
NICE = lambda f: (f.replace("_8A"," 8Å").replace("_10A"," 10Å").replace("_"," ")
                   .replace("rsa tien","exposure").replace("propka pKa","predicted pKa")
                   .replace("nz sasa","NZ exposure").replace("sasa sidechain","side-chain exp.")
                   .replace("dist centroid","dist. to centroid").replace("n lys","nearby lysines")
                   .replace("dist nearest pocket","dist. to pocket"))

fig = plt.figure(figsize=(st.DOUBLE, st.DOUBLE * 0.52))
gs = fig.add_gridspec(1, 4, left=0.085, right=0.980, top=0.905, bottom=0.150, wspace=0.62,
                      width_ratios=[1.30, 1.0, 1.0, 1.0])

# ---------------------------------------------------------------- a
ax = fig.add_subplot(gs[0, 0]); st.panel_label(ax, "a", dx=-0.52, dy=1.06)
S = T22[T22.group == "STP_ONLY"].set_index("feature")
N = T22[T22.group == "NASA_ONLY"].set_index("feature")
feats = list(S.d.sort_values().index)
for i, f in enumerate(feats):
    for df, c, mk, off in ((S, st.ORANGE, "o", 0.19), (N, st.GREY, "s", -0.19)):
        if f not in df.index: continue
        r = df.loc[f]
        ax.plot([r.d], [i + off], mk, color=c if r.q < 0.05 else st.LIGHTGREY,
                ms=3.4 if r.q < 0.05 else 2.4, mew=0)
ax.axvline(0, color=st.GREY, lw=0.8)
_A = []
for i, f in enumerate(feats):
    for _nm, df in (("anionic (STP)", S), ("neutral (NASA)", N)):
        if f not in df.index: continue
        r = df.loc[f]
        _A.append(dict(feature=f, label=NICE(f), probe=_nm, family=FAM(f), row=i,
                       d=float(r.d), q=float(r.q), n_proteins=int(r.n_proteins),
                       significant=bool(r.q < 0.05)))
fd.panel("figS8", "a", _A,
         note="each probe against unlabelled lysines, within protein; all 22 site features")
ax.set_yticks(range(len(feats))); ax.set_yticklabels([NICE(f) for f in feats], fontsize=4.4)
ax.set_xlabel("preference vs unlabelled (d)", fontsize=6.0)
ax.set_ylim(-2.6, len(feats) - 0.4); ax.set_xlim(-0.46, 0.50)
for j, (nm, c, mk) in enumerate([("anionic probe", st.ORANGE, "o"), ("neutral probe", st.GREY, "s")]):
    yy = -1.35 - j * 0.62
    ax.plot([-0.43], [yy], mk, color=c, ms=3.2, mew=0, clip_on=False)
    ax.text(-0.40, yy, nm, fontsize=4.9, color=c, va="center")
for s_ in ("top", "right"): ax.spines[s_].set_visible(False)

# ---------------------------------------------------------------- b
ax = fig.add_subplot(gs[0, 1]); st.panel_label(ax, "b", dx=-0.34, dy=1.06)
sub = M.sample(min(6000, len(M)), random_state=0)
ax.plot(sub.rsa_tien, sub.enclosure, "o", ms=1.0, mew=0, color=st.LIGHTGREY)
from scipy import stats as sstats
b = np.polyfit(M.rsa_tien, M.enclosure, 1)
# clip the fit to the plotted range: extended, the line leaves the 0-1 range of the measure
xs = np.linspace(M.rsa_tien.min(), M.rsa_tien.quantile(0.995), 50)
ax.plot(xs, np.polyval(b, xs), color=st.BLUE, lw=1.3)
rho = sstats.spearmanr(M.enclosure, M.rsa_tien)[0]
fd.panel("figS8", "b", sub[["uniprot", "position", "rsa_tien", "enclosure", "labelled"]],
         note=("the %d points drawn, a seed-0 subsample of %d lysines at pLDDT >= 70; "
               "Spearman rho %.4f and fit y = %.4f x + %.4f are computed on all %d"
               % (len(sub), len(M), rho, b[0], b[1], len(M))))
ax.set_xlabel("relative solvent accessibility", fontsize=6.0)
ax.set_ylabel("enclosure", fontsize=6.0)
ax.set_xlim(0, M.rsa_tien.quantile(0.999)); ax.set_ylim(0, 1.02)
ax.text(0.97, 0.95, "ρ = %.2f\n%.0f%% shared" % (rho, 100 * rho ** 2), transform=ax.transAxes,
        ha="right", va="top", fontsize=5.2, color=st.BLUE, fontweight="bold", linespacing=1.4)
for s_ in ("top", "right"): ax.spines[s_].set_visible(False)

# ---------------------------------------------------------------- c
ax = fig.add_subplot(gs[0, 2]); st.panel_label(ax, "c", dx=-0.34, dy=1.06)
# The pooled value was typed in as -0.032 until 2026-09-07. It is step 43b's pooled
# site-level enclosure effect and every input to it is already loaded here, so it is
# computed rather than transcribed: the panel and the prose cannot now drift apart.
_lab = M.loc[M.labelled.astype(bool), "enclosure"].values
_unl = M.loc[~M.labelled.astype(bool), "enclosure"].values
_sp = np.sqrt(((len(_lab) - 1) * _lab.var(ddof=1) + (len(_unl) - 1) * _unl.var(ddof=1))
              / (len(_lab) + len(_unl) - 2))
POOL_D = float((_lab.mean() - _unl.mean()) / _sp)
WITH_D = float(S43[S43.feature == "enclosure"].d.iloc[0])
fd.panel("figS8", "c",
         [dict(design="pooled across sites", feature="enclosure", d=POOL_D,
               n_labelled=int(len(_lab)), n_unlabelled=int(len(_unl)), n_proteins=""),
          dict(design="within protein", feature="enclosure", d=WITH_D,
               n_labelled="", n_unlabelled="",
               n_proteins=int(S43[S43.feature == "enclosure"].n_proteins.iloc[0]))],
         note="the same data, pooled and within protein; the sign reverses")
for i, (nm, v, c) in enumerate([("pooled\nacross sites", POOL_D, st.ORANGE),
                                ("within\nprotein", WITH_D, st.BLUE)]):
    ax.barh([i], [v], height=0.46, color=c)
    ax.text(v + (0.006 if v > 0 else -0.006), i, "%+.3f" % v, va="center",
            ha="left" if v > 0 else "right", fontsize=5.6, color=c, fontweight="bold")
ax.axvline(0, color=st.GREY, lw=0.9)
ax.set_yticks([0, 1]); ax.set_yticklabels(["pooled\nacross sites", "within\nprotein"],
                                          fontsize=5.2, linespacing=1.35)
ax.set_xlabel("enclosure, labelled vs unlabelled (d)", fontsize=6.0)
ax.set_xlim(-0.075, 0.115); ax.set_ylim(-0.7, 1.7)
for s_ in ("top", "right", "left"): ax.spines[s_].set_visible(False)
ax.tick_params(axis="y", length=0)

# ---------------------------------------------------------------- d
ax = fig.add_subplot(gs[0, 3]); st.panel_label(ax, "d", dx=-0.30, dy=1.06)
# OBS was typed in as 0.233. It is step 22's probe contrast on cationic_8A, already
# loaded above, so it is read from the table the panel is about.
_CON = T22[T22.group == "STP_ONLY vs NASA_ONLY"].set_index("feature")
OBS = float(_CON.loc["cationic_8A", "d"])
ION = 0.46   # the substantially-ionised endpoint: an assumed limit, not a measurement
for i, (nm, lo, hi, c) in enumerate([("NASA neutral\n(assumed)", OBS, OBS, st.BLUE),
                                     ("NASA partly\nionised", OBS, ION, st.ORANGE)]):
    if lo == hi:
        ax.plot([lo], [i], "o", color=c, ms=4.4, mew=0)
    else:
        ax.annotate("", xy=(hi, i), xytext=(lo, i),
                    arrowprops=dict(arrowstyle="-|>", color=c, lw=1.4))
        ax.plot([lo], [i], "o", color=c, ms=4.0, mew=0)
ax.axvline(0, color=st.GREY, lw=0.8)
ax.set_yticks([0, 1]); ax.set_yticklabels(["NASA neutral\n(assumed)", "NASA partly\nionised"],
                                          fontsize=5.2, linespacing=1.35)
ax.set_xlabel("true reagent effect,\ncationic 8Å (d)", fontsize=6.0, linespacing=1.4)
ax.set_xlim(0, 0.52); ax.set_ylim(-0.7, 1.7)
for s_ in ("top", "right", "left"): ax.spines[s_].set_visible(False)
ax.tick_params(axis="y", length=0)
fd.panel("figS8", "d",
         [dict(assumption="NASA neutral (assumed)", quantity="cationic_8A d", lo=OBS, hi=OBS,
               kind="measured probe contrast, step 22"),
          dict(assumption="NASA partly ionised", quantity="cationic_8A d", lo=OBS, hi=ION,
               kind="sensitivity range across assumed ionisation states, not a confidence interval")],
         note="a sensitivity range, not a confidence interval; the upper endpoint is assumed")

# ---------------------------------------------------------------- ledger
_NSITE = re.search(r"STP-only\s+(\d+)\s*\|\s*NASA-only\s+(\d+)\s+sites", R43)
assert _NSITE, "step43b report no longer states the per-probe site counts"
HEAD = {
    "a_n_features": int(len(feats)),
    "a_n_significant_anionic": int((S.q < 0.05).sum()),
    "a_n_significant_neutral": int((N.q < 0.05).sum()),
    "a_n_proteins_anionic": int(S.n_proteins.max()),
    "a_n_proteins_neutral": int(N.n_proteins.max()),
    "a_n_sites_anionic_only": int(_NSITE.group(1)),
    "a_n_sites_neutral_only": int(_NSITE.group(2)),
    "b_spearman_rho": round(float(rho), 3),
    "b_variance_shared_pct": round(100 * float(rho) ** 2, 1),
    "b_n_sites": int(len(M)),
    "c_enclosure_d_pooled": round(POOL_D, 3),
    "c_enclosure_d_within_protein": round(WITH_D, 3),
    "d_cationic_d_neutral_assumption": round(OBS, 3),
    "d_cationic_d_partly_ionised_bound": ION,
}

st.save(fig, "figS8")
fd.head("figS8", HEAD)
print("a: %d features per probe | STP sig %d, NASA sig %d"
      % (len(feats), (S.q < 0.05).sum(), (N.q < 0.05).sum()))
print("b: rho %.3f (%.0f%% shared) on %s sites" % (rho, 100*rho**2, f"{len(M):,}"))
print("c: pooled %+.3f vs within protein %+.3f" % (POOL_D, WITH_D))
