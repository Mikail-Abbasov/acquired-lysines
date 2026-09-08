#!/usr/bin/env python3
"""
Figure 7 — The probe is not a neutral reporter: reagent charge steers labelling.

Answers the exception Figure 6 leaves open. Figure 6d shows no compound property predicts which
lysines are engaged, and figS7 f shows charge was never among the properties tested. Here charge
is tested directly, on the reagent rather than the compound.

WHAT CHANGED FROM THE PLANNED VERSION (2026-09-03)
The architecture planned a second half arguing "proximity to a concavity, not occupancy of a
pocket", with "41.6% of lysines within 5 A of a cavity, only 10.3% large enough to matter".
Three problems: those numbers appear in no step report and do not reproduce (40.4% ungated,
48.3% at pLDDT>=70; no volume threshold gives 10.3%); step 19's `in_pocket` IS
`dist_nearest_pocket <= 5`, so the two halves of the sentence are one variable; and step 36 had
already found fpocket's cavities here are surface depressions, not binding sites. Panel e
replaces it with enclosure (step 43), measured from the AlphaFold structure and independent of
both fpocket and solvent accessibility.

THE DECOMPOSITION (panels c, d; step 44)
A57 said A49's charge effect is "at least partly" reagent. Bootstrap CIs over proteins split it
three ways rather than hedging: cationic attraction is the REAGENT (the neutral probe's interval
brackets zero and excludes the anionic probe's estimate); acidic avoidance, predicted pKa and
salt bridges are the SITE (no probe difference); net local charge is MIXED, and for that one
"at least partly" is still the right words.

NOT USED: the inference that cationic is all reagent because STP-vs-unlabelled is significant
while NASA-vs-unlabelled is not. A difference between significant and non-significant is not a
significant difference. Panels c and d show intervals so the reader can see the actual basis.
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
        if os.path.exists(p): return p
    raise FileNotFoundError(n)

R22 = open(find("step22_report.txt")).read()
T22 = pd.read_csv(find("step22_two_probe.csv"))
D44 = pd.read_csv(find("step44_charge_decomposition.csv"))
P43 = pd.read_csv(find("step43b_probe_contrast.csv"))
S43 = pd.read_csv(find("step43b_site_vs_unlabelled.csv"))

# step 22's own feature families, so the figure cannot drift from the analysis
CHARGE = ["anionic_8A","cationic_8A","net_charge_8A","propka_pKa","salt_bridge",
          "hbond_partners_NZ","n_lys_8A"]
POCKET = ["in_pocket","dist_nearest_pocket","pocket_volume","pocket_score","pocket_drugscore"]
SHAPE  = ["crowding_8A","atoms_10A","hydrophobic_8A","aromatic_8A","polar_8A","rsa_tien",
          "sasa_sidechain","nz_sasa","dist_centroid"]
FAM = {**{f:"charge" for f in CHARGE}, **{f:"pocket" for f in POCKET},
       **{f:"shape" for f in SHAPE}}
FCOL = {"charge": st.ORANGE, "pocket": st.BLUE, "shape": st.GREY, "context": st.LIGHTGREY}

NICE = lambda f: (f.replace("_8A"," 8Å").replace("_10A"," 10Å").replace("_"," ")
                   .replace("rsa tien","exposure").replace("propka pKa","predicted pKa")
                   .replace("nz sasa","NZ exposure").replace("sasa sidechain","side-chain exp.")
                   .replace("dist centroid","dist. to centroid").replace("n lys","nearby lysines")
                   .replace("dist nearest pocket","dist. to pocket")
                   .replace("hbond partners NZ","H-bonds at NZ"))

NSTP, NNASA, NBOTH = [int(x) for x in
    re.search(r"STP-only\s+(\d+)\s+\|\s+NASA-only\s+(\d+)\s+\|\s+both\s+(\d+)", R22).groups()]
CON = T22[T22.group == "STP_ONLY vs NASA_ONLY"].copy()
CON["fam"] = CON.feature.map(lambda f: FAM.get(f, "context"))
NSIG = int((CON.q < 0.05).sum())
NPROT = int(CON.n_proteins.iloc[0])
assert NSIG == 2, NSIG

HEAD = {"a_n_stp_only": NSTP, "a_n_nasa_only": NNASA, "a_n_both": NBOTH,
        "a_n_proteins_with_both": NPROT,
        "b_n_significant": NSIG, "b_n_features": int(len(CON)),
        "b_n_pocket_significant": int((CON[CON.fam == "pocket"].q < 0.05).sum()),
        "b_n_pocket": int((CON.fam == "pocket").sum()),
        "b_n_shape_significant": int((CON[CON.fam == "shape"].q < 0.05).sum()),
        "b_n_shape": int((CON.fam == "shape").sum())}

fd.panel("fig7", "b", CON[["feature", "fam", "d", "q", "n_proteins"]],
         note="all 22 site features, anionic-only vs neutral-only sites, within protein")
for _, _r in CON[CON.q < 0.05].iterrows():
    HEAD["b_%s_d" % _r.feature] = round(float(_r.d), 3)
    HEAD["b_%s_q" % _r.feature] = "%.1e" % float(_r.q)
# the next-largest difference is quoted in 2.7 as the one that does not clear the threshold
_ns = CON[CON.q >= 0.05].sort_values("d", ascending=False)
if len(_ns):
    HEAD["b_next_largest_feature"] = _ns.feature.iloc[0]
    HEAD["b_next_largest_d"] = round(float(_ns.d.iloc[0]), 3)

fd.panel("fig7", "cd", D44, note="bootstrap decomposition: each probe against unlabelled "
                                 "lysines, and the probe-minus-probe difference")

fd.panel("fig7", "e", pd.concat([S43.assign(contrast="labelled vs unlabelled"),
                                 P43.assign(contrast="anionic vs neutral probe")],
                                ignore_index=True),
         note="enclosure and pocket occupancy, site contrast against probe contrast")
for _lab, _df in (("site", S43), ("probe", P43)):
    _r = _df[_df.feature == "enclosure"]
    if len(_r):
        HEAD["e_enclosure_%s_d" % _lab] = round(float(_r.d.iloc[0]), 3)
        HEAD["e_enclosure_%s_q" % _lab] = ("%.1e" % float(_r.q.iloc[0])
                                           if float(_r.q.iloc[0]) < 0.01
                                           else round(float(_r.q.iloc[0]), 2))

fig = plt.figure(figsize=(st.DOUBLE, st.DOUBLE * 0.60))
gt = fig.add_gridspec(1, 2, left=0.075, right=0.980, top=0.945, bottom=0.585,
                      wspace=0.40, width_ratios=[1.05, 1.0])
gb = fig.add_gridspec(1, 3, left=0.075, right=0.980, top=0.455, bottom=0.085, wspace=0.52)

# ================================================================ a  the design
ax = fig.add_subplot(gt[0, 0]); st.panel_label(ax, "a", dx=-0.09, dy=1.08)
ax.set_xlim(0, 100); ax.set_ylim(0, 10); ax.set_axis_off()
for k, (nm, sub, c) in enumerate([("dTB-STP", "anionic leaving group", st.ORANGE),
                                  ("NASA", "neutral leaving group", st.GREY)]):
    ax.add_patch(Rectangle((2.0, 7.4 - k * 2.5), 26.0, 1.5, facecolor=c, lw=0))
    ax.text(15.0, 8.15 - k * 2.5, nm, fontsize=5.6, color=st.PAPER, fontweight="bold",
            ha="center", va="center")
    ax.text(2.0, 6.95 - k * 2.5, sub, fontsize=4.7, color=c, va="top")
ax.text(2.0, 3.35, "same proteome, same cells", fontsize=4.9, color=st.GREY)
ax.annotate("", xy=(38.0, 6.2), xytext=(30.0, 6.2),
            arrowprops=dict(arrowstyle="-|>", color=st.LIGHTGREY, lw=1.1))
BAR = [("STP-only", NSTP, st.ORANGE), ("both", NBOTH, st.LIGHTGREY), ("NASA-only", NNASA, st.GREY)]
tot = sum(b[1] for b in BAR); x0 = 41.0
for nm, n, c in BAR:
    w = 56.0 * n / tot
    ax.add_patch(Rectangle((x0, 5.5), w, 1.6, facecolor=c, lw=0))
    ax.text(x0 + w/2, 7.35, "%s" % f"{n:,}", fontsize=5.0, color=c, ha="center", fontweight="bold")
    ax.text(x0 + w/2, 5.15, nm, fontsize=4.7, color=c, ha="center", va="top")
    x0 += w
ax.annotate("", xy=(69.0, 3.55), xytext=(69.0, 4.55),
            arrowprops=dict(arrowstyle="-|>", color=st.GREY, lw=0.9))
ax.text(41.0, 2.95, "compared within protein", fontsize=4.9, color=st.GREY, va="top")
ax.text(41.0, 2.05, "%d proteins carry both" % NPROT, fontsize=4.9, color=st.BLUE,
        fontweight="bold", va="top")

# ================================================================ b  all 22 features
ax = fig.add_subplot(gt[0, 1]); st.panel_label(ax, "b", dx=-0.40, dy=1.08)
C = CON.sort_values("d")
y = np.arange(len(C))
for i, (_, r) in enumerate(C.iterrows()):
    sig = r.q < 0.05
    c = FCOL[r.fam] if sig else st.LIGHTGREY
    ax.plot([r.d], [i], "o", color=c, ms=4.2 if sig else 2.6, mew=0)
ax.axvline(0, color=st.GREY, lw=0.8)
ax.set_yticks(y)
ax.set_yticklabels([NICE(f) for f in C.feature], fontsize=4.5)
for t, f in zip(ax.get_yticklabels(), C.feature):
    if float(C.loc[C.feature == f, "q"].iloc[0]) < 0.05:
        t.set_color(st.ORANGE); t.set_fontweight("bold")
ax.set_xlabel("higher in the anionic probe's sites (d)", fontsize=6.0)
ax.set_xlim(-0.16, 0.30); ax.set_ylim(-1.6, len(C) - 0.4)
ax.text(0.02, 0.02, "%d of %d differ\nboth are charge" % (NSIG, len(C)), transform=ax.transAxes,
        fontsize=5.4, color=st.ORANGE, fontweight="bold", va="bottom", linespacing=1.4)
ax.text(0.98, 0.02, "0 of 5 pocket\n0 of 9 shape", transform=ax.transAxes, ha="right",
        va="bottom", fontsize=5.2, color=st.GREY, fontweight="bold", linespacing=1.4)
for s in ("top", "right"): ax.spines[s].set_visible(False)

# ================================================================ c  what each probe prefers
ORDER = ["cationic_8A", "net_charge_8A", "anionic_8A", "propka_pKa", "salt_bridge"]
ax = fig.add_subplot(gb[0, 0]); st.panel_label(ax, "c", dx=-0.30, dy=1.09)
CMPS = [("STP-only vs unlabelled", st.ORANGE, "o"), ("NASA-only vs unlabelled", st.GREY, "s")]
for i, f in enumerate(ORDER):
    for j, (cmp_, c, mk) in enumerate(CMPS):
        r = D44[(D44.feature == f) & (D44.comparison == cmp_)]
        if not len(r): continue
        r = r.iloc[0]; yy = i + (0.17 if j == 0 else -0.17)
        ax.plot([r.lo, r.hi], [yy, yy], color=c, lw=1.1, solid_capstyle="round")
        ax.plot([r.d], [yy], mk, color=c, ms=3.2, mew=0)
ax.axvline(0, color=st.GREY, lw=0.8)
ax.set_yticks(range(len(ORDER))); ax.set_yticklabels([NICE(f) for f in ORDER], fontsize=5.0)
ax.set_xlabel("preference vs unlabelled lysines (d)", fontsize=6.0)
# net charge's anionic CI reaches +0.45; a tighter limit clipped its marker off the panel.
_lo = min(D44[D44.comparison.str.contains("unlabelled")].lo) - 0.04
_hi = max(D44[D44.comparison.str.contains("unlabelled")].hi) + 0.04
ax.set_ylim(-1.15, len(ORDER) - 0.4); ax.set_xlim(_lo, _hi)
for j, (nm, c, mk) in enumerate([("anionic probe", st.ORANGE, "o"),
                                 ("neutral probe", st.GREY, "s")]):
    yy = -0.62 - j * 0.34
    ax.plot([_lo + 0.03], [yy], mk, color=c, ms=3.2, mew=0, clip_on=False)
    ax.text(_lo + 0.055, yy, nm, fontsize=5.0, color=c, va="center")
for s in ("top", "right"): ax.spines[s].set_visible(False)

# ================================================================ d  the difference
ax = fig.add_subplot(gb[0, 1]); st.panel_label(ax, "d", dx=-0.30, dy=1.09)
VERD = {"cationic_8A": ("reagent", st.ORANGE), "net_charge_8A": ("mixed", st.GREEN),
        "anionic_8A": ("site", st.BLUE), "propka_pKa": ("site", st.BLUE),
        "salt_bridge": ("site", st.BLUE)}
for i, f in enumerate(ORDER):
    r = D44[(D44.feature == f) & (D44.comparison == "STP-only vs NASA-only")]
    if not len(r): continue
    r = r.iloc[0]; lab, c = VERD[f]
    ax.plot([r.lo, r.hi], [i, i], color=c, lw=1.3, solid_capstyle="round")
    ax.plot([r.d], [i], "o", color=c, ms=3.8, mew=0)
    ax.text(r.hi + 0.018, i, lab, fontsize=4.8, color=c, va="center", fontweight="bold")
ax.axvline(0, color=st.GREY, lw=0.8)
ax.set_yticks(range(len(ORDER))); ax.set_yticklabels([NICE(f) for f in ORDER], fontsize=5.0)
ax.set_xlabel("anionic minus neutral probe (d)", fontsize=6.0)
ax.set_ylim(-0.7, len(ORDER) - 0.4); ax.set_xlim(-0.22, 0.52)
for s in ("top", "right"): ax.spines[s].set_visible(False)

# ================================================================ e  concavity is the site's
ax = fig.add_subplot(gb[0, 2]); st.panel_label(ax, "e", dx=-0.30, dy=1.09)
ROWS = [("labelled vs\nunlabelled", S43, st.BLUE), ("anionic vs\nneutral probe", P43, st.ORANGE)]
FE = ["enclosure", "enclosure_adj", "in_pocket"]
FL = {"enclosure": "enclosure", "enclosure_adj": "enclosure,\naccessibility removed",
      "in_pocket": "in a pocket"}
for i, f in enumerate(FE):
    for j, (nm, df, c) in enumerate(ROWS):
        r = df[df.feature == f]
        if not len(r): continue
        r = r.iloc[0]; yy = i + (0.18 if j == 0 else -0.18)
        sig = r.q < 0.05
        ax.plot([r.d], [yy], "o" if j == 0 else "s", color=c if sig else st.LIGHTGREY,
                ms=4.0 if sig else 3.0, mew=0)
ax.axvline(0, color=st.GREY, lw=0.8)
ax.set_yticks(range(len(FE))); ax.set_yticklabels([FL[f] for f in FE], fontsize=4.8,
                                                  linespacing=1.35)
ax.set_xlabel("effect size (d)", fontsize=6.0)
ax.set_xlim(-0.13, 0.20); ax.set_ylim(-1.15, len(FE) - 0.4)
for j, (nm, c, mk) in enumerate([("labelled vs unlabelled", st.BLUE, "o"),
                                 ("probe vs probe", st.ORANGE, "s")]):
    yy = -0.62 - j * 0.34
    ax.plot([-0.115], [yy], mk, color=c, ms=3.2, mew=0, clip_on=False)
    ax.text(-0.103, yy, nm, fontsize=5.0, color=c, va="center")
ax.text(0.98, 0.03, "pale = not significant", transform=ax.transAxes, ha="right", va="bottom",
        fontsize=4.8, color=st.GREY)
for s in ("top", "right"): ax.spines[s].set_visible(False)

st.save(fig, "fig7")
fd.head("fig7", HEAD)
print("a: STP-only %d | both %d | NASA-only %d | %d proteins" % (NSTP, NBOTH, NNASA, NPROT))
print("b: %d of %d features differ; families %s"
      % (NSIG, len(CON), CON[CON.q < 0.05].fam.tolist()))
print("c/d: %s" % {f: (VERD[f][0]) for f in ORDER})
print("e: enclosure labelled-vs-unlab d %.3f (q %.3g) | probe-vs-probe d %.3f (q %.2f)"
      % (S43[S43.feature=="enclosure"].d.iloc[0], S43[S43.feature=="enclosure"].q.iloc[0],
         P43[P43.feature=="enclosure"].d.iloc[0], P43[P43.feature=="enclosure"].q.iloc[0]))
