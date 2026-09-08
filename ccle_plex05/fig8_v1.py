#!/usr/bin/env python3
"""
fig8_v1.py — FIGURE 8. Acquired lysines sit in ordinary lysine environments.

THE CLAIM, worded to what the data support
Not "indistinguishable" -- the residual separation is real and survives every check. What the
data support is: acquired lysines occupy the same structural space as engaged lysines; 82% fall
inside the region containing 90% of engaged lysines against 89% for measured lysines that were
labelled but never engaged; they sit closer to engaged lysines than to labelled-but-unengaged
ones; and the residual excess is 1.2x the between-population yardstick once each protein votes
once. The barrier is detection, not chemistry.

NO SENTENCES INSIDE PANELS. Per the 2026-09-02 rule adopted for Figure 6 -- a figure is looked
at, not read. Everything explanatory belongs in the legend.

PANELS
  a  the anchor: why both sides must be measured from the same atom     step21 / step24
  b  the three populations, feature by feature                          step28c
  c  overlap -- do acquired lysines fall inside the engaged region      step28
  d  real crystal structures: does the new lysine move the site         step31
  e  AKT1 E17K -- drugged in vivo, invisible here. THE CLOSING IMAGE.   acquired_K_searched...

ROBUSTNESS (the same comparison under three anchorings, plus the artefact) is figS9.

INPUT   ligand_rebuild/step28c_medians.csv, step28_results.csv, step31_crystal.csv,
        ccle_plex05/acquired_K_searched_vs_detected.csv
OUTPUT  fig8.pdf, fig8.png
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import abbasov_style as st
import figure_data as fd
st.apply()
ORANGE, BLUE, GREY, LIGHT, INK = st.ORANGE, st.BLUE, st.GREEN, st.LIGHTGREY, st.GREY

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
LR = os.path.join(ROOT, "ligand_rebuild")
MED = pd.read_csv(os.path.join(LR, "step28c_medians.csv"))
RES = pd.read_csv(os.path.join(LR, "step28_results.csv"))
XT  = pd.read_csv(os.path.join(LR, "step31_crystal.csv"))
AK  = pd.read_csv(os.path.join(HERE, "acquired_K_searched_vs_detected.csv"))
AK  = AK[AK.gene == "AKT1"]

armA = RES[RES.arm.str.startswith("A")].iloc[0]

# P(a random acquired lysine sits closer to the engaged centre than a random engaged one).
# This was typed into panel c as "P = 0.443". It is computed by step28d and is now read from
# that report, so the artwork cannot drift from the analysis.
import re as _re
_r28d = open(os.path.join(LR, "step28d_report.txt")).read()
P_CLOSER = float(_re.search(r"=\s*([\d.]+)\s+\(0\.5 would mean", _r28d).group(1))

# AKT1 E17K: how many peptides did the search actually contain, and in how many plexes?
# The manuscript said "seven peptides in the 42 experiments". Counted from the search
# databases: the mutation is carried by one cell line, so it is searchable in one plex.
import glob as _glob
_akt = {}
for _f in sorted(_glob.glob(os.path.join(HERE, "plex*_varpep.fasta")) +
                 _glob.glob(os.path.join(HERE, "plex*_downpep.fasta"))):
    _plex = os.path.basename(_f).split("_")[0]
    _hdr = None
    for _ln in open(_f):
        if _ln.startswith(">"):
            _p = _ln[1:].strip().split("|")
            _hdr = _p if (len(_p) > 5 and _p[0] in ("VAR", "DOWN")
                          and _p[1] == "AKT1" and _p[5] == "17") else None
        elif _hdr:
            _akt.setdefault(_plex, {})[_ln.strip()] = _hdr[0]
            _hdr = None
AKT_PEPS = sorted({q for v in _akt.values() for q in v})
AKT_PLEXES = sorted(_akt)
AKT_VAR = sum(1 for v in _akt.values() for q, k in v.items() if k == "VAR")
AKT_DOWN = sum(1 for v in _akt.values() for q, k in v.items() if k == "DOWN")
# If the search databases are not present -- a deposit without plex*_{varpep,downpep}.fasta --
# these counts come out as zero, and zero here reads as "the site was never searched", which is
# the opposite of what the paper claims and is the number the caption quotes. Refuse to emit a
# count rather than emit a wrong one.
AKT_OK = bool(_akt)
if not AKT_OK:
    print("  WARNING fig8e: no plex*_varpep.fasta / plex*_downpep.fasta found, so the AKT1 E17K "
          "peptide count cannot be computed.\n"
          "          The deposited value is 3 peptides (2 variant, 1 downstream) in plex36. "
          "No count is written to the ledger.")

fig = plt.figure(figsize=(st.DOUBLE, st.DOUBLE * 0.60))
gs = fig.add_gridspec(2, 6, hspace=0.46, wspace=1.05,
                      left=0.075, right=0.985, top=0.95, bottom=0.105,
                      height_ratios=[1.0, 0.92])

# ---------------------------------------------------------------- a  the anchor
ax = fig.add_subplot(gs[0, 0:2]); st.panel_label(ax, "a", dx=-0.30, dy=1.09)
ax.set_title("One reference atom", fontsize=7.4, pad=6, color=INK)
pairs = [("before", 12, 6, "acquired\nCB", "measured\nNZ"),
         ("after", 12, 11, "acquired\nCB", "measured\nCB")]
xs = [0, 0.75, 2.0, 2.75]
vals = [12, 6, 12, 11]
cols = [ORANGE, LIGHT, ORANGE, BLUE]
ax.bar(xs, vals, 0.62, color=cols, edgecolor=INK, lw=0.5, zorder=3)
for x_, v in zip(xs, vals):
    ax.text(x_, v + 0.4, str(v), ha="center", va="bottom", fontsize=6.8, color=INK)
ax.plot([0, 0, 0.75, 0.75], [14.2, 15.0, 15.0, 14.2], color=ORANGE, lw=0.8)
ax.text(0.375, 15.3, "gap 6", ha="center", va="bottom", fontsize=6.4, color=ORANGE)
ax.plot([2.0, 2.0, 2.75, 2.75], [14.2, 15.0, 15.0, 14.2], color=BLUE, lw=0.8)
ax.text(2.375, 15.3, "gap 1", ha="center", va="bottom", fontsize=6.4, color=BLUE)
ax.set_xticks(xs)
ax.set_xticklabels(["acquired\nat CB", "measured\nat NZ", "acquired\nat CB", "measured\nat CB"],
                   fontsize=5.9, linespacing=1.3)
ax.set_ylim(0, 19)
ax.set_ylabel("Residues within 8 Å", fontsize=6.8)
ax.text(0.375, 16.6, "mismatched", ha="center", va="bottom", fontsize=6.3, color=ORANGE)
ax.text(2.375, 16.6, "matched", ha="center", va="bottom", fontsize=6.3, color=BLUE)

# ---------------------------------------------------------------- b  three populations
ax = fig.add_subplot(gs[0, 2:4]); st.panel_label(ax, "b", dx=-0.40, dy=1.09)
ax.set_title("Three populations, feature by feature", fontsize=7.4, pad=6, color=INK)
PRETTY = {"crowding_8A": "neighbours", "atoms_10A": "atom density",
          "hydrophobic_8A": "greasy", "aromatic_8A": "rings", "polar_8A": "polar",
          "anionic_8A": "acidic", "cationic_8A": "basic", "net_charge_8A": "net charge",
          "n_lys_8A": "other lysines"}
piv = MED.pivot(index="feature", columns="population", values="median")
sd = MED[MED.population == "acquired union"].set_index("feature")["sd"]
order = ["crowding_8A", "atoms_10A", "hydrophobic_8A", "aromatic_8A", "polar_8A",
         "anionic_8A", "cationic_8A", "net_charge_8A", "n_lys_8A"][::-1]
y = np.arange(len(order))
z_e = [(piv.loc[f, "engaged"] - piv.loc[f, "acquired union"]) / sd[f] for f in order]
z_n = [(piv.loc[f, "labelled, not engaged"] - piv.loc[f, "acquired union"]) / sd[f] for f in order]
ax.axvline(0, color=ORANGE, lw=1.1, zorder=1)
for yi, a_, b_ in zip(y, z_e, z_n):
    ax.plot([a_, b_], [yi, yi], color=LIGHT, lw=1.5, zorder=2)
ax.plot(z_e, y, "o", color=BLUE, ms=4.0, ls="none", zorder=4)
ax.plot(z_n, y, "o", color=GREY, ms=4.0, ls="none", zorder=4)
ax.set_yticks(y); ax.set_yticklabels([PRETTY[f] for f in order], fontsize=6.2)
lo = min(min(z_e), min(z_n)); hi = max(max(z_e), max(z_n))
pad = 0.10 * max(hi - lo, 0.2)
ax.set_xlim(lo - pad, hi + pad)
ax.set_ylim(-0.7, 9.9)
ax.set_xlabel("Difference from acquired lysines (SD)", fontsize=6.8)
span = hi - lo
ax.plot([lo], [9.3], "o", color=BLUE, ms=4.0)
ax.text(lo + span * 0.035, 9.3, "engaged", color=BLUE, fontsize=6.2, ha="left", va="center")
ax.plot([lo + span * 0.42], [9.3], "o", color=GREY, ms=4.0)
ax.text(lo + span * 0.455, 9.3, "labelled, not engaged", color=GREY, fontsize=6.2,
        ha="left", va="center")

HEAD = {}
fd.panel("fig8", "b", [dict(feature=_f, plain_name=PRETTY[_f],
                            median_acquired=float(piv.loc[_f, "acquired union"]),
                            median_engaged=float(piv.loc[_f, "engaged"]),
                            median_labelled_not_engaged=float(piv.loc[_f,
                                                        "labelled, not engaged"]),
                            sd_acquired=float(sd[_f]),
                            z_engaged=float((piv.loc[_f, "engaged"]
                                             - piv.loc[_f, "acquired union"]) / sd[_f]),
                            z_not_engaged=float((piv.loc[_f, "labelled, not engaged"]
                                                 - piv.loc[_f, "acquired union"]) / sd[_f]))
                       for _f in order[::-1]],
         note="medians of the nine environment features, three populations, and the engaged "
              "and unengaged offsets in SD of the acquired distribution")
HEAD["b_n_features"] = len(order)
HEAD["b_n_identical_medians"] = int(sum(
    1 for _f in order
    if piv.loc[_f, "acquired union"] == piv.loc[_f, "engaged"] == piv.loc[_f,
                                                             "labelled, not engaged"]))
for _f in order:
    HEAD["b_z_engaged_%s" % _f] = round(float((piv.loc[_f, "engaged"]
                                        - piv.loc[_f, "acquired union"]) / sd[_f]), 2)
HEAD["b_largest_abs_z"] = round(max(abs(_z) for _z in z_e), 2)

# ---------------------------------------------------------------- c  overlap, no threshold
ax = fig.add_subplot(gs[0, 4:6]); st.panel_label(ax, "c", dx=-0.34, dy=1.09)
ax.set_title("Distance from a typical engaged lysine", fontsize=7.4, pad=6, color=INK)
CV = pd.read_csv(os.path.join(LR, "step28d_curve.csv"))
ax.plot(CV.radius, 100 * CV.engaged, color=BLUE, lw=1.5)
ax.plot(CV.radius, 100 * CV.not_engaged, color=GREY, lw=1.3)
ax.plot(CV.radius, 100 * CV.union, color=ORANGE, lw=1.5)
ax.set_xlabel("Distance from the engaged centre\n(engaged SD)", fontsize=6.8)
ax.set_ylabel("Cumulative share of each population (%)", fontsize=6.8)
ax.set_xlim(0, 6.2); ax.set_ylim(0, 104)
ax.text(0.30, 96, "engaged", color=BLUE, fontsize=6.2, ha="left", va="center")
ax.text(0.30, 88, "labelled, not engaged", color=GREY, fontsize=6.2, ha="left", va="center")
ax.text(0.30, 80, "acquired", color=ORANGE, fontsize=6.2, ha="left", va="center")
ax.text(6.05, 8, "P = %.3f\n(0.5 = interchangeable)" % P_CLOSER, ha="right", va="bottom",
        fontsize=6.2, color=INK, linespacing=1.35)

fd.panel("fig8", "c", CV, note="cumulative share of each population within a given "
                                "distance of the engaged centre, in engaged SD")
HEAD["c_p_closer"] = P_CLOSER
HEAD["c_n_union_positions"] = int(armA.n_union)
HEAD["c_n_engaged"] = int(armA.n_eng)
HEAD["c_n_labelled_not_engaged"] = int(armA.n_not)
HEAD["c_separation_armA"] = round(float(armA.d_union_eng), 3)
HEAD["c_control_armA"] = round(float(armA.d_control), 3)
HEAD["c_overlap_union_pct"] = round(100 * float(armA.overlap_union), 0)
HEAD["c_overlap_not_engaged_pct"] = round(100 * float(armA.overlap_not), 0)

# ---------------------------------------------------------------- d  crystal validation
ax = fig.add_subplot(gs[1, 0:3]); st.panel_label(ax, "d", dx=-0.085, dy=1.10)
ax.set_title("Crystal structures: six acquired lysines against their own noise floor",
             fontsize=7.4, pad=6, color=INK)
SITES = [("PTPN11 76", "SHP2\nE76K"), ("NRAS 61", "NRAS\nQ61K"), ("AKT1 17", "AKT1\nE17K"),
         ("PIK3CA 542", "PI3Kα\nE542K"), ("PIK3CA 545", "PI3Kα\nE545K"),
         ("TP53 280", "p53\nR280K")]
rng = np.random.default_rng(3)
floor = {}
for site, _ in SITES:
    floor[site] = XT[(XT.site == site) & (XT.kind == "WT vs WT")].rmsd.median()
ax.axhline(1.0, color=BLUE, lw=1.2, zorder=2)
for i, (site, lab) in enumerate(SITES):
    f = floor[site]
    for kind, col, off in (("other sub vs WT", GREY, -0.19), ("K vs WT", ORANGE, 0.19)):
        v = XT[(XT.site == site) & (XT.kind == kind)].rmsd.values / f
        if not len(v):
            continue
        ax.plot(i + off + rng.uniform(-0.075, 0.075, len(v)), v, "o", color=col, ms=3.2,
                ls="none", zorder=4, markeredgewidth=0)
        ax.plot([i + off - 0.13, i + off + 0.13], [np.median(v)] * 2, color=INK, lw=1.5, zorder=5)
        ha_, dx_ = ("left", 0.20) if kind == "K vs WT" else ("right", -0.20)
        ax.text(i + off + dx_, np.median(v), "%.1f" % np.median(v), ha=ha_,
                va="center", fontsize=6.1, color=col)
ax.set_xticks(range(len(SITES)))
ax.set_xticklabels([s[1] for s in SITES], fontsize=6.2, linespacing=1.3)
for i, (site, _) in enumerate(SITES):
    nfl = int(((XT.site == site) & (XT.kind == "WT vs WT")).sum())
    ax.text(i, 5.85, "floor n=%d" % nfl, ha="center", va="center", fontsize=5.5, color=BLUE)
ax.set_ylim(0, 6.95); ax.set_xlim(-0.60, 5.60)
ax.set_ylabel("Local shift, relative to that protein's\nwild-type-vs-wild-type floor",
              fontsize=6.7)
ax.text(-0.55, 6.62, "acquired lysine", color=ORANGE, fontsize=6.2, ha="left", va="center")
ax.text(0.95, 6.62, "other substitutions at the same position", color=GREY,
        fontsize=6.2, ha="left", va="center")
ax.text(-0.55, 1.16, "wild-type floor", color=BLUE, fontsize=6.1, ha="left", va="bottom")

_drows = []
for _site, _lab in SITES:
    _f = floor[_site]
    _k = XT[(XT.site == _site) & (XT.kind == "K vs WT")].rmsd.values / _f
    _o = XT[(XT.site == _site) & (XT.kind == "other sub vs WT")].rmsd.values / _f
    _nfl = int(((XT.site == _site) & (XT.kind == "WT vs WT")).sum())
    _drows.append(dict(site=_site, label=_lab.replace(chr(10), " "),
                       wt_vs_wt_floor_angstrom=float(_f), n_floor_pairs=_nfl,
                       median_K_vs_WT_relative=float(np.median(_k)) if len(_k) else None,
                       n_K_pairs=len(_k),
                       median_other_sub_relative=float(np.median(_o)) if len(_o) else None,
                       n_other_sub_pairs=len(_o)))
    _key = _lab.replace(chr(10), "_")
    if len(_k):
        HEAD["d_%s_relative" % _key] = round(float(np.median(_k)), 2)
    if len(_o):
        HEAD["d_%s_other_sub_relative" % _key] = round(float(np.median(_o)), 2)
fd.panel("fig8", "d", _drows,
         note="local Ca displacement in a co-local shell, normalised to that protein's own "
              "wild-type-vs-wild-type floor")
HEAD["d_n_sites"] = len(SITES)

# ---------------------------------------------------------------- e  the drugged site
# Layout note: the assembly occupies a narrow horizontal band (v roughly -3 to +2) with the
# ligand's far end reaching v -7 and +5 at high u. Everything at v > 3 or v < -4 for u < 2 is
# empty, so the labels live there and nothing is placed over the structure.
ax = fig.add_subplot(gs[1, 3:6]); st.panel_label(ax, "e", dx=-0.22, dy=1.10)
ax.set_title("AKT1 E17K drugged: PDB 8UVY, measured here", fontsize=7.4, pad=6, color=INK)
AT = pd.read_csv(os.path.join(LR, "step32_akt1_atoms.csv"))
BD = pd.read_csv(os.path.join(LR, "step32_akt1_bonds.csv"))
u, v = AT.u.values, AT.v.values
for _, b in BD.iterrows():
    i, j = int(b.i), int(b.j)
    col = ORANGE if "ligand" in (AT.part[i], AT.part[j]) else INK
    ax.plot([u[i], u[j]], [v[i], v[j]], "-", color=col, lw=1.1, zorder=3,
            solid_capstyle="round")
zi = AT.index[AT.part == "zinc"][0]
nz = AT.index[(AT.part == "lys17") & (AT.name == "NZ")][0]
for _, r in AT[AT.part == "znligand"].iterrows():
    ax.plot([u[zi], r.u], [v[zi], r.v], ls=(0, (2, 2)), color=BLUE, lw=1.0, zorder=4)
lys = AT.part == "lys17"
ax.plot(u[lys], v[lys], "o", color=INK, ms=2.4, zorder=5)
lg = AT.part == "ligand"
ax.plot(u[lg], v[lg], "o", color=ORANGE, ms=2.4, zorder=5)
for _, r in AT[AT.part == "znligand"].iterrows():
    if r["name"].startswith("CYS"):
        ax.plot([r.u], [r.v], "o", color=BLUE, ms=4.4, zorder=6)
ax.plot([u[zi]], [v[zi]], "o", color=BLUE, ms=7.0, zorder=6)
ax.plot([u[nz]], [v[nz]], "o", color=INK, ms=4.6, zorder=7)
ax.text(u[zi] + 0.35, v[zi] + 1.0, "Zn", ha="left", va="bottom", fontsize=6.5, color=BLUE)

# Placed in the empty band BELOW the lysine chain. Two earlier positions (top-left, then a
# narrower block in the same place) both ran into the Cys310 and Zn atoms at u ~ -9.
ax.text(-15.6, -2.9,
        "Zn contacts\n"
        "Nζ 2.13 Å\n"
        "Cys296 2.33\n"
        "Cys310 2.22\n"
        "ligand 1.84",
        ha="left", va="top", fontsize=5.9, color=BLUE, linespacing=1.55)
ax.text(-13.2, 2.9, "Lys17", ha="left", va="center", fontsize=6.6, color=INK)
ax.text(-0.6, 6.1, "compound 3", ha="left", va="center", fontsize=6.4, color=ORANGE)
c25 = AT.index[(AT.part == "ligand") & (AT.name == "C25")]
if len(c25):
    k = c25[0]
    ax.annotate("covalent bond to Lys17 Nζ, 1.32 Å",
                xy=((u[nz] + u[k]) / 2, (v[nz] + v[k]) / 2), xytext=(3.2, -9.4),
                fontsize=6.4, color=ORANGE, ha="left", va="center",
                arrowprops=dict(arrowstyle="-", color=ORANGE, lw=0.7, shrinkB=3))
ax.text(0.5, -0.03, "the same site returns 0 identifications in every plex searched here",
        transform=ax.transAxes, ha="center", va="top", fontsize=6.4, color=INK)
ax.set_aspect("equal"); ax.axis("off")
ax.set_xlim(u.min() - 1.6, u.max() + 1.6); ax.set_ylim(v.min() - 4.6, v.max() + 2.6)

# NO FIGURE-LEVEL TITLE. 2026-09-05 (M.A.): the title belongs in the caption, not in the
# artwork. An in-figure title is pixels rather than text — it cannot be searched, indexed or
# copy-edited, it does not resize with the caption, and it duplicates the caption's first line
# on the page. The panel titles below stay; those are labels, not claims. The gridspec top
# margin was raised from 0.915 to 0.95 to reclaim the band the title occupied.

fd.panel("fig8", "e", AT, note="PDB 8UVY local shell as plotted: atom identities and "
                                "2-D projection coordinates")
HEAD["e_n_atoms_8uvy_shell"] = int(len(AT))
HEAD["e_akt1_e17k_spectra"] = int(AK[AK.change == "E17K"].detected.sum())
if AKT_OK:
    HEAD["e_akt1_e17k_peptides_searched"] = len(AKT_PEPS)
    HEAD["e_akt1_e17k_variant_peptides"] = AKT_VAR
    HEAD["e_akt1_e17k_downstream_peptides"] = AKT_DOWN
    HEAD["e_akt1_e17k_plexes_searchable"] = len(AKT_PLEXES)
    HEAD["e_akt1_e17k_peptide_sequences"] = AKT_PEPS
fd.head("fig8", HEAD)
st.save(fig, "fig8")
import shutil
for e in ("pdf", "png"):
    shutil.copy("fig8.%s" % e, os.path.join(ROOT, "results", "fig8.%s" % e))
print("fig8: 5 panels | P(closer) %.3f | sites in d: %d | 8UVY atoms %d"
      % (P_CLOSER, len(SITES), len(AT)))
print("  e AKT1 E17K: %d peptides (%d variant, %d downstream) in %d plex(es) %s -> %d spectra"
      % (len(AKT_PEPS), AKT_VAR, AKT_DOWN, len(AKT_PLEXES), AKT_PLEXES,
         int(AK[AK.change == "E17K"].detected.sum())))
print("     sequences: %s" % ", ".join(AKT_PEPS))
print("  b feature medians (acquired / engaged / labelled-not-engaged), z vs acquired SD:")
for _f in order[::-1]:
    _a = piv.loc[_f, "acquired union"]; _e = piv.loc[_f, "engaged"]
    _n = piv.loc[_f, "labelled, not engaged"]
    print("     %-16s %8.3g  %8.3g  %8.3g   | z_eng %+.3f  z_noteng %+.3f"
          % (PRETTY[_f], _a, _e, _n, (_e - _a) / sd[_f], (_n - _a) / sd[_f]))
print("  d crystallographic sites, median K-vs-WT relative to that protein's WT floor:")
for _site, _lab in SITES:
    _f = floor[_site]
    _k = XT[(XT.site == _site) & (XT.kind == "K vs WT")].rmsd.values / _f
    _o = XT[(XT.site == _site) & (XT.kind == "other sub vs WT")].rmsd.values / _f
    _nfl = len(XT[(XT.site == _site) & (XT.kind == "WT vs WT")])
    print("     %-14s K vs WT %s (n %d) | other subs %s (n %d) | floor pairs n = %d"
          % (_site.replace(chr(10), " "),
             ("%.2f" % float(np.median(_k))) if len(_k) else "-", len(_k),
             ("%.2f" % float(np.median(_o))) if len(_o) else "-", len(_o), _nfl))
