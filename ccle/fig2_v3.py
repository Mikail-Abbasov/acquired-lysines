#!/usr/bin/env python3
"""
fig2.py — Figure 2: acquired lysines are largely invisible to conventional proteomics. 3x3.

  a  the four peptide species an acquired lysine creates
  b  gate 1 — regions never covered (two search architectures agree)
  c  gate 2 — covered regions where the mutant is still missed
  d  three sites with every alternative explanation excluded
  e  peptide length predicts which regions are never covered
  f  loss of the wild-type allele inverts gate 1
  g  only peptide length separates identified species from missing ones
  h  detection does not track allele fraction
  i  two independent routes to the same sites

SCOPE. Every panel counts ACQUIRED-LYSINE sites only. The variant databases contain every
mutation carried by each cell line, so a protein with several mutations contributes a row per
position; 18% of rows in the detection table are other substitutions. Those were legitimately
searched and are reported separately in the Supplementary material, where they also serve as a
control: the same pipeline, spectra and FDR treatment applied to variant peptides that do not
create a lysine.
"""
import os, re, sys, glob
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import abbasov_style as st
import figure_data as fd

# Headline numbers, keyed. The prints at the foot of this file stay for a human reading the log;
# HEAD is what verify.py diffs against the manuscript.
HEAD = {}

# --- portability -------------------------------------------------------------------------
# The data search paths below were written as ~/Desktop/DepMap/... , which only resolves on
# the machine the analysis was built on. They are now anchored to this script's own location,
# so the build works wherever the folder lives (a collaborator's checkout, a container, a
# different user account). On the original machine the two are identical.
_TOP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

st.apply()

def find(n):
    for d in (".", "..", os.path.join(_TOP, 'ccle'),
              _TOP,
              os.path.join(_TOP, 'Datasets'),
              os.path.join(_TOP, 'ligand_rebuild'),
              os.path.join(_TOP, 'Datasets/pandepmap_out')):
        p = os.path.join(d, n)
        if os.path.exists(p): return p
    return None

MARGIN = 10          # the score-margin threshold section 2.3 defines


def acquired_lysine_only(df):
    """Keep rows whose stored `site` IS the acquired lysine.

    region_compare writes one row per mutated position, and a cell line carrying several
    mutations in one protein contributes a row for each. Those extra rows are real variant
    sites that were legitimately searched — but they are not acquired lysines, and counting
    them inflates every denominator in this figure. In the 42-plex table 1,902 of 10,372 rows
    (18%) are of this kind. They are reported separately (Supplementary), not discarded.
    """
    import re as _re3
    def _isk(ch, site):
        for t in str(ch).split(";"):
            m = _re3.match(r"^[A-Z](\d+)K$", t.strip())
            if m and int(m.group(1)) == site:
                return True
        return False
    return df[[_isk(c, s) for c, s in zip(df.change, df.site)]]


PEP  = find("region_compare_pep.tsv")
CONV = find("region_compare_conv.tsv")
HITS = find("hits62.tsv")
SITES= find("sites_covalent_residues.csv")
FLUX = find("residue_flux_by_line.csv")
notes = []

fig = plt.figure(figsize=(st.DOUBLE, st.DOUBLE * 0.60))
# A 12-column grid rather than 3, so panels can take the width their content needs. The bar
# charts (b, c, f, i) carry two or three bars each and were given a full third of the figure;
# they are narrowed to 3 columns, which leaves room for panel j without stretching it.
gs = fig.add_gridspec(2, 12, hspace=0.46, wspace=1.5,
                      left=0.055, right=0.985, top=0.958, bottom=0.055)

p_all = pd.read_csv(PEP, sep="\t") if PEP else None
p = acquired_lysine_only(p_all) if p_all is not None else None
if p_all is not None:
    print("region_compare rows: %d total, %d acquired lysine, %d other variant sites"
          % (len(p_all), len(p), len(p_all) - len(p)))

# ================================================================== six panels, 2 x 3
# The four alternative-explanation panels (allele loss, abundance, allele fraction) and the
# digestion panel move to the supplement: three nulls and a methods result are not three
# and a half main-text slots. Evidence depth, which had its own panel, is now a line inside
# the recovery panel, where it qualifies the number it belongs to.


# ------------------------------------------------------------------ a
ax = fig.add_subplot(gs[0, 0:5]); st.panel_label(ax, "a", dx=-0.07, dy=1.12)
ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.set_axis_off()
# backbone
ax.plot([3.6, 9.6], [9.2, 9.2], color=st.GREY, lw=1.1)
for x, lab, col in ((4.4, "K/R", st.GREY), (6.2, "E\u2192K", st.ORANGE), (8.6, "K/R", st.GREY)):
    ax.plot([x, x], [8.95, 9.45], color=col, lw=1.5 if col != st.GREY else 1.0)
    ax.text(x, 9.72, lab, ha="center", fontsize=5.2, color=col,
            fontweight="bold" if col != st.GREY else "normal")
# species: label sits directly ON its own row, bar to the right of it
# Four species from the five colours. The two mutant-carrying species are the SAME peptide
# digested two ways, so they share the teal and are told apart by fill, not by a sixth hue:
# solid = cleaved, hatched outline = missed cleavage. Wild type drops to the background grey.
rows = [(7.6, 4.4, 6.2, st.BLUE,      "mutant, cleaved",         True,  None),
        (6.0, 4.4, 8.6, st.BLUE,      "mutant, missed cleavage", True,  "///"),
        (4.4, 6.2, 8.6, st.GREEN,     "downstream",              False, None),
        (2.8, 4.4, 8.6, st.LIGHTGREY, "wild type",               False, None)]
for y, x0, x1, col, lab, blue, hat in rows:
    ax.add_patch(Rectangle((x0, y - 0.34), x1 - x0, 0.68,
                           facecolor=st.PAPER if hat else col,
                           edgecolor=col, lw=0.9 if hat else 0.5, hatch=hat))
    # labels are right-aligned and end just before the bars begin, so the longest one
    # ("mutant, missed cleavage") cannot run into its bar
    ax.text(4.15, y, lab, ha="right", va="center", fontsize=4.8,
            color=st.GREEN if col == st.GREEN else (st.BLUE if blue else st.GREY),
            fontweight="bold" if (blue or col == st.GREEN) else "normal")
ax.text(0.6, 1.55, "the first two species carry the mutated residue", fontsize=5.0,
        color=st.BLUE)
ax.text(0.6, 0.95, "the third carries none, but its tryptic boundary", fontsize=5.0,
        color=st.GREEN)
ax.text(0.6, 0.35, "exists only in the mutant protein", fontsize=5.0, color=st.GREEN)

# ------------------------------------------------------------------ b
ax = fig.add_subplot(gs[0, 6:8]); st.panel_label(ax, "b", dx=-0.42)
if PEP and CONV:
    c = acquired_lysine_only(pd.read_csv(CONV, sep="\t"))
    vals = [100 * c.wt_detected.sum() / len(c), 100 * p.wt_detected.sum() / len(p)]
    ax.bar([0, 1], vals, color=[st.LIGHTGREY, st.BLUE], width=0.55)
    for i, v in enumerate(vals):
        ax.text(i, v + 0.2, "%.2f%%" % v, ha="center", fontsize=6.2, fontweight="bold")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["conventional", "mutation-\naware"], fontsize=5.8)
    ax.set_ylabel("regions covered (%)"); ax.set_ylim(0, max(vals) * 1.45)
    ax.text(0.5, max(vals) * 1.28, "~93% never observed", ha="center", fontsize=5.8, color=st.GREY)
    HEAD["b_conv_covered_n"] = int(c.wt_detected.sum()); HEAD["b_conv_total"] = int(len(c))
    HEAD["b_conv_covered_pct"] = round(vals[0], 2)
    HEAD["b_pep_covered_n"] = int(p.wt_detected.sum()); HEAD["b_pep_total"] = int(len(p))
    HEAD["b_pep_covered_pct"] = round(vals[1], 2)
    HEAD["b_never_covered_n"] = int(len(p) - p.wt_detected.sum())
    HEAD["b_never_covered_pct"] = round(100 * (len(p) - p.wt_detected.sum()) / len(p), 2)
    fd.panel("fig2", "b", [
        dict(search="conventional 20,431-protein MANE", covered=int(c.wt_detected.sum()),
             observations=int(len(c)), pct=vals[0]),
        dict(search="mutation-aware 10.17M peptide", covered=int(p.wt_detected.sum()),
             observations=int(len(p)), pct=vals[1])], note="gate 1, both architectures")

# ------------------------------------------------------------------ c
ax = fig.add_subplot(gs[0, 9:12]); st.panel_label(ax, "c", dx=-0.26, dy=1.10)
if PEP:
    cov = p[p.wt_detected == 1]; unc = p[p.wt_detected == 0]
    rc = 100 * cov.mut_detected.sum() / len(cov); ru = 100 * unc.mut_detected.sum() / len(unc)
    ax.bar([0, 1], [ru, rc], color=[st.LIGHTGREY, st.BLUE], width=0.55)
    for i, (v, n, kk) in enumerate([(ru, len(unc), int(unc.mut_detected.sum())),
                                    (rc, len(cov), int(cov.mut_detected.sum()))]):
        ax.text(i, v + 0.1, "%.2f%%" % v, ha="center", fontsize=6.2, fontweight="bold")
        ax.text(i, -0.5, "%d/%s" % (kk, f"{n:,}"), ha="center", fontsize=5.4, color=st.GREY)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["not\ncovered", "covered"], fontsize=5.8)
    ax.set_ylabel("mutant recovered (%)"); ax.set_ylim(-0.9, rc * 1.32)
    ax.text(0.5, rc * 1.18, "%.0f\u00d7" % (rc / max(ru, 1e-9)), ha="center",
            fontsize=6.4, color=st.BLUE, fontweight="bold")
    HEAD["c_covered_found"] = int(cov.mut_detected.sum()); HEAD["c_covered_n"] = int(len(cov))
    HEAD["c_covered_pct"] = round(rc, 2)
    HEAD["c_uncovered_found"] = int(unc.mut_detected.sum()); HEAD["c_uncovered_n"] = int(len(unc))
    HEAD["c_uncovered_pct"] = round(ru, 2)
    HEAD["c_total_found"] = int(p.mut_detected.sum())
    HEAD["c_total_pct"] = round(100 * p.mut_detected.sum() / len(p), 2)

    # THE SAME DETECTIONS UNDER THE PAPER'S OWN SCORE MARGIN.
    # region_compare.py counts a detection at the search's reported 1% FDR and applies no margin
    # filter, so these counts are NOT the margin >= 10 population that section 2.3's error
    # analysis is defined over. Both are correct and they are different denominators, exactly as
    # Figure 3b (distinct peptides) and Figure 3c (identifications) are. Quoting one with the
    # other's counts is the mistake this block exists to make impossible.
    if HITS:
        _hm = acquired_lysine_only(pd.read_csv(HITS, sep="\t")).copy()
        _hm["margin"] = _hm.hyper - _hm.nextscore
        _hm["key"] = list(zip(_hm.gene, _hm.change, _hm.line, _hm.site))
        _bm = _hm.groupby("key").margin.max()
        _keep = set(_bm[_bm >= MARGIN].index)
        _pk = list(zip(p.gene, p.change, p.line, p.site))
        _sel = pd.Series([k in _keep for k in _pk], index=p.index) & (p.mut_detected == 1)
        _cv = p.wt_detected == 1
        HEAD["c_margin_threshold"] = MARGIN
        HEAD["c_total_found_margin"] = int(_sel.sum())
        HEAD["c_total_pct_margin"] = round(100 * int(_sel.sum()) / len(p), 2)
        HEAD["c_covered_found_margin"] = int((_sel & _cv).sum())
        HEAD["c_covered_pct_margin"] = round(100 * int((_sel & _cv).sum()) / int(_cv.sum()), 2)
        HEAD["c_uncovered_found_margin"] = int((_sel & ~_cv).sum())
        HEAD["c_uncovered_pct_margin"] = round(
            100 * int((_sel & ~_cv).sum()) / int((~_cv).sum()), 2)
        HEAD["c_genomic_sites"] = len({(g, st) for g, c_, l_, st in _bm.index})
        HEAD["c_genomic_sites_margin"] = len({(g, st) for g, c_, l_, st in _keep})
        HEAD["c_n_psms"] = int(len(_hm))
        HEAD["c_n_psms_margin"] = int((_hm.margin >= MARGIN).sum())
    fd.panel("fig2", "c", [
        dict(group="region not covered", recovered=int(unc.mut_detected.sum()),
             observations=int(len(unc)), pct=ru),
        dict(group="region covered", recovered=int(cov.mut_detected.sum()),
             observations=int(len(cov)), pct=rc)], note="gate 2, mutant-only route")

# ------------------------------------------------------------------ d
ax = fig.add_subplot(gs[1, 0:4]); st.panel_label(ax, "d", dx=-0.22)
# Driver status is the Hotspot flag carried in the DepMap variant table, not a gene list chosen
# here — which removes the question of who decided what counts as a driver.
if PEP and SITES and FLUX:
    import re as _re2
    dd = pd.read_csv(SITES, low_memory=False)
    nm2 = pd.read_csv(FLUX, low_memory=False)[["ModelID", "CellLine"]]
    kk2 = (dd[dd.to_aa == "K"].merge(nm2, on="ModelID")
           .rename(columns={"CellLine": "line", "Gene": "gene", "Position": "ksite"}))
    kk2 = kk2.groupby(["line", "gene", "ksite"], as_index=False).agg(
        AF=("AF", "max"), Hotspot=("Hotspot", "max"))
    def _ks(ch, fb):
        for t in str(ch).split(";"):
            m = _re2.match(r"^[A-Z](\d+)K$", t.strip())
            if m: return int(m.group(1))
        return fb
    cov = p[p.wt_detected == 1].copy()
    cov["ksite"] = [_ks(c, s) for c, s in zip(cov.change, cov.site)]
    cov = cov.merge(kk2, on=["line", "gene", "ksite"], how="left").dropna(subset=["AF"])
    cov["Hotspot"] = cov.Hotspot.fillna(False).astype(bool)
    miss = cov[cov.mut_detected == 0]; fnd = cov[cov.mut_detected == 1]
    hm = miss[miss.Hotspot]; hf = fnd[fnd.Hotspot]
    rng = np.random.default_rng(0)
    ax.scatter(miss.AF, rng.normal(0, 0.055, len(miss)), s=2.0, color=st.LIGHTGREY,
               linewidths=0)
    ax.scatter(fnd.AF, rng.normal(1, 0.055, len(fnd)), s=4.4, color=st.BLUE, linewidths=0)
    ax.scatter(hm.AF, rng.normal(0, 0.04, len(hm)), s=5, color=st.ORANGE,
               linewidths=0, zorder=6)
    ax.scatter(hf.AF, rng.normal(1, 0.04, len(hf)), s=5, color=st.ORANGE,
               linewidths=0, zorder=6)
    labs = sorted({"%s %s" % (r.gene, str(r.change).split(";")[0]) for r in hm.itertuples()})
    half = (len(labs) + 1) // 2
    # the labels need clear air ABOVE the found row, so the y-limit is extended rather than
    # the text being squeezed between the two rows
    ax.text(0.01, 0.99, "\n".join(labs[:half]), transform=ax.transAxes, fontsize=4.6,
            color=st.ORANGE, linespacing=1.55, va="top")
    ax.text(0.36, 0.99, "\n".join(labs[half:]), transform=ax.transAxes, fontsize=4.6,
            color=st.ORANGE, linespacing=1.55, va="top")
    ax.text(0.99, 0.99, "%d of %d driver\nhotspots missed" % (len(hm), len(hm) + len(hf)),
            transform=ax.transAxes, fontsize=5.4, color=st.ORANGE, ha="right", va="top",
            linespacing=1.4, fontweight="bold")
    ax.set_yticks([0, 1]); ax.set_yticklabels(["missed\n%d" % len(miss), "found\n%d" % len(fnd)],
                                              fontsize=5.6)
    ax.set_xlabel("allele fraction"); ax.set_xlim(0.1, 1.04); ax.set_ylim(-0.45, 3.1)

# ------------------------------------------------------------------ e
ax = fig.add_subplot(gs[1, 5:8]); st.panel_label(ax, "e", dx=-0.15)
# Restored. This is a GATE 1 result: where the region was never covered, the wild-type tryptic
# peptides are longer, and the two distributions diverge specifically in the tail beyond 20
# residues rather than shifting as a whole. That points at a subset of regions requiring
# peptides too long to observe — an argument for a protease that redraws cleavage boundaries,
# not for deeper fractionation.
if PEP and "shortest_wt" in p.columns:
    a = pd.to_numeric(p.loc[p.wt_detected == 1, "shortest_wt"], errors="coerce").dropna()
    b = pd.to_numeric(p.loc[p.wt_detected == 0, "shortest_wt"], errors="coerce").dropna()
    from scipy import stats as _st2
    bins = np.arange(4, 41, 2)
    ax.hist(a, bins=bins, density=True, color=st.BLUE, label="covered")
    ax.hist(b, bins=bins, density=True, histtype="step", color=st.ORANGE, lw=1.2,
            label="not covered")
    _, pv = _st2.mannwhitneyu(a, b)
    ax.set_xlabel("shortest tryptic peptide (residues)"); ax.set_ylabel("density")
    ax.legend(frameon=False, fontsize=5.6, handlelength=1.1)
    ax.text(0.97, 0.62, "mean %.1f vs %.1f\np = %.0e" % (a.mean(), b.mean(), pv),
            transform=ax.transAxes, ha="right", fontsize=5.6, color=st.GREY, linespacing=1.4)
    HEAD["e_covered_mean"] = round(float(a.mean()), 1); HEAD["e_covered_n"] = int(len(a))
    HEAD["e_uncovered_mean"] = round(float(b.mean()), 1); HEAD["e_uncovered_n"] = int(len(b))
    HEAD["e_mannwhitney_p"] = float("%.1e" % pv)
    fd.panel("fig2", "e", pd.DataFrame(dict(
        shortest_wt_peptide=pd.concat([a, b]).values,
        region=["covered"] * len(a) + ["never covered"] * len(b))),
        note="gate 1 explanation; one row per observation")

# ------------------------------------------------------------------ f
ax = fig.add_subplot(gs[1, 9:12]); st.panel_label(ax, "f", dx=-0.26, dy=1.10)
# Two independent routes to the same class of site. They are kept separate because their
# false-discovery behaviour differs: downstream peptides returned zero entrapment hits across
# 42 experiments (FDP 0.00%) against 0.58% for variant-containing peptides in the same searches,
# and their boundary control — the same region shifted one residue — was never identified.
DOWN = find("downstream_sites.csv")
if HITS and DOWN:
    h = acquired_lysine_only(pd.read_csv(HITS, sep="\t"))
    hv = h.groupby(["gene", "change", "line", "site"]).ngroups
    dn = pd.read_csv(DOWN)
    # YBX1 E4K was excluded here on a "+1.0087 Da, wrong molecular species" reading. That was
    # wrong: 1.0087 sits 0.0053 Da from a neutron, which at 5,213 Da is 1.0 ppm, while
    # deamidation is 0.0247 Da away. It is an isotope peak, and selecting the first isotope
    # rather than the monoisotopic happens in 7.7% of variant PSMs and 7.3% of background — the
    # same rate — and is expected at 5+ on a 5.2 kDa precursor. Reinstated. See A81.
    sites_d = dn.groupby(["gene", "change", "line"]).ngroups
    known = set(zip(h.gene, h.line))
    newd = dn[[(g, l) not in known for g, l in zip(dn.gene, dn.line)]]
    n_new = newd.groupby(["gene", "change", "line"]).ngroups
    n_conf = sites_d - n_new
    ax.bar([0], [hv], color=st.BLUE, width=0.6)
    ax.bar([1], [n_conf], color=st.GREEN, width=0.6)
    ax.bar([1], [n_new], bottom=[n_conf], color=st.LIGHTGREY, width=0.6)
    ax.text(0, hv + 1.5, "%d" % hv, ha="center", fontsize=6.4, fontweight="bold")
    ax.text(1, n_conf / 2, "%d" % n_conf, ha="center", va="center", fontsize=5.8,
            color=st.PAPER, fontweight="bold")
    ax.text(1, n_conf + n_new / 2, "%d" % n_new, ha="center", va="center", fontsize=5.8,
            color=st.GREY, fontweight="bold")
    ax.text(1, n_conf + n_new + 1.5, "%d" % sites_d, ha="center", fontsize=6.4,
            fontweight="bold")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["mutant\npeptide", "downstream\npeptide"], fontsize=6)
    ax.set_ylabel("acquired-lysine sites")
    ax.set_ylim(0, max(hv, sites_d) * 1.30)
    ax.text(0.99, 0.97, "solid: also found by\nthe other route", transform=ax.transAxes,
            ha="right", va="top", fontsize=5.0, color=st.GREY, linespacing=1.4)
    ax.text(0.5, -0.50, "of the 53 mutant-peptide sites, 41 rest on a single spectrum and 12 on\n"
            "two or more; that depth is in figure S3 e, and is what Figure 3 answers",
            transform=ax.transAxes, ha="center", va="top", fontsize=4.5,
            color=st.GREY, linespacing=1.5)
    ax.text(0.5, -0.34, "%d sites total, %d by both routes" % (hv + n_new, n_conf),
            transform=ax.transAxes, ha="center", fontsize=5.4, color=st.GREY)



# ------------------------------------------------------------------ j
# Which peptide species get identified. Three species, each comparing the peptides that were
# identified against those that were not. The x-axis is the difference in medians in standard
# deviations of that property across all searched peptides — a common scale, so three properties
# measured in different units can share one axis.
#
# The label's hydrophobic contribution was FITTED against observed retention times of 107
# identified peptides rather than assumed (best fit K = -3.9, rho 0.537).
# (The peptide-property panel that consumed peptide_properties.csv moved to figS2 in v3;
#  the loader and its series list were left behind here and are removed.)
st.save(fig, "fig2")
for _e in (".png", ".pdf"):          # review copies alongside fig4_v2 / figS_react_v2
    fig.savefig("../results/fig2_v3" + _e, dpi=400, bbox_inches="tight")
print("wrote fig2.pdf and fig2.png")
if PEP and CONV:
    _c = acquired_lysine_only(pd.read_csv(CONV, sep="\t"))
    print("b: conventional %.2f%% (%d/%d) | mutation-aware %.2f%% (%d/%d) regions covered"
          % (100 * _c.wt_detected.sum() / len(_c), int(_c.wt_detected.sum()), len(_c),
             100 * p.wt_detected.sum() / len(p), int(p.wt_detected.sum()), len(p)))
    print("b: conventional mutant detections %d (passed control: no variant sequences in that database)"
          % int(_c.mut_detected.sum()))
if PEP:
    _cov = p[p.wt_detected == 1]; _unc = p[p.wt_detected == 0]
    print("c: covered %.2f%% (%d/%s) | not covered %.2f%% (%d/%s) | overall %.2f%% (%d/%s)"
          % (100 * _cov.mut_detected.sum() / len(_cov), int(_cov.mut_detected.sum()), f"{len(_cov):,}",
             100 * _unc.mut_detected.sum() / len(_unc), int(_unc.mut_detected.sum()), f"{len(_unc):,}",
             100 * p.mut_detected.sum() / len(p), int(p.mut_detected.sum()), f"{len(p):,}"))
    if "shortest_wt" in p.columns:
        from scipy import stats as _st3
        _a = pd.to_numeric(p.loc[p.wt_detected == 1, "shortest_wt"], errors="coerce").dropna()
        _b = pd.to_numeric(p.loc[p.wt_detected == 0, "shortest_wt"], errors="coerce").dropna()
        _, _pv = _st3.mannwhitneyu(_a, _b)
        print("e: shortest tryptic peptide, covered mean %.1f (median %.0f, n=%s) vs never covered "
              "mean %.1f (median %.0f, n=%s), Mann-Whitney p = %.1e"
              % (_a.mean(), _a.median(), f"{len(_a):,}", _b.mean(), _b.median(), f"{len(_b):,}", _pv))
# headline numbers per panel, printed so no legend has to hold a number the script does not report
try:
    print("d: %d of %d driver hotspots missed | missed %d, found %d covered sites | %d gene-mutation labels across %d genes"
          % (len(hm), len(hm) + len(hf), len(miss), len(fnd), len(labs), len(set(hm.gene))))
    HEAD["d_hotspots_missed"] = int(len(hm)); HEAD["d_hotspots_total"] = int(len(hm) + len(hf))
    HEAD["d_covered_missed"] = int(len(miss)); HEAD["d_covered_found"] = int(len(fnd))
    HEAD["d_gene_mutation_labels"] = int(len(labs)); HEAD["d_genes"] = int(len(set(hm.gene)))
    fd.panel("fig2", "d", hm.assign(outcome="mutant missed")[["gene", "change", "line", "site"]]
             if hasattr(hm, "assign") else [], note="driver hotspots among covered sites")
except NameError:
    print("d: not built")
try:
    print("f: mutant-peptide route %d sites | downstream route %d (%d shared, %d route-unique) | union %d"
          % (hv, sites_d, n_conf, n_new, hv + n_new))
    HEAD["f_mutant_route_sites"] = int(hv); HEAD["f_downstream_sites"] = int(sites_d)
    HEAD["f_shared"] = int(n_conf); HEAD["f_downstream_unique"] = int(n_new)
    HEAD["f_union"] = int(hv + n_new)
    fd.panel("fig2", "f", [
        dict(route="mutant peptide", sites=int(hv)),
        dict(route="downstream peptide", sites=int(sites_d)),
        dict(route="shared by both", sites=int(n_conf)),
        dict(route="downstream only", sites=int(n_new)),
        dict(route="union", sites=int(hv + n_new))], note="two detection routes")
except NameError:
    print("f: not built")

# The scope line at the head of the run is itself a published quantity (8,470 / 8,055 / 372).
if PEP:
    HEAD["scope_acquiredK_observations"] = int(len(p))
    HEAD["scope_nonlysine_observations"] = int(len(p_all) - len(p))
    HEAD["scope_rows_total"] = int(len(p_all))
    if "line" in p.columns:
        HEAD["scope_cell_lines"] = int(p.line.nunique())
    if {"gene", "site"} <= set(p.columns):
        HEAD["scope_distinct_sites"] = int(p.groupby(["gene", "site"]).ngroups)

fd.head("fig2", HEAD)
for x in notes: print("  MISSING:", x)
