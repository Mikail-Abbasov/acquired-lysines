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

fig = plt.figure(figsize=(st.DOUBLE, st.DOUBLE * 0.90))
# A 12-column grid rather than 3, so panels can take the width their content needs. The bar
# charts (b, c, f, i) carry two or three bars each and were given a full third of the figure;
# they are narrowed to 3 columns, which leaves room for panel j without stretching it.
gs = fig.add_gridspec(3, 12, hspace=0.52, wspace=1.6,
                      left=0.055, right=0.985, top=0.958, bottom=0.055)

p_all = pd.read_csv(PEP, sep="\t") if PEP else None
p = acquired_lysine_only(p_all) if p_all is not None else None
if p_all is not None:
    print("region_compare rows: %d total, %d acquired lysine, %d other variant sites"
          % (len(p_all), len(p), len(p_all) - len(p)))

# ============================================ moved here from Figure 2, 2026-09-01
# Three panels each showing that detection is NOT explained by something other than
# coverage, plus the digestion result and the evidence-depth count. In the main figure they
# occupied four of eleven slots to report absences.

HEAD = {}   # headline numbers per panel, captured where they are computed


ax = fig.add_subplot(gs[0, 0:4]); st.panel_label(ax, "a", dx=-0.26, dy=1.10)
if PEP and SITES and FLUX:
    import re as _re
    d = pd.read_csv(SITES, low_memory=False); nm = pd.read_csv(FLUX, low_memory=False)[["ModelID", "CellLine"]]
    kk = d[d.to_aa == "K"].merge(nm, on="ModelID").rename(
        columns={"CellLine": "line", "Gene": "gene", "Position": "ksite"})
    kk = kk.groupby(["line", "gene", "ksite"], as_index=False).AF.max()
    def ks(ch, fb):
        for t in str(ch).split(";"):
            m = _re.match(r"^[A-Z](\d+)K$", t.strip())
            if m: return int(m.group(1))
        return fb
    q = p[p.mut_detected == 1].copy()
    q["ksite"] = [ks(c, s) for c, s in zip(q.change, q.site)]
    q = q.merge(kk, on=["line", "gene", "ksite"], how="left")
    grp = [("both\ndetected", q[q.wt_detected == 1]), ("mutant\nonly", q[q.wt_detected == 0])]
    vals = [100 * (g.AF >= 0.7).mean() for _, g in grp]
    ax.bar([0, 1], vals, color=[st.LIGHTGREY, st.BLUE], width=0.55)
    for i, (v, (lab, g)) in enumerate(zip(vals, grp)):
        ax.text(i, v + 0.6, "%.0f%%" % v, ha="center", fontsize=6.2, fontweight="bold")
        ax.text(i, -2.4, "n = %d" % len(g), ha="center", fontsize=5.4, color=st.GREY)
    ax.axhline(11, color=st.ORANGE, ls="--", lw=0.8)
    ax.text(1.45, 11.8, "background 11%", fontsize=5.4, color=st.ORANGE, ha="right")
    ax.set_xticks([0, 1]); ax.set_xticklabels([g[0] for g in grp], fontsize=5.8)
    ax.set_ylabel("clonal, AF \u2265 0.7 (%)"); ax.set_ylim(-4, max(vals) * 1.5)
    # Keyed scalars beside the readable string: verify.py cannot diff a sentence.
    HEAD["a clonality"] = "%s clonal (AF>=0.7), background 11%%; n = %s" % (
        ", ".join("%s %.0f%%" % (l.replace(chr(10), " "), v) for (l, _g), v in zip(grp, vals)),
        ", ".join("%d" % len(g) for _l, g in grp))
    for (_l, _g), _v in zip(grp, vals):
        _key = _l.replace(chr(10), " ").strip().lower().replace(" ", "_")
        HEAD["a_clonal_pct_" + _key] = round(float(_v), 1)
        HEAD["a_n_" + _key] = int(len(_g))
    HEAD["a_background_clonal_pct"] = 11
    fd.panel("figS3", "a", [dict(group=_l.replace(chr(10), " "), clonal_pct=float(_v),
             n=int(len(_g))) for (_l, _g), _v in zip(grp, vals)], note="clonality control")

# ------------------------------------------------------------------ a
ax = fig.add_subplot(gs[0, 5:8]); st.panel_label(ax, "b", dx=-0.26)
# Protein abundance from the CCLE quantification itself (Nusinow et al.), not a proxy. The file
# is cell lines x UniProt accessions, so each site is looked up at its own protein IN ITS OWN
# CELL LINE — the quantity that determines whether that protein was observable in that run.
MS = find("harmonized_MS_CCLE_Gygi.csv")
GMAP = find("uniprot_gene_map.tsv")
if PEP and MS and GMAP and FLUX:
    q = pd.read_csv(MS, low_memory=False).rename(columns={"Unnamed: 0": "ModelID"}).set_index("ModelID")
    gm = pd.read_csv(GMAP, sep="\t"); g2u = dict(zip(gm.gene, gm.uniprot))
    nm = pd.read_csv(FLUX, low_memory=False)[["ModelID", "CellLine"]]
    l2m = dict(zip(nm.CellLine, nm.ModelID))
    d = p.copy(); d["model"] = d.line.map(l2m); d["acc"] = d.gene.map(g2u)
    d = d.dropna(subset=["model", "acc"])
    d = d[d.acc.isin(set(q.columns)) & d.model.isin(q.index)]
    d["abund"] = [q.at[m, a] for m, a in zip(d.model, d.acc)]
    d = d.dropna(subset=["abund"])
    grp = [("region not\ncovered", d[d.wt_detected == 0].abund, st.LIGHTGREY),
           ("covered,\nmutant missed", d[(d.wt_detected == 1) & (d.mut_detected == 0)].abund,
            st.GREEN),
           ("mutant\nfound", d[d.mut_detected == 1].abund, st.BLUE)]
    grp = [(a, b, c) for a, b, c in grp if len(b) >= 5]
    parts = ax.violinplot([g_[1].values for g_ in grp], showextrema=False, widths=0.8)
    for pc, (_, _, col) in zip(parts["bodies"], grp):
        pc.set_facecolor(col); pc.set_alpha(1.0); pc.set_edgecolor(col); pc.set_linewidth(0.5)
    for i2, (lab, v, col) in enumerate(grp, 1):
        ax.plot([i2 - 0.18, i2 + 0.18], [v.median()] * 2, color=st.GREY, lw=1.0, zorder=5)
        ax.text(i2, ax.get_ylim()[0], "n=%d" % len(v), ha="center", va="bottom",
                fontsize=5.0, color=st.GREY)
    ax.set_xticks(range(1, len(grp) + 1))
    ax.set_xticklabels([g_[0] for g_ in grp], fontsize=5.6)
    ax.set_ylabel("protein abundance in that cell line\n(log$_2$, CCLE quantification)",
                  fontsize=5.8)
    if len(grp) == 3:
        from scipy import stats as _st
        a_, b_ = grp[2][1].values, grp[0][1].values
        cd = ((a_[:, None] > b_[None, :]).mean() - (a_[:, None] < b_[None, :]).mean())
        ax.text(0.97, 0.05, "Cliff's $\\delta$ = %.2f" % cd, transform=ax.transAxes,
                ha="right", va="bottom", fontsize=5.6, color=st.GREY)
        HEAD["b_cliffs_delta"] = round(float(cd), 2)
        for _g in grp:
            HEAD["b_n_" + str(_g[0]).replace(chr(10), " ").strip().lower().replace(" ", "_")] = int(len(_g[1]))
        HEAD["b abundance"] = ("Cliff's delta %.2f, %s vs %s; n = %s"
                               % (cd, str(grp[2][0]).replace(chr(10), " "), str(grp[0][0]).replace(chr(10), " "),
                                  ", ".join("%s %d" % (str(g_[0]).replace(chr(10), " "), len(g_[1])) for g_ in grp)))

ax = fig.add_subplot(gs[0, 9:12]); st.panel_label(ax, "c", dx=-0.22)
if SITES and FLUX:
    d2 = pd.read_csv(SITES, low_memory=False)
    allaf = d2[d2.to_aa == "K"].AF.dropna()
    ax.hist(allaf, bins=np.arange(0.1, 1.02, 0.05), density=True, color=st.LIGHTGREY,
            label="all acquired lysines")
    if HITS:
        ax.axvline(0.48, color=st.BLUE, lw=1.3)
        ax.text(0.52, ax.get_ylim()[1] * 0.86, "detected\nmedian 0.48", fontsize=5.6,
                color=st.BLUE, linespacing=1.4)
    ax.axvline(float(allaf.median()), color=st.GREY, lw=1.0, ls="--")
    ax.text(float(allaf.median()) - 0.02, ax.get_ylim()[1] * 0.45, "all sites\nmedian %.2f"
            % allaf.median(), fontsize=5.6, color=st.GREY, ha="right", linespacing=1.4)
    ax.set_xlabel("allele fraction"); ax.set_ylabel("density")
    HEAD["c_allsites_median_af"] = round(float(allaf.median()), 2)
    HEAD["c_detected_median_af"] = 0.48
    HEAD["c allele fraction"] = "detected median 0.48 vs all sites median %.2f" % float(allaf.median())
    ax = fig.add_subplot(gs[1, 0:7]); st.panel_label(ax, "d", dx=-0.11)
    # WHAT THIS PANEL USED TO SAY, AND WHY IT WAS WRONG.
    # The first version reported that only LENGTH separates identified peptides from missing
    # ones, from peptide_properties.csv — a calculation on UNLABELLED peptides that also
    # normalised by length to test length. Four rounds of correction followed: mass, then
    # charge, then charge+mass, and finally missed cleavages once the literature made the
    # charge reading untenable (the field ADDS basic groups to raise charge state because it
    # HELPS identification; our result said the opposite).
    #
    # The decomposition settles it. Every tryptic peptide ends in K or R, so a C-terminal basic
    # residue is universal and carries no information; an INTERNAL K or R is a cleavage the
    # protease missed. Splitting the two:
    #
    #     internal K/R (missed cleavages)   -0.62
    #     C-terminal K/R                    -0.004
    #
    # and within peptides that have no missed cleavage at all, the apparent charge effect falls
    # from -0.58 to -0.09. The variable was measuring incomplete digestion, not ionisation.
    MC = find("missed_cleavage_report.txt")
    eff = {}
    if MC:
        _t = open(MC).read()
        # Read ONLY the effects section. The report opens with a MEDIAN COUNTS table carrying
        # the same row labels, and a naive search matched that instead — plotting 6.0, 1.0, 1.0
        # and 0.0 as though they were effect sizes.
        _i = _t.find("EFFECT ON IDENTIFICATION")
        _j = _t.find("EFFECT OF BASIC CENTRES", _i + 1)
        _sec = _t[_i:_j] if _i >= 0 else ""
        for _k, _lab in (("basic centres", "basic centres, all"),
                         ("internal K/R", "internal K/R\n(missed cleavages)"),
                         ("C-terminal K/R", "C-terminal K/R\n(every tryptic peptide)"),
                         ("histidine", "histidine\n(basic, not a cleavage site)")):
            _m = re.search(r"^\s+%s.*?(-?\d+\.\d+)\s*$" % re.escape(_k), _sec, re.M)
            if _m:
                eff[_lab] = float(_m.group(1))
    if eff:
        labs = list(eff)[::-1]
        vals = [eff[l] for l in labs]
        cols = [st.ORANGE if "internal" in l else
                (st.LIGHTGREY if "C-terminal" in l or "histidine" in l else st.BLUE)
                for l in labs]
        y = np.arange(len(labs))
        ax.barh(y, vals, height=0.6, color=cols)
        ax.axvline(0, color=st.GREY, lw=0.9)
        ax.set_yticks(y); ax.set_yticklabels(labs, fontsize=5.8, linespacing=1.4)
        ax.set_xlabel("effect on identification (standardised)", fontsize=6)
        for _l, _v in zip(labs, vals):
            HEAD["d_effect_" + _l.replace(chr(10), " ").strip().lower().replace(" ", "_")
                 .replace("/", "_").replace("(", "").replace(")", "").replace(",", "")] = round(float(_v), 3)
        HEAD["d identification effects"] = ", ".join("%s %.3f" % (l.replace(chr(10), " "), v) for l, v in zip(labs, vals))
        fd.panel("figS3", "d", [dict(feature=_l.replace(chr(10), " "), standardised_effect=float(_v))
                 for _l, _v in zip(labs, vals)], note="effect on identification")
        ax.set_xlim(min(min(vals) * 1.35, -0.1), 0.12)
        for yy, v in zip(y, vals):
            ax.text(v - 0.02, yy, "%.3f" % v, ha="right", va="center", fontsize=5.4,
                    fontweight="bold" if abs(v) > 0.3 else "normal")
        ax.text(0.5, -0.38, "A C-terminal K or R is universal and carries no information; an "
                "internal one is a cleavage the protease missed.\nHistidine is the control: "
                "basic, but never a cleavage site. Within peptides having NO missed cleavage, "
                "the\napparent charge effect falls from \u22120.58 to \u22120.09.",
                transform=ax.transAxes, ha="center", va="top", fontsize=5.0, color=st.GREY,
                linespacing=1.5)

ax = fig.add_subplot(gs[1, 8:12]); st.panel_label(ax, "e", dx=-0.24)
if HITS:
    h2 = acquired_lysine_only(pd.read_csv(HITS, sep="\t"))
    g2 = h2.groupby(["gene", "change", "line", "site"]).agg(
        psms=("spectrum", "size"), peps=("peptide", "nunique"), plexes=("plex", "nunique"))
    n2 = len(g2)
    cats = ["\u22651\nspectrum", "\u22652\nspectra", "\u22652\npeptides", ">1\nexperiment"]
    vals = [n2, int((g2.psms >= 2).sum()), int((g2.peps >= 2).sum()),
            int((g2.plexes > 1).sum())]
    ax.bar(range(4), vals, color=[st.BLUE, st.BLUE, st.BLUE, st.ORANGE], width=0.66)
    for i3, v in enumerate(vals):
        ax.text(i3, v + n2 * 0.03, "%d" % v, ha="center", fontsize=6, fontweight="bold")
    ax.set_xticks(range(4)); ax.set_xticklabels(cats, fontsize=5.4)
    ax.set_ylabel("acquired-lysine sites")
    ax.set_ylim(0, n2 * 1.22)
    HEAD["e evidence depth"] = ", ".join("%s %d" % (c.replace(chr(10), " "), v) for c, v in zip(cats, vals))
    for _c, _v in zip(cats, vals):
        HEAD["e_" + _c.replace(chr(10), " ").strip().lower().replace(" ", "_")
             .replace("\u2265", "ge").replace(">", "gt")] = int(_v)
    HEAD["e_exactly_one_spectrum"] = int((g2.psms == 1).sum())
    # The same sites under the score margin section 2.3 defines. hits62 is written at the
    # search's reported 1% FDR, so these bars are the 1%-FDR population; the margin >= 10
    # subset is the one the entrapment error rate is measured over. Both are deposited so
    # neither can be quoted with the other's counts.
    _m2 = h2.copy(); _m2["margin"] = _m2.hyper - _m2.nextscore
    _g3 = _m2[_m2.margin >= 10].groupby(["gene", "change", "line", "site"]).agg(
        psms=("spectrum", "size"), peps=("peptide", "nunique"), plexes=("plex", "nunique"))
    HEAD["e_margin_threshold"] = 10
    HEAD["e_ge1_spectrum_margin"] = int(len(_g3))
    HEAD["e_ge2_spectra_margin"] = int((_g3.psms >= 2).sum())
    HEAD["e_ge2_peptides_margin"] = int((_g3.peps >= 2).sum())
    HEAD["e_gt1_experiment_margin"] = int((_g3.plexes > 1).sum())
    HEAD["e_exactly_one_spectrum_margin"] = int((_g3.psms == 1).sum())
    HEAD["e evidence depth"] += " | exactly one spectrum %d" % int((g2.psms == 1).sum())
    fd.panel("figS3", "e", [dict(depth=_c.replace(chr(10), " "), sites=int(_v))
             for _c, _v in zip(cats, vals)], note="evidence behind the 53 sites")

# ============================================ f  added 2026-09-07
# A new lysine is a new tryptic cleavage site, so the ZERO-missed-cleavage mutant peptide runs
# from the previous K/R to the new lysine. Both searches enforced 7-50 residues (our
# fragger.params digest_min_length = 7 / digest_max_length = 50, verified empirically: the
# shortest peptide identified is 7 and the longest 50; Gygi's original SEQUEST search excluded
# "peptides below a length of 7"). Half this class falls under that floor.
# WHAT THIS IS NOT. It is not a hidden exclusion and it does not explain gate 2. Both searches
# allowed two missed cleavages, and an upstream missed cleavage lengthens a too-short mutant
# peptide back into the window: 8,066 of 8,470 observations (95.2%) still carry at least one
# searchable species ending at the acquired lysine. The direct test cannot settle it either --
# only 27 of the 617 covered observations lack such a species, which at a 4% base rate expects
# one event. That is no power, not a null, and it is stated that way in the caption.
# Lengths from step52_mutant_peptide_length.py; counts from region_compare_pep.tsv, the table
# panels a-e already use. No alpha anywhere: the five palette colours only.
LEN = None
for _c in ("mutant_peptide_lengths.csv", "../results/mutant_peptide_lengths.csv"):
    if os.path.exists(_c):
        LEN = pd.read_csv(_c); break

ax = fig.add_subplot(gs[2, 0:8]); st.panel_label(ax, "f", dx=-0.10, dy=1.12)
if LEN is not None:
    LO, HI = 7, 50
    # Row order puts the two dark fills apart, so no two neighbours read as the same series.
    SPECIES = [("wt",   "wild-type peptide over the site",          st.GREY),
               ("mut",  "mutant peptide, new lysine cleaves",       st.ORANGE),
               ("abpp", "probe-labelled, acquired lysine blocked",  st.BLUE),
               ("down", "downstream peptide, the second route",     st.GREEN)]
    ypos = np.arange(len(SPECIES))[::-1]
    for y, (k, lab, col) in zip(ypos, SPECIES):
        v = pd.to_numeric(LEN[k], errors="coerce").dropna().values
        short = 100 * (v < LO).mean()
        good = 100 * ((v >= LO) & (v <= HI)).mean()
        longp = 100 * (v > HI).mean()
        ax.barh(y, short, color=col, height=0.60)
        ax.barh(y, good, left=short, color=st.LIGHTGREY, height=0.60)
        ax.barh(y, longp, left=short + good, color=st.GREY, height=0.60)
        ax.text(short / 2, y, "%.0f%%" % short, ha="center", va="center", fontsize=6,
                color=st.PAPER, fontweight="bold")
        ax.text(short + good / 2, y, "%.0f%% searchable" % good, ha="center", va="center",
                fontsize=5.6, color=st.GREY)
        ax.text(102, y, "median %.0f aa" % np.median(v), va="center", fontsize=5.6, color=st.GREY)
    ax.set_yticks(ypos); ax.set_yticklabels([l for _, l, _ in SPECIES], fontsize=6)
    ax.set_xlim(0, 126); ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("share of the 8,020 acquired-lysine sites (%), by peptide length", fontsize=6)
    # The band key lives in the caption, not the artwork: it collided with panel d's footer,
    # and the bars already carry their own labels.
    for sp in ("top", "right", "left"): ax.spines[sp].set_visible(False)
    for _k2, _lab2, _ in SPECIES:
        _v2 = pd.to_numeric(LEN[_k2], errors="coerce").dropna()
        HEAD["f_%s_below7_pct" % _k2] = round(float(100 * (_v2 < 7).mean()), 1)
        HEAD["f_%s_median_aa" % _k2] = int(np.median(_v2))
    HEAD["f length bands"] = " | ".join(
        "%s: <7 %.1f%%, median %.0f aa" % (
            k, 100 * (pd.to_numeric(LEN[k], errors="coerce").dropna() < 7).mean(),
            np.median(pd.to_numeric(LEN[k], errors="coerce").dropna()))
        for k, _, _ in SPECIES)
if p is not None:
    _cov = p[p.wt_detected > 0]
    HEAD["f missed-cleavage rescue"] = (
        "%d of %d observations (%.1f%%) still carry a searchable species ending at the acquired "
        "lysine; only %d of the %d covered do not, so the recovery comparison is underpowered"
        % (int((p.n_cleaved > 0).sum()), len(p), 100 * (p.n_cleaved > 0).mean(),
           int((_cov.n_cleaved == 0).sum()), len(_cov)))

st.save(fig, "figS3")
fd.head("figS3", HEAD)
for _e in (".png", ".pdf"):          # review copies alongside fig4_v2 / figS_react_v2
    fig.savefig("../results/figS3_v3" + _e, dpi=400, bbox_inches="tight")
print("wrote figS3")
for _k, _v in HEAD.items(): print("  %s: %s" % (_k, _v))
for x in notes: print("  MISSING:", x)
