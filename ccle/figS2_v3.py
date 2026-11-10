#!/usr/bin/env python3
"""
Supplementary figure for Figure 2, controls and supporting detail for the detection analysis.

a,b  THE NON-LYSINE CONTROL. 1,902 of the 10,372 region_compare rows are variant sites that do
     NOT create a lysine, they come from multi-mutation records, one row per position. They
     went through the identical pipeline: same database construction, same search, same FDR,
     same filters.

     The result is not a null. Acquired lysines are COVERED more often (7.3% vs 4.6%) and
     RECOVERED less often once covered (4.4% vs 6.9%). A non-lysine variant does not create a
     tryptic site, so its mutant peptide shares boundaries with the wild-type peptide that
     established coverage; an acquired lysine produces a different peptide. The control
     therefore demonstrates the mechanism rather than controlling it away.

c    CONVENTIONAL vs MUTATION-AWARE, in full. Figure 2b reports two coverage numbers; this
     shows the distributions behind them.

f    Score margin of every detection, with the retention threshold marked.
g    Spectrum quality: are variant identifications supported by worse spectra than ordinary
     ones? If so the detections would be suspect regardless of FDR.
h    An annotated spectrum for one named site.

     THE DOWNSTREAM BOUNDARY CONTROL IS NOT A PANEL. Across 42 experiments it gave 25
     downstream identifications, ZERO boundary controls and ZERO entrapment hits, three
     numbers, two of them zero, which is a sentence rather than a figure. It belongs in the
     text.

d,e  MISSED CLEAVAGE. The rates measured the way the field reports them, and the comparison
     between acquired-lysine peptides and background. Reported without a benchmark: published
     values for Lys-C/trypsin range from ~5% in an optimised protocol to >40% of unique peptides
     across public databases, so a single comparator would mislead.
"""
import os, sys, re
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import abbasov_style as st
import figure_data as fd

# --- portability -------------------------------------------------------------------------
# The data search paths below were written as ~/Desktop/DepMap/... , which only resolves on
# the machine the analysis was built on. They are now anchored to this script's own location,
# so the build works wherever the folder lives (a collaborator's checkout, a container, a
# different user account). On the original machine the two are identical.
_TOP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# One body type size for this figure; every other size is a ratio to it.
BODY = 7.0
st.apply(BODY)

HEAD = {}       # headline numbers per panel, captured where they are computed
MARGIN = 10     # score-margin threshold -- one definition for the axis line, the label and the count

def find(n):
    for d in (".", "..", os.path.join(_TOP, 'ccle_plex05')):
        p = os.path.join(d, n)
        if os.path.exists(p): return p
    return None

def is_kgain(ch, site):
    for t in str(ch).split(";"):
        m = re.match(r"^[A-Z](\d+)K$", t.strip())
        if m and int(m.group(1)) == site:
            return True
    return False

P = pd.read_csv(find("region_compare_pep.tsv"), sep="\t")
P["kgain"] = [is_kgain(c, s) for c, s in zip(P.change, P.site)]
K = P[P.kgain]; N = P[~P.kgain]

fig = plt.figure(figsize=(st.DOUBLE, st.SINGLE * 2.76))
gs = fig.add_gridspec(5, 12, hspace=0.46, wspace=2.3,
                      height_ratios=[1, 1, 1, 0.18, 0.85],
                      left=0.065, right=0.985, top=0.972, bottom=0.070)

# ---------------------------------------------------------------- a: the two gates, both classes
ax = fig.add_subplot(gs[0, 0:6]); st.head(fig, ax, "a")
g1 = [100 * K.wt_detected.mean(), 100 * N.wt_detected.mean()]
cov_k = K[K.wt_detected == 1]; cov_n = N[N.wt_detected == 1]
g2 = [100 * cov_k.mut_detected.mean() if len(cov_k) else np.nan,
      100 * cov_n.mut_detected.mean() if len(cov_n) else np.nan]
x = np.arange(2); w = st.grp(2)
ax.bar(x - w/2, g1, w, color=st.BLUE, label="gate 1: region covered")
ax.bar(x + w/2, g2, w, color=st.GREEN, label="gate 2: mutant recovered\nwhere covered")
for i in range(2):
    ax.text(i - w/2, g1[i] + 0.15, "%.1f%%" % g1[i], ha="center", fontsize=st.sz(5.6))
    if not np.isnan(g2[i]):
        ax.text(i + w/2, g2[i] + 0.15, "%.1f%%" % g2[i], ha="center", fontsize=st.sz(5.6))
ax.set_xticks(x)
ax.set_xticklabels(["acquired lysine\n(%s sites)" % f"{len(K):,}",
                    "other variant\n(%s sites)" % f"{len(N):,}"], fontsize=st.sz(5.6))
ax.set_ylabel("% of sites (note: the two bars\nhave different denominators)", fontsize=st.sz(5.6),
              linespacing=1.5)
# 1.35 -> 1.58: the key sits in the upper right and at the raised body the "6.9%" bar
# value reached it. The headroom is cheaper than moving the key onto the bars.
ax.set_ylim(0, max(g1 + [v for v in g2 if not np.isnan(v)]) * 1.58)
ax.legend(frameon=False, fontsize=st.sz(5.0), loc="upper right", handlelength=1.0, labelspacing=0.5)
HEAD["a_n_acquired_lysine_observations"] = int(len(K))
HEAD["a_n_other_variant_observations"] = int(len(N))
HEAD["a_gate1_covered_pct_acquired_lysine"] = round(float(g1[0]), 1)
HEAD["a_gate1_covered_pct_other_variant"] = round(float(g1[1]), 1)
HEAD["a_gate2_recovered_pct_acquired_lysine"] = round(float(g2[0]), 1)
HEAD["a_gate2_recovered_pct_other_variant"] = round(float(g2[1]), 1)
HEAD["a_n_covered_acquired_lysine"] = int(len(cov_k))
HEAD["a_n_covered_other_variant"] = int(len(cov_n))
fd.panel("figS2", "a", [
    dict(cls="acquired lysine", gate="gate 1: region covered", pct=float(g1[0]),
         n_events=int(K.wt_detected.sum()), n_denominator=int(len(K))),
    dict(cls="other variant", gate="gate 1: region covered", pct=float(g1[1]),
         n_events=int(N.wt_detected.sum()), n_denominator=int(len(N))),
    dict(cls="acquired lysine", gate="gate 2: mutant recovered where covered", pct=float(g2[0]),
         n_events=int(cov_k.mut_detected.sum()), n_denominator=int(len(cov_k))),
    dict(cls="other variant", gate="gate 2: mutant recovered where covered", pct=float(g2[1]),
         n_events=int(cov_n.mut_detected.sum()), n_denominator=int(len(cov_n))),
    ], note="one row per bar; the two gates have different denominators")
# The gates do NOT behave the same, and that is the result: acquired lysines are covered more
# often but recovered LESS often once covered. A non-lysine variant does not create a tryptic
# site, so its mutant peptide has the same boundaries as the wild-type one; an acquired lysine
# produces a different peptide. This control demonstrates the mechanism the figure is about.

# ---------------------------------------------------------------- b: peptide length, both classes
ax = fig.add_subplot(gs[0, 7:12]); st.head(fig, ax, "b")
bins = np.arange(4, 42, 2)
_brows = []
for d, lab, col in ((K, "acquired lysine", st.BLUE), (N, "other variant", st.ORANGE)):
    v = pd.to_numeric(d.shortest_wt, errors="coerce").dropna()
    ax.hist(v.clip(upper=40), bins=bins, density=True, histtype="step", lw=1.3,
            color=col, label="%s (median %.0f)" % (lab, v.median()))
    HEAD["b_median_length_%s_aa" % lab.replace(" ", "_")] = int(round(float(v.median())))
    HEAD["b_n_%s_with_length" % lab.replace(" ", "_")] = int(len(v))
    _dn, _ = np.histogram(v.clip(upper=40), bins=bins, density=True)
    _ct, _ = np.histogram(v.clip(upper=40), bins=bins)
    _brows += [dict(cls=lab, bin_left=float(bins[_i]), bin_right=float(bins[_i + 1]),
                    density=float(_dn[_i]), n=int(_ct[_i])) for _i in range(len(bins) - 1)]
fd.panel("figS2", "b", _brows,
         note="one row per histogram bin per class; shortest wild-type tryptic peptide, clipped at 40 aa")
ax.set_xlabel("shortest tryptic peptide (residues)", fontsize=st.sz(6))
ax.set_ylabel("density", fontsize=st.sz(6))
ax.legend(frameon=False, fontsize=st.sz(5.0), loc="upper right", handlelength=1.1)

# Panel c of the previous version repeated Figure 2b exactly — the same two bars, 7.56%
# conventional against 7.28% mutation-aware, from the same file. One of the two had to go, and
# the main text keeps its gate-1 framing, so it goes from here. The remaining panels move up
# one letter.

# ---------------------------------------------------------------- c,e: missed cleavage
rep = find("missed_cleavage_analysis_report.txt")
psm, seq = {}, {}
nseq, npsm = {}, {}                 # the counts behind the two rates, for the panel deposit
acq = bg = np.nan
n_acq = n_bg = np.nan               # the denominators behind panel d
if rep:
    txt = open(rep).read()
    for r in ("K", "R"):
        m = re.search(r"^\s+%s\s+([\d,]+)\s+([\d.]+)%%\s+([\d,]+)\s+([\d.]+)%%" % r, txt, re.M)
        if m:
            seq[r] = float(m.group(2)); psm[r] = float(m.group(4))
            nseq[r] = int(m.group(1).replace(",", "")); npsm[r] = int(m.group(3).replace(",", ""))
    mv = re.search(r"PSM-weighted is the comparable pair: ([\d.]+)% vs ([\d.]+)%", txt)
    if mv:
        acq, bg = float(mv.group(1)), float(mv.group(2))
    for _ln in txt.splitlines():                    # the PSM counts behind those two rates
        _pp = _ln.split()
        if len(_pp) > 1 and _pp[0] == "background" and _pp[1][:1].isdigit():
            n_bg = int(_pp[1].replace(",", ""))
        if len(_pp) > 2 and _pp[:2] == ["acquired", "lysine"] and _pp[2][:1].isdigit():
            n_acq = int(_pp[2].replace(",", ""))

ax = fig.add_subplot(gs[1, 0:5]); st.head(fig, ax, "c")
if psm:
    x = np.arange(2); w = st.grp(2)
    ax.bar(x - w/2, [psm["K"], psm["R"]], w, color=st.BLUE, label="PSM-weighted")
    ax.bar(x + w/2, [seq["K"], seq["R"]], w, color=st.LIGHTGREY, label="distinct sequences")
    for i, r in enumerate(("K", "R")):
        ax.text(i - w/2, psm[r] + 0.7, "%.1f%%" % psm[r], ha="center", fontsize=st.sz(5.8), color=st.BLUE)
        ax.text(i + w/2, seq[r] + 0.7, "%.1f%%" % seq[r], ha="center", fontsize=st.sz(5.8),
                color=st.GREY)
    ax.set_xticks(x)
    ax.set_xticklabels(["first basic residue is K", "first basic residue is R"], fontsize=st.sz(5.8))
    ax.set_ylabel("peptides carrying a missed\ncleavage at that site (%)", fontsize=st.sz(6))
    ax.set_ylim(0, 44)
    ax.legend(frameon=False, fontsize=st.sz(5.4), loc="upper left", handlelength=1.0)
    for _r in ("K", "R"):
        HEAD["c_missed_cleavage_first_basic_%s_psm_pct" % _r.lower()] = float(psm[_r])
        HEAD["c_missed_cleavage_first_basic_%s_seq_pct" % _r.lower()] = float(seq[_r])
    fd.panel("figS2", "c", [
        dict(first_basic_residue=_r, weighting=_w, missed_cleavage_pct=float(_v), n=_n)
        for _r in ("K", "R")
        for _w, _v, _n in (("PSM-weighted", psm[_r], npsm.get(_r)),
                           ("distinct sequences", seq[_r], nseq.get(_r)))],
        note="missed / (missed + cleaved) among identified peptides, split by first basic residue")

ax = fig.add_subplot(gs[1, 6:12]); st.head(fig, ax, "d")
if not np.isnan(acq):
    ax.bar([0, 1], [bg, acq], st.BARW, color=[st.LIGHTGREY, st.BLUE])
    for i, v in enumerate([bg, acq]):
        ax.text(i, v + 0.6, "%.1f%%" % v, ha="center", fontsize=st.sz(6.2))
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["all identified peptides",
                        "peptides carrying an acquired lysine"], fontsize=st.sz(5.8))
    ax.set_ylabel("peptides carrying a missed\ncleavage (%)", fontsize=st.sz(6))
    ax.set_ylim(0, max(acq, bg) * 1.42)
    ax.text(0.5, 0.90, "acquired-lysine peptides are digested like any other",
            transform=ax.transAxes, ha="center", fontsize=st.sz(5.6), color=st.GREY)
    HEAD["d_missed_cleavage_acquired_lysine_pct"] = float(acq)
    HEAD["d_missed_cleavage_background_pct"] = float(bg)
    fd.panel("figS2", "d", [
        dict(group="all identified peptides", missed_cleavage_pct=float(bg), n_psms=n_bg),
        dict(group="peptides carrying an acquired lysine", missed_cleavage_pct=float(acq),
             n_psms=n_acq)], note="PSM-weighted, the comparable pair")

# the note belongs directly under panels d and e. 0.415 put it across f and g; 0.560 pushed
# it into the x tick labels of c and d at 83%, far worse than the 21% clip of panel f's
# letter that it was meant to fix. The gap between these two rows is not tall enough for a
# two-line note, so this is the least-bad position until the note is shortened.
# The 30-word definition that sat here is in the legend now. There was never a good place
# for it: the gap between these two rows is not tall enough for a two-line note, and both
# positions tried put it through the x tick labels of c and d. It is a definition, which is
# what a legend is for.

# ---------------------------------------------------------------- d: score margin
ax = fig.add_subplot(gs[2, 0:5]); st.head(fig, ax, "e")
H = find("hits62.tsv")
if H:
    h = pd.read_csv(H, sep="\t")
    # is_kgain, not a regex on the first term only: a record like HSPA8 T273A;R77K at site 77
    # carries the K-gain in its SECOND term, and a first-term test drops it. That is the
    # difference between 52 and the 53 sites Figure 2 reports.
    h = h[[is_kgain(c, s_) for c, s_ in zip(h.change, h.site)]]
    d = pd.to_numeric(h.delta, errors="coerce").dropna()
    ax.hist(d.clip(upper=45), bins=np.arange(0, 47, 3), color=st.BLUE)
    ax.axvline(MARGIN, color=st.ORANGE, lw=1.2, ls="--")
    ax.text(0.97, 0.97, "margin \u2265 %d retained" % MARGIN,
            transform=ax.transAxes, ha="right",
            va="top", fontsize=st.sz(5.2), color=st.ORANGE)
    ax.set_xlabel("hyperscore \u2212 nextscore", fontsize=st.sz(6))
    ax.set_ylabel("PSMs", fontsize=st.sz(6))
    nsite = h.groupby(["gene", "change", "line", "site"]).ngroups
    ax.text(0.97, 0.80, "%d of %d PSMs clear the threshold,\nacross %d sites"
            % (int((d >= MARGIN).sum()), len(d), nsite), transform=ax.transAxes, ha="right",
            va="top", fontsize=st.sz(5.2), color=st.GREY, linespacing=1.4)
    HEAD["e_margin_threshold"] = int(MARGIN)
    HEAD["e_n_psms_total"] = int(len(d))
    HEAD["e_n_psms_above_margin"] = int((d >= MARGIN).sum())
    HEAD["e_n_sites"] = int(nsite)
    HEAD["e_median_margin"] = round(float(d.median()), 2)
    fd.panel("figS2", "e",
             h.assign(margin=pd.to_numeric(h.delta, errors="coerce"))[
                 ["plex", "gene", "change", "line", "site", "peptide", "spectrum",
                  "hyper", "nextscore", "margin"]],
             note="one row per acquired-lysine PSM; margin = hyperscore - nextscore, the binned quantity")

# ---------------------------------------------------------------- e: spectrum quality
ax = fig.add_subplot(gs[2, 6:12]); st.head(fig, ax, "f")
# Are variant identifications supported by worse spectra than ordinary ones? If so, the
# detections would be suspect regardless of FDR.
import glob as _g
vh, bh = [], []
for f in sorted(_g.glob("out_plex*_pep/*/psm.tsv"))[:12]:
    try:
        q = pd.read_csv(f, sep="\t", low_memory=False,
                        usecols=lambda c: c in ("Hyperscore", "Protein"))
    except Exception:
        continue
    if "Hyperscore" not in q.columns:
        continue
    pr = q.get("Protein", pd.Series([""] * len(q))).astype(str)
    vh.append(q.loc[pr.str.startswith("VAR|"), "Hyperscore"])
    bg_ = q.loc[~pr.str.contains(r"VAR\||DOWN\||BOUNDARY\||rev_", regex=True, na=False),
                "Hyperscore"]
    bh.append(bg_.sample(min(len(bg_), 20000), random_state=0))
if vh and bh:
    V = pd.concat(vh).dropna(); B = pd.concat(bh).dropna()
    bins = np.linspace(0, max(V.max(), B.quantile(0.99)), 26)
    ax.hist(B, bins=bins, density=True, color=st.LIGHTGREY, label="background (%s)"
            % f"{len(B):,}")
    ax.hist(V, bins=bins, density=True, histtype="step", lw=1.3, color=st.BLUE,
            label="variant (%d)" % len(V))
    ax.axvline(B.median(), color=st.GREY, lw=0.9, ls="--")
    ax.axvline(V.median(), color=st.BLUE, lw=0.9, ls="--")
    ax.set_xlabel("hyperscore", fontsize=st.sz(6)); ax.set_ylabel("density", fontsize=st.sz(6))
    ax.legend(frameon=False, fontsize=st.sz(5.0), loc="upper right", handlelength=1.1)
    ax.text(0.97, 0.55, "medians %.1f vs %.1f" % (V.median(), B.median()),
            transform=ax.transAxes, ha="right", fontsize=st.sz(5.2), color=st.GREY)
    HEAD["f_median_hyperscore_variant"] = round(float(V.median()), 1)
    HEAD["f_median_hyperscore_background"] = round(float(B.median()), 1)
    HEAD["f_n_variant_psms"] = int(len(V))
    HEAD["f_n_background_psms"] = int(len(B))
    _frows = []
    for _lab, _vv in (("variant", V), ("background", B)):
        _dn, _ = np.histogram(_vv, bins=bins, density=True)
        _ct, _ = np.histogram(_vv, bins=bins)
        _frows += [dict(series=_lab, bin_left=float(bins[_i]), bin_right=float(bins[_i + 1]),
                        density=float(_dn[_i]), n=int(_ct[_i])) for _i in range(len(bins) - 1)]
    fd.panel("figS2", "f", _frows,
             note="one row per histogram bin per series; background sampled at 20,000 PSMs per search output")

# ---------------------------------------------------------------- f: an annotated spectrum
# Score distributions and false-discovery estimates are arguments about a population. For at
# least one site a reviewer should see the fragment ions, and specifically the ions containing
# the mutated residue — those are what distinguish the mutant peptide from its wild-type
# counterpart. Everything else is inference.
try:
    import figS_spectrum as SP
    axl = fig.add_subplot(gs[3:5, 0:12]); axl.set_axis_off()
    st.head(fig, axl, "g")
    sub = gs[3:5, 0:12].subgridspec(2, 1, height_ratios=[0.22, 1.05], hspace=0.02)
    SP.draw_spectrum(fig, sub[0], sub[1])
    # the panel's numbers are read back from the module that draws it, not retyped here
    _bseen = {_r.resi for _r in SP.M.itertuples() if _r.series == "b"}
    _yseen = {_r.resi for _r in SP.M.itertuples() if _r.series == "y"}
    HEAD["g_peptide"] = SP.PEPTIDE
    HEAD["g_scan"] = int(SP.SCAN)
    HEAD["g_n_peaks"] = int(len(SP.mz))
    HEAD["g_base_peak_intensity"] = float(SP.inten.max())
    HEAD["g_n_ions_matched"] = int(len(SP.M))
    HEAD["g_n_theoretical_ions"] = int(len(SP.frags))
    HEAD["g_match_tolerance_da"] = float(SP.TOL_DA)
    HEAD["g_n_backbone_bonds_observed"] = int(len(_bseen | _yseen))
    HEAD["g_n_backbone_bonds"] = int(len(SP.PEPTIDE) - 1)
    fd.panel("figS2", "g", SP.M,
             note="matched fragment ions for %s, scan %d" % (SP.PEPTIDE, SP.SCAN))
except Exception as e:
    ax = fig.add_subplot(gs[3:5, 0:12]); ax.set_axis_off()
    ax.text(.5, .5, "spectrum panel unavailable: %s" % str(e)[:70], ha="center",
            transform=ax.transAxes, fontsize=st.sz(6))

# headline numbers per panel, printed so no legend holds a number the script does not report
try:
    print("c: missed cleavage at the site, first basic K %.1f%% PSM / %.1f%% seq | first basic R %.1f%% / %.1f%%"
          % (psm["K"], seq["K"], psm["R"], seq["R"]))
except Exception as _e: print("c: not built (%s)" % _e)
try:
    print("d: missed cleavage overall, acquired-lysine peptides %.1f%% vs background %.1f%%" % (acq, bg))
except Exception as _e: print("d: not built (%s)" % _e)
try:
    print("e: score margin, %d of %d PSMs clear margin>=10, across %d sites"
          % (int((d >= 10).sum()), len(d), nsite))
except Exception as _e: print("e: not built (%s)" % _e)
try:
    print("f: hyperscore, variant median %.1f vs background median %.1f (n %d vs %d)"
          % (V.median(), B.median(), len(V), len(B)))
except Exception as _e: print("f: not built (%s)" % _e)
st.save(fig, "figS2")
fd.head("figS2", HEAD)
for _e in (".png", ".pdf"):          # review copies alongside fig4_v2 / figS_react_v2
    fig.savefig("../results/figS2_v3" + _e, dpi=400, bbox_inches="tight")
print("acquired lysine %s sites, other variant %s | gate1 %.1f%% vs %.1f%% | gate2 %.1f%% vs %.1f%%"
      % (f"{len(K):,}", f"{len(N):,}", g1[0], g1[1], g2[0], g2[1]))
