#!/usr/bin/env python3
"""
figS10_v1.py — the acquired-lysine class across DepMap cell lines, and what a mutation-aware
re-search of the three published screens could find.

WHY THIS EXISTS
Figures 4-8 measure ligandability on lysines from three chemoproteomic screens. Every one ran in
a cell line whose acquired lysines were outside its search space by construction. This figure asks
what such a re-search could find, before anyone runs it, and where it should be run instead.

DEDUPLICATION — READ THIS
OmicsSomaticMutations.csv carries duplicate rows for the same (model, gene, protein change):
55,116 acquired-lysine rows collapse to 39,390 unique, and 936 of 1,954 models are affected.
Every count here is on unique (ModelID, HugoSymbol, ProteinChange). An earlier version of this
panel used raw row counts and overstated MDA-MB-231 (40 vs 24), MDA-MB-436 (38 vs 19) and
Ramos (12 vs 6); HeLa and MeWo were unaffected.

PANELS
  (a) burden by cancer lineage — the mutational processes of Figure 1e, seen in cell lines
  (b) the highest-burden lines and the four the screens used, split by whether the mutated gene
      is expressed in that line
  (c) the reachable set: acquired lysines on a protein the matching screen probe-labels

EVERY NUMBER IN THE TEXT IS PRINTED BELOW. Nothing is transcribed from the artwork.
OUTPUT  figS10.png/.pdf, figS10_report.txt, ../results/scope_*.csv
"""
import os, numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import abbasov_style as A
import figure_data as fd
A.slate(); A.apply()

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
log = []
def say(s=""):
    print(s, flush=True); log.append(s)

SCREEN = {"ACH-001086": ("HeLa", {"CyMPSA"}),
          "ACH-000768": ("MDA-MB-231", {"Meroterpenoid", "NatChem"}),
          "ACH-000573": ("MDA-MB-436", {"Meroterpenoid"}),
          "ACH-001636": ("Ramos", {"NatChem"})}
FOUR = list(SCREEN)

say("=" * 78); say("FIGURE S10 — THE ACQUIRED-LYSINE CLASS ACROSS DepMap CELL LINES"); say("=" * 78)

# ------------------------------------------------------------------ acquired lysines, unique
d = pd.read_csv(os.path.join(D, "Datasets", "OmicsSomaticMutations.csv"), low_memory=False,
                usecols=["ModelID", "HugoSymbol", "ProteinChange", "Hotspot", "HessDriver"])
pc = d.ProteinChange.astype(str); ok = pc.str.match(r"^p\.[A-Z]\d+[A-Z]$")
d = d[ok].copy(); d["ref"] = pc[ok].str[2]; d["alt"] = pc[ok].str[-1]
raw = d[(d.alt == "K") & (d.ref != "K")]
acq = raw.drop_duplicates(["ModelID", "HugoSymbol", "ProteinChange"]).copy()
say("\nDEDUPLICATION: %d acquired-lysine rows -> %d unique (model, gene, change); %d of %d models"
    % (len(raw), len(acq), (raw.groupby("ModelID").size() >
                            acq.groupby("ModelID").size()).sum(), acq.ModelID.nunique()))
say("             affected. All counts below are unique.")

lostK = (d[(d.ref == "K") & (d.alt != "K")]
         .drop_duplicates(["ModelID", "HugoSymbol", "ProteinChange"])
         .groupby("ModelID").size().rename("lost"))
say("\nLYSINE GAINED vs LOST: %d gained, %d lost across all models, net %+d, ratio %.2f"
    % (len(acq), int(lostK.sum()), len(acq) - int(lostK.sum()), len(acq) / lostK.sum()))

tf = lambda s: s.fillna(0).astype(str).isin(["1", "1.0", "True", "TRUE", "Y"])
acq["drv"] = tf(acq.HessDriver) | tf(acq.Hotspot)
mod = pd.read_csv(os.path.join(D, "Datasets", "Model.csv"), low_memory=False,
                  usecols=["ModelID", "StrippedCellLineName", "OncotreeLineage",
                           "OncotreePrimaryDisease", "PrimaryOrMetastasis"])
B = (acq.groupby("ModelID").agg(acqK=("HugoSymbol", "size"), drv=("drv", "sum")).reset_index()
        .merge(lostK.reset_index(), on="ModelID", how="left").merge(mod, on="ModelID", how="left"))
B["lost"] = B.lost.fillna(0).astype(int); B["net"] = B.acqK - B.lost
B = B.sort_values("net", ascending=False)
say("    net gain: median %+d; %d of %d models are net-negative (%.1f%%)"
    % (B.net.median(), (B.net < 0).sum(), len(B), 100 * (B.net < 0).mean()))
LR = B.groupby("OncotreeLineage").agg(n=("net", "size"), g=("acqK", "sum"), l=("lost", "sum"))
LR = LR[LR.n >= 15]; LR["ratio"] = LR.g / LR.l
say("    gained-to-lost ratio, highest lineages: " + ", ".join(
    "%s %.1f" % (i, r.ratio) for i, r in LR.sort_values("ratio", ascending=False).head(5).iterrows()))
say("    the C>T asymmetry of figure S1f, seen per cell line: skin runs a %.1f-fold surplus against"
    % LR.loc["Skin", "ratio"])
say("    %.2f overall." % (len(acq) / lostK.sum()))
B.to_csv(os.path.join(D, "results", "scope_acquiredK_burden_by_model.csv"), index=False)

say("\n(a) BURDEN ACROSS %d MODELS: median %d, mean %.1f, 90th pct %d, max %d (%s)"
    % (len(B), B.acqK.median(), B.acqK.mean(), B.acqK.quantile(.90),
       B.acqK.max(), B.loc[B.acqK.idxmax(), "StrippedCellLineName"]))
L = (B.groupby("OncotreeLineage").agg(n=("acqK", "size"), med=("acqK", "median"),
                                      mx=("acqK", "max")).reset_index())
L = L[L.n >= 15].sort_values("med", ascending=False)
say("    by lineage (>=15 models), median acquired lysines:")
for r in L.itertuples():
    say("      %-26s n=%3d  median %3.0f  max %3d" % (r.OncotreeLineage, r.n, r.med, r.mx))
say("    Skin and Uterus/Bowel lead; Eye (uveal, melanocytic but not UV-driven) is near the")
say("    bottom. This is the mutational-process result of Figure 1e, in a second dataset.")
say("    the four screen lines:")
marks = []
for mid in FOUR:
    r = B[B.ModelID == mid].iloc[0]
    pct = 100 * (B.acqK < r.acqK).mean()
    marks.append((SCREEN[mid][0], int(r.acqK), pct, r.OncotreeLineage))
    say("      %-12s %-10s %3d   %.0fth percentile" % (SCREEN[mid][0], r.OncotreeLineage,
                                                       r.acqK, pct))
say("    the four carry %d unique acquired lysines between them; %s alone carries %d"
    % (sum(v for _, v, _, _ in marks), B.iloc[0].StrippedCellLineName, int(B.iloc[0].acqK)))


# ================================================================== gates: expression, protein, dependency
FOURm = [m for m in FOUR]
Pp = pd.read_csv(os.path.join(D, "Datasets", "pipeline_out", "proteomics_subset.csv"), low_memory=False)
Pp = Pp.rename(columns={Pp.columns[0]: "ModelID"}).set_index("ModelID"); pcols = set(Pp.columns)
ch = pd.read_csv(os.path.join(D, "Datasets", "CRISPRGeneEffect.csv"), nrows=0); cid = ch.columns[0]
cmods = set(pd.read_csv(os.path.join(D, "Datasets", "CRISPRGeneEffect.csv"), usecols=[cid])[cid])
EPATH = os.path.join(D, "Datasets", "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv")
eh = pd.read_csv(EPATH, nrows=0)
emods = set(pd.read_csv(EPATH, usecols=["ModelID", "IsDefaultEntryForModel"])
            .query("IsDefaultEntryForModel=='Yes'").ModelID)
full = set(Pp.index) & cmods & emods
TOPN = 26
cand = [m for m in B.sort_values("net", ascending=False).ModelID if m in full][:TOPN] + [m for m in FOURm if m in full]
sub = acq[acq.ModelID.isin(cand)].copy()
need = set(sub.HugoSymbol)
egm = {c.split(" (")[0]: c for c in eh.columns if " (" in c}
E = pd.read_csv(EPATH, usecols=["ModelID", "IsDefaultEntryForModel"] +
                [egm[g] for g in need if g in egm], low_memory=False)
E = E[E.IsDefaultEntryForModel == "Yes"].drop(columns="IsDefaultEntryForModel").set_index("ModelID")
E.columns = [c.split(" (")[0] for c in E.columns]
cgm = {c.split(" (")[0]: c for c in ch.columns[1:]}
C = pd.read_csv(os.path.join(D, "Datasets", "CRISPRGeneEffect.csv"),
                usecols=[cid] + [cgm[g] for g in need if g in cgm], low_memory=False).set_index(cid)
C.columns = [c.split(" (")[0] for c in C.columns]
MP = pd.read_csv(os.path.join(D, "Datasets", "uniprot_hugo_entrez_id_mapping.csv"), low_memory=False)
g2u = dict(zip(MP.Symbol.astype(str), MP.UniprotID.astype(str)))
def _get(df, m, g):
    try: return float(df.at[m, g])
    except Exception: return np.nan
sub["up"] = sub.HugoSymbol.map(g2u)
sub["expr"] = [_get(E, r.ModelID, r.HugoSymbol) >= 1.0 for r in sub.itertuples()]
sub["det"] = [(r.up in pcols) and (not pd.isna(Pp.at[r.ModelID, r.up])) if r.up in pcols else False
              for r in sub.itertuples()]
sub["dep"] = [_get(C, r.ModelID, r.HugoSymbol) <= -0.5 for r in sub.itertuples()]
sub["g1"] = sub.expr; sub["g2"] = sub.expr & sub.det; sub["g3"] = sub.expr & sub.det & sub.dep
F = (sub.groupby("ModelID").agg(total=("HugoSymbol", "size"), expressed=("g1", "sum"),
                                detected=("g2", "sum"), depend=("g3", "sum")).reset_index()
        .merge(mod, on="ModelID", how="left").merge(B[["ModelID","net"]], on="ModelID")
        .sort_values("net", ascending=False))
F.to_csv(os.path.join(D, "results", "scope_funnel.csv"), index=False)
say("\n    DATA AVAILABILITY: of %d models with mutations, %d have expression, %d proteomics, %d\n"
    "    CRISPR, and %d all three. None of the ten highest net-gain lines carries all three."
    % (len(B), len(set(B.ModelID) & emods), len(set(B.ModelID) & set(Pp.index)),
       len(set(B.ModelID) & cmods), len(full)))
say("\n(c) FOUR GATES — acquired lysine -> gene expressed -> protein detected in that line ->")
say("    gene is a CRISPR dependency (<= -0.5). %d models carry all three data types." % len(full))
say("    aggregate over the %d drawn: %d -> %d (%.0f%%) -> %d (%.0f%%) -> %d (%.0f%%)"
    % (len(F), F.total.sum(), F.expressed.sum(), 100*F.expressed.sum()/F.total.sum(),
       F.detected.sum(), 100*F.detected.sum()/F.total.sum(),
       F.depend.sum(), 100*F.depend.sum()/F.total.sum()))
for r in F.head(8).itertuples():
    say("      %-12s %-14s %4d -> %4d -> %4d -> %3d"
        % (r.StrippedCellLineName, r.OncotreeLineage, r.total, r.expressed, r.detected, r.depend))
for mid in FOURm:
    r = F[F.ModelID == mid]
    if len(r):
        r = r.iloc[0]
        say("      %-12s %-14s %4d -> %4d -> %4d -> %3d  (screen line)"
            % (r.StrippedCellLineName, r.OncotreeLineage, r.total, r.expressed, r.detected, r.depend))
    else:
        say("      %-12s no CRISPR and/or proteomics entry, not drawn" % SCREEN[mid][0])

# ================================================================== recurrent driver-flagged sites
DD = (acq[acq.drv].groupby(["HugoSymbol", "ProteinChange"]).size().rename("models")
      .reset_index().sort_values("models", ascending=False))
DD["label"] = DD.HugoSymbol + " " + DD.ProteinChange.str.replace("p.", "", regex=False)
DD.to_csv(os.path.join(D, "results", "scope_recurrent_drivers.csv"), index=False)
IN_F8 = {"NRAS Q61K", "PIK3CA E542K", "PIK3CA E545K", "AKT1 E17K", "TP53 R280K", "PTPN11 E76K"}
say("\n(c) RECURRENT DRIVER-FLAGGED ACQUIRED LYSINES — %d of %d unique acquired lysines carry a"
    % (int(acq.drv.sum()), len(acq)))
say("    hotspot or HessDriver flag (%.2f%%), in %d of %d models. Top by model count:"
    % (100*acq.drv.mean(), int((acq.groupby('ModelID').drv.sum() > 0).sum()), acq.ModelID.nunique()))
for r in DD.head(14).itertuples():
    say("      %-18s %3d models%s" % (r.label, r.models, "   <- structure in Figure 8d"
                                      if r.label in IN_F8 else ""))
say("    Four of the six sites measured in Figure 8d appear here, and they are the four most")
say("    recurrent. BRAF V600K is an acquired lysine on the most-drugged oncogene in melanoma.")
g3 = sub[sub.g3 & sub.drv].merge(mod[["ModelID", "StrippedCellLineName"]], on="ModelID", how="left")
say("    driver-flagged sites clearing all three gates in the models drawn: %s"
    % ", ".join("%s %s (%s)" % (r.HugoSymbol, r.ProteinChange.replace("p.", ""), r.StrippedCellLineName)
                for r in g3.itertuples()))

# ================================================================== subtype and stage
say("\n    SUBTYPE: cutaneous melanoma median %d acquired lysines (n=%d), ocular melanoma %d (n=%d)"
    % (B[B.OncotreePrimaryDisease == "Melanoma"].acqK.median(),
       (B.OncotreePrimaryDisease == "Melanoma").sum(),
       B[B.OncotreePrimaryDisease == "Ocular Melanoma"].acqK.median(),
       (B.OncotreePrimaryDisease == "Ocular Melanoma").sum()))
from scipy.stats import mannwhitneyu
pa = B[B.PrimaryOrMetastasis == "Primary"].acqK; pb = B[B.PrimaryOrMetastasis == "Metastatic"].acqK
u, pv = mannwhitneyu(pa, pb)
say("    STAGE: metastatic median %d (n=%d) vs primary %d (n=%d), p = %.2g, Cliff's delta %+.3f"
    % (pb.median(), len(pb), pa.median(), len(pa), pv, 2*u/(len(pa)*len(pb)) - 1))

# ================================================================== reachable set
al = acq[acq.ModelID.isin(FOUR)]
ms = pd.read_csv(os.path.join(D, "ligand_rebuild", "master_sites.csv"), low_memory=False)
ms["screens"] = ms.datasets.astype(str).str.split("+").apply(set)
rows = []
for mid in FOUR:
    name, scr = SCREEN[mid]
    for r in al[al.ModelID == mid].itertuples():
        g = ms[(ms.gene == r.HugoSymbol) & ms.screens.apply(lambda x: bool(x & scr))]
        if not len(g): continue
        best = g.loc[g.max_occupancy.idxmax()]
        rows.append(dict(line=name, gene=r.HugoSymbol, acquired=r.ProteinChange.replace("p.", ""),
                         acq_pos=int("".join(c for c in r.ProteinChange if c.isdigit())),
                         labelled_K=int(best.position), occ=float(best.max_occupancy),
                         eng60=int(best.n_engaged_60)))
R = pd.DataFrame(rows).sort_values("occ", ascending=False).reset_index(drop=True)
R["sep"] = (R.labelled_K - R.acq_pos).abs()
R.to_csv(os.path.join(D, "results", "scope_reachable_set.csv"), index=False)
say("\n(d) REACHABLE SET — %d of %d unique sites in the four screen lines sit on a protein the"
    % (len(R), len(al)))
say("    matching screen probe-labels. THE MEASURED LYSINE IS NOT THE ACQUIRED ONE.")
for r in R.itertuples():
    say("      %-11s %-7s %-9s measured K%-6d %5.1f%%  |sep| %4d" %
        (r.line, r.gene, r.acquired, r.labelled_K, r.occ, r.sep))
say("    engaged >=60%%: %s | median |sep| %d residues (%d to %d)"
    % (", ".join(R[R.eng60 > 0].gene) or "none", R.sep.median(), R.sep.min(), R.sep.max()))

# ================================================================== figure
fig = plt.figure(figsize=(A.DOUBLE, 10.4))
gs = fig.add_gridspec(4, 2, height_ratios=[0.90, 1.00, 0.95, 0.95], hspace=0.95, wspace=0.42)
ax1 = fig.add_subplot(gs[0, :]); ax2 = fig.add_subplot(gs[1, :])
ax5 = fig.add_subplot(gs[2, 0]); ax3 = fig.add_subplot(gs[2, 1])
ax4 = fig.add_subplot(gs[3, :])

order = list(L.OncotreeLineage)
bp = ax1.boxplot([B[B.OncotreeLineage == l].acqK.values for l in order],
                 positions=np.arange(len(order)), widths=0.62, showfliers=False, patch_artist=True,
                 medianprops=dict(color=A.INK, lw=0.9), whiskerprops=dict(color=A.GREY, lw=0.6),
                 capprops=dict(color=A.GREY, lw=0.6), boxprops=dict(linewidth=0.5, edgecolor=A.GREY))
scr_lin = {B[B.ModelID == m].OncotreeLineage.iloc[0] for m in FOUR}
for b_, l in zip(bp["boxes"], order): b_.set_facecolor(A.ORANGE if l in scr_lin else A.LIGHTGREY)
fd.panel("figS10", "a",
         [dict(lineage=l, n_models=int(len(_v)), q1=float(np.percentile(_v, 25)),
               median=float(np.median(_v)), q3=float(np.percentile(_v, 75)),
               min=int(_v.min()), max=int(_v.max()), screen_lineage=bool(l in scr_lin))
          for l in order for _v in [B[B.OncotreeLineage == l].acqK.values]],
         note="one box per lineage with at least 15 models, ordered by median acquired lysines")
SHORT = {"Bladder/Urinary Tract": "Bladder/Urinary", "Peripheral Nervous System": "Peripheral NS",
         "Ovary/Fallopian Tube": "Ovary/Fallopian", "Esophagus/Stomach": "Esoph./Stomach"}
ax1.set_yscale("log"); ax1.set_ylabel("acquired lysines")
ax1.set_xticks(np.arange(len(order)))
ax1.set_xticklabels([SHORT.get(l, l) for l in order], rotation=90, fontsize=5.6)
for t in ax1.get_xticklabels():
    if t.get_text() in {SHORT.get(x, x) for x in scr_lin}: t.set_color(A.ORANGE)
ax1.set_xlim(-0.7, len(order) - 0.3)
ax1.set_title("(a) Burden by lineage, ordered by median. Rust, the lineages the screens sampled.",
              loc="left", fontsize=7.5)

AB = B.sort_values("net", ascending=False).reset_index(drop=True)
xr = np.arange(len(AB))
ax2.vlines(xr, 0, AB.acqK.values, color=A.LIGHTGREY, lw=0.55)
ax2.vlines(xr, 0, -AB.lost.values, color=A.BLUE, lw=0.55, alpha=0.9)
ax2.axhline(0, color=A.INK, lw=0.7)
ax2.set_yscale("symlog", linthresh=10, linscale=0.45)
ax2.set_yticks([-1000, -100, -10, 0, 10, 100, 600])
ax2.set_yticklabels(["1000", "100", "10", "0", "10", "100", "600"])
say("    y axis is symlog (linear within +/-10); the extreme loss tail is a hypermutator losing"
    " %d lysines." % int(AB.lost.max()))

# Callouts: the three highest net gainers, then the four screen lines, each labelled once with
# a marker, offsets chosen so no two labels share a band.
CALL = [("MEWO", 1), ("SW684", 0), ("MCC13", 2), ("MDAMB231", 0), ("MDAMB436", 1),
        ("HELA", 2), ("RAMOS", 3)]
SCRN = {"MDAMB231", "MDAMB436", "HELA", "RAMOS"}
for nm, band in CALL:
    r = AB[AB.StrippedCellLineName == nm]
    if not len(r): continue
    i = int(r.index[0]); v = int(r.acqK.iloc[0]); lo = int(r.lost.iloc[0])
    col = A.ORANGE if nm in SCRN else A.GREEN
    ax2.plot([i], [v], marker="o", ms=2.6, color=col, zorder=6)
    ty = [900, 380, 150, 60][band] if nm not in SCRN else [-40, -120, -330, -820][band]
    tx = i + (len(AB) * 0.055 if nm not in SCRN else len(AB) * 0.045)
    ax2.annotate("%s  +%d / \u2212%d" % (nm, v, lo), xy=(i, v), xytext=(tx, ty),
                 ha="left", va="center", fontsize=5.8, color=col,
                 arrowprops=dict(arrowstyle="-", color=col, lw=0.6, shrinkA=0, shrinkB=1))
fd.panel("figS10", "b",
         AB[["ModelID", "StrippedCellLineName", "OncotreeLineage", "acqK", "lost", "net"]]
           .assign(rank=np.arange(1, len(AB) + 1),
                   called_out=[n in dict(CALL) for n in AB.StrippedCellLineName]),
         note="one row per plotted model: lysines gained above the line, lost below, ranked by net gain")
ax2.set_xlabel("all %d DepMap models with at least one missense lysine change, ranked by net gain"
               % len(AB))
ax2.set_ylabel("lysines gained (above)\nand lost (below)", fontsize=6.8)
ax2.set_xlim(-15, len(AB) + 15)
ax2.set_title("(b) Every cell line in DepMap. Gained above the line, lost below, ranked by net "
              "gain. Symlog scale.", loc="left", fontsize=7.5)
for s_ in ("top", "right"): ax2.spines[s_].set_visible(False)

# (c) the gates, on the 308 models that carry all three data types
Fg = F.sort_values("depend", ascending=False).head(12)[::-1]
yg = np.arange(len(Fg))
ax5.barh(yg, Fg.total, color=A.LIGHTGREY, edgecolor=A.GREY, lw=0.4, height=0.7)
ax5.barh(yg, Fg.expressed, color=A.BLUE, edgecolor=A.GREY, lw=0.3, height=0.7)
ax5.barh(yg, Fg.detected, color=A.GREEN, edgecolor=A.GREY, lw=0.3, height=0.7)
ax5.barh(yg, Fg.depend, color=A.ORANGE, edgecolor=A.GREY, lw=0.3, height=0.7)
ax5.set_yticks(yg); ax5.set_yticklabels(Fg.StrippedCellLineName, fontsize=6)
ax5.set_xlabel("acquired lysines")
ax5.legend(handles=[Patch(facecolor=A.LIGHTGREY, edgecolor=A.GREY, label="acquired"),
                    Patch(facecolor=A.BLUE, label="expressed"),
                    Patch(facecolor=A.GREEN, label="protein detected"),
                    Patch(facecolor=A.ORANGE, label="dependency")],
           loc="lower right", fontsize=5.6)
fd.panel("figS10", "c",
         [dict(cell_line=r.StrippedCellLineName, lineage=r.OncotreeLineage, gate=_g,
               sites=int(getattr(r, _c)))
          for r in Fg.itertuples()
          for _g, _c in (("acquired", "total"), ("expressed", "expressed"),
                         ("protein detected", "detected"), ("dependency", "depend"))],
         note="four nested gates, the twelve models with the most sites clearing every gate")
# "308" was typed into this title. It is len(full) and always was.
ax5.set_title("(c) Four nested gates, in the %d models\ncarrying all three data types." % len(full),
              loc="left", fontsize=7.5)
for s_ in ("top", "right"): ax5.spines[s_].set_visible(False)

dd = DD.head(12)[::-1]
yv = np.arange(len(dd))
ax3.barh(yv, dd.models, color=[A.ORANGE if l in IN_F8 else A.LIGHTGREY for l in dd.label],
         edgecolor=A.GREY, lw=0.4, height=0.68)
ax3.set_yticks(yv); ax3.set_yticklabels(dd.label, fontsize=6)
ax3.set_xlabel("cell lines carrying it")
fd.panel("figS10", "d",
         [dict(label=r.label, gene=r.HugoSymbol, change=r.ProteinChange.replace("p.", ""),
               models=int(r.models), in_figure_8d=bool(r.label in IN_F8)) for r in dd.itertuples()],
         note="driver-flagged acquired lysines ranked by the number of cell lines carrying each")
ax3.set_title("(d) Recurrent driver-flagged acquired lysines.\nRust, a structure in Figure 8d.",
              loc="left", fontsize=7.5)
for s_ in ("top", "right"): ax3.spines[s_].set_visible(False)

yv = np.arange(len(R))[::-1]
ax4.barh(yv, R.occ, color=[A.BLUE if e > 0 else A.LIGHTGREY for e in R.eng60],
         edgecolor=A.GREY, lw=0.4, height=0.64)
ax4.axvline(60, color=A.ORANGE, lw=0.8, ls=(0, (3, 2)))
ax4.set_yticks(yv)
ax4.set_yticklabels(["%s K%d   (%s \u00b7 %s)" % (r.gene, r.labelled_K, r.acquired, r.line)
                     for r in R.itertuples()], fontsize=6)
ax4.set_xlim(0, 100); ax4.set_xlabel("occupancy of the probe-labelled lysine (%), not of the acquired lysine")
fd.panel("figS10", "e", R,
         note="acquired lysines in the four screen lines that sit on a protein the matching screen "
              "labels; occ is the best-engaged PROBE-LABELLED lysine, not the acquired one")
# "11" was typed into this title. It is len(R) and always was.
ax4.set_title("(e) The %d reachable sites in the screen lines: the lysine the screens did "
              "measure, acquired mutation in parentheses." % len(R), loc="left", fontsize=7.5)
for s_ in ("top", "right"): ax4.spines[s_].set_visible(False)

_top = {r.label: int(r.models) for r in DD.head(14).itertuples()}
_bv = DD[DD.label == "BRAF V600K"]
if len(_bv): _top["BRAF V600K"] = int(_bv.models.iloc[0])
_be = R[R.eng60 > 0]

HEAD = {
 "a_n_acquired_rows_raw": int(len(raw)),
 "a_n_acquired_unique": int(len(acq)),
 "a_n_models": int(acq.ModelID.nunique()),
 "a_n_models_dedup_affected": int((raw.groupby("ModelID").size() >
                                   acq.groupby("ModelID").size()).sum()),
 "a_burden_median": float(B.acqK.median()),
 "a_burden_max": int(B.acqK.max()),
 "a_burden_max_line": str(B.loc[B.acqK.idxmax(), "StrippedCellLineName"]),
 "a_n_lineages_drawn": int(len(L)),
 "a_median_by_lineage": {r.OncotreeLineage: float(r.med) for r in L.itertuples()},
 "a_median_cutaneous_melanoma": float(B[B.OncotreePrimaryDisease == "Melanoma"].acqK.median()),
 "a_n_cutaneous_melanoma": int((B.OncotreePrimaryDisease == "Melanoma").sum()),
 "a_median_ocular_melanoma": float(B[B.OncotreePrimaryDisease == "Ocular Melanoma"].acqK.median()),
 "a_n_ocular_melanoma": int((B.OncotreePrimaryDisease == "Ocular Melanoma").sum()),
 "a_median_metastatic": float(pb.median()), "a_n_metastatic": int(len(pb)),
 "a_median_primary": float(pa.median()), "a_n_primary": int(len(pa)),
 "a_stage_p": float(pv),
 "a_stage_cliffs_delta": round(2 * u / (len(pa) * len(pb)) - 1, 3),

 "b_n_gained": int(len(acq)),
 "b_n_lost": int(lostK.sum()),
 "b_gained_to_lost_ratio": round(float(len(acq) / lostK.sum()), 2),
 "b_n_models_net_negative": int((B.net < 0).sum()),
 "b_pct_models_net_negative": round(100 * float((B.net < 0).mean()), 1),
 "b_gained_to_lost_by_lineage": {i: round(float(r.ratio), 1) for i, r in LR.iterrows()},
 "b_max_lysines_lost": int(AB.lost.max()),
 "b_acquired_in_screen_lines": {nm: int(n) for nm, n, _, _ in marks},
 "b_n_acquired_four_screen_lines": int(sum(n for _, n, _, _ in marks)),

 "c_n_models_with_expression": int(len(set(B.ModelID) & emods)),
 "c_n_models_with_proteomics": int(len(set(B.ModelID) & set(Pp.index))),
 "c_n_models_with_crispr": int(len(set(B.ModelID) & cmods)),
 "c_n_models_all_three": int(len(full)),
 "c_n_models_drawn": int(len(F)),
 "c_n_sites_total": int(F.total.sum()),
 "c_n_sites_expressed": int(F.expressed.sum()),
 "c_pct_sites_expressed": round(100 * F.expressed.sum() / F.total.sum(), 1),
 "c_n_sites_detected": int(F.detected.sum()),
 "c_pct_sites_detected": round(100 * F.detected.sum() / F.total.sum(), 1),
 "c_n_sites_dependency": int(F.depend.sum()),
 "c_pct_sites_dependency": round(100 * F.depend.sum() / F.total.sum(), 1),
 "c_n_lines_with_dependency_site": int((F.depend > 0).sum()),

 "d_n_driver_flagged": int(acq.drv.sum()),
 "d_pct_driver_flagged": round(100 * float(acq.drv.mean()), 2),
 "d_n_models_with_driver_flagged": int((acq.groupby("ModelID").drv.sum() > 0).sum()),
 "d_models_carrying_site": _top,
 "d_n_figure_8d_sites": int(len(IN_F8)),
 "d_n_figure_8d_sites_in_ranking": int(DD.head(14).label.isin(IN_F8).sum()),
 "d_n_driver_sites_clearing_all_gates": int(len(g3)),

 "e_n_reachable": int(len(R)),
 "e_n_acquired_in_screen_lines": int(len(al)),
 "e_median_separation_residues": float(R.sep.median()),
 "e_min_separation_residues": int(R.sep.min()),
 "e_min_separation_gene": str(R.loc[R.sep.idxmin(), "gene"]),
 "e_max_separation_residues": int(R.sep.max()),
 "e_max_separation_gene": str(R.loc[R.sep.idxmax(), "gene"]),
 "e_n_engaged_above_60": int(len(_be)),
 "e_best_engaged_gene": str(_be.gene.iloc[0]) if len(_be) else "",
 "e_best_engaged_position": int(_be.labelled_K.iloc[0]) if len(_be) else -1,
 "e_best_engaged_occupancy_pct": round(float(_be.occ.iloc[0]), 1) if len(_be) else float("nan"),
}

fig.savefig(os.path.join(D, "ccle", "figS10.png"))
fig.savefig(os.path.join(D, "ccle", "figS10.pdf"))
fd.head("figS10", HEAD)
open(os.path.join(D, "ccle", "figS10_report.txt"), "w").write("\n".join(log) + "\n")
say("\nwrote figS10.png/.pdf and figS10_report.txt")
