#!/usr/bin/env python3
"""
step2_occupancy.py — put all three lysine ligandability datasets on one scale.

BUILT FROM SCRATCH. Reads only the three published supplementary workbooks. Reuses no
previously generated intermediate. Asperugins is excluded — it is cysteine data.

THE THREE SOURCES
  Nature Chemistry 2021   14,101 lysines, reported as a RATIO (control / treated),
                          328 compound x concentration columns.
                          LIMITATION: Tab 3 records how many isoTOP-ABPP datasets each
                          compound was run in per cell line (Ramos, MDA-MB-231), but Tab 2
                          carries ONE column per compound x concentration. The values are
                          therefore already aggregated across cell lines and replicates.
                          NatChem CANNOT be stratified by cell line; its stratum is
                          protein-level only. CyMPSA and Meroterpenoid can be stratified
                          by experiment and cell line.
  JACS CyMPSA             6 experiments, RAW TMT ABUNDANCES, 18 channels each:
                          3 DMSO + 5 compounds x 3 replicates.
                          2 probes (dTB-STP, Probe 1) x 3 compound sets (Ac-1..Ac-15), HeLa
  JACS Meroterpenoids     4 experiments, RAW TMT ABUNDANCES, same 18-channel design,
                          MDA-MB-231 and MDA-MB-436, named natural products

CONVERSION TO PERCENT OCCUPANCY — occupancy is the fraction of that lysine blocked by
the compound, so that all three datasets answer the same question.

  ratio data (Nat Chem):     occupancy% = 100 * (1 - 1/R)
                             R = 4  ->  75%   (the published "liganded" cut)
  abundance data (both JACS): occupancy% = 100 * (1 - mean(compound reps) / mean(DMSO reps))

DECISIONS MADE HERE (agreed before writing)
  - replicate channels for the same compound at the same concentration are AVERAGED first,
    then converted, giving one occupancy per compound per site
  - enantiomers count as TWO INDEPENDENT COMPOUNDS ((+)- and (-)-guadial B are separate),
    because a site engaged by one and not its mirror image is engaged by BINDING, not by
    intrinsic reactivity — which is the central claim being tested
  - occupancy is clipped to [0, 100]. Negative values (compound channel above DMSO) are
    noise or enrichment artefacts, not negative blocking; clipping to 0 is conservative
    and keeps the scale interpretable
  - nothing is thresholded here. Cases and controls are defined in step 5, so the
    threshold choice stays visible and separable

OUTPUT
  occ_long.csv   one row per site x experiment x compound: the atomic measurement
  step2_report.txt  what was read, what was dropped, and why

Nothing downstream should read the source workbooks again.
"""
import re, sys, os
import numpy as np, pandas as pd
import warnings
warnings.filterwarnings("ignore")

SRC = os.environ.get("SRC_DIR", ".")

def find(*patterns):
    """Locate a workbook by any of several published filenames."""
    import glob
    for p in patterns:
        hits = glob.glob(os.path.join(SRC, p))
        if hits:
            return hits[0]
    raise SystemExit("cannot find any of %s in %s" % (patterns, SRC))

NATCHEM = find("NatChem*data*sites*.xlsx", "NatChemdatasites.xlsx")
CYMPSA  = find("CyMPSA*rawdata*sites*.xlsx", "CyMPSArawdatasites.xlsx")
MERO    = find("Meroterpenoid*rawdata*sites*.xlsx", "Meroterpenoidsrawdatasites.xlsx")
print("sources:\n  %s\n  %s\n  %s" % (NATCHEM, CYMPSA, MERO))

log = []
def say(m):
    print(m)
    log.append(m)


def occ_from_ratio(r):
    r = pd.to_numeric(r, errors="coerce")
    return (100.0 * (1.0 - 1.0 / r)).clip(0, 100).where(r > 0)


def occ_from_abundance(cmpd_mean, dmso_mean):
    with np.errstate(divide="ignore", invalid="ignore"):
        o = 100.0 * (1.0 - cmpd_mean / dmso_mean)
    return pd.Series(o).clip(0, 100).where(pd.Series(dmso_mean) > 0)


# --------------------------------------------------------------- Nature Chemistry
say("=" * 72)
say("NATURE CHEMISTRY 2021  (ratio data)")
d = pd.read_excel(NATCHEM, "Tab 2 | related to 2b")
ann = list(d.columns)[:7]
val = list(d.columns)[7:]
say("  %d lysine sites, %d compound x concentration columns" % (len(d), len(val)))

recs = []
for c in val:
    o = occ_from_ratio(d[c])
    m = o.notna()
    if not m.any():
        continue
    recs.append(pd.DataFrame({
        "dataset": "NatChem", "experiment": "NatChem|isoTOP-ABPP",
        "probe": "P1", "cell_line": "aggregated",
        "uniprot": d.loc[m, "uniprot"].values,
        "position": pd.to_numeric(d.loc[m, "probe P1-labeled site"], errors="coerce").values,
        "peptide": d.loc[m, "peptide sequence (K* represents probe P1-labeled sites)"].values,
        "peptide_marked": d.loc[m, "peptide sequence (K* represents probe P1-labeled sites)"].values,
        "cell_fraction": np.nan,
        "compound": c, "occupancy_pct": o[m].values}))
nat = pd.concat(recs, ignore_index=True)
say("  -> %d site x compound measurements" % len(nat))
say("     occupancy: median %.1f%%, %d at >=75%% (published cut, ratio>=4)"
    % (nat.occupancy_pct.median(), (nat.occupancy_pct >= 75).sum()))


# --------------------------------------------------------------- TMT abundance sheets
def parse_tmt(path, dataset, sheet, cell_line, probe):
    d = pd.read_excel(path, sheet)
    abund = [c for c in d.columns if "Abundance" in str(c)]
    groups = {}
    for c in abund:
        # "Abundances (Normalized): 128N, Ac-1"  /  "...: 126, Sample, DMSO"
        tail = str(c).split(":", 1)[1]
        parts = [p.strip() for p in tail.split(",")]
        label = parts[-1]                       # compound or DMSO
        # All meroterpenoids were screened at 100 uM; Tab 2 suffixes the concentration
        # while Tabs 4/5 do not, so "(+)-guadial B (100 uM)" and "(+)-guadial B" are the
        # same compound. Strip the suffix or the same compound is counted twice, which
        # would inflate the "number of compounds tried" that defines a confident negative.
        label = re.sub(r"\s*\(100\s*uM\)\s*$", "", label).strip()
        groups.setdefault(label, []).append(c)
    dmso_cols = [v for k, v in groups.items() if k.upper() == "DMSO"]
    if not dmso_cols:
        return None, "no DMSO channel"
    dmso = d[dmso_cols[0]].apply(pd.to_numeric, errors="coerce").mean(axis=1)

    accession = next(c for c in d.columns if "accession" in str(c).lower())
    site = next(c for c in d.columns if "lysine residue" in str(c).lower())
    pep = next((c for c in d.columns if "Annotated peptide" in str(c)), None)
    # the marked-lysine form (TK*VFQLK / TK[307]VFQLK) is what lets step 3 verify that
    # the reported position really is a lysine. The "Annotated ... with flanking" column
    # carries no marker, so retaining only that leaves JACS numbering unverifiable.
    pepmark = next((c for c in d.columns
                    if "probe-label" in str(c).lower() and "sequence" in str(c).lower()), None)
    frac = next((c for c in d.columns if "Cell fraction" in str(c)), None)

    out = []
    for label, cols in groups.items():
        if label.upper() == "DMSO":
            continue
        cm = d[cols].apply(pd.to_numeric, errors="coerce").mean(axis=1)   # average replicates
        o = occ_from_abundance(cm.values, dmso.values)
        m = o.notna().values
        if not m.any():
            continue
        out.append(pd.DataFrame({
            "dataset": dataset, "experiment": f"{dataset}|{sheet}",
            "probe": probe, "cell_line": cell_line,
            "uniprot": d.loc[m, accession].values,
            "position": pd.to_numeric(d.loc[m, site], errors="coerce").values,
            "peptide": d.loc[m, pep].values if pep else np.nan,
            "peptide_marked": d.loc[m, pepmark].values if pepmark else np.nan,
            "cell_fraction": d.loc[m, frac].values if frac else np.nan,
            "compound": label, "occupancy_pct": o[m].values}))
    if not out:
        return None, "no compound channels"
    return pd.concat(out, ignore_index=True), None


say("")
say("=" * 72)
say("JACS CyMPSA  (raw TMT abundances -> occupancy)")
cy = []
for sheet in pd.ExcelFile(CYMPSA).sheet_names:
    if not sheet.lower().startswith("tab"):
        continue
    probe = "dTB-STP" if "dTB-STP" in sheet else ("Probe 1" if "Probe 1" in sheet else "?")
    r, err = parse_tmt(CYMPSA, "CyMPSA", sheet, "HeLa", probe)
    if r is None:
        say("  %-34s SKIPPED (%s)" % (sheet[:34], err)); continue
    say("  %-34s %6d measurements, %2d compounds" % (sheet[:34], len(r), r.compound.nunique()))
    cy.append(r)
cy = pd.concat(cy, ignore_index=True)
say("  -> %d measurements, %d distinct compounds, %d experiments"
    % (len(cy), cy.compound.nunique(), cy.experiment.nunique()))

say("")
say("=" * 72)
say("JACS MEROTERPENOIDS  (raw TMT abundances -> occupancy)")
me = []
for sheet in pd.ExcelFile(MERO).sheet_names:
    if not sheet.lower().startswith("tab"):
        continue
    cl = "MDA-MB-231" if "231" in sheet else ("MDA-MB-436" if "436" in sheet else "?")
    r, err = parse_tmt(MERO, "Meroterpenoid", sheet, cl, "probe")
    if r is None:
        say("  %-34s SKIPPED (%s)" % (sheet[:34], err)); continue
    say("  %-34s %6d measurements, %2d compounds" % (sheet[:34], len(r), r.compound.nunique()))
    me.append(r)
me = pd.concat(me, ignore_index=True)
say("  -> %d measurements, %d distinct compounds (enantiomers counted separately)"
    % (len(me), me.compound.nunique()))
ena = sorted(c for c in me.compound.unique() if c.startswith(("(+)", "(\u2013)", "(-)")))
say("     enantiomer-labelled compounds: %s" % ", ".join(ena))


# --------------------------------------------------------------- combine
allo = pd.concat([nat, cy, me], ignore_index=True)
allo = allo[allo.position.notna()].copy()
allo["position"] = allo.position.astype(int)
allo["uniprot"] = allo.uniprot.astype(str).str.strip()

say("")
say("=" * 72)
say("COMBINED")
say("  %d site x compound measurements" % len(allo))
for ds, g in allo.groupby("dataset"):
    say("    %-14s %7d measurements | %6d distinct sites | %3d compounds | %d experiments"
        % (ds, len(g), g.groupby(["uniprot", "position"]).ngroups,
           g.compound.nunique(), g.experiment.nunique()))
say("  distinct (uniprot, position) overall: %d"
    % allo.groupby(["uniprot", "position"]).ngroups)
say("  occupancy distribution: median %.1f%%  90th %.1f%%  99th %.1f%%"
    % (allo.occupancy_pct.median(), allo.occupancy_pct.quantile(.9),
       allo.occupancy_pct.quantile(.99)))

allo.to_csv("occ_long.csv", index=False)
open("step2_report.txt", "w").write("\n".join(log) + "\n")
say("")
say("wrote occ_long.csv (%d rows) and step2_report.txt" % len(allo))
say("NOTHING IS THRESHOLDED YET — cases/controls are defined in step 5.")
