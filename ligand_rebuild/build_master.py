#!/usr/bin/env python3
"""
build_master.py — one searchable master table: every site-compound measurement, with the site
annotation and the compound annotation carried alongside.

WHY A LONG TABLE RATHER THAN A MATRIX
occ_long.csv already holds every measurement, but on its own it answers only "what occupancy did
compound X give at site Y". Joining the atlas onto it makes every row self-describing, so the
file can be filtered on any column — gene, position, compound, occupancy, pathway, pocket,
cell line, dataset — without a second lookup. That is what makes it searchable the way a
spreadsheet is.

The cost is redundancy: site annotation repeats for every compound tested there. At 1.6 M rows
that is a large file, which is why it is CSV and gzipped rather than xlsx (Excel stops at
1,048,576 rows and would silently truncate).

WHAT EACH ROW IS
One measurement: one compound, at one site, in one experiment. A site tested against 15
compounds has 15 rows. A site detected as two peptide species has rows for each.

WHAT IS AND IS NOT IN IT — this governs every question the file can answer
  - ZEROS ARE REAL. 890,456 of 1.64 M measurements are exactly 0% occupancy: tested, not
    engaged. Absence of a row is NOT zero.
  - COMPLETE WITHIN 10 OF 11 EXPERIMENTS. In CyMPSA and the meroterpenoids every detected site
    was tested against every compound in that experiment. In the single NatChem experiment
    (328 compounds) only 24% of site-compound pairs were measured, so a missing pair there
    means NOT MEASURED, not "not engaged".
  - ONLY DETECTED LYSINES. A lysine that never produced a quantifiable peptide has no row at
    all. The file describes 32,586 OBSERVED lysines, not every lysine in these proteins.
    "Not in the file" must never be read as "not ligandable".
  - REACTIVITY IS NOT LIGANDABILITY. Occupancy measures blocking, and a molecule can block by
    reacting rather than by binding. Ribosomal proteins engage at 19% against 9.3% elsewhere,
    which survives depth matching, but their structural profile matches reactivity. Low hit
    rate can mean selectivity OR a weak electrophile; the two are not separable here.

OUTPUT
  master_site_compound.csv.gz   every measurement, fully annotated  (~1.6 M rows)
  master_sites.csv              one row per site, the atlas plus peptide and sample columns
  master_compounds.csv          one row per compound, with hit rate
  MASTER_README.md              what each column means and what the file cannot answer
"""
import os, sys, gzip
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")

def find(n):
    for d in (".", "..", os.path.expanduser("~/Desktop/DepMap/ligand_rebuild")):
        p = os.path.join(d, n)
        if os.path.exists(p): return p
    return None

print("reading ...")
OCC = pd.read_csv(find("occ_long.csv"), low_memory=False)
ATL = pd.read_csv(find("lysine_ligandability_atlas_v2.csv"), comment="#", low_memory=False)
print("  occ_long %s rows | atlas %s sites" % (f"{len(OCC):,}", f"{len(ATL):,}"))

# ---------------------------------------------------------------- compound table + hit rate
print("building compound table ...")
C = OCC.groupby("compound").agg(
        datasets=("dataset", lambda s: ";".join(sorted(set(s)))),
        n_experiments=("experiment", "nunique"),
        n_sites_tested=("position", "size"),
        n_sites_engaged_60=("occupancy_pct", lambda s: int((s >= 60).sum())),
        n_sites_engaged_80=("occupancy_pct", lambda s: int((s >= 80).sum())),
        n_sites_engaged_90=("occupancy_pct", lambda s: int((s >= 90).sum())),
        median_occupancy=("occupancy_pct", "median"),
        max_occupancy=("occupancy_pct", "max")).reset_index()
# hit rate: the single most useful derived column, and it reorders the table completely.
# A compound engaging 486 of 4,974 sites is not "strong", it is unselective.
for t in (60, 80, 90):
    C["hit_rate_%d" % t] = 100 * C["n_sites_engaged_%d" % t] / C.n_sites_tested
C = C.sort_values("n_sites_engaged_60", ascending=False)

# SMILES, if the atlas carries them for the best-compound column
if "best_compound_SMILES" in ATL.columns and "best_compound" in ATL.columns:
    sm = (ATL.dropna(subset=["best_compound_SMILES"])
             .drop_duplicates("best_compound")[["best_compound", "best_compound_SMILES"]]
             .rename(columns={"best_compound": "compound", "best_compound_SMILES": "SMILES"}))
    C = C.merge(sm, on="compound", how="left")
print("  %d compounds, hit rate 60%%: median %.2f%%, range %.2f-%.2f%%"
      % (len(C), C.hit_rate_60.median(), C.hit_rate_60.min(), C.hit_rate_60.max()))

# ---------------------------------------------------------------- site table
print("building site table ...")
# peptide, cell line, fraction and probe come from occ_long, which the atlas dropped
meta = (OCC.sort_values("occupancy_pct", ascending=False)
           .groupby(["uniprot", "position"])
           .agg(peptide=("peptide", "first"),
                cell_lines=("cell_line", lambda s: ";".join(sorted(set(map(str, s))))),
                cell_fractions=("cell_fraction", lambda s: ";".join(sorted(set(map(str, s))))),
                probes=("probe", lambda s: ";".join(sorted(set(map(str, s))))),
                experiments=("experiment", "nunique"),
                compounds_engaged_60=("compound", lambda s: "")).reset_index())
eng = (OCC[OCC.occupancy_pct >= 60].groupby(["uniprot", "position"]).compound
          .agg(lambda s: ";".join(sorted(set(s)))).rename("compounds_engaged_60").reset_index())
meta = meta.drop(columns=["compounds_engaged_60"]).merge(eng, on=["uniprot", "position"],
                                                         how="left")
S = ATL.merge(meta, on=["uniprot", "position"], how="left")
S["sample_type"] = np.where(S.cell_fractions.astype(str).str.contains("soluble|lysate",
                                                                      case=False, na=False),
                            "lysate", "live cell / other")
print("  %s sites, %d columns" % (f"{len(S):,}", len(S.columns)))

# ---------------------------------------------------------------- the master long table
print("joining ...")
site_cols = [c for c in S.columns if c not in ("peptide",)]
M = OCC.merge(S[site_cols], on=["uniprot", "position"], how="left", suffixes=("", "_site"))
comp_cols = ["compound", "hit_rate_60", "n_sites_tested", "n_sites_engaged_60"]
M = M.merge(C[comp_cols].rename(columns={
                "n_sites_tested": "compound_n_sites_tested",
                "n_sites_engaged_60": "compound_n_sites_engaged_60",
                "hit_rate_60": "compound_hit_rate_60"}), on="compound", how="left")

# rank of this site inside the compound's own footprint — the column that separates a ligand
# from a reactive lysine, as far as this data can
M["rank_within_compound"] = M.groupby("compound").occupancy_pct.rank(ascending=False,
                                                                    method="min")
M["engaged_60"] = (M.occupancy_pct >= 60).astype(int)

front = ["gene", "uniprot", "position", "compound", "occupancy_pct", "engaged_60",
         "dataset", "experiment", "probe", "cell_line", "cell_fraction", "peptide",
         "rank_within_compound", "compound_hit_rate_60"]
M = M[[c for c in front if c in M.columns] + [c for c in M.columns if c not in front]]
M = M.sort_values(["gene", "position", "occupancy_pct"], ascending=[True, True, False])

print("writing ...")
M.to_csv("master_site_compound.csv.gz", index=False, compression="gzip")
S.to_csv("master_sites.csv", index=False)
C.to_csv("master_compounds.csv", index=False)
print("  master_site_compound.csv.gz  %s rows x %d columns" % (f"{len(M):,}", len(M.columns)))
print("  master_sites.csv             %s rows x %d columns" % (f"{len(S):,}", len(S.columns)))
print("  master_compounds.csv         %s rows x %d columns" % (f"{len(C):,}", len(C.columns)))
