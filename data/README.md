# Site-level tables

These fifteen tables hold the site-level results behind the paper. `panel_data/` holds what each figure panel plots;
the tables here hold the rows those panels are computed from, so the sites can be searched, filtered and joined directly.

| file | rows | columns | one row is | read by |
|---|---:|---:|---|---|
| `union_sites_gnomad_v2.tsv` | 295,411 | 11 | the acquired-lysine site catalog; one row per site in the union of DepMap, TCGA MC3 and COSMIC | ccle/fig1_v3.py, ccle/figS1_v3.py; written by ccle/step28_gnomad_rerun.py |
| `region_compare_conv.tsv` | 10,372 | 15 | per-site search outcome against the conventional MANE database; one row per variant site per cell line | ccle/fig2_v3.py, ccle/figS3_v3.py |
| `region_compare_pep.tsv` | 10,372 | 15 | per-site search outcome against the mutation-aware peptide database; same rows as the conventional table | ccle/fig2_v3.py, ccle/figS2_v3.py, ccle/figS3_v3.py, ccle/peptide_properties.py |
| `acquired_K_searched_vs_detected.csv` | 8,182 | 15 | per-site searched and detected flags for acquired lysines; one row per site per cell line | ccle/fig8_v2.py |
| `hits62.tsv` | 99 | 14 | per-spectrum detail for the acquired-lysine detections | ccle/fig2_v3.py, ccle/figS2_v3.py, ccle/figS3_v3.py, ccle/analyse_downstream.py |
| `fig3_margins.csv` | 1,489 | 4 | per-peptide score margins from the entrapment searches | ccle/fig3_v3.py, ccle/figS4_v3.py |
| `fdp_p2.csv` | 79,365 | 12 | per-peptide q-values and entrapment false-discovery proportions for plex 2 | ccle/figS4_v3.py |
| `residue_flux_by_line.csv` | 1,968 | 59 | per-cell-line counts of residues gained and lost by missense mutation | ccle/fig2_v3.py, ccle/figS3_v3.py |
| `features_final_full.csv` | 45,552 | 28 | per-site structural features for the lysines measured in the three chemoproteomic screens, before any confidence gate | ccle/test_behaviour.py |
| `step43_enclosure.csv` | 246,900 | 3 | per-site enclosure computed on the AlphaFold models | ccle/figS8_v1.py, ligand_rebuild/step43b_enclosure_tests.py |
| `step27_modnz_union.csv` | 155,008 | 23 | environment features at modeled side-chain amine positions for union sites | ligand_rebuild/step28_three_ways.py, ligand_rebuild/step28b_sensitivity.py |
| `step27_modnz_measured.csv` | 30,805 | 23 | the same features at measured lysines, with the error of the modeled amine position | ccle/figS9_v1.py, ligand_rebuild/step28_three_ways.py, ligand_rebuild/step28b_sensitivity.py |
| `step10_roc_splits.csv` | 201 | 4 | ROC band of the ligandability composite across protein-held-out splits | ccle/fig5_v2.py, ccle/figS6_v1.py |
| `step10_roc_pooled.csv` | 201 | 2 | pooled ROC curve of the ligandability composite | ccle/fig5_v2.py, ccle/figS6_v1.py |
| `step37_within_NatChem.csv` | 84 | 7 | compound descriptor against site-preference correlations within the lysine-targeting compound screen | ccle/fig6_v2.py, ccle/figS7_v1.py |

Row counts exclude the header line. SHA-256 checksums are in `SHA256SUMS.txt`.

## The acquired-lysine site catalog (`union_sites_gnomad_v2.tsv`)

295,411 sites across 20,365 genes. A site is a missense position whose substituted residue is lysine.

| column | meaning |
|---|---|
| `gene` | gene symbol |
| `ref_aa` | reference residue |
| `position` | residue number of the substituted position |
| `alt_aa` | substituted residue, `K` in every row |
| `in_depmap`, `in_tcga`, `in_cosmic` | 1 if the site is present in DepMap Public 26Q1, TCGA MC3 v0.2.8-PUBLIC, COSMIC Cancer Mutation Census v103 |
| `match_status` | `in_dbnsfp_with_AF`, `in_dbnsfp_no_gnomad` or `not_in_dbnsfp` (not assessed) |
| `gnomad_AF` | gnomAD allele frequency taken from dbNSFP 5.4a |
| `gnomad_noncancer_AF` | allele frequency in the gnomAD non-cancer subset, from the same record |
| `band` | the frequency band drawn in Figure 1C |

Sites by source: DepMap only 22,717; TCGA only 11,604; COSMIC only 167,321; DepMap and TCGA 147; DepMap and COSMIC 6,547; TCGA and COSMIC 83,937; all three 3,138.

## Licence and upstream terms

The tables are this study's own derived results and are released under CC BY 4.0 (`LICENSE-DATA`), subject to the upstream terms listed there.

COSMIC. The catalog includes 260,943 sites listed in the COSMIC Cancer Mutation Census v103 (Wellcome Sanger Institute), 167,321 of them listed by COSMIC alone. They are disclosed here to enable publication of the accompanying paper. COSMIC data remain subject to COSMIC's own licence terms, which CC BY 4.0 does not replace; commercial use requires a licence from the Wellcome Sanger Institute.

## Not included

Six larger intermediates named in `external_inputs.json` are not redistributed: `fig3_deltamass.csv` (97 MB), `proteomics_subset.csv` (64 MB), `depmap_variants_annotated.parquet` (48 MB), `mc3_bysample.tsv` (35 MB), `mc3_nt.tsv` (31 MB) and `sites_covalent_residues.csv` (29 MB). The variant and decoy search databases are not redistributed either.
