# Absent by construction, not by chemistry: the lysines cancer creates

Code and deposited data for the manuscript. Everything the paper reports is either produced by a
script here or read from a public source named in `external_inputs.json`.

## What this deposit is for

The Methods claim that *"every number in this paper is printed by the script that produces the
corresponding panel rather than transcribed from the artwork."* That claim is checkable here, in
both directions, and the checks are part of the deposit rather than something you have to take on
trust.

## Two tiers

**Tier 1 — redraw any figure from deposited data.** `panel_data/` holds one tidy CSV per figure
panel plus an `<figure>_head.json` ledger of the numbers that figure's caption and text quote.
This tier needs nothing but the files in this deposit.

**Tier 2 — regenerate the panel data from public inputs.** The analysis scripts in `ccle_plex05/`
and `ligand_rebuild/` rebuild those CSVs from the sources in `external_inputs.json`. This tier
needs those downloads, and for the proteomic re-search it needs the search engine as well
(MSFragger 4.4.1 / FragPipe v24.0 / Philosopher).

Some intermediates could not be traced to a producing script by static analysis — the filename is
assembled at run time or written through a config constant — so they ship as data under `data/`
and are labelled `deposited intermediate` in `external_inputs.json`. Tier 1 works for those;
tier 2 does not, and the deposit says so rather than implying otherwise.

## Layout

```
ccle_plex05/          the proteomic re-search, Figures 1-3 and S1-S4, S10
ccle_plex05/manuscript/  the manuscript build: prose, captions, references, assembly
ligand_rebuild/       the chemoproteomic ligandability analysis, Figures 4-8 and S5-S9
panel_data/           one CSV per panel + one JSON ledger per figure  (tier 1)
data/                 analysis intermediates that ship as data         (tier 1)
repo/                 the manifest, the release builder, and the checks
external_inputs.json  every input at the edge, with its source and version
```

## The checks

Run from `repo/`.

| script | what it answers |
|---|---|
| `build_manifest.py` | which script draws each figure, what it reads, what it writes, which numbers it prints |
| `build_release.py` | assembles this deposit; **fails** if any input is undeclared; reports output-filename collisions |
| `verify.py` | ledger → prose: does every number a script computes still appear in the manuscript? |
| `reverse_check.py` | prose → ledger: is there a number in the text that no script produces? |
| `../ccle_plex05/manuscript/seams_check.py` | citations, figure pointers, note numbering, recurring quantities |
| `acceptance.py` | does the assembled deposit actually reproduce the deposited ledgers? |

`verify.py` and `reverse_check.py` are not symmetric. A ledger value missing from the prose is
usually fine — scripts compute more than a paper quotes. A prose number that no ledger produces is
the dangerous direction, and `reverse_check.py` exists because that direction cannot be spotted by
reading.

Current state: 18 of 18 figures carry a ledger; 692 of 717 ledger values (97%) appear in the text;
`seams_check.py` reports 0 issues.

## Does the deposit actually reproduce the paper? (`acceptance.py`)

`acceptance.py` runs each figure script inside the assembled tree and diffs the ledger it writes
against the ledger deposited with the paper. Not "the code is here" — "the code here gives the
same answers". Result on the deposited tree:

| verdict | n | figures |
|---|---|---|
| **MATCH** — every deposited value reproduced | 11 | 2, 4, 5, 7, 8, S1, S5, S6, S7, S8, S9 |
| **PARTIAL** — reproduced values agree, some panels skipped | 4 | 1, 3, S2, S3 |
| **NOT RUN** — needs an input not in the deposit | 3 | 6, S4, S10 |

Nothing DIFFERS. Every gap is a named input that is public but not redistributed here:

| figure | needs | why not included |
|---|---|---|
| 1 (panel e), S10 | `OmicsSomaticMutations.csv` | DepMap 26Q1 — free download, not ours to mirror |
| 3 (panel d) | `SpecQuality/sq_out/` | 364 MB of per-spectrum features for one panel |
| S2, S3 | `out_plex*/psm.tsv` | raw search output, several GB |
| 6 | `SMILES.xlsx` | compound structures from three papers' SI |
| S4 | `plex*_varpair.txt` | 319 MB of entrapment pairings |

**Tier 1 covers all of them**: every panel above can be redrawn from `panel_data/` with no
download at all. The gaps are only in re-deriving that data from scratch.

Two of these were found by this test rather than by reading, and both were silent:

- Figure 1g reported `NaN` for the germline comparison because its two `.out` inputs were invisible
  to the manifest's file-extension list.
- Figure 8e reported **0 peptides** for AKT1 E17K, when the answer is 3. Zero there reads as "the
  site was never searched", which is the opposite of the paper's claim. The script now refuses to
  emit a count it cannot compute, and the search databases ship so that it can.

## Known output-filename collisions

Static analysis finds seven filenames written by two different shipping scripts. They are listed
by `build_release.py` on every run. None affects the deposited panel data, but if you re-run the
pipeline the later script wins, so the authoritative producer is named here:

| file | authoritative producer | note |
|---|---|---|
| `step29_report.txt` | `step29_ceiling_fast.py` | the deposited report is the fast variant; figure S6 parses it |
| `step21_report.txt`, `step21_methods.txt` | `step21_union_characterise.py` | |
| `step28_report.txt` | directory-scoped: each writes into its own folder | |
| `step32_report.txt` | `step32_akt1_site.py` | |
| `features.csv`, `step6_report.txt` | `step6_features_v2.py` | v1 ships only as a producer of other files |

This class of bug is not hypothetical: `step31_crystal_validation.py` used to overwrite
`step31_compound_site_matching.py`'s report, and one arm of that analysis was lost until it was
re-run on 7 September 2026. The matching script's report has since been renamed, and
`build_release.py` now checks for collisions on every build.

## Reproducing a figure

```bash
cd ccle_plex05
python3 fig2_v3.py          # writes fig2.png/.pdf, panel CSVs and fig2_head.json
```

Each figure script prints its own headline numbers and deposits them. Shared plotting style lives
in `abbasov_style.py`; the deposit helper is `figure_data.py`.

## Licence

Code: MIT (`LICENSE`). Data: CC BY 4.0 (`LICENSE-DATA`). Several inputs carry upstream terms that
this licence does not supersede — COSMIC in particular restricts commercial use. See
`LICENSE-DATA`.

## Citation

See `CITATION.cff`. The archived release has its own DOI; cite that alongside the paper.
