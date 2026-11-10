# A mutation-aware search of 375 cancer cell-line proteomes recovers 53 of 8,470 acquired-lysine observations

Code and deposited data for the manuscript. Everything the paper reports is either produced by a
script here or read from a public source named in `external_inputs.json`.

## What this deposit is for

The Methods claim that *"every number in this paper is printed by the script that produces the
corresponding panel rather than transcribed from the artwork."* That claim is checkable here, in
both directions, and the checks are part of the deposit rather than something you have to take on
trust.

## Two tiers

**Tier 1. Redraw any figure from deposited data.** `panel_data/` holds one tidy CSV per figure
panel plus an `<figure>_head.json` ledger of the numbers that figure's caption and text quote.
This tier needs nothing but the files in this deposit.

**Tier 2. Regenerate the panel data from public inputs.** The analysis scripts in `ccle/`
and `ligand_rebuild/` rebuild those CSVs from the sources in `external_inputs.json`. This tier
needs those downloads, and for the proteomic re-search it needs the search engine as well
(MSFragger 4.4.1 / FragPipe v24.0 / Philosopher).

Twenty intermediates could not be traced to a producing script by static analysis, because the
filename is assembled at run time or written through a config constant. They are labelled
`deposited intermediate` in `external_inputs.json`. They are this project's own outputs rather
than public downloads, they total 343 MB, and they are not redistributed. Tier 1 covers every
figure that reads one; tier 2 does not, and the deposit says so rather than implying otherwise.

## What this repository holds

The code, `panel_data/` and `external_inputs.json`. It does not carry the 343 MB of deposited
intermediates or the 120 MB of variant and decoy search databases; `.gitignore` says so, and the
table of gaps below says what that costs.

Tier 1 therefore runs from a clone of this repository with no download at all. Tier 2 runs for
every step whose inputs are public and named in `external_inputs.json`, and not for the steps
that read an undeposited intermediate.

## Layout

```
ccle/                 the proteomic re-search, Figures 1-3 and S1-S4, S10
ccle/manuscript/      the manuscript build: prose, captions, references, assembly
ligand_rebuild/       the chemoproteomic ligandability analysis, Figures 4-8 and S5-S9
panel_data/           one CSV per panel + one JSON ledger per figure  (tier 1)
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
| `../ccle/manuscript/seams_check.py` | citations, figure pointers, note numbering, recurring quantities |
| `acceptance.py` | does the assembled deposit actually reproduce the deposited ledgers? |

`verify.py` and `reverse_check.py` are not symmetric. A ledger value missing from the prose is
usually fine: scripts compute more than a paper quotes. A prose number that no ledger produces is
the dangerous direction, and `reverse_check.py` exists because that direction cannot be spotted by
reading.

Current state: 19 ledgers cover 18 figures and the graphical abstract; 713 of 736 ledger values
(97%) appear in the text;
`seams_check.py` reports 0 issues.

## Does the deposit actually reproduce the paper? (`acceptance.py`)

`acceptance.py` runs each figure script inside the assembled tree and diffs the ledger it writes
against the ledger deposited with the paper. Not "the code is here" but "the code here gives the
same answers". Measured on an assembled tree that also held the undeposited intermediates, so
these verdicts are the ceiling rather than what a clone of this repository reaches:

| verdict | n | figures |
|---|---|---|
| **MATCH**, every deposited value reproduced | 11 | 2, 4, 5, 7, 8, S1, S5, S6, S7, S8, S9 |
| **PARTIAL**, reproduced values agree, some panels skipped | 4 | 1, 3, S2, S3 |
| **NOT RUN**, needs an input not in the deposit | 3 | 6, S4, S10 |

Nothing DIFFERS. Every gap is a named input that is public but not redistributed here:

| figure | needs | why not included |
|---|---|---|
| 1 (panel e), S10 | `OmicsSomaticMutations.csv` | DepMap 26Q1; free download, not ours to mirror |
| 3 (panel d) | `SpecQuality/sq_out/` | 364 MB of per-spectrum features for one panel |
| S2, S3 | `out_plex*/psm.tsv` | raw search output, several GB |
| 6 | `SMILES.xlsx` | compound structures from three papers' SI |
| S4 | `plex*_varpair.txt` | 319 MB of entrapment pairings |

One class of gap is different in kind. The twenty `deposited intermediate` entries in
`external_inputs.json` are this project's own outputs, not public downloads, so where a step reads
one there is no download that closes the gap and tier 2 cannot be run from this deposit at all.

**Tier 1 covers all of them**: every panel above can be redrawn from `panel_data/` with no
download at all. The gaps are only in re-deriving that data from scratch.

Two of these were found by this test rather than by reading, and both were silent:

- Figure 1g reported `NaN` for the germline comparison because its two `.out` inputs were invisible
  to the manifest's file-extension list.
- Figure 8e reported **0 peptides** for AKT1 E17K, when the answer is 3. Zero there reads as "the
  site was never searched", which is the opposite of the paper's claim. The script now refuses to
  emit a count it cannot compute rather than print a zero it cannot stand behind.

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

## Published numbering is not the numbering in these filenames

The paper has seven figures and this repository has eight sets of scripts and ledgers, because
source figures 5 and 6 are merged into one published Figure 5. `panel_data/figure_numbering.json`
is written by the manuscript build and is the authority.

| published | script | ledgers and panel CSVs |
|---|---|---|
| Graphical abstract | `ccle/ga_v1.py` | reads the ledgers below, deposits none |
| Figure 1 to 4 | `ccle/fig1_v3.py` to `fig4_v3.py` | `fig1_*` to `fig4_*` |
| Figure 5 | `ccle/fig5pub_v1.py` | `fig5_*` and `fig6_*` |
| Figure 6 | `ccle/fig7_v2.py` | `fig7_*` |
| Figure 7 | `ccle/fig8_v2.py` | `fig8_*` |
| Figure S1 to S10 | `ccle/figS1_v3.py` to `figS10_v1.py` | `figS1_*` to `figS10_*` |

`fig5pub_v1.py` draws no panels of its own. It reads the six surviving panel blocks out of
`fig5_v2.py` and `fig6_v2.py` at run time and composes them, so those two scripts ship as well and
are what the ledgers are keyed on. Six of their twelve panels are not in the published figure;
their text is in Extended Results in the Supporting Information.

Anything in this file written before 9 November 2026 that names a figure by number, including the
Figure 8e note above, uses the source numbering.

## Reproducing a figure

```bash
cd ccle
python3 fig2_v3.py          # writes fig2.png/.pdf, panel CSVs and fig2_head.json
```

Each figure script prints its own headline numbers and deposits them. Shared plotting style lives
in `abbasov_style.py`; the deposit helper is `figure_data.py`.

## Licence

Code: MIT (`LICENSE`). Data: CC BY 4.0 (`LICENSE-DATA`). Several inputs carry upstream terms that
this licence does not supersede; COSMIC in particular restricts commercial use. See
`LICENSE-DATA`.

## Citation

See `CITATION.cff`, and cite the paper.
