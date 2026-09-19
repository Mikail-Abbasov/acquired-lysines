# -*- coding: utf-8 -*-
"""Main-text content, one entry per section, assembled by build.py into a single manuscript.

Cite references as [[key]] or [[key1,key2]]; the builder numbers them globally in order of first
appearance. `fig` names the figure image (in ccle/) and `cap` its caption paragraphs.
`status` is "draft" for written sections and "todo" for placeholders carried from the outline.
"""

# ---------------------------------------------------------------- BMC required front/back matter
# Every item below is a stated requirement or encouragement in the Bioorganic & Medicinal Chemistry
# guide for authors, read 6 September 2026. Where the guide is explicit the requirement is quoted.

HIGHLIGHTS = dict(status="draft", spec=(
 "REQUIRED FORMAT: 3 to 5 bullet points, each a maximum of 85 characters INCLUDING SPACES, "
 "submitted as a separate editable file with \"highlights\" in the file name. Encouraged, not "
 "mandatory. All five are within the limit; character counts are appended in brackets by the build, measured rather than asserted, and are not part of the bullet."), body=[
 # Bullets 1 and 2 were merged on 7 September to free the fifth slot for the probe result
 # without losing anything: the p-value belongs in the abstract and §2.1, not in a highlight.
 "Cancer mutation creates 295,411 new lysines, 1.68-fold enriched on driver hotspots",
 "Conventional search finds none; mutation-aware 0.63%, 0.35% past the margin",
 "Reported 1% FDR is 4.70% for variant peptides without a score margin",
 "Structure predicts ligandability at AUC 0.617 against a \u2265 0.72 ceiling",
 "Probe charge steers labelling: the only 2 of 22 site features that differ",
])

GRAPHICAL_ABSTRACT = dict(status="todo", spec=(
 "REQUIRED AT SUBMISSION \u2014 the guide states \"You are required to provide a graphical abstract\". "
 "531 \u00d7 1328 pixels (h \u00d7 w) or proportionally more, readable at 5 \u00d7 13 cm, submitted as a separate "
 "file. Preferred formats TIFF, EPS, PDF or MS Office. Not yet made. The natural content is the "
 "two-gate cartoon from Figure 2a beside the AKT1 E17K contrast from Figure 8e \u2014 a lysine that "
 "can be drugged and cannot be seen."))

DECLARATIONS = [
 ("Acknowledgements", "todo",
  "Separate section immediately before the reference list \u2014 the guide is explicit that it must not "
  "sit on the title page or as a title footnote. Include anyone who helped with language, writing or "
  "proof-reading."),
 ("Author contributions", "todo",
  "REQUIRED. The corresponding author must acknowledge co-author contributions using CRediT roles: "
  "conceptualization, data curation, formal analysis, funding acquisition, investigation, "
  "methodology, project administration, resources, software, supervision, validation, "
  "visualization, writing \u2013 original draft, writing \u2013 review and editing."),
 ("Declaration of competing interests", "todo",
  "REQUIRED, and uploaded separately as a .doc/.docx via the declarations tool. Discloses employment, "
  "consultancies, stock, honoraria, expert testimony, patents, grants and any journal affiliation. "
  "Authors with nothing to declare select that option explicitly."),
 ("Funding", "todo",
  "REQUIRED, in the guide's standard form: \"Funding: This work was supported by the National "
  "Institutes of Health [grant numbers xxxx, yyyy]\". The role of sponsors in study design, "
  "collection, analysis, interpretation, writing and the decision to submit must be stated, or their "
  "non-involvement stated."),
 ("Declaration of generative AI and AI-assisted technologies in the manuscript preparation process", "todo",
  "REQUIRED IF APPLICABLE, as a new section immediately before the reference list, in the guide's "
  "form: \"During the preparation of this work the author(s) used [TOOL] in order to [REASON]. After "
  "using this tool/service, the author(s) reviewed and edited the content as needed and take(s) full "
  "responsibility for the content of the published article.\" The exemption for basic grammar, "
  "spelling and reference tools does not cover drafting assistance. AI tools must not be listed as "
  "authors. THE WORDING IS THE PI'S DECISION."),
 ("Data and code availability", "draft", None),
]

FRONT = dict(
 # Settled 7 September 2026. Every word traces to a sentence already in the manuscript:
 # "absent by construction" is §2.2 — a reference database cannot contain these peptides, which is
 # why the conventional search returns zero and why that zero is a passed control rather than a
 # failure; "not by chemistry" is the Conclusion verbatim, resting on §2.8's nine environment
 # features. The earlier working title claimed "ligandable", which §2.8 explicitly refuses — those
 # sites have been tested against nothing. That was critic finding 1, and this closes it.
 title="Absent by construction, not by chemistry: the lysines cancer creates",
 title_note="Settled 7 September 2026. Avoid any title containing \"drugable\" or \"ligandable\" — the paper measures neither for this class.",
 authors="[PLACEHOLDER] author list and affiliations",
 venue="Bioorganic & Medicinal Chemistry — Tetrahedron Prize special issue honouring B. F. Cravatt (invited Research Article).",
)

KEYWORDS = dict(status="draft", body=[
 "chemoproteomics", "lysine", "covalent", "ligandability", "proteomics", "variant", "AlphaFold",
], note=("Seven, the maximum the guide allows, all single words as the guide prefers. The earlier "
         "eight-keyword list used multi-word entries the guide discourages."))

ABSTRACT = dict(status="drafted", body=[
 "[NOT RENDERED. The abstract is drafted in drafted.py as P_ABSTRACT and build.py reads it "
 "from there. This entry is kept only so the section table stays complete; edit drafted.py.]",
])

INTRO = dict(status="todo", body=[
 "[TO BE DRAFTED — target ~700 words, four moves. (1) Covalent lysine chemistry is mature enough to ask the question; cite the three screens. (2) Mutation makes new lysines and nobody has counted them as a target class; E→K is a known common substitution, so cite the priors and do not claim novelty there; name the counter-position (passenger-like) and flag that Figure 1d answers it. (3) The reason nobody has is a measurement problem: a new lysine changes peptide boundaries, and variant-peptide FDR has mostly been tested incorrectly. (4) What this paper does: a reanalysis of published data, no new data, every script published.]",
])

RESULTS = [
 dict(n="2.1", status="draft", heading="Cancer mutation creates a large, driver-enriched class of new lysines",
      fig="fig1.png", figno=1,
      body=["[TO BE DRAFTED — ~800 words. Must establish: the class is large and multi-source (295,411 sites, 20,365 genes; 93,769 multi-source); E→K is the commonest route but not because lysine is special (rank 1 of 363 observed pairs; trinucleotide-adjusted only 1.12-fold; lysine 4th by residues gained, 5th by net flux, cysteine 2nd); the sites are somatic (91.12% of assessed absent or below 1e-5, with 6.73% never assessed counted separately); they are driver-enriched (OR 1.68, p = 1.3e-10) — and the Hess flag is conservative by construction because that method exists to strip passenger hotspots; identifiable processes write them (UV 23.5/9.6, MSI 20.6/6.3, POLE 15.2/6.0, APOBEC 8.3/4.3, with smoking and CpG running the other way); burden is lineage-specific with an internal control (melanoma 39 per tumour, uveal melanoma 0, all 32 lineages shown); and the excess is real against two independent germline references (6.72% against 3.29%, non-cancer control 3.41%).]"],
      cap=["[Caption to be finalised from the figure legends document.]"]),
 dict(n="2.2", status="draft", heading="The class is largely invisible to conventional proteomics",
      fig="fig2.png", figno=2,
      body=["[TO BE DRAFTED — ~900 words. Gate 1: 617 of 8,470 regions ever covered (7.28%); the 6.99%/6.79% figures are the cross-architecture agreement check over all 10,372 variant sites and must be labelled as that, not as the gate-1 result. Gate 2: 53 of 8,470 (0.63%), 27 of 617 (4.38%); conventional search returns zero, a passed control. Drivers: of the ten driver hotspots whose region was covered, none yielded a mutant peptide — and NRAS Q61K was recovered in HCC15 only by the mutant-only route, where the wild-type region was never observed. Two routes: 53 + 15 (7 shared) = 61 sites. State plainly that this is what has been looked for and found, not evidence of absence.]"],
      cap=["[Caption to be finalised from the figure legends document.]"]),
 dict(n="2.3", status="draft", heading="Reported error rates do not hold for variant peptides",
      fig="fig3.png", figno=3,
      body=["[TO BE DRAFTED — ~700 words. Burden-matched entrapment at r = 1: FDP 4.6% at margin 0 falling to 0.83% at margin ≥ 10, 1,033 variant identifications against 6 entrapment hits; pooled 0.58% across 42 plexes; decoys underestimate by 1.7-fold on plex 2. Do not claim burden-matched entrapment as a new method — cite Wen et al. and describe it as correct application. Figure S4h is not comparable with the 0.58% and the text must say so.]"],
      cap=["[Caption to be finalised from the figure legends document.]"]),
 dict(n="2.4", status="draft", heading="Reactivity and ligandability separate only where the structure is trustworthy",
      fig="fig4.png", figno=4,
      body=["[TO BE DRAFTED — ~750 words. 245,684 lysines, 6,201 proteins, 30,805 labelled (12.5%); 14 site features, 5 charge and 9 shape, with pLDDT excluded as the gating variable. The sentence the section exists to write: under the confidence gate, shape features stop explaining reactivity (crowding +0.197 → −0.042, atom density +0.244 → +0.003, polar +0.183 → −0.002 at pLDDT ≥ 90) while holding in ligandability (crowding +0.143 → +0.100, polar +0.185 → +0.197), and that divergence is what separates the two properties. Report the gate ladder (shape tracking 0.97 → 0.90 → 0.77 → 0.35; p 0.487 / 0.081 / 0.019 / 0.008) and both permutation arms (+0.008 ungated, +0.638 gated). Say that the separation emerges only at the strictest gate.]"],
      cap=["[Caption to be finalised from the figure legends document.]"]),
 dict(n="2.5", status="draft", heading="Ligandability is reproducible, and poorly predicted by the site features tested",
      fig="fig5.png", figno=5,
      body=["[TO BE DRAFTED — ~850 words. AUC 0.617 against a measured reproducibility ceiling of at least 0.72; train 0.619, optimism +0.002. Report the composite BOTH ways — 0.617 including model confidence as a feature and 0.611 [0.590, 0.632] without it — which forecloses the objection that the model is partly predicting how well a site was modelled, and name this feature set distinctly from Figure 4's. The ceiling is flat across five confidence bands. No useful operating point: 11% of engaged sites recovered at 5% false positives, 19% at 10%. The shortfall is not a feature-count problem: 14 features carry about five independent properties. Breadth across 34 chemotypes. Enantiomers as a bounded observation, n = 2 pairs, leading with the agreement rather than the difference.]"],
      cap=["[Caption to be finalised from the figure legends document.]"]),
 dict(n="2.6", status="draft", heading="Engagement is a property of the lysine, over the descriptors tested",
      fig="fig6.png", figno=6,
      body=["[TO BE DRAFTED — ~700 words. 13 of 14 site features consistent within every one of 84 compounds; 0 of 84 descriptor correlations within NatChem alone and 0 of 72 in the 173-compound within-dataset test, against 3.6 expected; the pooled 21 of 84 was provenance, with molecular weight all but naming the screen (306/386/480 Da, p = 1.3e-9). Keep n = 84 and n = 173 apart. The heading is bounded because figure S7f shows the tested descriptors contain nothing encoding charge and every shape term is achiral by construction — and Figure 7 and Figure 5f both show molecular properties that do matter. Never word this as \"no property of the molecule matters\".]"],
      cap=["[Caption to be finalised from the figure legends document.]"]),
 dict(n="2.7", status="draft", heading="The probe is not a neutral reporter: reagent charge steers labelling",
      fig="fig7.png", figno=7, body=None, cap=None),   # filled from sec_2_7.py
 dict(n="2.8", status="draft", heading="Acquired lysines are structurally ordinary",
      fig="fig8.png", figno=8, body=None, cap=None),   # filled from sec_2_8.py
]

DISCUSSION = dict(status="todo", body=[
 "[TO BE DRAFTED — ~1,200 words, five moves. (1) What the paper establishes, in one paragraph. (2) The closest prior work, named: the three screens mapped wild-type lysines; Craven et al. drugged one acquired lysine by target-based design; neither asked whether the class is measurable. Rebut the passenger reading of E→K explicitly. (3) What follows for probe design — Figure 7 is the ABPP design principle, and it stops at the principle; the five-probe panel is not named here. (4) What follows for measurement — enrichment removes the abundance gate, a margin filter is required, entrapment must be burden-matched; keep this short. (5) Limitations in one honest block: cell lines throughout, screens not saturated so every ligandability estimate is a lower bound, cross-dataset generality the weakest evidence, descriptor coverage excludes charge and chirality, environment features only in Figure 8. Do not convert cell-line engagement into disease relevance anywhere.]",
])

CONCLUSION = dict(status="todo", body=[
 "[TO BE DRAFTED — ~120 words. The class is large, driver-enriched and effectively unmeasured; the obstacle is detection; the chemistry to engage such sites already exists. Confirm whether this journal expects a separate Conclusion or folds it into the Discussion.]",
])

METHODS = dict(status="todo", body=[
 "[TO BE DRAFTED. Sections required, with what each must state: variant catalogues and the union (DepMap 26Q1, TCGA MC3, COSMIC CMC v103; and that TCGA-OV and TCGA-LAML are whole-genome-amplified, so a PASS-only filter deletes them). Germline comparison (dbNSFP 5.4a, MANE, gnomAD v4.1 joint and v2.1.1 non-cancer; expected diploid count as the sum of 2 × allele frequency; not a count of distinct variants). Proteomics reanalysis (42 TMT plexes, 375 cell lines; both search architectures; the two-gate definition). Entrapment and FDR (FDRBench v1.1.1, r = 1, burden matching, the score-margin definition, and which estimator supports which claim). Site assignment (97.82% two-route concordance, 44,736 of 45,552; why assignments are called provisional). Occupancy and engagement (100 × (1 − compound/DMSO); engaged at ≥ 60%; ratio ≥ 4 equals 75% engagement and is a reporting convention). Structural features (AlphaFold models; the Cβ anchor; TWO feature sets named distinctly — Figure 4's 14 exclude pLDDT as the gating variable, Figure 5's 14 include model confidence, and the composite is reported both ways; TWO gate rungs — Figures 4b, c display pLDDT ≥ 90, Figure 8 is gated at ≥ 70). Cell lines (Figures 2–3 are 375 CCLE lines; Figures 4–8 are the chemoproteomic screens: CyMPSA in HeLa, meroterpenoids in MDA-MB-231 and MDA-MB-436, NatChem to be confirmed). Software and code (every script published, with versions).]",
])

DATA_AVAIL = dict(status="drafted", body=[
 # Drafted 2026-09-07 alongside the code deposit. Every version string here is also carried, in
 # machine-readable form, in the deposit's external_inputs.json, which build_release.py refuses
 # to assemble if any input is missing from it. The two must be changed together.
 "**No new data were generated.** Every source is public and is listed with its version below; "
 "each is also recorded in the deposit\u2019s external_inputs.json, which the release build "
 "refuses to assemble if any input is undeclared.",

 "**Variant catalogues.** DepMap Public 26Q1 (depmap.org/portal/download), comprising "
 "OmicsSomaticMutations, Model, CRISPRGeneEffect, CRISPRInferredCommonEssentials and "
 "OmicsExpressionTPMLogp1HumanProteinCodingGenes. TCGA MC3 v0.2.8-PUBLIC "
 "(gdc.cancer.gov/about-data/publications/mc3-2017). COSMIC Cancer Mutation Census v103 "
 "(cancer.sanger.ac.uk/cmc); COSMIC requires a licence for commercial use, which this deposit "
 "does not confer.",

 "**Population frequency and reference sequence.** dbNSFP 5.4a; gnomAD v4.1 joint and v2.1.1 "
 "non-cancer exomes (gnomad.broadinstitute.org); MANE Select v1.2 (ncbi.nlm.nih.gov/refseq/MANE); "
 "UniProt and AlphaFold DB (uniprot.org, alphafold.ebi.ac.uk).",

 "**Proteomes.** The 42 TMT plexes reanalysed here are the published CCLE proteomes of Nusinow "
 "et al. [[nusinow2020]]; raw files were obtained from that study\u2019s repository and "
 "re-searched with MSFragger 4.4.1, FragPipe v24.0, Percolator and Philosopher. Entrapment sets "
 "were built with FDRBench v1.1.1.",

 "**Chemoproteomic screens.** Occupancy and compound structures were taken from the supporting "
 "information of the three published screens [[abbasov2021,bracken2024,ryan2024]]. No new "
 "chemoproteomic experiments were performed.",

 "**Code and deposited data.** All analysis and figure scripts, the panel-level data behind every "
 "figure, and the checks described below are deposited at [PLACEHOLDER repository] and archived "
 "at [PLACEHOLDER DOI]. The deposit has two tiers: every published panel can be redrawn from the "
 "deposited CSVs alone, and those CSVs can be regenerated from the public sources above. Each "
 "figure script also writes a machine-readable ledger of the numbers its panels contribute to "
 "the text; two scripts check those ledgers against the manuscript in both directions \u2014 "
 "numbers computed but no longer quoted, and numbers quoted that no script produces. A small "
 "number of intermediates could not be tied to a producing script by static analysis and are "
 "deposited as data rather than regenerated; they are labelled as such in the deposit.",
])
