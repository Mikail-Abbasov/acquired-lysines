# -*- coding: utf-8 -*-
"""Main-text content, one entry per section, assembled by build.py into a single manuscript.

Cite references as [[key]] or [[key1,key2]]; the builder numbers them globally in order of first
appearance. `fig` names the figure image (in ccle_plex05/) and `cap` its caption paragraphs.
`status` is "draft" for written sections and "todo" for placeholders carried from the outline.
"""

# ---------------------------------------------------------------- BMC required front/back matter
# Every item below is a stated requirement or encouragement in the Bioorganic & Medicinal Chemistry
# guide for authors, read 6 September 2026. Where the guide is explicit the requirement is quoted.

# ---------------------------------------------------------------- Resource availability
# Required by Cell Press for all research articles, with exactly these three subsections.
# Added 19 October 2026; the manuscript had none. It is main text and counts against the 7,000.
# The two bracketed items are the only things that cannot be written from what is on record.
RESOURCE_AVAIL = dict(status="draft", body=[
 "**Lead contact.** Requests for further information should be directed to the lead contact, "
 "Mikail E. Abbasov (mikail-abbasov@ou.edu).",
 "**Materials availability.** This study generated no new unique reagents. The lysine-directed "
 "probes and compounds reanalyzed here were reported in the three source studies and are "
 "available from their corresponding authors.",
 "**Data and code availability.** This study generated no new primary mass-spectrometry data. "
 "The proteomic data reanalyzed here are the CCLE proteomes at MassIVE MSV000085836, and the "
 "three chemoproteomic screens are available from their original publications. Variant calls are "
 "from DepMap Public 26Q1, the TCGA MC3 controlled call set and COSMIC Cancer Mutation Census "
 "v103; sequences are MANE v1.2 and GENCODE v38; population frequencies are gnomAD v4.1 and the "
 "gnomAD v2.1.1 non-cancer exome subset; structures are AlphaFold DB. All derived site tables, "
 "the mutation-aware database construction, the entrapment analysis and every figure script are "
 "deposited at github.com/Mikail-Abbasov/acquired-lysines and are sufficient to reproduce every "
 "number reported here. Any "
 "additional information required to reanalyze the data reported in this paper is available from "
 "the lead contact on request.",
])

HIGHLIGHTS = dict(status="drafted", spec=(
 "REQUIRED, and read from the Cell Chemical Biology final-submission page on 20 October 2026: "
 "\"You may include up to four highlights. Avoid the use of abbreviations. The length of each "
 "highlight must not exceed 85 characters (including spaces).\" They are uploaded in a single "
 "Word document SEPARATE from the main text, together with the eTOC blurb and the Significance "
 "paragraph, under the Editorial Manager option that includes \"highlights and eTOC blurb\". "
 "Because that document is separate from the main text, none of its three components counts "
 "against the 7,000. Character counts are appended in brackets by the build, measured rather "
 "than asserted, and are not part of the bullet."), body=[
 # Corrected 20 October 2026 against the final-submission page. Two earlier readings were wrong:
 # the limit is FOUR, not three to five, and it is mandatory rather than encouraged; and the
 # 85-character limit is Cell Press's own, not a carry-over from the BMC guide.
 # Abbreviations are now spelled out, since the page asks for them to be avoided: FDR became
 # "false discovery" and AUC became "area under the curve".
 # The probe bullet was the fifth and is the one dropped, because bullets 1 to 4 trace the
 # paper's argument end to end (the class exists, it is unseen, the error rate does not hold,
 # and prediction cannot route around it) while the probe result is a separate finding about
 # reagents. Swapping it back in for bullet 4 is a one-line change if MEA prefers it:
 #   "Probes differing in the acyl leaving group split on 2 of 22 features, both charge"
 "Cancer mutation creates 295,411 acquired lysines, driver-hotspot odds ratio 1.68",
 "Mutation-aware search sees 0.63% of 8,470 observations, 0.35% past the margin",
 "Reported 1% false discovery is 4.70% for variant peptides without a score margin",
 "Structure predicts ligandability at an area under the curve of 0.617, ceiling 0.72",
])

GRAPHICAL_ABSTRACT = dict(status="drafted", fig="graphical_abstract.png", spec=(
 "REQUIRED AT SUBMISSION: the guide states \"You are required to provide a graphical abstract\". "
 "531 \u00d7 1328 pixels (h \u00d7 w) or proportionally more, readable at 5 \u00d7 13 cm, submitted as a "
 "separate file. Preferred formats TIFF, EPS, PDF or MS Office. Drawn at exactly that aspect, "
 "13.28 \u00d7 5.31 cm at 600 dpi, by ccle_plex05/graphical_abstract_v2.py. Confirm the pixel "
 "minimum and the accepted formats with the editor, since the Elsevier author pages cannot be "
 "fetched from here.\n\n"
 "FOUR PICTURES, each carrying one number, because the earlier version was a slide: four columns "
 "of heading, subtitle, number and two description lines, which a reader outside the field meets "
 "as a wall of type. (1) a protein chain with one residue rewritten to lysine, Figure 1a; (2) one "
 "hundred acquired lysines as dots with 82 of them drawn inside the region the screens already "
 "engage, Figures 4a, 8b and 8c; (3) the same chain digested twice, a new cysteine leaving one "
 "peptide and a new lysine splitting it into two, Figure 2a; (4) 8,470 cells, one per "
 "acquired-lysine observation, with the 53 recovered picked out, Figures 2b and 2c.\n\n"
 "PANEL 4 IS THE ARGUMENT FOR DRAWING THIS IN THE PIPELINE. The grid is exactly 110 \u00d7 77 = "
 "8,470 cells and exactly 53 are filled, both read from the ledgers, and the rendered file is "
 "re-counted by connected components to confirm 53 blocks survive to the artwork. A generative "
 "tool draws an approximate field, and an approximate field is a fabricated one.\n\n"
 "BOUNDS HELD IN THE ARTWORK. Panel 2 says the sites occupy the same structural space as engaged "
 "lysines, not that they have been engaged, and the footer says so. Panel 3 is a mechanism drawn "
 "from the four peptide species of Figure 2a, not a measurement. Nothing asserts that the barrier "
 "is detection rather than chemistry, which is still open for approval.\n\n"
 "Every number is read from results/panel_data/*_head.json, so the artwork cannot drift from the "
 "panels. The conventional search's zero, which no ledger held, is computed and deposited by the "
 "same script. The specialist results, the composite AUC against its ceiling and the two of "
 "twenty-two probe features, sit in the footer rather than in a panel."))

DECLARATIONS = [
 # Acknowledgements and the generative-AI declaration were removed on the PI's instruction,
 # 19 September 2026. The guide records the AI section as required if applicable, with the
 # exemption for grammar and reference tools not covering drafting assistance.
 ("Author contributions", "drafted", [
  "**CRediT authorship contribution statement.** "
  "**Rajib Molla:** Data curation, Investigation, Writing \u2013 review & editing. "
  "**Sona Tiwari:** Data curation, Investigation, Writing \u2013 review & editing. "
  "**Ayan Mukherjee:** Visualization, Writing \u2013 review & editing. "
  "**Md Khaled Hossain:** Investigation, Writing \u2013 review & editing. "
  "**Mikail E. Abbasov:** Conceptualization, Methodology, Software, Formal analysis, Validation, "
  "Visualization, Writing \u2013 original draft, Supervision, "
  "Project administration, Funding acquisition.",
  # The attribution rationale was printed under the CRediT statement. It is addressed to the
  # co-authors, to be confirmed before submission, and it is not a statement for a reader.
  # Kept here so the basis of each role stays on the record.
  # Rajib Molla assembled the acquired-lysine site catalogue and its reported mechanisms;
  # Sona Tiwari supplied the primary literature and the patient-level frequencies behind
  # table S1; Ayan Mukherjee is attributed Visualization for the chemical structures in
  # figures 6a and 7e, and the role should be dropped if they do not appear;
  # Md Khaled Hossain located references and reviewed the draft.
  ]),
 ("Declaration of competing interests", "drafted", ["None."]),
 ("Funding", "drafted", [
  "This work was supported by the National Institute of General Medical Sciences of the National "
  "Institutes of Health under award R35GM143027, and by a University of Oklahoma Health Sciences "
  "Center College of Pharmacy Seed Grant, PHR20102."]),
 ("Data and code availability", "draft", None),
]


FRONT = dict(
 # Reset twice on 18 September 2026. The 7 September title, "Absent by construction, not by
 # chemistry: the lysines cancer creates", went for two reasons: the thesis-and-colon shape read
 # as a Perspective, and "not by chemistry" was the Conclusion sentence §2.8 forbids.
 # The first replacement stated that reference-database proteomics cannot detect the class. That
 # is true and trivial: a reference database cannot contain a peptide whose boundaries mutation
 # redrew, so the conventional zero is a passed control and needs no experiment. Putting the
 # control on the cover is what made it read as obvious.
 # The title now carries the result one gate further in, which no reader can call obvious, and
 # states it as the positive count rather than as a negative. Three earlier attempts failed in
 # three different ways: one was trivially true, one called a plex a proteome, and one said the
 # search "still misses" the class, which is false because 53 detections were made. A title
 # that prints both numbers lets the reader do the division and cannot overshoot.
 # A mutation-aware database of 10.17 million peptide entries, searched over the 42 TMT plexes
 # that carry the 375 CCLE cell-line proteomes (nine lines per plex plus a bridge channel),
 # still recovers the mutant peptide in only 0.63% of 8,470 site-in-
 # cell-line observations, on 52 distinct genomic sites, with 30 of the 53 clearing the score
 # margin (§2.2). The region is covered in 7.28% of cases, and where it is covered the mutant
 # is missed 95.6% of the time, so the barrier survives the fix everyone would name first.
 # The title gives the count and its denominator rather than the percentage: 0.63% is a rate per
 # site-in-cell-line observation, and alone on a cover it would be read against §2.1's 295,411
 # distinct sites, which is the denominator trap the paper polices everywhere else. Printing
 # 53 and 8,470 together removes the ambiguity and names the unit.
 # The earlier working title claimed "ligandable", which §2.8 explicitly refuses; those sites
 # have been tested against nothing. That was critic finding 1 and it stays closed.
  # Replaced 19 October 2026 with the title agreed in the 13 October front-matter pass and
 # confirmed by MEA. 12 words, 96 characters, inside the Cell Press guidance of 10 to 12 words
 # and 145 characters. The 21 September title named the result as a bare count; this one names
 # what the paper measures, which is where the coverage and error-rate limits are.
 title="Proteome-wide reanalysis identifies coverage and error-rate limits on acquired lysines in cancer",
 running_title="Coverage and error-rate limits on acquired lysines",   # 50 characters
 title_old_2026_09_21="A mutation-aware search of 375 cancer cell-line proteomes recovers 53 of 8,470 acquired-lysine observations"
       "acquired-lysine observations",
 title_note="Reset 18 September 2026, twice. Avoid any title containing \"drugable\" or \"ligandable\"; the paper measures neither for this class. Avoid any title asserting that the obstacle is not the chemistry, for the same reason §2.8 refuses it in the text.",
 authors="Rajib Molla\u00b9, Sona Tiwari\u00b9, Ayan Mukherjee\u00b9, Md Khaled Hossain\u00b9, Mikail E. Abbasov\u00b9,\u00b2,\u00b3,\u2074",
 affiliations=("\u00b9Department of Pharmaceutical Sciences, College of Pharmacy, University of "
               "Oklahoma Health Sciences Center, Oklahoma City, OK 73117. \u00b2Department of "
               "Biochemistry and Physiology, College of Medicine, University of Oklahoma Health "
               "Sciences Center, Oklahoma City, OK 73104. \u00b3Stephenson Cancer Center, "
               "University of Oklahoma Health Sciences Center, Oklahoma City, OK 73104. "
               "\u2074Harold Hamm Diabetes Center, University of Oklahoma Health Sciences Center, "
               "Oklahoma City, OK 73104."),
 venue="Cell Chemical Biology (Resource). Retargeted 9 October 2026 from the Bioorganic & Medicinal Chemistry Tetrahedron Prize special issue honoring Benjamin F. Cravatt; that commitment still needs withdrawing.",
)

KEYWORDS = dict(status="drafted", body=[
 # Set by MEA, 9 November 2026. Nine of the ten Cell Press allows, comma-separated.
 "acquired lysines", "cancer mutation", "mass-spectrometry-based chemoproteomics",
 "activity-based protein profiling", "covalent modification", "ligandability",
 "false discovery rate", "variant peptide", "AlphaFold",
], note="Nine of the ten Cell Press allows, comma-separated.")

ABSTRACT = dict(status="drafted", body=[
 "[NOT RENDERED. The abstract is drafted in drafted.py as P_ABSTRACT and build.py reads it "
 "from there. This entry is kept only so the section table stays complete; edit drafted.py.]",
])

INTRO = dict(status="drafted", body=[
 "[superseded by spine.py, which carries this section's prose; see git history]",
])

RESULTS = [
 dict(n="2.1", status="drafted", heading="Cancer mutation creates a large, driver-enriched class of new lysines",
      fig="fig1.png", figno=1,
      body=["[superseded by spine.py, which carries this section's prose; see git history]"],
      cap=["[Caption to be finalised from the figure legends document.]"]),
 dict(n="2.2", status="drafted", heading="The class is mostly missed by a search built to find it",
      fig="fig2.png", figno=2,
      body=["[superseded by spine.py, which carries this section's prose; see git history]"],
      cap=["[Caption to be finalised from the figure legends document.]"]),
 dict(n="2.3", status="drafted", heading="Reported error rates do not hold for variant peptides",
      fig="fig3.png", figno=3,
      body=["[superseded by spine.py, which carries this section's prose; see git history]"],
      cap=["[Caption to be finalised from the figure legends document.]"]),
 dict(n="2.4", status="drafted", heading="Reactivity and ligandability separate only where the structure is trustworthy",
      fig="fig4.png", figno=4,
      body=["[superseded by spine.py, which carries this section's prose; see git history]"],
      cap=["[Caption to be finalised from the figure legends document.]"]),
 dict(n="2.5", status="drafted", heading="Ligandability is reproducible, and poorly predicted by the site features tested",
      fig="fig5.png", figno=5,
      body=["[superseded by spine.py, which carries this section's prose; see git history]"],
      cap=["[Caption to be finalised from the figure legends document.]"]),
 dict(n="2.6", status="drafted", heading="Engaged lysines share site properties across compounds, and no tested descriptor predicts which",
      fig="fig6.png", figno=6,
      body=["[superseded by spine.py, which carries this section's prose; see git history]"],
      cap=["[Caption to be finalised from the figure legends document.]"]),
 dict(n="2.7", status="drafted", heading="The probe is not a neutral reporter: reagent charge steers labelling",
      fig="fig7.png", figno=7, body=None, cap=None),   # filled from sec_2_7.py
 dict(n="2.8", status="drafted", heading="Acquired lysines are structurally ordinary",
      fig="fig8.png", figno=8, body=None, cap=None),   # filled from sec_2_8.py
]

DISCUSSION = dict(status="drafted", body=[
 "[superseded by spine.py, which carries this section's prose; see git history]",
])

CONCLUSION = dict(status="drafted", body=[
 "[superseded by spine.py, which carries this section's prose; see git history]",
])

METHODS = dict(status="drafted", body=[
 "[superseded by spine.py, which carries this section's prose; see git history]",
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

 # Deposit shape settled 19 September 2026 and narrowed 9 November: one public repository, no
 # archived release. The 343 MB of deposited intermediates and the 120 MB of search databases
 # are redistributed nowhere, which is stated here rather than implied.
 "**Code and deposited data.** All analysis and figure scripts, the panel-level data behind every "
 "figure, and the checks described below are available at "
 "github.com/Mikail-Abbasov/acquired-lysines. The "
 "code carries the MIT licence and the deposited data CC BY 4.0. The deposit has two tiers. "
 "Every published panel can be redrawn from the deposited panel data alone, with no download. "
 "Those files can in turn be regenerated from the public sources above, with two exceptions that "
 "are stated rather than implied: twenty intermediates that are this study\u2019s own outputs "
 "rather than public downloads, and the variant and decoy search databases. Neither is "
 "redistributed, both are named in external_inputs.json, and the second tier is therefore "
 "incomplete for the steps that read them while the first tier covers every figure affected. "
 "Each figure script also writes a machine-readable ledger of the numbers its panels contribute "
 "to the text; two scripts check those ledgers against the manuscript in both directions, "
 "numbers computed but no longer quoted, and numbers quoted that no script produces.",
])

# ================================================================ the separate upload
# Read from the Cell Chemical Biology final-submission page, 20 October 2026. Highlights, the
# eTOC blurb and the Significance paragraph "should appear together in a single Word document,
# separate from the main manuscript text", uploaded under the Editorial Manager option that
# includes "highlights and eTOC blurb". Because that document is separate from the main text,
# none of it counts against the 7,000. This is a journal-mandated third file, which is the one
# exception to the one-manuscript-one-SI rule.
ETOC = dict(status="drafted", spec=(
 "REQUIRED. 50 words or fewer, two sentences, third person, no \"we\" or \"I\". Should showcase "
 "the context and significance of the findings and be readable by a non-specialist. Cell Press "
 "house style names the authors; the name goes in once the author list is final."), body=(
 "Covalent drugs need a residue that exists only in diseased cells, and cancer mutation writes "
 "lysines by the hundred thousand. This reanalysis counts the class, finds that conventional "
 "proteomics recovers under 1% of it, and measures the error rate on those detections, "
 "placing the obstacle in detection rather than chemistry."))

SIGNIFICANCE = dict(status="drafted", spec=(
 "REQUIRED. A single paragraph of no more than 300 words, explaining the relevance and "
 "significance of the study to the broad readership. No reference citations. Uploaded in the "
 "same separate document as the Highlights and the eTOC blurb, so it does not count against the "
 "7,000."), body=(
 "A covalent drug is written around the residue it modifies, so the most valuable residue is one "
 "that exists in the diseased cell and nowhere else. Cancer supplies such residues by mutation, "
 "and the argument has already been made for cysteine, where a mutation-created thiol became a "
 "drug target and then a drugged one. Lysine has not been counted the same way, although "
 "mutation creates more of them. This work counts the class and then asks the question that "
 "decides whether it can be worked on: has anyone been able to see it? Re-searching 42 "
 "published proteomes with a mutation-aware database recovers the mutant peptide at well under "
 "1% of the opportunities, and more than nine in ten are never observed in any "
 "form. The reason is specific rather than general. A new lysine is a new cleavage site, so the "
 "peptide a search must find is not the peptide it would have looked for, and the error control "
 "that would license the few identifications it does make has been calibrated over whole "
 "searches rather than within this atypical subset of one. Measured directly, with entrapment "
 "sequences built to be indistinguishable from their targets by mass and length, the realized "
 "error rate at a reported 1% filter is nearly fivefold higher, and a score margin restores it. "
 "Set against what three chemoproteomic screens show about which lysines compounds engage, "
 "nothing in the chemistry of these sites explains their absence from the record: they sit "
 "within the ordinary range of the lysines this chemistry already reaches. The barrier is "
 "detection, which is the tractable kind of problem. For anyone designing lysine-directed "
 "chemistry, the practical consequence is that this target class has to be measured before it "
 "can be prioritized, and the measurements that would do it are named here."))
