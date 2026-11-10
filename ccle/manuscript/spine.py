# -*- coding: utf-8 -*-
"""Condensed main text for the Cell Chemical Biology Resource.

WHY THIS FILE EXISTS
Cell Chem Biol Resource: under 7,000 words counting main text and main figure legends,
excluding STAR Methods, supplemental legends and references; at most 7 main figures.
The draft measured 16,950 countable words over 8 figures, so 59% has to come out.

The decision (19 Oct 2026) was to keep every result and move prose to the Supporting
Information rather than to narrow the paper. So:
  - Figure 8 leaves the main text and becomes figure S9.
  - Each main section keeps a condensed spine here; the full original prose and the full
    original legend move to the SI verbatim as an Extended Results note, so nothing is lost
    and the cut is auditable against drafted.py.

SPINE[n]    replaces drafted.P_2_<n>   for the main text
LEGEND[fig] replaces drafted.CAP_2_<n> for the main text
A section absent from these dicts falls through to its full original text, so the build
always runs and the word gate reports honestly how far over the limit it still is.

House rules observed here: US English, numerals and symbols, claim-first
paragraphs, and every bound and statement of evidentiary status kept. The detail that moves
to the SI is panel arithmetic and method, never a caveat.
"""

SPINE = {}
LEGEND = {}

# EXTENDED[section] holds paragraphs that leave the COUNTED main text and are reproduced in
# Extended Results instead. Nothing is lost: Cell Press excludes supplemental item legends and
# the supplement itself from the 7,000, so an audit that a main-text reader does not need in
# order to follow the argument costs nothing there and is still fully verifiable.
# The rule applied: the main text carries the figure's spine and every bound on it. A secondary
# audit, a control whose evidence is entirely in the supplement, and a bookkeeping note about
# denominators are moved. A bound on a spine claim never moves.
EXTENDED = {}

# ---------------------------------------------------------------- 2.1 / Figure 1
# Split used throughout: the BODY carries the argument and every bound, the LEGEND carries the
# panel arithmetic, the SI carries the method. No number appears in both body and legend, which
# is where most of the reduction comes from. Bounds are never what moves.
#
# One paragraph per panel and one legend entry per panel, so panels.py can move a panel out
# without having to cut a sentence in half.
SPINE["2.1"] = [
# opener, no panel
"Cancer makes a question available that a proteome-wide survey cannot: how many lysines does the "
"disease itself create? Whether that is a class worth building chemistry around depends first on "
"how large it is.",
# panel A
"The class is large, and it is not the product of one dataset. Pooling missense sites from DepMap "
"[[tsherniak2017]], the TCGA MC3 call set [[ellrott2018]] and COSMIC [[tate2019]] gives 295,411 "
"sites that create a lysine, across 20,365 genes (Figure 1A). The three sets are nested rather than "
"independent, so the union is a catalogue of what has been reported and not three replicate surveys.",
# panel B
"Glutamate to lysine is the commonest substitution in cancer, and that is a statement about "
"glutamate and about C>T mutation rather than about lysine. E to K ranks first of the 363 "
"substitution pairs observed, but against a trinucleotide-context expectation its enrichment is only "
"1.12-fold, against more than three for each of the CpG hotspots R to H, R to C and R to Q "
"(Figure 1B). Lysine is fourth of twenty residues by sites gained. What distinguishes it is the "
"asymmetry of the routes: lysine is gained overwhelmingly by C>T and never lost by C>T, so a "
"C>T-dominated process creates it and does not destroy it. The class is abundant for a mechanistic "
"and unremarkable reason.",
# panel C
"The sites are somatic rather than common polymorphism. More than nine in ten assessed sites are "
"absent from gnomAD v4.1 [[karczewski2020]] or below an allele frequency of 10\u207b\u2075 "
"(Figure 1C), and the sites never assessed are counted separately rather than folded into a "
"denominator that would flatter the result. Frequency cannot exclude a private germline variant; "
"only matched normals can, and these data have none.",
# panel D
"Acquired lysines fall on driver hotspots more often than other substitutions do, at an odds ratio "
"of 1.68 (Figure 1D). The flag is the one carried in the DepMap variant table rather than a gene "
"list chosen here, and it is conservative by construction, since the method exists to strip "
"recurrent passenger sites from hotspot lists [[hess2019]]. Hotspot recurrence is a statement about "
"where mutation is selected, not about whether the resulting protein is expressed or reachable.",
# panel E
"The excess tracks mutational processes that make lysine, not mutation in general. Split by assigned "
"signature, the four processes that generate the relevant base changes are all enriched among "
"acquired-lysine driver hotspots, and the two that do not make lysine run the other way; both "
"directions are labeled (Figure 1E).",
# panel F
"Burden is lineage-specific, and the ranking carries its own control: cutaneous melanoma leads all "
"32 lineages, and uveal melanoma carries zero (Figure 1F). Both are melanocytic and only the first "
"is ultraviolet-driven, so a property of the lineage rather than of the mutational process would "
"make both high.",
# panel G
"The excess survives a germline reference chosen to remove the obvious objection. Acquired lysines "
"are twice as large a share of missense in the somatic union as in the germline (Figure 1G). Because "
"gnomAD contains cancer-affected individuals, the same pass was run over its non-cancer exome "
"subset, and the two agree to within 0.12 percentage points. They are not independent cohorts, both "
"being gnomAD with heavy overlap, so what the second controls is specifically the presence of "
"cancer-affected individuals and not reference bias in general. Repeating the pass for every residue "
"on the identical file, filter and frequency field ranks lysine third of twenty and cysteine fourth, "
"so the comparison the field draws with acquired cysteine can be made on one calculation instead of "
"two, and on that calculation lysine is not the weaker case. The ratio is not a measure of "
"selection: germline shares are depressed by purifying selection, which is the likeliest reading of "
"phenylalanine and tyrosine, the only two residues ranking above lysine.",
# closer, no panel
"Together these exclude the cheapest alternative to the premise of this paper, that acquired lysines "
"are simply what mutation does everywhere. They are created at twice the germline rate, they are "
"enriched on driver hotspots, and the processes that make them are identifiable. None of that says "
"the resulting proteins are expressed, detectable or reachable by a molecule, and the rest of the "
"paper is about how far each of those holds.",
]

LEGEND[1] = [
# header, no panel
"Sites are missense positions in the union of DepMap, TCGA and COSMIC; a site is an acquired lysine "
"if the substituted residue is lysine. The set and frequency panels are union-wide, the "
"driver-hotspot panel is DepMap only because the flag is carried in that table, and the germline "
"panel compares the union against two germline references. Counts are of sites, not of tumors or of "
"individuals, except where the unit is stated.",

"(A) The class, assembled from three variant sources. 295,411 acquired-lysine sites across 20,365 "
"genes, from 4,395,403 missense sites; 93,769 (31.7%) supported by more than one source, which is "
"the shaded territory. The three small overlap regions are unlabeled and are listed in Table S1.",

"(B) Codon routes into and out of lysine. E to K ranks first of 363 observed substitution pairs at "
"181,252 sites, 4.12% of all missense; context-corrected enrichment 1.12-fold, against 3.59, 3.38 "
"and 3.16 for R to H, R to C and R to Q. Lysine is gained 85.5% by C>T and lost 0% by C>T. Rust "
"marks the C>T routes, which are G:C to A:T changes, half of all cancer missense mutation. No "
"C>T route leaves lysine because AAA carries no G or C, and the only G to A available to AAG "
"gives AAA, which is still lysine. Arrow width increases with the share shown beside it.",

"(C) gnomAD v4.1 allele-frequency bands over all 295,411 sites: absent 154,188; below "
"10\u207b\u2075 96,892; 10\u207b\u2075 to 10\u207b\u2074 19,194; 10\u207b\u2074 to "
"10\u207b\u00b3 3,297; above 10\u207b\u00b3 1,971; not assessed 19,869. 91.12% of assessed "
"sites are absent or below 10\u207b\u2075, and the 6.73% never assessed is excluded from that "
"denominator.",

"(D) 1.024% of 19,239 acquired lysines carry the Hess driver-hotspot flag against 0.6138% of "
"311,166 other substitutions; odds ratio 1.68, Fisher exact p = 1.3 \u00d7 10\u207b\u00b9\u2070, "
"two-sided.",

"(E) Share of driver hotspots by assigned signature, acquired lysines (n = 315) against other "
"substitutions (n = 3,465): ultraviolet 23.5% against 9.6%, microsatellite instability 20.6% "
"against 6.3%, POLE 15.2% against 6.0%, APOBEC 8.3% against 4.3%, smoking 18.7% against 26.1%, CpG "
"deamination 10.5% against 28.2%.",

"(F) Median acquired lysines per tumor across 32 TCGA lineages. Cutaneous melanoma 39, thirteen "
"times the pan-TCGA median; uveal melanoma 0. All 32 are shown because a truncated ranking cannot "
"carry that control.",

"(G) Somatic union 6.72% (295,411 of 4,395,403) against germline 3.29%, 2.04-fold. The germline "
"value is an expected diploid count per individual, 471.1 acquired lysines among 14,342 missense "
"variants, and not a count of distinct variants. The gnomAD v2.1.1 non-cancer exome subset gives "
"508.0 per individual, 3.41%. Across all twenty residues lysine ranks third at 2.05 and cysteine "
"fourth at 1.76, in a distribution from 0.43 to 2.63, median 1.07.",

]


# ---------------------------------------------------------------- 2.2 / Figure 2
# Panels after the cut: a gate 1, b gate 2, c driver hotspots, d the second route.
# Kept verbatim in substance: the two thresholds and that neither is quoted with the other's
# counts; the conventional zero as a passed control; the NRAS Q61K caveat, without which panel C
# reads as contradicting the census; "identification quality rather than existence"; that the
# non-lysine control demonstrates the mechanism rather than controlling it away; that the two
# routes are not pooled; and the closing gate-1-against-gate-2 bound.
SPINE["2.2"] = [
"Whether anyone has been looking at these lysines is answerable on the 42 TMT plexes of the CCLE "
"proteomes [[nusinow2020]], which cover 375 cell lines and are among the deepest shotgun datasets "
"available. The unit is the site-in-cell-line observation rather than the distinct genomic site, "
"because each line is an independent chance to detect.",

"Gate 1: the region is almost never covered. A conventional search and a mutation-aware search "
"over 10.17 million peptide entries observe the wild-type region at the same rate, agreeing to "
"within 0.3 percentage points, so coverage is a property of the sample rather than of the "
"database (Figure 2B). More than nine in ten of these sites are never covered in any of the 42 "
"plexes, and 26 of them are detected only through the mutant peptide itself (Figure 2C).",

"Gate 2: where the region is covered, the mutant peptide is still usually missed. Fifty-three "
"detections across the two routes are 0.63% of the class, and 30 of them also clear the score "
"margin at which section 2.3 measures the realized error rate, 0.35% (Figure 2C). The two "
"thresholds are not interchangeable and neither should be quoted with the other's counts: the "
"larger is the census of what the searches returned, the smaller the subset over which the "
"false-discovery proportion is 0.83% rather than 4.70% (section 2.3). The acquired-lysine class "
"contributes zero entrapment peptides against 38 variant peptides at that margin, so its own "
"proportion is bounded at 7.89% by the rule of three rather than measured at zero. "
"Under the conventional "
"search the mutant "
"count is zero, which is a passed control rather than a failure, since that database contains no "
"variant sequences and can speak to coverage and to nothing else.",

"The misses include the driver hotspots. Every driver hotspot among the covered sites falls on "
"the missed side, across five genes (Figure 2D). Driver status is the DepMap hotspot flag rather "
"than a gene list chosen here [[hess2019]], and allele fraction is on the axis so that the misses "
"can be seen not to be confined to subclonal sites. The panel counts hotspots among the covered "
"sites, which is the population it plots: NRAS Q61K was recovered in one line by the mutant-only "
"route, where the wild-type region was never observed, so that observation is not in this "
"population.",

"The mutant peptide has the opposite geometry, and the search absorbs it rather than being "
"stopped by it. A new lysine is a new cleavage site, so the fully cleaved mutant peptide runs "
"about half the length of the wild-type peptide over the same site, and half the class falls "
"below the seven-residue minimum that this search and the original one both enforced (figure "
"S3f). That is not a second hidden exclusion. Both searches allowed two missed cleavages, and an "
"upstream missed cleavage returns a too-short peptide to the window, so 95% of observations still "
"carry a searchable species ending at the acquired lysine. What the geometry costs is "
"identification quality rather than existence, and whether it costs detections cannot be settled "
"here.",

"The non-lysine control demonstrates the mechanism rather than controlling it away: variant sites "
"that do not create a lysine are covered less and recovered more, the direction expected if "
"peptide boundaries are what matters, and the median lengths differ by one residue so it is not a "
"length artefact (figure S2).",

"A second, independent route reaches the same class of site. The tryptic peptide released "
"downstream of the new cleavage site carries no mutated residue but exists only if the cleavage "
"happened (Figure 2F). The two routes are reported separately rather than pooled because their "
"error properties differ, which is the subject of Figure 3.",

"What the evidence behind these counts will bear is limited, and Figure 3 exists because of it: "
"four in five of the mutant-peptide sites rest on a single spectrum, and none on more than one "
"experiment (figure S3e). Three alternative explanations for the gate-1 failure were tested and "
"none accounts for it, and what predicts identification is incomplete digestion rather than "
"charge (figure S3a to S3d). Gate 1 is a statement about tryptic geometry and gate 2 about search "
"sensitivity, and neither licenses a claim about what the proteome contains.",
]

LEGEND[2] = [
"42 TMT plexes across 375 cell lines. The unit is a site-in-cell-line observation, not a distinct "
"genomic site: 8,470 acquired-lysine observations covering 8,055 distinct sites across 372 lines, "
"with 1,902 non-lysine observations carried through the identical pipeline as the control (figure "
"S2). Every rate is per observation. Detection is counted at two gates that fail for different "
"reasons and are never combined: gate 1 is whether the wild-type region was observed at all, "
"gate 2 whether the mutant peptide was recovered where it was.",

"(B) Both bars are the same 8,470 observations. A conventional 20,431-protein search covers 640 "
"(7.56%); the mutation-aware peptide-level search over 10.17 million entries covers 617 (7.28%). "
"7,853 sites, 92.7%, are never covered in any of the 42 plexes.",

"(C) The mutant peptide was detected at 27 of the 617 covered sites (4.38%) and at 26 of the "
"7,853 uncovered ones (0.33%) by the mutant-only route, 53 sites in all, 0.63% of the class, on "
"52 distinct genomic sites, at the search's reported 1% false-discovery rate, which is where the "
"panel is drawn. Thirty of those 53 also clear the score margin of Figure 3B: 16 of 617 covered "
"and 14 of 7,853 uncovered, 0.35% of the class on 29 distinct genomic sites. Detection is by a "
"cleaved species, a missed-cleavage species, or both.",

"(D) The 617 covered sites arrayed by allele fraction and split into the 590 whose mutant peptide "
"was missed and the 27 where it was found. All ten driver hotspots in that set fall on the missed "
"row, eight gene and mutation labels across five genes, listed in Extended Results. "
"mutant-only route, where the wild-type region was never observed, so that observation is one of "
"the 26 in figure S3a and not in this population.",

"(F) The downstream peptide gives 15 sites, 7 shared with the first route and 8 new, 61 in the "
"union. It returned zero entrapment hits across 42 experiments (figure S4d) against 0.58% for "
"variant-containing peptides in the same searches, which is why the two are kept apart. Of the 53 "
"mutant-peptide sites, 41 rest on a single spectrum and 12 on two or more; of the 30 past the "
"margin, 24 and six (figure S3e).",

"*This figure states what has been looked for and found, in these cell lines, with this "
"workflow. It is not evidence that the peptides are absent from the sample.*",
]


# ---------------------------------------------------------------- 2.3 / Figure 3
# Original panels: a design, b margin dose-response, c per-plex consistency (cut), d spectrum
# quality, e precursor mass error (cut), f paired score separation. Written with the ORIGINAL
# letters; the build cuts C and E and reletters a,b,d,f to a,b,c,d.
# Every bound kept: that target-decoy is calibrated over a whole search; that neither prior
# platform reports a realized proportion within the variant class and that what follows is a
# property of THIS search; "upper-bound test rather than a decoy in a different costume";
# "counted directly rather than estimated"; that borrowing the whole-class figure is not the
# conservative choice; "bounds ... by the rule of three rather than measuring it at zero"; the two
# denominators; and that a single spectrum past the margin remains a single spectrum.
SPINE["2.3"] = [
"Figure 2 recovered a small number of mutant peptides, and a result carried by single-spectrum "
"identifications has to survive an error-rate audit before anything is built on it. The usual "
"answer, that the search was filtered at 1% false discovery, does not settle it: target-decoy "
"control is calibrated over a whole search [[elias2007]], and variant peptides are a tiny, "
"atypical subset of one, so a rate that holds globally need not hold within them.",

"Chemoproteomic work on mutation-created residues has answered this in two ways. DrugMap "
"searches canonical and cell-line-specific mutant sequences in one pass, takes consensus across "
"four search engines, and filters to 1% peptide false discovery at run level against reversed "
"canonical decoys [[takahashi2024]]. The chemoproteogenomic platform searches in two stages, "
"removing spectra confidently assigned to canonical sequences before the variant database is "
"searched [[desai2024]]. Both control the error. Neither reports a realized false-discovery "
"proportion measured within the variant class, which is what follows, and what follows is a "
"property of this search and this database, not of theirs.",

"The test is burden-matched entrapment, an application of a published method rather than a new "
"one [[wen2025]]. Each variant peptide is given a partner built from its own residues, "
"permuted, with the C-terminal lysine or arginine kept in place (Figure 3A). The pair therefore "
"shares an exact precursor mass and an identical length, and only fragment ions can separate "
"them, which is what makes this an upper-bound test rather than a decoy in a different costume. "
"An entrapment sequence cannot be in the sample, so every entrapment identification is a false "
"one, counted directly rather than estimated.",

"At the reported filter the realized error rate is roughly fivefold higher than the nominal "
"one, and a score margin recovers the intended rate (Figure 3B). The margin is the gap between "
"a spectrum's best interpretation and its next best, and requiring it to be large removes "
"exactly the identifications that a mass-matched impostor could have produced. Every error rate "
"quoted in this paper, and every identification quoted from the entrapment analysis, is taken "
"at a margin of 10 or more; Figure 2's detection census is reported at the search's own 1% "
"filter and is given at both thresholds where it appears.",

"The margin does work rather than simply discarding the weakest identifications of every kind. "
"Entrapment peptides are identified, and they score like a population that cannot be present "
"should (Figure 3F), and it is that separation the margin exploits. Neither rate is a per-site "
"confidence, and a single spectrum past the margin remains a single spectrum, which is why "
"Figure 2's detections are treated as a lower bound on coverage rather than a catalogue of "
"sites.",
]

# Moved out of the COUNTED main text and reproduced in Extended Results: the class-restricted
# entrapment audit, the two-denominator bookkeeping note, and the spectrum-quality audit.
# Each is an audit of the audit rather than the figure's spine, and each stays fully
# verifiable where it now sits. The spine and every bound on it remain above.
EXTENDED["2.3"] = [
"That rate is measurable on the acquired-lysine class alone, and measured there it is not "
"lower. A rate estimated on every variant peptide and then applied to one kind of substitution "
"is an inference across populations, and it does not have to hold, since a substitution to "
"lysine creates a tryptic cleavage site and the peptide the search must find is not the peptide "
"it would have looked for. The restriction is legitimate here because the entrapment design is "
"pairwise, so keeping only lysine-creating substitutions preserves the one-to-one match instead "
"of slicing an unmatched comparison (note N18). Restricted that way the class runs higher than "
"the rest of the variant class, and although the counts cannot resolve the difference, the "
"direction is the one that matters, because it means borrowing the whole-class figure is not "
"the conservative choice. Past the margin the class contributes no entrapment peptides at all, "
"which bounds its false-discovery proportion by the rule of three rather than measuring it at "
"zero. Seven entrapment peptides across 42 plexes is the entire supply of counted errors, so no "
"single class of substitution can be resolved to better than a few percent on its own.",

"Two denominators are in use and they are not interchangeable. Figure 3B counts distinct "
"peptides, one best spectrum each, which is the unit the entrapment pairing is defined over; "
"Extended Results counts identifications. Both are correct, and neither number should be quoted "
"with the other's counts. Decoys do not see any of this, which is why the margin is applied "
"rather than the reported q-value trusted (Extended Results).",

"The surviving identifications are ordinary spectra, not marginal ones. Judged on properties of "
"the spectrum alone, with no reference to the search engine or the database, the "
"acquired-lysine hit spectra are at the background peak count and a little above its "
"signal-to-noise (Figure 3D). The one feature that invites suspicion is a tail at plus one "
"dalton, the signature of the first isotope being selected instead of the monoisotopic peak, "
"and it occurs at the same rate in variant and in background spectra, so it is instrument "
"behaviour and not a property of variant peptides.",
]

LEGEND[3] = [
"Burden-matched entrapment at r = 1 (Methods). Every entrapment identification is by construction "
"false. Two denominators appear and are not interchangeable: panel B counts distinct peptides, "
"one best spectrum each; Extended Results counts identifications.",

"(A) The design. Entrapment sequences are matched to the variant set on the properties that drive "
"identification, so the comparison is like for like: 589,588 variant-pair rows across eight "
"plexes, entrapment sets of 147,397 each, all of them sharing residues with their target "
"(figure S4a to S4c).",

"(B) Margin dose-response over 1,489 rows, counted as distinct peptides. The false-discovery "
"proportion is 4.70% at margin 0, 1,419 variant peptides against 70 entrapment, and falls to "
"3.91% at a margin of 2, 2.36% at 5 and 0.83% at 10, where 832 variant peptides survive against "
"seven entrapment. Past a margin of 15 only one entrapment peptide remains, so the rate there "
"reflects a single count rather than a cleaner set. A 1% reported false-discovery rate is not a "
"1% realized one for this class of peptide unless a margin is applied. Restricted to the "
"acquired-lysine class the proportion is bounded at 7.9% by the rule of three, with the counts "
"and the comparison against the rest of the variant class in Extended Results.",

"(D) Median spectral-quality features for the 38 acquired-lysine identifications expressed "
"relative to 802,067 background spectra from 20 of the 77 runs in the same searches, so 1.0 "
"is the background. Peak count 0.91-fold, signal-to-noise 1.19-fold, from medians of 631 "
"against 695 peaks and 0.221 against 0.185. These features score a spectrum on its own properties, "
"with no reference to the search engine, the database or the false-discovery rate. The plus one "
"dalton tail accounts for 7.7% of variant and 7.3% of background spectra. SQS$_{GM}$ is omitted "
"because the MGF export zeroed its intensity term, which would collapse a geometric mean.",

"(F) Score-margin distributions for variant and entrapment identifications, each drawn as a "
"percentage of its own total because 1,419 cannot be compared with 70 on counts. Entrapment "
"peptides are identified, and they score like a population that cannot be present should, which "
"is what makes the margin in panel B do work rather than simply discard the weakest "
"identifications of every kind. 7 of the 70 entrapment peptides clear the margin of 10 and 90% "
"sit below it.",

"*Pooled false-discovery proportion across all plexes is 0.58% as identifications, six among "
"1,039 and appearing in four of 42 plexes (Extended Results). Panel b's 0.83% is seven among 839 "
"distinct peptides, in five plexes. The difference is one peptide at a margin of 10.53, which "
"clears the threshold in the recomputed margins behind panel B and not in the per-plex summary "
"behind panel C. Neither number should be quoted with the other's counts.*",
]


# ---------------------------------------------------------------- 2.4 / Figure 4
# Original panels: a design, B and C the paired display ungated and gated, d the permutation
# ladder, e occupancy dose-response (cut), f the gate step by step (cut).
# Bounds kept: that the gate is applied to both groups at once; "one property measured twice";
# that the spread shows it is not range restriction AND that two shape features are exceptions
# where it partly is; that panel D measures whether the two move together rather than whether
# either is stable; that both arms must be reported because either alone misleads; and the whole
# closing paragraph, which is the conditional the section rests on.
SPINE["2.4"] = [
"Two questions are routinely run together in work on covalent lysine chemistry: what makes a "
"lysine reactive enough to be labelled by a probe at all, and what makes a labelled lysine "
"blockable by a competing molecule [[hacker2017]]. One contrast cannot answer both, so two were "
"run against two different control groups (Figure 4A). Reactivity compares probe-labelled lysines "
"against every other lysine on the same protein; ligandability compares engaged lysines against "
"all labelled ones. Both run through the same ladder of structural-confidence gates, always "
"applied to both groups at once, since gating one side would manufacture a difference rather than "
"reveal one.",

"Ungated, the two contrasts are indistinguishable. The mean distance between a feature's two "
"effects is the same for charge and for shape features, and within each class the two contrasts "
"track each other almost perfectly (Figure 4B). An exact permutation over all 5,005 ways of "
"splitting fifteen features six against nine gives p = 0.412. Read on its own, this panel says "
"that reactivity and ligandability are one property measured twice.",

"At a confidence of 90 or above the two contrasts come apart, and only across the shape features "
"(Figure 4C). The mean pair distance for shape nearly doubles while charge is essentially "
"unchanged, the shape class correlation falls from 0.97 to 0.35 while charge holds at 0.98 to "
"0.99, and the permutation value rises to p = 0.003.",

"The packing features are where it happens, and the spread shows it is not range restriction: in "
"each case the feature's own standard deviation holds or widens under the gate, so the reactivity "
"effect is not shrinking because the spread shrank (Extended Results). Two shape features are "
"exceptions and are reported as such. Exposure and distance to the centroid lose 12% and 29% of "
"their spread under the gate, so for those two the collapse is partly range restriction rather "
"than artefact removal.",

"What Figure 4E measures is whether the two contrasts move together, not whether either is "
"stable. Several charge effects change substantially under the gate, but they change in step on "
"both sides. The shape features do not: their reactivity arms collapse toward zero or reverse "
"sign while their ligandability arms hold direction. That asymmetry, not the size of any single "
"effect, is the result.",

"Both arms of the test have to be reported together, because either alone misleads. Across the "
"ladder the shape coupling falls monotonically while charge holds throughout (Figure 4E). The "
"separation therefore does not exist before the gate is applied, and it is not manufactured by "
"the gate keeping only the strongest sites, which is the control in figure S5g. Nor does it "
"depend on the 60% cut existing, since the shift from the lowest occupancy tier rises "
"monotonically across all ten tiers with no discontinuity where the cut falls (Extended Results).",

"The claim is conditional. Where predicted structure is trustworthy, shape distinguishes "
"ligandability from reactivity and charge does not; ungated, neither does, because the apparent "
"shape contribution to reactivity is carried by sites whose geometry cannot be measured. This "
"does not say that shape is irrelevant to reactivity in general; it says that the evidence for it "
"in these data comes from the sites where the evidence is weakest. It also inherits the gate's "
"own bias: AlphaFold confidence is not random with respect to biology, and a claim that holds "
"only above a confidence of 90 is a claim about ordered protein.",
]

# The per-feature packing numbers and the occupancy dose-response are the arithmetic behind two
# claims the body states with their bounds. They move; the bounds stay above.
EXTENDED["2.4"] = [
"Packing features, feature by feature, reactivity against ligandability under the gate: crowding "
"runs +0.197 to −0.042 in reactivity while holding at +0.143 to +0.100 in ligandability; atom "
"density +0.244 to +0.003 against +0.157 to +0.109; polar contacts +0.183 to −0.002 against "
"+0.185 to +0.197. The feature's own standard deviation holds or widens in each case (crowding "
"4.56 to 4.70, polar 1.66 to 1.81, aromatic 1.08 to 1.25). Exposure and distance to the centroid "
"are the two exceptions, losing 12% and 29% of their spread.",

"Charge effects under the gate, for comparison with the shape features: anionic contacts "
"strengthen from −0.307 to −0.454 and salt bridges from −0.221 to −0.327, while "
"the ligandability arms run −0.210 to −0.336 and −0.183 to −0.286, so both sides "
"move in step.",
]

LEGEND[4] = [
"Reactivity and ligandability are measured on the same lysines and with the same estimator, one "
"difference per protein tested against zero, so the two contrasts are directly comparable. The "
"figure compares site sets; individual residue numbers within them are provisional (Methods). "
"Every structural gate is applied to both groups of a contrast at once. Effects are "
"within-protein standardised differences, measured in cell lines.",

"(A) Two comparisons against two control groups, run before and after a structural-confidence "
"gate so the gate's effect is visible rather than assumed. 245,684 lysines across 6,201 proteins, "
"30,805 of them labelled (12.5%). Fifteen features, six charge and nine shape, the identical set "
"on both ladders [[tien2013,olsson2011]]; AlphaFold terms are excluded because confidence is the "
"gating variable rather than a site property.",

"(B, C) The paired display, ungated and at a confidence of 90 or above. Each feature's two "
"effects are drawn on one line, a filled circle for the effect on being labelled and an open "
"diamond for the effect on being blocked, so the movement under the gate reads directly. Average "
"distance between the pair: ungated charge 0.074, shape 0.075; gated charge 0.068, shape 0.126.",

"(E) The permutation test, over the gate ladder and against the null. Two stacked axes. Agreement "
"is the correlation, across the features of one kind, between a feature's effect on reactivity "
"and its effect on ligandability. The upper axis is that agreement as the gate tightens: charge "
"stays at 0.98 to 0.99 throughout while shape falls 0.97, 0.90, 0.77, 0.35 at gates of none, 50, "
"70 and 90. Printed above each gate is the exact one-sided permutation p for the gap between the "
"two kinds at that gate: 0.412, 0.034, 0.009 and 0.003. The lower is the null itself, all 5,005 ways of "
"splitting these fifteen features six against nine, with both observed values marked: +0.009 "
"ungated and +0.639 gated. Both are stated because either alone misleads, and figure S5g is the "
"control against the gate keeping only the strongest sites.",
]

# ---------------------------------------------------------------- 2.7 / Figure 7
# Panel b is cut, but the two-features-are-charge result is the section's spine claim, so it
# stays in the body and cites Extended Results instead of a panel. Its q values and the full
# twenty-two-feature test move with the panel. The power bound (eight of twenty-two significant
# for each probe on its own) is a bound on a spine claim and does not move.
SPINE["2.7"] = [
# opener, no panel
"Charge was never among the descriptors tested for ligandability (Figure S7F), so the question is "
"open. It is put to the reagent rather than to the compound, because in a competitive experiment "
"the probe decides which lysines are visible to be competed for at all.",
# panel A
"Two probes were compared that share a desthiobiotin-hexanoyl core and differ only in the acyl "
"leaving group (Figure 7A). One carries a 4-sulfo-2,3,5,6-tetrafluorophenyl ester, anionic at "
"physiological pH [[hacker2017]]; the other an N-(cyanomethyl)-N-(4-methylbenzenesulfonyl)amide, "
"neutral overall [[ryan2024]]. Both were applied to the same HeLa proteome, and the discriminating "
"contrast is drawn within the proteins carrying sites of both kinds, so a probe that happened to "
"sample a different set of proteins cannot create a difference.",
# the spine claim; panel B itself is in Extended Results
"Two site features of twenty-two differ between the probes, and both are charge: lysines labelled "
"only by the anionic probe sit among more cationic residues, and in more positive local environments "
"overall. No pocket feature and no shape feature differs (Extended Results). That pattern is "
"what a "
"steering effect should look like and what a general difference in proteome sampling should not. "
"Neither probe is individually quiet: measured against unlabelled lysines, eight of the twenty-two "
"features are significant for each (Figure S8A), so two of twenty-two is a statement about the "
"contrast between the probes rather than about the power to detect one.",
# panels C, d
"Bootstrap intervals over proteins separate the reagent's contribution from the site's "
"(Figure 7C, D). Cationic attraction is the reagent: the neutral probe's interval brackets zero and "
"excludes the anionic probe's estimate. Acidic avoidance, predicted pKa and salt bridges are the "
"site, with the two probes indistinguishable on all three. Net local charge is mixed, and it is the "
"one feature for which a partial attribution is the honest reading. Reagent attribution is not "
"inferred from one probe reaching significance while the other does not, because a difference "
"between significant and non-significant is not itself a significant difference, which is why the "
"panels show intervals rather than asterisks.",
# no panel: the interpretive correction
"This changes how an earlier reading of the same data has to be worded. Labelled lysines sit in "
"cationic, acid-poor environments with fewer salt bridges and lower predicted pKa, and the natural "
"reading of that is nucleophilicity. The decomposition says the reading is right for the pKa and "
"salt-bridge terms and wrong for cationic attraction, which tracks the reagent. Both effects are "
"real; they are not the same effect.",
# panel E
"Concavity, by contrast, belongs to the site. Labelled lysines sit in more enclosed environments "
"than unlabelled ones, and the two probes do not differ on it, so this is a property of where "
"reactive lysines are and not of what reaches them (Figure 7E). Two bounds belong with that result. "
"Enclosure and relative solvent accessibility share about 74% of their variance across the 167,057 "
"sites, so the concavity effect is largely, though not wholly, a solvent-accessibility effect "
"(Figure S8B); and the effect is reported within protein, because pooling across sites reverses its "
"sign (Figure S8C).",
# no panel: the enolization bound
"One assumption bounds the charge result rather than threatening it. The neutral probe carries "
"hydrogens α to an electron-withdrawing acylsulfonamide and may enolize at physiological pH, "
"acquiring partial anionic character; its pKa is not known. If it does, the two probes are less "
"different than assumed, and the measured reagent effect is an underestimate: it rises from the "
"neutral-assumption value to +0.46 if the probe is substantially ionized, which is a sensitivity "
"range and not a confidence interval (Figure S8D). The assumption cannot explain the result away.",
# closer, no panel
"A chemoproteomic probe is not a neutral window onto lysine reactivity: part of what a "
"lysine-directed experiment measures is where its reagent is drawn. Coverage and reactivity are "
"entangled at the point of capture, so a probe panel is a sampling variable to be varied "
"deliberately rather than held fixed by convention.",
]

LEGEND[7] = [
"Two probes sharing a desthiobiotin-hexanoyl core and differing only in the acyl leaving group, "
"applied to the same HeLa proteome. Twenty-two site features measured on 167,057 lysines at "
"pLDDT ≥ 70 [[jumper2021,varadi2022]], including relative solvent accessibility [[tien2013]], "
"predicted pKa [[olsson2011]] and pocket descriptors [[leguilloux2009]]. Labelling is expressed as "
"% occupancy.",
"(A) The design. dTB-STP, a 4-sulfo-2,3,5,6-tetrafluorophenyl ester and anionic, against a neutral "
"N-(cyanomethyl)-N-(4-methylbenzenesulfonyl)amide probe. 6,294 sites labelled only by the anionic "
"probe, 5,851 only by the neutral one, 12,298 by both; the discriminating contrast is drawn within "
"the 383 proteins carrying both.",
"(C, D) What each probe prefers, and which part of the charge effect belongs to which. Effects are "
"standardized differences across proteins; intervals are 95% percentile bootstrap confidence "
"intervals, resampling proteins (Methods). Cationic residues within 8 Å, d = +0.233, and net charge within 8 Å, d = +0.180, are the "
"two features that separate the probes; the next largest, aromatic contacts at d = 0.118, does not. "
"Zero of five pocket features and zero of nine shape features differ. The q values and the full "
"twenty-two-feature test are in Figure 7B.",
"(E) Concavity belongs to the site, not the reagent. Enclosure d = +0.073 (q = 1.1 × 10⁻⁵) for "
"labelled against unlabelled lysines, and d = 0.019 (q = 0.71) between the two probes. See also "
"Figures S8B and S8c.",
]

EXTENDED["2.7"] = [
"Site environment of labelled lysines, the effect sizes behind the nucleophilicity reading, "
"measured at pLDDT ≥ 70: fewer acidic neighbors (d = −0.46), more cationic ones (+0.40), fewer salt "
"bridges (−0.32) and lower predicted pKa (−0.30).",
]

# ---------------------------------------------------------------- 2.8 / Figure 8
# Panels B and D are cut. The anchoring-artifact paragraph is the section's honesty passage and
# stays almost intact; the crystallographic audit is a secondary audit whose detail and caveats
# travel with panel D, leaving one sentence of result in the body. Every bound in the
# what-follows paragraph stays, since each one is a limit on the section's own claim.
SPINE["2.8"] = [
# opener, no panel
"The acquired lysines of Figure 1 have been tested against nothing. Do they sit inside the same "
"population as the lysines these rules were measured on, or does mutation create lysines in "
"environments unlike the ones the atlas has sampled?",
# panel A
"The comparison had to be anchored before it could be made. Our first version of it measured union "
"positions at Cβ and measured lysines at Nζ, and returned a separation that would have been the "
"strongest structural result in the paper. It is an artifact: Cβ and Nζ lie about 6.4 Å apart in an "
"extended lysine side chain, so the two populations were being measured from different atoms. "
"Figure 8A shows the artifact beside its correction, and three legitimate anchorings agree with "
"each other (Figures S9A and S9b); the mismatched arm is retained there as a negative control, "
"because it is the analysis we first ran. Where side chains had to be built, all-trans rotamers "
"reproduce the amine positions of 30,805 measured lysines to a median of 4.26 Å, and all-trans is "
"the weakest possible construction, so that is an upper bound on the error rather than an estimate "
"of it (Figure S9E).",
# panels B (moved) and c
"On a common anchor the three populations are difficult to tell apart. Four of the nine features "
"have identical medians across acquired, engaged and unengaged lysines, and a fifth is identical "
"between acquired and engaged sites (Figure 8B). Where the populations do differ they differ by "
"less than one standard deviation of the acquired distribution. Reported without a threshold, the "
"probability that a randomly chosen acquired lysine lies closer to the center of the engaged "
"population than a randomly chosen engaged lysine does is 0.443, against 0.5 for populations that "
"cannot be distinguished at all (Figure 8C). The reading does not depend on where a boundary is "
"drawn: at the 50, 75, 90 and 95% engaged regions acquired lysines sit inside the engaged region "
"less often than the labelled-but-unengaged control does, at every cut and under either "
"anchoring and by 3.9 to 11.5 percentage points, so acquired lysines fall below both "
"comparators rather than between them (Figure S9F).",
# panel D is in Extended Results; the result itself stays
"Six of these sites have crystal structures, and they behave the same way. Local Cα displacement in "
"a co-local shell, normalized to each protein's own wild-type floor, puts five of the six at that "
"floor; SHP2 E76K is the exception, and it is the one site carrying a substitution series, where "
"other substitutions at the same position perturb the environment as much (Extended Results).",
# no panel: what follows and what does not
"Acquired lysines are not structurally unusual, which removes the cheapest reason to expect the "
"engagement rules measured in Figures 4 to 7 not to reach them. It does not establish that they can "
"be engaged. Eight of the nine features compared here are among the fourteen from which the "
"composite model of Figure 5 was built, every one except the count of neighboring lysines, and that "
"model reaches an AUC of 0.617 against a measured reproducibility ceiling of at least 0.72. "
"Features that predict engagement that poorly cannot be turned around and used to infer "
"engageability from similarity. The engaged and unengaged comparators come from the three "
"chemoproteomic screens and describe those cells, not a tumor and not a patient. What this "
"comparison supports is narrower and more useful: nothing about the chemistry of these sites "
"explains why they are absent from the record. The barrier is detection.",
# panel E
"One acquired lysine makes that concrete. AKT1 E17K, the transforming "
"pleckstrin-homology-domain mutation identified in breast, colorectal and ovarian cancer "
"[[carpten2007]], has since been engaged selectively in cells through lysine targeting and "
"neo-zinc chelation [[craven2025]], and it returns zero spectra here (Extended Results). That "
"negative "
"rests on a single experiment: of the five DepMap models carrying the mutation, one is in the "
"375-line proteomic set, so its plex is the only one of the 42 in which the mutation was present to "
"be found. It is one informative observation, not a survey. It is worth stating anyway, because the "
"one acquired lysine that has been drugged is the one this workflow did not see.",
]

LEGEND[8] = [
"Nine environment features at 156,781 union positions with an AlphaFold confidence of at least 70 "
"[[jumper2021,varadi2022]], against the engaged and unengaged lysines measured in the three screens "
"[[abbasov2021,bracken2024,ryan2024]] (Methods). Engagement is expressed as % occupancy.",
"(A) One reference atom, on both sides. The first version of this comparison measured union "
"positions at Cβ and measured lysines at Nζ, and returned a separation of 2.256 against a control "
"of 0.330, a 6.84-fold ratio that is an artifact of the mismatch. The three legitimate anchorings, "
"both sides at Cβ, both modelled at Nζ, and both at Cβ under a symmetric convention, give 0.571, "
"0.535 and 0.583 against controls of 0.410, 0.331 and 0.407. The panel shows the artifact and the "
"correction.",
"(C) The distance between acquired and engaged lysines, against permutation nulls. Reported "
"threshold-free: the probability that a randomly chosen acquired lysine sits closer to the engaged "
"center than a randomly chosen engaged lysine does, against 0.5 for indistinguishable populations. "
"The populations overlap heavily.",
"(E) AKT1 E17K, engaged selectively in cells by others, and undetected here. Site geometry from PDB "
"8UVY, 54 atoms in the local shell. Searched as the three peptides the databases carried, two "
"spanning the new lysine and one immediately downstream. The one model carrying the mutation that "
"is in the proteomic set is KU-19-19 (RRID:CVCL_1344). The point is the contrast between what is "
"chemically reachable and what this workflow can see.",
"The engaged and unengaged comparators come from the three chemoproteomic screens (the fragment "
"screen in Ramos and MDA-MB-231 [[abbasov2021]]; the meroterpenoid screen in MDA-MB-231 and "
"MDA-MB-436 [[bracken2024]]; the CyMPSA screen in HeLa [[ryan2024]]). Nothing in this "
"figure is a statement about a tumor or a patient.",
]

EXTENDED["2.8"] = [
"Boundary independence, the numbers behind the threshold-free reading: acquired lysines fall inside "
"the 50, 75, 90 and 95% engaged regions at 42, 69, 82 and 88% under one anchoring and 45, 67, 83 "
"and 89% under another, against a labelled-but-unengaged control running 47, 73, 89 and 95% and 50, "
"78, 92 and 96%.",
]

# ------------------------------------------------- 2.5 + 2.6 / the merged Figure 5
# The merge takes the FIRST panel-less paragraph of 2.6 as the opener, then the panel
# paragraphs in MERGED_ORDER (6a, 6b, 5a, 5c, 6d, 5e), then the LAST panel-less paragraph of
# 2.5 as the closer. Any other panel-less paragraph is carried to Extended Results, so the two
# spines below are written with exactly one panel-less paragraph each, in those positions.
# Citations are written in ORIGINAL numbering; the single global pass maps them.
SPINE["2.6"] = [
# opener for the merged section, panel-less, must come first
"If engagement is a property of the site rather than of the molecule, the same lysines should be "
"preferred by compounds that have little else in common. Figure 4 measured which site properties "
"separate engaged lysines from merely labelled ones; what follows is how well they predict "
"engagement, and whether any property of the compound predicts which lysines it prefers.",
# panel 6a -> merged a
"Across the three screens [[abbasov2021,bracken2024,ryan2024]], 0.49% of site × compound "
"measurements reach the 60% occupancy that defines engagement here, and the distribution is smooth "
"through 60%, so the threshold is a line drawn through a continuum rather than a boundary in the "
"data (Figure 6A). For each compound separately, the sites it engaged are compared with the sites "
"it was quantified at and left below the cut, and each site feature is standardized within that "
"compound, so what is compared is compounds' preferences rather than their site lists, and a "
"compound tested on a different set of lysines cannot manufacture a difference.",
# panel 6b -> merged b
"Thirteen of fourteen site features point the same way in every compound: engaged lysines are "
"consistently more buried, more crowded, more hydrophobically packed and lower in predicted pKa "
"(Figure 6B). The single exception is the count of nearby lysines within 8 Å, which is a property "
"of the neighborhood rather than of the lysine in question. The three charge features clear "
"significance but are the three smallest of the thirteen that do, and two of them do not strengthen "
"with occupancy (Figure S7B).",
# panel 6d -> merged e
"No tested compound descriptor predicts which lysines a compound prefers. Correlating every "
"compound descriptor against every site-feature preference gives 21 of 84 pairs clearing q < 0.05 "
"when the three screens are pooled, against 4.2 expected, but 0 of 84 within one screen alone and "
"0 of 72 in the larger within-dataset slope analysis (Figure 6D). The pooled positives are "
"provenance rather than chemistry: molecular weight all but names the screen a compound came from, "
"each screen is a different cell line in a different sample state, and the descriptors that "
"separate the three screens are the descriptors that produced the pooled positives "
"(Extended Results).",
]

SPINE["2.5"] = [
# panel 5a -> merged c, first of two
"A logistic model was fitted on fourteen features, thirteen describing the lysine's environment and "
"one describing how confidently the residue was modelled, and cross-validated by protein rather "
"than by site, because lysines on one protein are correlated and a site-level split leaks "
"information between training and test. Held out, the composite reaches an AUC of 0.617, and its "
"optimism against the training fit is +0.002, so the number is not depressed by overfitting "
"(Figure 5A). Refitting on the thirteen environment features alone gives 0.611, which forecloses "
"the objection that the model is partly predicting how well a site was modelled rather than what "
"the site is.",
# panel 5a -> merged c, second of two: the ceiling and the bound on the shortfall
"Whether 0.617 is disappointing depends on what a perfect model could reach, and that limit is set "
"by the assay. Engagement is the outcome of a measurement: if the same lysine is assayed twice and "
"the two assays disagree, no predictor can agree with both. Splitting each site's compound "
"measurements into two independent halves, taking the label from one half and the score from the "
"other, gives the AUC a perfect structural model would reach if structure captured everything the "
"assay captures. At the deepest split the data support, that ceiling is 0.723 and the composite "
"falls 0.106 short of it (Figure 5A). The shortfall is the smallest the data permit, on two counts. "
"Each half carries fewer compounds than the full screen and is therefore noisier than the real "
"label, so the ceiling is a lower bound and is still rising at the deepest split; and 0.723 is the "
"lowest of the three screens, the other two sitting between 0.90 and 0.95 in every band where they "
"are populated. The comparison is drawn against the weakest ceiling in the data, and is quoted as "
"≥ 0.72.",
# panel 5c -> merged d
"No operating point on the curve would change a decision. At a 5% false-positive rate the composite "
"recovers 11% of engaged sites, and at 10% it recovers 19% (Figure 5C). Nothing on this curve would "
"alter which lysines a screen went after, which is the practical content of the result: "
"ligandability at these sites, in these cell lines, has to be measured.",
# panel 5e -> merged f
"What separates ligandable sites is breadth of engagement rather than a structural signature. "
"Grouping sites by how many of the 34 tested chemotypes engage them at ≥ 60% occupancy, the "
"distance from a never-engaged site rises monotonically across the four groups, from 0.49 to 2.65, "
"against nulls of 0.12 to 0.49 from label shuffles (Figure 5E). The direction is secure across the "
"gradient; the magnitude at the far end rests on 82 sites, and adjusting for how often each site "
"was tested does not remove it (Figure S6E).",
# closer for the merged section, panel-less, must come last
"Both halves of the heading are bounded. Reproducible is a measurement rather than an assumption, "
"and the same split-half analysis supplies the ceiling. Poorly predicted is bounded to these "
"features, these compounds and these cell lines. The descriptors tested contain nothing that "
"encodes charge, the largest correlation any of them reaches with a charge-related preference being "
"0.193 against a median of 0.552 for molecular weight, and every shape term among them is achiral "
"by construction, since reflecting a conformer leaves the principal-moment ratios unchanged "
"[[sauer2003]] (Figure S7F). Engagement is a property of the lysine over the descriptors tested, "
"and the analysis is silent on exactly the two axes this paper goes on to test, charge in Figure 7 "
"and chirality in the enantiomer pairs of Extended Results, both of which show molecular properties "
"that do matter. It must not be read as a claim that no property of the molecule matters.",
]

LEGEND[5] = [
"Ligandability measured in three chemoproteomic screens (the fragment screen in Ramos and "
"MDA-MB-231 [[abbasov2021]]; the meroterpenoid screen in MDA-MB-231 and MDA-MB-436 "
"[[bracken2024]]; the CyMPSA screen in HeLa [[ryan2024]]), with engagement defined as "
"≥ 60% occupancy, 100 × (1 − compound / DMSO). The composite is a logistic model on fourteen site "
"features, cross-validated by protein (Methods); the site features of panel B are a different set, "
"and the two are never pooled. Intervals on model AUCs are the 2.5th to 97.5th percentile across 20 "
"held-out splits, and on ceilings across half-splits (Methods).",
"(A) The measured model against the reproducibility ceiling. AUC 0.617 [0.594, 0.640] against a "
"ceiling of 0.723 [0.706, 0.744] derived from agreement between independent halves of the same "
"screen, a gap of 0.106; train AUC 0.619. 30,805 lysines on 6,201 proteins, 3,047 of them (9.9%) "
"engaged by at least one compound. The ceiling drawn is the fragment screen's at twenty compounds per half, the "
"deepest split the data support. That gap is reproducible signal the structure does not reach. It is not overfitting, since optimism is +0.002 AUC. The ceiling is still rising at the deepest split, so 0.723 is a lower bound. See also Figures S6B and S6D.",
"(C) What the model buys at any operating point. Receiver-operating curve for the composite, every "
"site predicted once, with the band spanning 95% of 20 held-out splits. There is no operating point where this model is useful. At 5% false positives it recovers 11% of engaged sites, and at 10% it recovers 19%.",
"(E) Breadth across chemotypes. Sites grouped by how many of the 34 tested chemotypes engage them "
"at ≥ 60% occupancy, against multivariate distance from a never-engaged site: one compound, "
"n = 1,529, 0.491; two or three, n = 624, 0.913; four or five, n = 142, 1.215; six or more, n = 82, "
"2.650. Gray is the null over 2,000 label shuffles (0.120, 0.183, 0.369, 0.487); all four bins "
"clear it at p < 0.0005, the resolution limit of 2,000 shuffles. What marks a site is "
"not which chemistry engages it but how many can. See also Figure S6E.",
"Engagement measured in cell lines. The figure compares site sets; individual residue numbers "
"within them are provisional (Methods).",
]

LEGEND[6] = [
"(A) What is measured, and what is done with it. A probe labels lysines; a compound competes with "
"the probe, and occupancy is the share of probe labelling it competes away. All 1,642,330 "
"site × compound measurements in the three screens are shown, 8,019 of them (0.49%) reaching 60%. "
"The median compound quantifies 6,332 lysines and engages 14. Right: the design for one compound, "
"its engaged sites against the sites it was quantified at and left below the cut.",
"(B) Thirteen of fourteen site features point the same way in every compound. Mean standardized "
"shift between engaged and below-cut sites with 95% intervals, across n = 84 compounds. Rust marks "
"the single exception, nearby lysines within 8 Å (+0.038, q = 0.097). See also Figure S7B.",
"(D) No compound property predicts which lysines a compound prefers. Every "
"descriptor × site-feature-preference correlation, one point per pair, with rust marking pairs "
"clearing q < 0.05. Pooled across the three screens (82 compounds), 21 of 84 pairs clear it against "
"4.2 expected; within the fragment screen alone (62 compounds), 0 of 84; in the within-dataset slope analysis "
"(173 compounds), 0 of 72 against 3.6 expected. Descriptors were computed with RDKit [[rdkit]]. "
"This is a statement about the descriptors that were tested, which Figure S7F and note N8 list.",
]

# Carried out of the counted legends by the second pass, 19 Oct. Nothing here is a bound: the
# bound each one supports stays in the main text or in its legend.
EXTENDED["2.2"] = [
"The ten driver hotspots among the 617 covered sites, all of them on the missed side: BRAF V600K "
"(IGR1, MDST8), NRAS Q61K (HT1080, SKNAS), PIK3CA E542K (CAL51), PTEN I67K (LS180), and TP53 "
"M133K (HEL), R280K (MDAMB231), E258K (MEWO) and E285K (RPMI8226).",
]

EXTENDED["2.3"] = EXTENDED.get("2.3", []) + [
"The entrapment audit restricted to the acquired-lysine class: 66 peptides compete with five "
"entrapment at margin 0, 7.04% against 4.58% for the rest of the variant class, a difference the "
"counts cannot resolve (Fisher exact, odds ratio 1.58, p = 0.38). Past the margin the class "
"contributes zero entrapment against the 38 variant peptides quoted in the main text, which is "
"what bounds its proportion at 7.9% by the rule of three.",
]

# ================================================================ back matter
# Cell Press's components list has no Conclusion section, so the Conclusion is folded into the
# last Discussion paragraph rather than kept as a heading. House rule observed here: no bold on
# a topic sentence, which the original back matter used throughout.
INTRO = [
"Covalent chemistry directed at lysine has reached the point where the limiting question is no "
"longer whether a lysine can be engaged but which one to engage. Reactivity and ligandability "
"vary enormously between sites [[hacker2017]], and three subsequent screens mapped thousands of "
"ligandable lysines with chemically distinct probes [[abbasov2021,bracken2024,ryan2024]]. "
"Newer aminophilic chemistry has enlarged that set again, quantifying over 7,000 modifiable "
"lysines on its own [[zhou2025]]. None of that work supplies a principled way to choose one. Selectivity in a covalent campaign "
"is usually argued from the chemistry of the molecule; the alternative is to choose a residue "
"that exists only in the diseased cell.",

"Cancer supplies such residues. A residue created by missense mutation is present in tumor cells "
"and absent from every normal cell in the same patient, which is the logic that made KRAS G12C "
"tractable [[ostrem2013]]. That argument has been generalized for cysteine: acquired cysteines "
"have been catalogued across cancer genomes and pursued chemoproteomically as a target class "
"[[desai2024]], and the experimentally observed cysteines are kept as a curated community "
"repository carrying functional, druggability and genetic-variant annotation [[boatner2023]]. "
"Lysine has not been counted the same way, although mutation creates more of "
"them, 295,411 acquired lysines against 211,386 acquired cysteines in the union of three variant "
"catalogues. On net flux the order reverses, because cysteine is rarely lost and lysine is not; "
"neither residue is exceptional, and the point is only that one of the two has been catalogued "
"as a target class and the other has not.",

"That E to K is common is not a discovery of this work, and the reason is unremarkable: "
"glutamate is abundant and the dominant mutational process in cancer is C>T. The question is not "
"whether these lysines are numerous but whether they are worth targeting, and the obvious "
"counter-position is that they are passengers, residues mutation produces incidentally on "
"proteins under no selection.",

"A second question is why nobody has counted them. The answer is a measurement problem with two "
"parts, both specific to this class rather than general to proteomics. Lysine is a tryptic "
"cleavage site, so creating one redraws the peptide boundaries: the mutant peptide is not the "
"wild-type peptide with a letter changed but a different species with a different mass, length "
"and charge, which a reference database cannot contain and a conventional search cannot find. "
"And the error control that would license the identifications it does make has mostly been "
"tested in the wrong place, since target-decoy rates are calibrated over a whole search "
"[[elias2007]] and variant peptides are a small, atypical subset of one. Strategies that control "
"error separately for the variant space exist [[desai2024]], and entrapment methods that measure "
"the realized rate within a subset now exist [[wen2025]]. What is rarely done is to turn that "
"measurement on a variant search and report what comes back.",

# REWRITTEN 9 November. The paragraph it replaces opened "This paper is a reanalysis. It
# generates no new data", then listed what the work does not report and closed on numbers
# being "taken on trust". Three sentences of what the paper is not, in the position where a
# reader decides whether to keep reading. The scope is still here and still bounded, as the
# negation pair: the data are public, the measurement is not.
"This work turns it on one. The data are public; the measurement is not. Re-searching 42 "
"published TMT plexes of the CCLE proteomes [[nusinow2020]] against a mutation-aware database "
"measures how much of the acquired-lysine class a search can see, and burden-matched entrapment "
"measures the realised error on what it returns. Three published chemoproteomic screens answer "
"the chemical question on the same residues, across four cell lines and 34 chemotypes. Every "
"script, intermediate and panel-level table is deposited, so the census extends to a new search "
"engine, a larger variant catalogue or the next screen rather than having to be rebuilt.",
]

DISCUSSION = [
"Mutation-created lysines are numerous, largely absent from common germline variation, enriched "
"on driver hotspots, written by identifiable mutational processes, and almost entirely missing "
"from the proteomic record. Every one of those is measured above, and the error rate that would "
"license even the few identifications recovered does not hold at the reported filter.",

"The closest prior work approached this from two directions, and neither reaches the question "
"asked here. The three chemoproteomic screens reanalyzed in this paper mapped ligandable lysines "
"across several proteomes [[abbasov2021,bracken2024,ryan2024]], building on the reactivity and "
"ligandability survey established earlier [[hacker2017]], but they searched reference "
"databases, so a mutation-created lysine was outside their search space by construction. From "
"the other direction, AKT1 E17K has been engaged selectively in cells through lysine targeting "
"and neo-zinc chelation [[craven2025]]: an acquired lysine drugged by target-based design, one "
"site at a time, without any claim about the class it belongs to. The measurability question has "
"been asked and addressed for cysteine, where a chemoproteogenomic pipeline with variant-aware "
"databases and staged error control catalogued acquired cysteines directly [[desai2024]]. "
"Coverage itself has been audited there too, by pairing predicted side-chain accessibility with "
"a meta-analysis of eighteen public cysteine profiling datasets, which found the sampled "
"residues biased by warhead chemistry and concluded that current approaches reach only a small "
"fraction of the ligandable population [[white2023]]. Lysine "
"has had no equivalent, and it is the harder case for a specific reason: an acquired cysteine "
"leaves the tryptic peptide it sits in intact, while an acquired lysine creates a cleavage site "
"and produces a different peptide altogether.",

"One of those obstacles is partly dissolved by the experiment that would be used to look for it. "
"An acquired lysine creates a tryptic cleavage site in a shotgun experiment but not in a "
"lysine-chemoproteomic one, because the probe acylates the ε-amine and trypsin does not cleave "
"after an acylated lysine, so a labelled acquired lysine leaves the peptide boundaries near "
"wild-type while the peptide is streptavidin-enriched.",

"The passenger reading of this class deserves a direct answer. Glutamate to lysine is the "
"commonest substitution in cancer for reasons that carry no biological weight, and it would be "
"reasonable to conclude from that alone that these are incidental residues. One observation "
"argues otherwise: acquired lysines carry the driver-hotspot flag at an odds ratio of 1.68 "
"against other substitutions, and that flag is conservative by construction, since the method "
"behind it exists to strip recurrent passenger sites from hotspot lists [[hess2019]]. The "
"signature composition of those hotspots is not independent support and should not be offered as "
"such: it compares lysine-creating hotspots against other hotspots, and since lysine is made "
"85.5% by C>T and never by T>C, any lysine-creating set would show the same enrichment whether "
"or not it were under selection. The test that would separate the two, acquired-lysine hotspots "
"against acquired-lysine non-hotspots, is not in this paper. What the odds ratio establishes is "
"that these sites are selected more often than other substitutions; it does not establish that "
"the affected proteins are expressed, detectable, or reachable by a molecule.",

"Which cells carry the class is not arbitrary, and it can be chosen. The four lines the three "
"screens ran in carry 56 unique acquired lysines between them, against a maximum of 623 across "
"1,954 DepMap models, and ranked on the sites a screen could actually reach the leading lines "
"carry 59 to 83 times what HeLa does. Which cells those are, which sites recur across "
"models, how far burden survives gating to expression, protein detection and dependency, and "
"what is independently known about the recurrent driver-flagged sites are set out in Extended "
"Results with the per-line and per-site tables.",

"A probe panel is a sampling variable rather than a fixed choice, and that is as far as these "
"data reach. Establishing which panel would maximize coverage requires prospective experiments a "
"reanalysis cannot contain, and a two-probe comparison supports the principle without licensing "
"a design.",

"What would actually recover this class is not what intuition suggests. Deeper or "
"enrichment-based proteomics addresses abundance, and abundance does not look like the limiting "
"gate: it separates sites whose mutant peptide was found from those whose region was never "
"covered by Cliff's δ = 0.31, on the 1,769 of 8,470 observations that carry an abundance value "
"and only 40 found sites. A subset of the failures is instead associated with long tryptic "
"peptides, the two distributions diverging in the tail beyond twenty residues rather than "
"shifting as a whole, and no effect size is available for that comparison, so neither test is "
"decisive. What can be said is that changing which proteins are sampled leaves the peptide "
"boundaries untouched, and that three measured results point at the same lever. Half this class "
"has no fully cleaved mutant peptide above the seven-residue floor the searches enforce, so its "
"only searchable species carries a missed cleavage. An internal missed cleavage is the strongest "
"negative predictor of identification measured here, at a standardized effect of −0.621. And "
"the cleavage that fails is specifically at arginine, 29.4% of identified peptides whose first "
"basic residue is arginine carrying a missed cleavage against 15.0% for lysine. What those three "
"identify is completeness at the lysine and arginine sites trypsin already targets, not a "
"different specificity: a protease cutting at additional residues would shorten peptides "
"further, and short is where this class already fails. Digesting the same proteome more "
"completely is a prediction this paper does not test.",

"The limitations are specific, and several of them bound results this paper otherwise leans on. "
"Everything measured in cells is measured in cell lines, 375 in the proteomic re-search and four "
"across the three chemistry screens, and no result is a statement about a tumor or a patient. "
"Figure 1 is a census of reported variants and a population frequency comparison rather than a "
"measurement in any cell. "
"The screens are not saturated, so every ligandability figure is a lower bound: a site recorded "
"as unengaged was unengaged by the compounds tested. The structural separation of reactivity "
"from ligandability emerges progressively as the confidence gate tightens and is strongest above "
"pLDDT 90, which makes it a claim about ordered protein, and AlphaFold confidence is not random "
"with respect to biology. The compound descriptors tested encode neither charge nor chirality, "
"so the finding that none of them predicts which lysines a compound prefers is silent on exactly "
"the two axes this paper shows to matter. The nine environment features of the structural "
"comparison describe a site's surroundings within 8 to 10 Å; the lysine's own Nζ chemistry is "
"not compared, and at most acquired positions that Nζ is a built all-trans rotamer carrying a "
"median error of 4.26 Å against a reference that is itself predicted, which is why the "
"six-structure check matters, and that check crosses a method boundary at two PI3Kα pairs and "
"rests on a single wild-type pair at NRAS. The cross-dataset version of the breadth result is "
"the weakest evidence in the paper, 25 sites in 21 dimensions with two of six comparisons short "
"of their null, and is reported as a limit rather than as support. And the detections themselves "
"are thin: of 53, none rests on more than one experiment, which is why the error-rate audit "
"exists and why the coverage counts are a lower bound rather than a catalogue of sites.",

"Taken together the results describe a class that is large, selected more often than other "
"substitutions, chemically unremarkable and effectively unmeasured, and they locate the obstacle "
"precisely. It is not that acquired lysines are structurally unusual: nine environment features "
"place them within the ordinary range of the lysines chemoproteomics already engages. Nor is the "
"chemistry unimaginable, since one such lysine has been engaged selectively by a bespoke "
"campaign against a single site [[craven2025]], which establishes possibility there rather than "
"tractability for the class. It is that the same site returns zero spectra in the one plex of "
"the 42 whose cell line carries it, and that a class which cannot be seen cannot be prioritized "
"or argued about with data. Screening does not require proteomic visibility, as the three "
"reanalyzed screens show; knowing which sites exist, in which cells, and with what confidence "
"does. The measurement problem is the tractable one.",
]

# Moved out of the counted Discussion: the prospectus for the next experiment. None of it bounds
# a claim this paper makes, and all of it is already drawn in figure S10 and table S1.
EXTENDED_BACK = [
"Which cells carry the class, and what gating does to it. The four cell lines of the three "
"chemoproteomic screens carry 56 unique acquired lysines between them, against a median of 10 "
"and a maximum of 623 across 1,954 DepMap models: MDA-MB-231 24, MDA-MB-436 19, HeLa seven and "
"Ramos six (figure S10b). Eleven of the 56 fall on a protein the matching screen probe-labels at "
"one or more lysines, and one of them, SF3A3, already yields a lysine engaged above 60% (figure "
"S10e). The lysine measured is not the acquired one and usually sits far from it, a median of "
"271 residues away, so what the eleven establish is that the protein is reachable by a probe "
"rather than that the acquired site is. A re-search of those three screens is worth running, and "
"it would test eleven sites rather than a class.",

"Burden by lineage runs skin, uterus, bowel and lung at the top and eye second from last; at "
"subtype resolution cutaneous melanoma carries a median of 37 against ocular melanoma's three, "
"which is the cutaneous-versus-uveal contrast recovered in cell lines rather than tumors (figure "
"S10a). The asymmetry that makes the class is visible one cell line at a time: across 1,954 "
"models 39,390 lysines are gained against 25,531 lost, and the ratio runs from 6.4 in skin to "
"1.3 in esophagus and stomach (figure S10b). A quarter of models are net-negative. Burden is not "
"the same as opportunity either. Passing an acquired lysine through three gates, the gene "
"expressed in that line, the protein quantified in that line's CCLE proteome, and the gene a "
"CRISPR dependency, leaves 46%, then 31%, then 6% of 2,236 sites (figure S10c). Six percent is "
"not a small number in absolute terms: 132 sites in 27 cell lines are on a protein that is "
"present, detectable and required. The gates also expose a gap in the reference data: only 308 "
"of 1,954 models carry expression, proteomic and CRISPR data together, and none of the ten "
"highest net-gain lines does. The cells that carry most of the class are the least characterized.",

"The sites that recur are the sites this paper already measured. Of the 39,390 unique acquired "
"lysines, 248 carry a hotspot or driver flag, and ranked by how many cell lines carry them the "
"leaders are NRAS Q61K in 36, PIK3CA E545K in 33, PIK3CA E542K in 15, XPO1 E571K in nine, TP53 "
"E286K in seven, TP53 E285K in six and AKT1 E17K in five (figure S10d). Four of the six sites "
"whose structures are measured in the structural comparison are on that list: three of them are "
"its three most recurrent, and AKT1 E17K is seventh, which was not how they were chosen. BRAF "
"V600K, in three models, is an acquired lysine on the most-drugged oncogene in melanoma. Four "
"driver-flagged sites clear all three gates in the models examined, in five site-in-line "
"instances: MTOR E1799K in HEC59, PIK3CA E545K in MCF7, BRAF V600K in MDST8 and in IGR1, and "
"PIK3CD E1021K in REH. PIK3CD already carries an approved drug against another site. What is "
"independently known about these sites is collected in table S1, which is provenance rather than "
"a survey. Twelve of the thirteen are carried by TCGA patients as well as by cell lines (PIK3CA "
"E545K by 263 of the 10,224 participants in the MC3 call set, E542K by 169, NRAS Q61K by 68, "
"AKT1 E17K by 54), and nine carry a functional consequence measured by someone else, among them "
"impaired GTP hydrolysis at NRAS Q61K and accelerated leukemogenesis at XPO1 E571K. The "
"thirteenth is instructive rather than absent: BRAF V600K returns no MC3 patients because valine "
"reaches lysine only through two adjacent nucleotide changes and that call set does not merge "
"adjacent variants, so the site is invisible to it by construction while DepMap reports it in "
"three lines. These are the prospective screen's first targets, and naming them, with the gaps "
"marked, is what a catalogue is for.",

"A choice of cells follows. Lysine chemoproteomics has been run in lines that barely carry the "
"class, while the chemoproteogenomic work on acquired cysteine ran in high-burden lines "
"[[desai2024]]. Ranked on the sites a screen could actually reach, acquired lysines on a gene "
"the line expresses, the leaders are MAMEL46 with 248, MCC13 with 178 and HEC251 with 176, "
"against three in HeLa. A prospective screen in any of them would put the same chemistry against "
"59 to 83 times the reachable sites.",
]

# ---------------------------------------------------------------- late condensation, 19 Oct
# Two paragraphs of 2.3 restated the Introduction's own argument about where target-decoy control
# is calibrated, and the per-platform detail duplicates the Discussion's prior-work paragraph.
# Applied here by index rather than by string, because the source strings wrap mid-phrase.
SPINE["2.3"][0] = (
"Figure 2 recovered a small number of mutant peptides, and a result carried by single-spectrum "
"identifications has to survive an error-rate audit before anything is built on it. A global 1% "
"filter does not settle it, for the reason given above: a rate that holds over a whole search need "
"not hold within a tiny, atypical subset of one [[elias2007]].")

SPINE["2.3"][1] = (
"Chemoproteomic work on mutation-created residues has answered this in two ways, and both control "
"the error: one searches canonical and cell-line-specific mutant sequences in a single pass and "
"filters globally [[takahashi2024]], the other in two stages, removing spectra confidently "
"assigned to canonical sequences before the variant database is searched [[desai2024]]. Neither "
"reports a realized false-discovery proportion measured within the variant class, which is what "
"follows, and what follows is a property of this search and this database, not of theirs.")

# The passenger counter-position is stated twice, here and in the Discussion, where the answer to
# it sits. Keep it once, in the Discussion.
INTRO[2] = (
"That E to K is common is not a discovery of this work, and the reason is unremarkable: glutamate "
"is abundant and the dominant mutational process in cancer is C>T. Whether these lysines are worth "
"targeting is a separate question from whether they are numerous.")

# ---------------------------------------------------------------- method out of the counted text
# Verified 19 Oct from the Resource format on cell.com: "Detailed description of experimental
# procedures and analytical/statistical methods must be included in the STAR Methods section of
# the article", and STAR Methods text is excluded from the 7,000. So procedure description sitting
# in a legend or a body paragraph is a compliance defect as well as a cost. What stays is what
# Cell Press requires of a legend: the test, the measure, and every error bar defined.
def _sub(seq, i, old, new):
    assert old in seq[i], "method-relocation miss: %s" % old[:60]
    seq[i] = seq[i].replace(old, new)

# Figure 7 scope: the feature inventory and its three tool citations are method.
_sub(LEGEND[7], 0,
     "Twenty-two site features measured on 167,057 lysines at pLDDT ≥ 70 "
     "[[jumper2021,varadi2022]], including relative solvent accessibility [[tien2013]], "
     "predicted pKa [[olsson2011]] and pocket descriptors [[leguilloux2009]].",
     "Twenty-two site features measured on 167,057 lysines at pLDDT ≥ 70 "
     "[[jumper2021,varadi2022]], listed with their sources in the Methods.")

# Figure 3A: how the entrapment partner is built is method; the row counts are the panel.
_sub(LEGEND[3], 1,
     "(A) The design. Entrapment sequences are matched to the variant set on the properties that "
     "drive identification, so the comparison is like for like: ",
     "(A) The design, matched to the variant set on the properties that drive identification, so "
     "the comparison is like for like (Methods): ")

# Figure 4: the gating rule is already stated in the body, and why AlphaFold terms are excluded
# is method.
_sub(LEGEND[4], 0,
     " Every structural gate is applied to both groups of a contrast at once.", "")
_sub(LEGEND[4], 1,
     " [[tien2013,olsson2011]]; AlphaFold terms are excluded because confidence is the gating "
     "variable rather than a site property.",
     " [[tien2013,olsson2011]], with the AlphaFold terms excluded (Methods).")

# Figure 6D: the toolchain is method.
_sub(LEGEND[6], 2, " Descriptors were computed with RDKit [[rdkit]].", "")

# 2.5: why the split is by protein is a statistical-method justification.
_sub(SPINE["2.5"], 0,
     "and cross-validated by protein rather than by site, because lysines on one protein are "
     "correlated and a site-level split leaks information between training and test.",
     "and cross-validated by protein rather than by site, since a site-level split leaks "
     "information between training and test (Methods).")

# What moved, written out where the journal requires it.
METHODS_RELOCATED = (
"Site features and their sources. Site features were measured on AlphaFold Monomer v2.0 models, "
"and the twenty-two-feature set comprises relative solvent accessibility on the Tien scale "
"[[tien2013]], predicted pKa from PROPKA3 [[olsson2011]], pocket descriptors from fpocket "
"[[leguilloux2009]] and the neighborhood composition and shape terms listed in the data "
"deposit. AlphaFold per-residue confidence is excluded from any feature set used inside a "
"confidence-gated contrast, because there it is the gating variable rather than a property of "
"the site; every structural gate is applied to both groups of a contrast at once, since gating "
"one side would manufacture a difference rather than reveal one. Compound descriptors were "
"computed with RDKit [[rdkit]]. Entrapment partners were built from each variant peptide's own "
"residues, permuted with the C-terminal lysine or arginine held in place, so that a pair shares "
"an exact precursor mass and an identical length and only fragment ions can separate them. The "
"composite classifier was cross-validated by protein rather than by site, because lysines on one "
"protein are correlated and a site-level split leaks information between training and test.")

# ---------------------------------------------------------------- replicate statements
# Cell Press requires every legend to give the number of biological and technical replicates.
# This paper generated no data, so the honest statement names the replicate structure of the
# original studies and says plainly that the intervals here are computational resampling rather
# than replication. Added as a closing entry per figure.
# Figure 5 is the merged item: its statement is appended INTO the existing closing entry, because
# a new trailing panel-less entry would displace that closer into Extended Results.
# The per-figure closer that used to sit at the foot of every legend is gone at the PI's
# instruction, and so are the four interval sentences it had carried. Those were the only
# statement of interval construction left in four legends, which is why they were kept on the
# first pass, but Methods already specifies every one of them and in more detail: percentile
# ranges, 20 protein-held-out splits, 50 random half-splits, the exact C(15,6) = 5,005
# enumeration for Figure 4E, and 2,000 label shuffles for the breadth nulls. The captions no
# longer stand alone on interval construction. Figure 6 is unaffected: it states its own
# intervals inline, "95% percentile bootstrap confidence intervals, resampling proteins".

# ---------------------------------------------------------------- claims-audit fix, 20 Oct
# The Discussion's passenger paragraph said E to K is common "for reasons that carry no biological
# weight" and then rebutted that reading with a number, the 1.68 odds ratio. But the quantitative
# basis for the reading it was rebutting had left the main text: I cut it from the Discussion as a
# restatement of Results, and the Results version left too, because panel 1b is cut. The paragraph
# was quantifying one side of its own contrast and asserting the other. Restored here, where the
# claim now lives.
_old = ("Glutamate to lysine is the commonest substitution in cancer for reasons that carry no "
        "biological weight, and it would be reasonable to conclude from that alone that these are "
        "incidental residues.")
_new = ("Glutamate to lysine is the commonest substitution in cancer for reasons that carry no "
        "biological weight: glutamate is abundant, the dominant mutational process is C>T, and "
        "against a trinucleotide-context expectation the substitution is enriched only 1.12-fold, "
        "against more than threefold for each of the CpG hotspots. It would be reasonable to "
        "conclude from that alone that these are incidental residues.")
_hits = [i for i, p in enumerate(DISCUSSION) if _old in p]
assert len(_hits) == 1, "passenger paragraph not matched (%d)" % len(_hits)
DISCUSSION[_hits[0]] = DISCUSSION[_hits[0]].replace(_old, _new)

# ---------------------------------------------------------------- terminology pass, 20 Oct
# Engagement is the ligand-residue event; ligandability is the property of the residue. Of 77
# uses in the manuscript, 71 were already right. These six used the event word for the property,
# and two of them leaned on a coinage, "engageability", that the paper does not define anywhere.
# "Engagement is the outcome of a measurement" is deliberately NOT changed: there the event
# reading is the point, and the sentence exists to deny that it is a fixed property.
def _fix(seq, i, pairs):
    for old, new in pairs:
        assert old in seq[i], "terminology miss: %s" % old[:60]
        seq[i] = seq[i].replace(old, new)

_fix(SPINE["2.6"], 0, [
    ("If engagement is a property of the site rather than of the molecule",
     "If ligandability is a property of the site rather than of the molecule"),
    ("what follows is how well they predict engagement",
     "what follows is how well they predict ligandability"),
])
_fix(SPINE["2.5"], len(SPINE["2.5"]) - 1, [
    ("Engagement is a property of the lysine over the descriptors tested",
     "Ligandability is a property of the lysine over the descriptors tested"),
])
for _i, _p in enumerate(SPINE["2.8"]):
    if "the engagement rules measured in Figures 4 to 7" in _p:
        _fix(SPINE["2.8"], _i, [("the engagement rules measured in Figures 4 to 7",
                                 "the ligandability rules measured in Figures 4 to 7")])
    if "infer engageability from similarity" in _p:
        _fix(SPINE["2.8"], _i, [
            ("Features that predict engagement that poorly cannot be turned around and used to "
             "infer engageability from similarity.",
             "A feature set that predicts ligandability this weakly cannot support the reverse inference, from similarity to ligandability.")])

# ---------------------------------------------------------------- panels restored, 20 Oct
# Three citations were pointed at Extended Results when their panels were cut. The panels are
# back in the main figures, so the citations have to follow them.
def _recite(seq, old, new):
    hits = [i for i, p in enumerate(seq) if old in p]
    assert len(hits) == 1, "recite miss (%d): %s" % (len(hits), old[:60])
    seq[hits[0]] = seq[hits[0]].replace(old, new)

_recite(SPINE["2.7"], "No pocket feature and no shape feature differs (Extended Results).",
                      "No pocket feature and no shape feature differs (Figure 7B).")
_recite(SPINE["2.8"], "other substitutions at the same position perturb the environment as much "
                      "(Extended Results).",
                      "other substitutions at the same position perturb the environment as much "
                      "(Figure 8D).")
_recite(SPINE["2.8"], "and it returns zero spectra here (Extended Results)",
                      "and it returns zero spectra here (Figure 8E)")

# ---------------------------------------------------------------- citations for restored panels
# A restored panel needs calling out in the text, and the numbers stay in its legend.
def _append_to(seq, anchor, text):
    hits = [i for i, p in enumerate(seq) if anchor in p]
    assert len(hits) == 1, "anchor miss (%d): %s" % (len(hits), anchor[:60])
    seq[hits[0]] = seq[hits[0]].rstrip() + " " + text

# 2e: why a region is never covered at all
_append_to(SPINE["2.2"], "More than nine in ten of these sites are never covered",
    "Where the region is never covered, the shortest tryptic peptide spanning the site is the "
    "longer one, and the two distributions diverge in the tail rather than shifting as a whole "
    "(Figure 2E).")
# 2a: the four species, which is the setup for the whole detection problem
_append_to(SPINE["2.2"], "The mutant peptide has the opposite geometry",
    "Four peptide species can carry an acquired-lysine site and only two of them contain it "
    "(Figure 2A).")

# 3c, 3d, 3e: three properties of the identifications that refer to nothing in the search
SPINE["2.3"].insert(4,
    "Three properties of the identifications themselves agree with that reading, and none of "
    "them refers to the search engine or the database. The entrapment rate is spread across "
    "plexes rather than carried by one experiment (Figure 3C), the spectra score close to the "
    "background on peak count and signal-to-noise (Figure 3D), and the precursor mass errors sit "
    "in the same window as the background, including the first-isotope tail, which appears at the "
    "same rate in both (Figure 3E).")

# 4d: what the gate costs, which is a bound and belongs in the main text
SPINE["2.4"].insert(4,
    "Gating costs sites, and the price is worth stating. 30,805 labelled lysines fall to roughly "
    "17,000 at a confidence of 90 or above (Figure 4D), and some of the discarded sites are "
    "engaged, so this is a real cost rather than a free cleanup.")

# 4f and 7d: citations that pointed elsewhere while their panels were out
_recite(SPINE["2.4"], "rises monotonically across all ten tiers with no discontinuity where the "
                      "cut falls (Extended Results).",
                      "rises monotonically across all ten tiers with no discontinuity where the "
                      "cut falls (Figure 4F).")
_recite(SPINE["2.7"], "separate the reagent's contribution from the site's (Figure 7C, D).",
                      "separate the reagent's contribution from the site's (Figure 7C and "
                      "Figure 7D).")
