# -*- coding: utf-8 -*-
"""Supporting Information: the nine supplementary figures with publishable captions.

Captions are the publishable form of the entries in FIGURE_LEGENDS_all.md; internal [VERIFY]
notes, decision records and drafting commentary are stripped and live in the section log instead.
"""

SI_FIGS = [
 dict(no="S1", fig="figS1.png", parent="Figure 1", title="The substitution matrix in full.", cap=[
  "**(A) The complete 20 × 20 substitution matrix.** Direct counts from the union, normalised over the whole matrix rather than per row, so cells are comparable across rows.",
  "**(B) Against a trinucleotide-context expectation, E→K is barely enriched.** Observed over expected for the twelve commonest substitutions: E→K = 1.12, against three CpG hotspots that dominate the panel (R→H 3.59, R→C 3.38, R→Q 3.16).",
  "**(C) gnomAD frequency bands as shares.** Absent 52.2%, below 10⁻⁵ 32.8%, 10⁻⁵–10⁻⁴ 6.5%, 10⁻⁴–10⁻³ 1.1%, above 10⁻³ 0.7%, not assessed 6.7%. 91.1% of assessed sites are absent or below 10⁻⁵.",
  "**(D) Which residues are gained.** Sites gained per residue across the union. Lysine is fourth at 295,411; cysteine eleventh at 211,386.",
  "**(E) Net flux, gains minus losses.** Lysine gains 295,411 and loses 164,237, netting 131,174 and ranking fifth. Cysteine nets 144,683 and ranks second, because cysteine is rarely lost (66,703).",
  "**(F) Gains are one substitution class; losses are every other one.** Lysine gain is 85.5% C>T and 0% T>C; lysine loss is 0% C>T and spread across the other five. All missense sits between, at 51.1% C>T.",
  "**(G) The asymmetry is not inherited from the background spectrum.** The signature control, per lineage.",

  "**(H) Every residue, not only lysine, and cysteine on the identical calculation.** Somatic "
  "share of the missense burden divided by germline share, for all twenty residues, ordered by the "
  "ratio. Computed by step54_germline_by_residue.py on the same pass over dbNSFP 5.4a that produced "
  "Figure 1G: same file, same MANE restriction, same gnomAD v4.1 joint frequency field, same "
  "expected-diploid-count definition. It returns 471.1 acquired lysines per individual, which is "
  "Figure 1G\u2019s value, and that agreement is the check that generalising the calculation did "
  "not change it. Rust, lysine at 2.05, third of twenty; teal, cysteine at 1.76, fourth. The "
  "distribution runs from 0.43 (alanine) to 2.63 (phenylalanine) with a median of 1.07, and eight "
  "residues fall below 1. The panel exists because panel F shows all 32 lineages for the same "
  "reason: a single ratio cannot carry its own control. **It is not a measure of selection in "
  "cancer.** Germline shares are lowered by purifying selection, so residues whose introduction is "
  "most disruptive are depleted from inherited variation and score high here independently of any "
  "cancer-specific process; that is the likeliest explanation for phenylalanine and tyrosine, "
  "the only two residues above lysine. Selection is Figure 1D\u2019s odds ratio of 1.68; the context-corrected "
  "enrichment for E\u2192K is panel B\u2019s 1.12-fold.",
 ]),
 dict(no="S2", fig="figS2.png", parent="Figure 2", title="The non-lysine control and what limits detection.", cap=[
  "**(A) The two gates, both classes.** 1,902 variant sites that do not create a lysine, through the identical pipeline. Acquired lysines (n = 8,470): gate 1 7.3% covered, gate 2 4.4% recovered. Other variants (n = 1,902): 4.6% and 6.9%. Acquired lysines are covered more and recovered less, because a non-lysine variant shares peptide boundaries with the wild-type peptide that established coverage while an acquired lysine produces a different peptide.",
  "**(B) Peptide length, both classes.** Medians 13 against 14 residues, so the difference in panel A is not a length artefact.",
  "Missed cleavage in panels C and D is missed / (missed + cleaved) among IDENTIFIED peptides, as the field reports it. A cleavage site that produced no observable peptide appears in neither term, so this is not a digestion efficiency in the absolute sense.",
  "**(C) Cleavage after arginine fails about twice as often as cleavage after lysine.** Peptides carrying a missed cleavage at the site, split by whether the first basic residue is K or R, PSM-weighted and by distinct sequence: K 15.0% / 23.6%; R 29.4% / 37.0%.",
  "**(D) Acquired-lysine peptides are digested like any other.** Missed cleavage across all identified peptides: 23.8% for peptides carrying an acquired lysine against 21.6% background.",
  "*Panels C and D report missed / (missed + cleaved) among identified peptides, as the field reports it. A cleavage site that produced no observable peptide appears in neither term, so each is a lower bound on how often the protease fails rather than a digestion efficiency in the absolute sense.*",
  "**(E) Score margin.** Hyperscore minus nextscore for acquired-lysine PSMs: 53 of the 89 PSMs clear margin ≥ 10, and they fall on 30 of the 53 sites plotted.",
  "**(F) Spectrum quality.** Hyperscore densities, variant against background: medians 33.6 against 31.2.",
  "**(G) An annotated spectrum.** GAPDH E106K in WM2664, scan 72116: 635 peaks, base peak 8378, 41 of 76 theoretical ions matched within 0.6 Da, and 18 of 19 backbone bonds observed, with y1–y3 directly evidencing the C-terminal lysine.",
 ]),
 dict(no="S3", fig="figS3.png", parent="Figure 2", title="Four explanations that are not the answer.", cap=[
  "**(A) The detections are not carried by subclonal sites.** Of the 53 sites where the mutant peptide was detected, the 27 whose wild-type region was also covered are 11% clonal (allele fraction ≥ 0.7) and the 26 seen only as the mutant are 27%, against an 11% background.",
  "**(B) Protein abundance does not explain which regions are covered.** Abundance from the CCLE quantification rather than a proxy; 1,769 of 8,470 observations have a value and the panel is drawn on those. Cliff's δ = 0.31 between found and never-covered, a real shift, far too small to account for roughly 93% of regions never being observed.",
  "**(C) Allele fraction does not either.** Detected sites have a median allele fraction of 0.48 against 0.44 across all acquired-lysine sites.",
  "**(D) What predicts identification is incomplete digestion, not charge.** Standardised effects: internal K/R (missed cleavages) −0.621; basic centres, all −0.584; histidine (basic, but never a cleavage site) −0.220; C-terminal K/R −0.004. Within peptides having no missed cleavage the apparent charge effect falls from −0.58 to −0.09.",
  "**(E) The evidence behind the 53 sites, at full depth.** 53 sites rest on at least one spectrum, 12 on two or more, 4 on two or more distinct peptides, and none on more than one experiment. Bars are drawn at the search’s reported 1% false-discovery rate, which is the population this panel bounds; of these, 30 sites clear the score margin ≥ 10 of Figure 3B, with 24 on a single spectrum, six on two or more and two on two or more distinct peptides.",

  "**(F) The mutant peptide is short, and the missed-cleavage allowance absorbs it.** Peptide "
  "length for the four species that can carry an acquired-lysine site, over the 8,020 sites that "
  "map onto a reference sequence, at perfect digestion. Coloured, shorter than seven residues; "
  "pale, seven to fifty and therefore searchable; dark, longer than fifty. The window is not a "
  "convention chosen here; both searches enforced it: digest_min_length = 7 and "
  "digest_max_length = 50 in the MSFragger parameters used for both architectures, confirmed in "
  "the output, where the shortest peptide identified is seven residues and the longest fifty; and "
  "the original CCLE search excluded \u201cpeptides below a length of 7\u201d [[nusinow2020]]. "
  "Because a new lysine is a new cleavage site, the fully cleaved mutant peptide has a median "
  "length of six residues against thirteen for the wild-type peptide over the same site, and "
  "52.3% of it falls below the floor. Acylating the acquired lysine removes that cleavage site "
  "and restores the geometry: median fourteen residues, 16.8% below the floor, the least "
  "truncated of the four. **This panel is description, not explanation.** Both searches allowed "
  "two missed cleavages, so an upstream missed cleavage lengthens a too-short mutant peptide back "
  "into the window, and 8,066 of 8,470 observations (95.2%) still carry a searchable species "
  "ending at the acquired lysine. The direct test does not resolve whether the geometry costs "
  "detections: only 27 of the 617 covered observations lack a fully cleaved species, and at that "
  "population\u2019s 4.4% recovery rate one event is expected, so the comparison, 1 of 27 "
  "against 26 of 590, has no power rather than a null result. What the panel does "
  "establish is that for half this class the only searchable species carries a missed cleavage, "
  "and panel D shows internal missed cleavages are the strongest negative predictor of "
  "identification measured here. Computed by step52_mutant_peptide_length.py; counts from the "
  "same per-observation table as panels A\u2013e.",
 ]),
 dict(no="S4", fig="figS4.png", parent="Figure 3", title="The entrapment design in full.", cap=[
  "**(A) The four-way design.** 589,588 variant-pair rows read from 8 plexes; entrapment sets of 147,397 each. 100% of entrapment sequences share residues with their target.",
  "**(B) Length matching.** Peptide-length densities for the two sets, identical by construction.",
  "**(C) What the shuffle preserves.** C-terminal residue composition. Each entrapment peptide is its target's own residues permuted with the C-terminal K or R kept, so the pair shares an exact precursor mass and only fragment ions distinguish them. 100% of pairs share the same residues in a different order, so they have an identical precursor mass and only fragment ions discriminate. The R:K split survives because the C-terminal residue is held fixed.",
  "**(D) Boundary control.** Identifications across 42 experiments, all against those surviving margin ≥ 10: downstream target 58 → 25; boundary control (the same region shifted one residue, tryptic in neither protein) 2 → 0; entrapment for the downstream route 0 → 0.",
  "**(E) Per-plex FDP.** 42 plexes ordered by identifications, 0.00% to 9.52%, pooled 0.58%.",
  "**(F) Nominal against observed FDR.** On plex 2 (79,365 identifications) the nominal q-value from decoys reaches 0.60% where the observed proportion from entrapment reaches 0.99%, decoys underestimate by 1.7-fold.",
  "**(G) All four sets.** Score-margin densities. Entrapment peptides are identified, and score like a population that cannot be present should.",
  "**(H) A published estimator, on the whole search.** The Wen et al. estimator over the entire plex-2 search (79,365 identifications, not the variant subset), so it is not comparable with the 0.58% above.",
 ]),
 dict(no="S5", fig="figS5.png", parent="Figure 4", title="What the confidence gate does.", cap=[
  "**(A) The tiering null.** Five series. The steepest slope per tier is 0.013, against 0.083 for ligandability, 6.3 times smaller, and the same feature appears on both sides. **(B) Every feature, both gates.** 15 features. **(C) Leave-one-out robustness.** 15 features, the same set as Figure 4, each dropped in turn.",
  "**(D) Dimensionality.** Principal components explain 46.3, 28.8, 18.5, 6.4 and 0.0% of variance for reactivity and 48.9, 14.5, 7.8, 7.3, 4.9 and 4.6% for ligandability; reactivity is close to three-dimensional, ligandability is not. 5 charge features reduce to 3 axes and 14 ligandability features to 5. Each set is redundant and neither collapses to one axis, but they are different sets on different populations, so this is not 3 against 5.",
  "**(E) Mutual adjustment.** Each charge feature on its own against the same feature adjusted for the others. Salt bridge was anionic-neighbor density in another form, and does not survive the adjustment.",
  "**(F) What the gate removes, shown rather than asserted.** Two real engaged sites, each rendered whole-protein and as the lysine with everything within 8 Å, from the same AlphaFold models the features are measured on. PFKL K677, pLDDT 95, 43 compounds at 100% occupancy, kept. LYSMD4 K72 (UniProt Q6MZR2, which the deposited site table carries under its TrEMBL clone name DKFZp779J2370), pLDDT 28, 3 compounds at 95% occupancy, discarded. The same 8 Å sphere: in one of them there is no structure to measure.",
  "**(G) The gate is not a strength filter**, the control that excludes the objection that gating simply keeps the strongest sites.",
 ]),
 dict(no="S6", fig="figS6.png", parent="Figure 5", title="The ceiling, and the limit.", cap=[
  "**(A) The ceiling per dataset**, against how deeply each screen is split. Intervals are the 2.5th to 97.5th percentile across 50 random half-splits, the same construction as Figure 5A. The two shallow screens reach 0.881 to 0.938, but at engaged rates of 0.6 to 2.3%, where a few dozen positives make an AUC unstable. The fragment screen is the only screen deep enough to measure a ceiling at a realistic engaged rate, and its ceiling is still rising at 20 compounds per half, so 0.723 is a lower bound. The deeper the split, the less noisy the label.",
  "**(B) Single features against the composite.** Each of the 14 features alone, against the composite, held out by protein. The composite reaches 0.617 against 0.598 for the best single feature, distance to the protein centroid, a gain of +0.019. The 14 features are not independent.",
  "**(C) Every coefficient.** The 14 composite coefficients, averaged over the 20 protein-held-out splits, with each feature’s own AUC on the twin axis. An AUC below 0.5 is an inverted feature rather than a weak one: distance to the centroid reads 0.402, which is 0.598 with the sign flipped.",
  "**(D) Train against held-out.** Train 0.619, held-out 0.617, optimism +0.002. Intervals are the 2.5th to 97.5th percentile across the 20 protein-held-out splits, the same construction as Figure 5A and 5C. The model is not overfit. It is fitting all there is to fit, and that is the point.",
  "**(E) The breadth gradient, adjusted for testing depth**: 21 breadth features, the control for the objection that broadly engaged sites are simply the sites tested most often. Holding depth constant strengthens most of the gradient and flips 3 features outright.",
  "**(F) Two chemotypes agree.** 13 chemotype features. None of the 13 differs significantly between the two chemotypes, and the same sites are engaged by both (26 of 8,211 shared, 5.5 expected, OR 6.05), so one general composite is defensible. The three most discordant features are named in the panel.",
  "**(G) Across datasets, a limit, not a result.** This is the only cross-dataset evidence in the paper and the weakest evidence in it, and Figure 5F depends on it. The group carrying the cross-dataset claim is 25 sites in 21 dimensions, and 2 of 6 pairs fall short of their null over 2,000 label shuffles. Figure 5F repeats this test along the breadth gradient within one dataset, and this panel is what happens when it is asked across two. Gray is the null, teal clears it, rust does not.",
 ]),
 dict(no="S7", fig="figS7.png", parent="Figure 6", title="The confound, and what no descriptor can see.", cap=[
  "**(A) The provenance confound named.** Pooled across the three screens, 21 of 84 descriptor × site-feature pairs clear q < 0.05 over 82 compounds, and every one of them falls on molecular weight (8), volume (7) or logP (6). The three pure-shape descriptors, NPR1, NPR2 and rod-likeness, produce none, so what survives pooling is the size and lipophilicity axis rather than shape. Within the fragment screen alone the count is 0 of 84 (62 compounds), and across the within-dataset slopes 0 of 72 (173 compounds). This is the same analysis as Figure 6D, at the same inclusion rule of at least fifteen engaged sites per compound; the panel names which descriptors carry that figure’s pooled column, and Figure 6F shows those are the descriptors that separate the three screens.",

  "**(B) Per-compound slopes across occupancy.**",
  "**(C) Three threshold-free continuous measures** (mean occupancy 14 features, maximum occupancy 12, normalised sum 15), so the Figure 6 result does not depend on a 60% cut existing at all.",
  "**(D) The decile profile.** **(E) Where group sizes collapse**, which is why two of the five cuts in Figure 6E are not strictly monotonic.",
  "**(F) What no tested descriptor can see.** Charge: maximum |ρ| across tested descriptors is 0.193, against a median |ρ| of 0.552 for molecular weight. Chirality: NPR1 and NPR2 are ratios of principal moments of inertia and achiral by construction, so reflecting a conformer changes NPR1 by 0.000 while re-embedding an enantiomer changes it by up to 0.281, against a between-compound standard deviation of 0.123, so that figure is conformer sampling rather than chirality. This panel bounds Figure 6's claim rather than supporting it.",
 ]),
 dict(no="S8", fig="figS8.png", parent="Figure 7", title="Each probe alone, and the pooling reversal.", cap=[
  "**(A) Each probe against unlabelled lysines.** 22 features per probe; 8 significant for the anionic probe and 8 for the neutral one, so the two-of-twenty-two difference in Figure 7B is a statement about the contrast rather than about power.",
  "**(B) Enclosure is not a restatement of solvent accessibility.** Enclosure against relative solvent accessibility across 167,057 sites: Spearman ρ = −0.862, so the two share about 74% of their variance and roughly a quarter is enclosure's own.",
  "**(C) Pooling reverses the sign of the enclosure effect.** Pooled across sites, labelled lysines look less enclosed (d = −0.032); within protein they are more (d = +0.073). The within-protein value is the one Figure 7E reports.",
  "**(D) What the neutral-probe assumption bounds.** If the neutral probe is in fact partly ionised, the true reagent effect on cationic environment is larger, not smaller: +0.233 under the neutral assumption, running to +0.46 if the probe is substantially ionised. The assumption makes the reported effect a lower bound. **This is a sensitivity range across assumed ionisation states, not a confidence interval**: it carries no sampling statement, and the two endpoints are the neutral and substantially-ionised limits rather than percentiles.",
 ]),
 dict(no="S9", fig="figS9.png", parent="Figure 8", title="The anchor tested rather than defended.", cap=[
  "**(A) Distances, four ways**, and **(B) each arm against its own control.** Three legitimate anchorings (both sides at Cβ, both sides modelled at Nζ, and both at Cβ under a symmetric convention) give acquired-versus-engaged distances of 0.571, 0.535 and 0.583 against controls of 0.410, 0.331 and 0.407. The mismatched anchoring, union positions at Cβ against measured lysines at Nζ, gives 2.256 against 0.330, a 6.84-fold ratio. The anchor does not change the answer; the mismatch does. The mismatched arm is retained as a negative control because it is the analysis this project originally ran.",
  "**(C) Overlap, four ways.**",
  "**(D) One site per protein.** The correction for proteins contributing many lysines; the ratio shrinks but the direction holds.",
  "**(E) The built side chain, checked against real lysines.** All-trans rotamers were constructed at every position and compared with the amine positions in the measured set: median error 4.26 Å across 30,805 measured lysines, against a Cβ-to-Nζ separation of 6.4 Å. All-trans is the weakest possible construction, so this is an upper bound on the construction error rather than an estimate of it, and the reference amine positions are themselves predicted, since pLDDT is a backbone measure that does not certify side-chain placement.",
  "**(F) Overlap at every cut.** Four quantile cuts (50, 75, 90, 95), both legitimate arms, both comparison populations. Acquired lysines are at 42 / 69 / 82 / 88% inside the 50 / 75 / 90 / 95% engaged region under one anchoring and 45 / 67 / 83 / 89% under the other; labelled-but-unengaged lysines are at 47 / 73 / 89 / 95% and 50 / 78 / 92 / 96%. The reading does not depend on where the cut is drawn.",
 ]),
]

SI_NOTES_HEADER = (
 "Audited against the finalised figures on 6 September 2026. Five topics are now carried by a "
 "panel and are struck rather than written; two rest on an argument that was withdrawn; one needs "
 "rewriting. The remainder stand."
)

SI_NOTES = [
 ("N1",  "keep",     "Site validation: two-route agreement, the peptide correction and the marker rescue.", [
  "Every chemoproteomic site used in Figures 4 to 7 was assigned by two independent routes and the "
  "two assignments reconciled before any analysis ran.",
  "**Of the 45,552 usable sites, 44,736 agree between the two routes, 98.21%.** The remaining 816 "
  "were recovered rather than discarded: 470 were corrected from the peptide and 346 resolved by "
  "the marked-lysine method below. The verdict is retained as a column, so any analysis here can "
  "be restricted to the two-route agreement set alone.",
  "**The marker was validated before it was used.** On 4,000 sites whose position was already "
  "confirmed two other ways and which also carried a marked peptide, the marker pointed to the "
  "same residue in 3,999 of them, 99.97%, and disagreed once.",
  "**What the rescue recovered.** 529 rows were unusable before it: 442 ambiguous, 49 matching "
  "several positions and 38 with no positional evidence. The marker resolved 346 of them, 281 in "
  "CyMPSA and 65 in the fragment screen; 121 had neither a marker nor a structure and 62 sat on a peptide that "
  "was not unique. All 346 came from the ambiguous class, in which the peptide matched the sequence "
  "exactly once, carried more than one lysine, and the number the source reported named none of them.",
  "**How far the resolved positions moved.** No resolved position coincides with the number "
  "originally reported, and none can, for the reason just given. The signed displacement has a "
  "median of +1 residue, but 166 displacements are negative and 180 positive, so that figure is "
  "cancellation rather than agreement and is not the one to quote. The magnitude is: the resolved "
  "lysine sits a median of 34 residues from the reported number, interquartile range 12 to 111, "
  "maximum 3,234. As raw geography, the reported number falls inside the matched peptide but on a "
  "residue that is not a lysine in 36 rows, elsewhere in the same protein in 276, and outside the "
  "protein\u2019s modelled residue range altogether in 34; beta-actin, 375 residues long, carries a "
  "reported position of 913.",
  "**Most of that displacement is a numbering frame, not a wrong site.** Sorting the 346 by the "
  "strongest explanation available for each: in 64 the reported number lands exactly on the marked "
  "lysine in an annotated alternative isoform of the same protein, so the source coordinate was "
  "never wrong, only written in another frame; in 31 it falls inside the matched peptide on a "
  "non-lysine, a median of one residue away, and these are almost the whole of the 21 rows that sit "
  "within a single residue; in 27 it falls outside the protein entirely; the remaining 224 fall "
  "elsewhere in the protein with no isoform that accounts for them. Of the 52 proteins contributing "
  "more than one rescued row, 38, covering 127 rows, are displaced by one identical signed offset "
  "throughout, which is what a frame difference looks like and not what a randomly misassigned site "
  "looks like. Capicua is displaced by +909 residues in all three of its rows, and its shorter "
  "annotated isoform is exactly 909 residues shorter than the modelled one.",
  "**Why a wrong reported number cannot reach the analysis.** The marker never reads the reported "
  "number as a coordinate. It locates the peptide in the AlphaFold sequence, requires that the "
  "peptide occur exactly once, and takes the lysine at the marked offset; the reported number serves "
  "only to join back to the source row. Geometry is therefore computed at a residue fixed by the "
  "peptide, not by the number. These 346 rows are 0.76% of the 45,552 usable sites, they retain the "
  "verdict MARKER, and the engagement-rate model carries the verdict as a covariate, so every result "
  "can be reproduced with them excluded. step6c_audit_shift.py regenerates every figure in these two "
  "paragraphs from features_final.csv and UniProt.",
  "**A denominator warning, because both numbers are real.** 44,736 of the 45,735 raw feature rows "
  "is 97.82%, and 44,736 of the 45,552 usable sites is 98.21%. The usable-site denominator is the "
  "one this work quotes, since rows that never became usable sites are not sites the two routes "
  "could have agreed about."]),
 ("N2",  "keep",     "Occupancy scale construction across the three datasets.", [
  "Occupancy is the share of probe labelling that a compound competes away. The three screens "
  "measure it on two different instruments, and the scale was built so that one number means the "
  "same thing in all of them.",
  "For the TMT competition datasets, CyMPSA and the meroterpenoids, occupancy is 100 \u00d7 "
  "(1 \u2212 mean(compound channels) / mean(DMSO channels)). For the isoTOP-ABPP ratio dataset, "
  "the fragment screen, where R is the DMSO-to-compound ratio, occupancy is 100 \u00d7 (1 \u2212 1/R). The "
  "published criterion of the fragment screen of liganded at a ratio of 4 or more is therefore a 75% occupancy "
  "cut, not the 60% cut the other two screens use; a ratio of 2.5 is 60%. This work uses 60% "
  "throughout, and figure S7c shows the Figure 6 result does not depend on a 60% cut existing at "
  "all.",
  "**Three scale caveats travel with every row as explicit columns rather than as prose.** The "
  "number of conditions per site differs enormously between screens, a fragment-screen median of 41 "
  "against 5 to 15 for the others, so a per-site maximum taken over more conditions is larger for "
  "purely statistical reasons; n_conditions is a covariate in the models rather than an "
  "afterthought. isoTOP quantifies isotopologue pairs at MS1 while the TMTpro screens use isobaric "
  "reporters, which compress ratios toward one and therefore understate occupancy. A residual "
  "platform offset survives matching on n_conditions, so datasets are held out as "
  "cross-validation folds and raw occupancies are never pooled as though they were one "
  "measurement.",
  "Every row also carries a single-target-residue flag, set where the peptide contains exactly one "
  "residue of the probed type, so any analysis here can be repeated on the unambiguous subset as "
  "a sensitivity check."]),
 ("N3",  "keep",     "Screening depth by dataset, and why every ligandability estimate here is a lower bound.", [
  "A site that no compound engaged may not be ligandable, or may simply not have been tried hard "
  "enough. The two are separated by depth, and depth is wildly uneven across the three screens.",
  "**Depth per site**, the denominator behind every negative: CyMPSA and the meroterpenoids test "
  "every site at exactly 5 compounds. The fragment screen runs from 1 to 328, median 41, 90th percentile 214. "
  "A negative from 5 compounds and a negative from 41 are not the same statement, which is why "
  "n_conditions is a covariate throughout rather than a footnote.",
  "**Discovery curves**, compounds added in random order over 30 permutations, median reported. "
  "CyMPSA reaches 77 engaged sites at 1 compound, 307 at 5 and 721 at 15, with the last compound "
  "still adding 19 sites, 2.6% of the total. The meroterpenoids reach 27, 95 and 143 at 1, 5 and "
  "10 compounds, the last still adding 5, 3.5%. The fragment screen reaches 5, 51, 166, 1,611 and 2,465 at 1, "
  "5, 10, 162 and 324 compounds, and its last compound adds 1 site, 0.0%.",
  "**Two of the three screens are plainly unsaturated and the third is not evidence of "
  "saturation.** CyMPSA and the meroterpenoids are still gaining sites at their last compound, so "
  "for them not engaged means not yet engaged by what was tried. The fragment screen's curve has flattened, "
  "but it flattened over 328 compounds applied to a population whose per-site depth runs from 1 "
  "upward, so the flattening describes the compounds that were tried rather than the sites that "
  "remain. Every ligandability figure in this work is therefore a lower bound.",
  "**Per-compound reach**, for scale: CyMPSA has 15 compounds engaging a median of 64 sites "
  "(16 to 304), the meroterpenoids 10 compounds engaging a median of 26 (10 to 47), and the fragment screen "
  "328 compounds engaging a median of 6 (0 to 467)."]),
 ("N4",  "keep",     "Cross-screen concordance of site ligandability, and the compound-treatment "
                     "condition.", [
  "Compounds were applied to intact cells in the CyMPSA and meroterpenoid screens and to cell "
  "lysate in the fragment screen. In all three, the probe was applied to lysate after compound "
  "treatment, so the readout is identical and the screens differ only in where the compound "
  "encountered the proteome.",

  "No compound in these data was assayed in both conditions, so no controlled comparison of the "
  "condition is available. What can be measured is whether ligandability agrees between screens. "
  "Among sites assayed in both members of a pair, ligandability in one screen predicts "
  "ligandability in the other in every pair: the fragment and CyMPSA screens, 6,811 shared sites, odds ratio "
  "1.60, Fisher exact p = 5.5 \u00d7 10\u207b\u2075; the fragment and meroterpenoid screens, 3,712 shared "
  "sites, 3.12, p = 2.4 \u00d7 10\u207b\u2076; CyMPSA and the meroterpenoids, 8,211 shared sites, "
  "6.05, p = 2.0 \u00d7 10\u207b\u00b9\u00b9.",

  "The largest odds ratio belongs to the pair that shares a condition. That ordering is consistent "
  "with a condition effect and equally consistent with library similarity (CyMPSA and the "
  "meroterpenoids are 15 and 11 elaborated natural products, the third is 186 fragments), and "
  "across three screens the two explanations are collinear. No condition effect is claimed from the "
  "ordering, and none is excluded. The ligandable base rates carry the same warning: the fragment screen calls "
  "17.5% of its sites ligandable against 2.9% and 1.5%, which is fragment promiscuity and compound "
  "count rather than condition.",

  "These odds ratios are therefore a lower bound on site-intrinsic transfer, measured across "
  "condition, cell line and compound library at once, and not an estimate of any one of them.",

  "Lysate-phase measurement is the standard on which ligandability atlases are built. The "
  "chemoproteogenomic platform closest to this work profiles cell lysate and makes its druggability "
  "claim from it [[desai2024]]. DrugMap compiles cysteine ligandability across 416 cancer cell "
  "lines from lysate, and reports ligandable newly encoded cysteines among them [[takahashi2024]]. "
  "SLC-ABPP reports that the same workflow runs on intact cells, lysates and tissue with minimal "
  "change and without loss of quantitative cysteine depth [[kuljanin2021]].",

  "Computed by step50_condition_concordance.py; a site is ligandable in a screen at \u2265 60% "
  "occupancy by at least one of that screen\u2019s compounds.",
 ]),
 ("N5",  "keep",     "Pocket membership survives the confidence gate, and the pocket axis is not the charge axis.", [
  "Figure 4 shows every shape feature collapsing once the structure is required to be "
  "well-modelled while the charge features strengthen. Pocket membership is a shape feature by "
  "any ordinary reading, so it should collapse too. It does not.",
  "**Overall.** Of 30,805 probe-labelled lysines 44.7% sit in a pocket, median distance 5.53 \u00c5, "
  "against 39.8% and 6.47 \u00c5 for the 214,879 unlabelled ones.",
  "**Across the gate**, within protein, reported as effect over the feature\u2019s own standard "
  "deviation so that a shrinking effect can be checked against a shrinking spread. Being in a "
  "pocket runs 0.187, 0.111, 0.107 and 0.111 at gates of none, 50, 70 and 90, on a standard "
  "deviation flat at 0.49 to 0.50. Pocket volume **strengthens**, 0.047 to 0.064. Distance to the "
  "nearest pocket collapses, \u22120.316 to \u22120.066, but its standard deviation collapses with "
  "it, 8.65 to 4.14 \u00c5, so that decline is spread rather than signal. Pocket score and "
  "druggability score are flat and never significant.",
  "**Read against Figure 4**, where crowding reverses from +0.197 to −0.042, atom count collapses "
  "from +0.244 to +0.003 and polar contacts reverse from +0.183 to −0.002 on flat spreads while anionic "
  "contacts strengthen from −0.307 to −0.454. Pocket membership behaves like the charge features "
  "rather than like the shape features. Reactive lysines are at cavities, and probe capture is "
  "already enriching for sites a molecule could occupy.",
  "**The pocket axis is not the charge axis.** In a principal-component analysis over 167,035 "
  "lysines at pLDDT \u2265 70 and 17 features, PC1 carries 40.1% of the variance and is a burial "
  "axis: side-chain solvent accessibility \u22120.362, relative accessibility \u22120.350, "
  "N\u03b6 accessibility \u22120.338, atom count +0.369 and crowding +0.368. The pocket terms ride "
  "on it, distance \u22120.135 and in-pocket +0.118, while the charge terms are at zero, net charge "
  "\u22120.032 and predicted pKa \u22120.020. PC2 carries 13.9% and is the charge axis: predicted "
  "pKa \u22120.471, anionic contacts \u22120.406, salt bridges \u22120.384, with the pocket terms "
  "back at zero, distance \u22120.066 and volume +0.024. The two properties travel on different "
  "components, so the pocket result is not the charge result under another name."]),
 ("N6",  "withdraw", "Cavity volumes and the covalent-proximity argument. WITHDRAWN: the figures it rested on (41.6% of lysines within 5 Å of a cavity, 10.3% of cavities large enough) appear in no step report and do not reproduce, and step 19's in_pocket IS dist_nearest_pocket ≤ 5, so the two halves of the sentence were one variable. Figure 7E replaces the argument with enclosure, measured independently of both fpocket and solvent accessibility."),
 ("N7",  "keep",     "Promiscuity control: the breadth gradient sharpens when promiscuous compounds are excluded.", [
  "The objection is that sites engaged by many compounds are simply the sites a few very "
  "promiscuous compounds hit, so the breadth gradient would measure the compounds rather than the "
  "lysines.",
  "**The library is genuinely uneven.** Across 353 compounds the median engages 7 sites, the 90th "
  "percentile 58 and the most promiscuous 467. The ten most promiscuous compounds account for "
  "**34.0%** of all site-compound engagements. The test is therefore to delete them and rerun.",
  "**The gradient survives every deletion and steepens.** On the 3,513 lysines measured in all "
  "three datasets, mean distance from a never-engaged site rises monotonically with breadth in "
  "every run. With all compounds the tiers run 0.264, 0.627, 0.825, 1.129 and 2.616 against nulls "
  "of 0.177 to 0.546. Excluding the five most promiscuous compounds they run 0.247, 0.783, 0.835, "
  "1.378 and 2.853; excluding the ten, 0.247, 0.921, 1.030, 1.199 and 3.201. The most broadly "
  "engaged tier moves from 2.616 to 3.201 as the promiscuous compounds come out, which is the "
  "opposite of what the objection predicts.",
  "**What weakens instead is the singly-engaged tier**, from p = 0.058 with all compounds to "
  "p = 0.185 with ten removed, on a falling n. Sites engaged by exactly one compound are the tier "
  "promiscuity was most likely to manufacture, and they are the tier that does not survive. The "
  "gradient is carried by the sites that several distinct compounds reach."]),
 ("N8",  "keep",     "Compound descriptors and 3D shape, the nulls: the full descriptor list behind figure S7, which is what makes figure S7f's bound checkable.", [
  "Figure 6 is a null result over the descriptors that were tested, so the descriptors have to be "
  "listed and their blind spots measured rather than asserted.",
  "**The ten tested descriptors** are molecular weight, logP, topological polar surface area, "
  "hydrogen-bond donors, hydrogen-bond acceptors, ring count, aromatic ring count, rotatable "
  "bonds, heavy-atom count and the fraction of sp3 carbons. The shape terms are NPR1 and NPR2, "
  "ratios of principal moments of inertia.",
  "**Charge: the claim is a contrast, not an absence.** Over the 212 compounds with a parsable "
  "structure, net charge at pH 7.4 is 0 for 171, +1 for 35, +2 for 2, +3 for 2 and \u22121 for 2. "
  "Two descriptors are weakly but significantly correlated with that charge, ring count at "
  "|\u03c1| = 0.193 (p = 0.0047) and hydrogen-bond acceptors at 0.186 (p = 0.0066). The largest "
  "charge correlation in the set is therefore 0.193 against a **median size correlation of "
  "0.552**, a ratio of 0.35. The defensible statement is that the tested descriptors carry about a "
  "third as much charge information as size information, not that none of them encodes charge.",
  "**Handedness: the shape terms are blind to it by construction.** Reflecting a conformer changes "
  "NPR1 and NPR2 by exactly zero. Re-embedding an enantiomer changes them by up to 0.281 and "
  "0.251 across 25 chiral compounds, against a between-compound standard deviation of 0.121 and "
  "0.090. That difference is conformer sampling rather than chirality, and an earlier check "
  "measured it and would have been read as a chirality effect. The descriptors respond to "
  "handedness not at all, and their own sampling noise is already a third of the spread they "
  "resolve between compounds.",
  "This is what bounds Figure 6. Charge and chirality are the two axes the descriptor set sees "
  "least well, and Figure 7 and Figure 5F are where each is shown to matter."]),
 ("N9",  "withdraw", "Volume and shape against cavity. WITHDRAWN with N6, same cavity argument."),
 ("N10", "cut",      "Compound-count tiers. CUT: carried by figure S7e, which shows where the group sizes collapse."),
 ("N11", "cut",      "Occupancy threshold, 40–80%. CUT: carried by figure S7c, three threshold-free continuous measures, and by Figure 6C."),
 ("N12", "cut",      "Every feature against continuous occupancy. CUT: carried by figure S7c."),
 ("N13", "cut",      "Mutant-only detections and clonality. CUT: carried by figure S3a, which reports the 27/26 split and the 11% background."),
 ("N14", "cut",      "Spectrum quality. CUT: carried by Figure 3D and figure S2f."),
 ("N15", "rewrite",  "Crystallographic validation candidates, rewritten.", [
  "Figure 8D asks whether an acquired lysine moves the structure around it more than two crystals "
  "of the same wild-type protein differ from each other. Superposition is on residues beyond "
  "12 \u00c5 of the site and the reported displacement is over C\u03b1 atoms within 8 \u00c5 of "
  "it, so wild type against wild type is the noise floor rather than zero.",
  "**Six sites were defined in advance and all six were measured; none was excluded.** They are "
  "PTPN11 76, TP53 280, AKT1 17, PIK3CA 542, PIK3CA 545 and NRAS 61. Where a substitution series "
  "exists it was carried through: PTPN11 76 has E76K, E76A, E76D and E76Q against five wild-type "
  "structures, which is what makes the 3.53 of Figure 8D interpretable against 3.21 for other "
  "substitutions at the same position. An earlier version of this note described two validated "
  "sites from ten candidates with eight excluded; that is the state before 4 September 2026 and no "
  "longer holds.",
  "**The NRAS floor rests on a single pair, and deliberately so.** The wild-type set was restricted "
  "to one nucleotide state so that every pair is state-matched. Residue 61 sits in switch II, and "
  "an earlier run that mixed GDP and GTP structures inflated that floor to 4.25 \u00c5 and "
  "flattered the mutant. The panel is annotated with that n of 1 rather than leaving it to be "
  "inferred. Two of the six comparisons also cross a method boundary the others do not: the "
  "PIK3CA mutant structures 8GUA and 8GUD are cryo-EM while their wild-type partners are crystal "
  "structures.",
  "**Global rearrangement is reported separately from the local question.** Eleven pairs align "
  "locally but not globally, ten of them involving 6CRF, which is deposited under its own title as "
  "SHP2 E76K in the open conformation. E76K is an activating mutation that releases the "
  "autoinhibited state, so a large whole-protein difference alongside an ordinary local one is "
  "that effect. It argues for the figure rather than against it: the site itself is unremarkable "
  "even when the protein as a whole has moved.",
  "**The limit, stated because it is real.** Crystal structures differ for crystallographic "
  "reasons, among them space group, resolution, bound ligands and temperature, and this measure "
  "cannot separate those from the effect of the mutation. That is exactly why wild type against "
  "wild type is the comparator."]),
 ("N16", "keep",     "PgxSAVy answers a different question from the one this work asks.", [
  "PgxSAVy is the obvious tool to reach for here, and it does not do the job this work needs, "
  "which is worth saying plainly rather than leaving as an omission.",
  "**What it does.** PgxSAVy rescores individual variant peptide-spectrum matches and combines "
  "match quality, variant-event and search-result features into a variant ambiguity score, then "
  "labels each match confident, semi-confident or doubtful. On a reanalysis of roughly 2.8 million "
  "spectra it labelled 3,028 of 12,705 variant matches confident, 1,409 semi-confident and 8,268 "
  "doubtful.",
  "**Why that is not what is needed.** A per-match confidence label is not a false-discovery "
  "proportion for a class. The question here is what fraction of the variant identifications in a "
  "defined population are wrong, and a classifier answers it only if its calibration on this data "
  "is known, which nothing establishes. Burden-matched entrapment measures the proportion "
  "directly, on this search and this database, which is why Figure 3 uses it.",
  "**Two further mismatches.** PgxSAVy evaluates matches one at a time and does not conflate them "
  "into peptides, whereas the unit that matters here is the distinct peptide and the "
  "identification, kept apart throughout Figure 3. And its scoring was developed and calibrated on "
  "a different search pipeline, so its thresholds are not transferable to this one without a "
  "recalibration nobody has done.",
  "The two approaches are complementary rather than competing. PgxSAVy would improve confidence in "
  "individual variant matches; it would not tell this work what its realised error rate is."]),
 ("N17", "keep",     "Entrapment prior art: the survey this design answers to. Cited, not claimed.", [
  "Entrapment is not new and this work claims no novelty for it. What the recent literature "
  "establishes is that entrapment is easy to run and easy to run wrongly, which is why the design "
  "is spelled out in Figure 3 rather than named.",
  "**The survey finding.** Wen et al. identify three validation methods in common use: one "
  "invalid, one that provides only a lower bound, and one valid but under-powered. Reviewing "
  "published entrapment experiments, they report that **only three of the studies they tabulate "
  "use entrapment estimation correctly, and all three were data-dependent acquisition**. The two "
  "common errors are treating the lower-bound estimator as though it established error control, "
  "and using the sample-entrapment method.",
  "**Where this work sits.** The design here is burden-matched at a ratio of one, with each "
  "entrapment sequence built from its own target\u2019s residues so that the pair shares an exact "
  "precursor mass and an identical length. That is an upper-bound construction rather than a lower "
  "bound, it is applied to data-dependent acquisition, and the proportion is measured within the "
  "variant class rather than assumed from the whole search.",
  "**What is cited and what is not.** The survey and its taxonomy are cited. No claim is made here "
  "that burden-matched entrapment is a new method; the contribution is its correct application to "
  "a variant population, and note N18 states the bound that follows for the acquired-lysine class "
  "specifically."]),
 ("N18", "keep",     "Class-restricted entrapment: the false-discovery proportion measured on the acquired-lysine class itself rather than borrowed from the whole variant class.",
  ["Figure 3B measures the realised false-discovery proportion over every variant peptide in the "
   "search: 4.70% at margin 0, falling to 0.83% at margin ≥ 10. The main text then applies "
   "that number to acquired lysines. That is an inference from one population to a subset of it, "
   "and it is not automatic. A substitution to lysine introduces a tryptic cleavage site, so the "
   "peptide the search must match differs systematically from the peptide it would have looked "
   "for at the same position without the mutation; there is a mechanism by which the class could "
   "carry its own error rate.",

   "**Why the subset can be taken.** The entrapment design is pairwise. Each target variant "
   "peptide is paired at build time with one entrapment peptide of matched burden "
   "(peptide_pair_index), and the pair inherits the same class tag from the same variant record: "
   "KGAIN for a substitution that creates a lysine, KNOVEL for a frameshift or novel stretch that "
   "contains one, NOK otherwise. Restricting to KGAIN therefore preserves the r = 1 burden match "
   "inside the subset; it is not a post-hoc slice of an unmatched comparison. step51_kgain_"
   "entrapment.py verifies this rather than assuming it, and reports zero class-mismatched pairs "
   "across all 42 plexes, over 817,808 target/entrapment pairs.",

   "**Validation of the restriction.** Run without the class filter the script reproduces Figure "
   "3B exactly: 1,419 variant peptides against 70 entrapment at margin 0 (4.70%) and 832 "
   "against 7 at margin ≥ 10 (0.83%), which is the check that it is slicing the "
   "published population and not a different one. The unit throughout is one best spectrum per "
   "distinct peptide sequence per plex, and the estimator is e / (t + e), both as in Figure 3B.",

   "**Result.** At margin 0, 66 acquired-lysine peptides compete with five entrapment: 7.04% "
   "(Wilson 95% CI 3.05–15.45%), against 4.58% for the rest of the variant class (1,353 "
   "against 65). Fisher's exact test does not resolve the difference (odds ratio 1.58, p = 0.38). "
   "At margin ≥ 10 the class contributes zero entrapment peptides against 38 variant "
   "peptides, so its false-discovery proportion is bounded at 7.89% by the rule of three, or "
   "9.18% by the Wilson interval, rather than measured at zero. KNOVEL, which is a different "
   "class and is not the subject of this paper, is worse at every margin, on counts too small to "
   "interpret (eight peptides at margin 0).",

   "**What this settles and what it does not.** Seven entrapment peptides across 42 plexes at "
   "margin ≥ 10 is the entire supply of counted errors, so no single class of substitution "
   "can be resolved on its own to better than a few percent; the class-specific run cannot "
   "replace the whole-class estimate, and 0.83% remains the number quoted in the main text. What "
   "it establishes is the direction of the borrowing. Before the margin filter the acquired-"
   "lysine point estimate is higher than the rest of the class, not lower, so quoting the "
   "whole-class rate is not the conservative choice; after the margin filter the class "
   "contributes no counted errors at all, so the filtered set that the paper reports is not "
   "carrying a concealed excess. Both statements are reported here because only the second one "
   "flatters the result.",

   "Computed by step51_kgain_entrapment.py from the per-plex varpair definitions and the "
   "mutation-aware search PSM tables; outputs in kgain_margins.csv and kgain_entrapment_sweep.csv."]),
]

SI_FIGS.append(dict(
  no="S10", fig="figS10.png", parent="Figure 8",
  title="The acquired-lysine class across cell lines, and what a re-search of the three screens "
        "could find.",
  cap=[
  "Counts are unique (model, gene, protein change) throughout. DepMap 26Q1 carries duplicate rows: "
  "55,116 acquired-lysine rows collapse to 39,390 unique across 1,954 models, and 936 models are "
  "affected.",

  "**(A) Burden by lineage recovers the mutational processes of Figure 1E in a second dataset.** "
  "Acquired lysines per cell line, log scale, boxes at the quartiles with medians marked, lineages "
  "with at least 15 models, ordered by median. Skin leads at 35, then uterus 22, bowel 20 and lung "
  "19; eye sits second from last at 4. At subtype resolution the same contrast is sharper still: "
  "cutaneous melanoma 37 (n = 130) against ocular melanoma 3 (n = 15), which is the "
  "cutaneous-versus-uveal comparison of Figure 1F in cell lines rather than tumours. Metastatic "
  "models carry more than primary ones (median 13 against 9, n = 576 and 1,089, p = 7.7 "
  "\u00d7 10\u207b\u2078, Cliff\u2019s \u03b4 = \u22120.160). Rust, the three lineages the "
  "reanalysed screens sampled: cervix (HeLa), breast (MDA-MB-231, MDA-MB-436), lymphoid (Ramos).",

  "**(B) Every cell line in DepMap.** All 1,954 models carrying at least one missense lysine "
  "change, ranked by net gain; grey above the line, lysines gained, teal below, lysines lost. "
  "Symlog scale, linear within \u00b1 10. Across all models 39,390 lysines are gained and 25,531 "
  "lost, a ratio of 1.54, and 494 models (25.3%) are net-negative. The ratio is strongly "
  "lineage-dependent: 6.4 gained per lost in skin, 2.8 in bladder and soft tissue, 1.3 in "
  "oesophagus and stomach, the C>T asymmetry of figure S1f seen one cell line at a time, "
  "since lysine is gained 85.5% by C>T and never lost by it. Labelled: the three highest net "
  "gainers and the four cell lines the reanalysed screens were run in (rust). The extreme loss "
  "tail is a hypermutator losing 999 lysines.",

  "**(C) Four nested gates, and the lines that carry the class are the least characterised.** Of "
  "the 1,954 models, 1,630 have expression, 1,204 CRISPR and 372 CCLE proteomic data; only 308 "
  "have all three, and none of the ten highest net-gain lines does; MeWo has no CRISPR "
  "entry, SW684 none of the three. Bars are the twelve of those 308 with the most sites clearing "
  "every gate. Grey, unique acquired lysines; teal, gene at log\u2082(TPM+1) \u2265 1 in that "
  "line; charcoal, protein also quantified in that line\u2019s CCLE proteome; rust, gene also a "
  "CRISPR dependency at gene effect \u2264 \u22120.5. Aggregate over the 28 models examined: "
  "2,236 \u2192 1,035 expressed (46%) \u2192 686 also detected (31%) \u2192 132 also a "
  "dependency (6%).",

  "**(D) The recurrent driver-flagged acquired lysines are the sites Figure 8D already measured.** "
  "248 of 39,390 unique acquired lysines carry a hotspot or HessDriver flag (0.63%), in 210 of "
  "1,954 models. Ranked by the number of cell lines carrying each: NRAS Q61K 36, PIK3CA E545K 33, "
  "PIK3CA E542K 15, XPO1 E571K 9, TP53 E286K 7, TP53 E285K 6, AKT1 E17K 5, the sites in five or "
  "more models. Rust marks the four with a structure in Figure 8D: three of them are the three "
  "most recurrent, and the fourth, AKT1 E17K, ranks seventh. BRAF V600K, three models, "
  "is an acquired lysine on the most-drugged oncogene in melanoma. Driver-flagged sites clearing "
  "all three gates of panel C in the models drawn: MTOR E1799K (HEC59), PIK3CA E545K (MCF7), "
  "BRAF V600K (MDST8 and IGR1), PIK3CD E1021K (REH).",

  "**(E) Eleven acquired lysines in the four screen lines sit on a protein the matching screen "
  "probe-labels, and the measured lysine is never the acquired one.** A site is reachable when the "
  "screen run in that line labels at least one lysine on that protein, the precondition for "
  "measuring occupancy. Each bar is the highest occupancy at the best-engaged **probe-labelled** "
  "lysine on that protein in that screen, named on the axis with the acquired mutation in "
  "parentheses; the acquired lysine was outside the search space of all three screens and has no "
  "occupancy value. Labelled and acquired positions are a median of 271 residues apart (7 for "
  "FCHO1 to 2,391 for BIRC6), so a bar reports that the protein is reachable by a probe, not that "
  "the acquired site is. Filled bar, the one protein whose labelled lysine is engaged at "
  "\u2265 60% (SF3A3 K264, 72.8%). Cell lines: SF3A3, NARS2, SFXN4 HeLa; CAST Ramos; BIRC6 "
  "MDA-MB-436; the rest MDA-MB-231. TP53 R280K is one of the six sites measured in Figure 8D. "
  "Dashed line, the 60% engagement convention.",

  "Computed by figS10_v1.py from DepMap 26Q1 somatic mutations, expression, CRISPR gene effect and "
  "model annotation, the CCLE proteomes, and the per-site occupancy table behind Figures 4 to 8; "
  "every number above is printed by that script. Sites are matched to screens by cell line: HeLa "
  "to CyMPSA, MDA-MB-231 to the meroterpenoid and fragment screens, MDA-MB-436 to the meroterpenoid "
  "screen, Ramos to the fragment screen.",
  ]))

# ================================================================ Supplementary tables
# Table S1 added 7 September 2026. It is provenance for sites this paper measures, not a survey:
# every row is an acquired lysine that appears in Figure 8D, in figure S10d's recurrence ranking,
# or in figure S10c's three-gate set, and the third column records what OTHERS have reported the
# mutation does. The paper measures none of that itself, which is why the table is here and not in
# the Results, and why the Discussion cites it in one sentence.
#
# CITATION STATUS IS PART OF THE TABLE. A row carries a reference only where the primary paper has
# been verified against PubMed. Rows marked [ref pending] have a mechanism drawn from a trainee's
# literature summary whose primary citation has not yet been checked; they must not be submitted in
# that state. Recording the gap in the table itself is what stops an unverified string reaching the
# bibliography by inertia.
SI_TABLES = [dict(
  no="S1",
  title="The recurrent acquired lysines: how many patients carry them, and what they are "
        "reported to do.",
  head=("Site", "TCGA patients", "Where this paper measures it", "Reported functional consequence",
        "Evidence", "Primary source"),
  rows=[
   ("NRAS Q61K", "Fig 8d structure; 36 models (fig S10d)",
    "Slows GTP hydrolysis rather than abolishing it, and the effect belongs to the loss of "
    "Gln61 rather than to the lysine that replaces it. Every codon-61 substitution tested in "
    "HRAS, whose catalytic domain is identical to that of NRAS through residue 86, hydrolysed "
    "GTP eight- to tenfold more slowly than wild type, while transforming potency across the "
    "same set varied by more than a thousandfold, so the rate loss alone does not account for "
    "transformation [[der1986]]; the tenfold figure and the decoupling were reproduced "
    "independently, with crystal structures attributing the rate loss to displacement of the "
    "water molecule bridging Tyr32 to the \u03b3-phosphate [[buhrman2007]]. GAP-stimulated "
    "hydrolysis, by contrast, is not uniformly impaired at codon 61: NMR profiling found Q61L "
    "unexpectedly susceptible to it, in direct contrast with G12V [[smith2013]]. For NRAS "
    "itself, a crystal structure of the GTP-bound Q61K mutant resolves a ligand-induced pocket "
    "beside switch II [[gebregiworgis2024]]. In vivo, NRasQ61K on an INK4a-deficient background "
    "produces metastasising melanoma in mice [[ackermann2005]], and the allele drives aberrant "
    "melanocyte signalling, survival and invasiveness through Rac1 [[li2012]].",
    "Established driver",
    "[[der1986]], [[buhrman2007]], [[smith2013]], [[gebregiworgis2024]], [[ackermann2005]], "
    "[[li2012]]"),
   ("PIK3CA E545K", "Fig 8d structure; 33 models (fig S10d)",
    "Helical-domain hotspot. Oncogenic in vivo in the chick chorioallantoic membrane, with "
    "increased angiogenesis and Akt activation [[bader2006]]; separately, releases p85\u03b2 from "
    "p110\u03b1 to translocate to the nucleus, where it stabilises EZH1/EZH2 [[hao2022]].",
    "Established driver", "[[bader2006]], [[hao2022]]"),
   ("PIK3CA E542K", "Fig 8d structure; 15 models (fig S10d)",
    "As E545K: oncogenic in vivo [[bader2006]] and acting through the same p85\u03b2 "
    "nuclear-translocation route rather than canonical PI3K\u2013AKT alone [[hao2022]].",
    "Established driver", "[[bader2006]], [[hao2022]]"),
   ("XPO1 E571K", "9 models (fig S10d); not among the sites clearing all three gates in the models drawn (fig S10c)",
    "Alters cargo recognition at the NES-binding groove. Across 1,286 CLL patients, XPO1 E571 "
    "mutations associate with high-risk features; in mice, E571K with TCL1 accelerates "
    "leukemogenesis. Selinexor and eltanexor binding is unaffected [[walker2021]].",
    "Established driver", "[[walker2021]]"),
   ("AKT1 E17K", "Fig 8d structure; 5 models (fig S10d); the site of Fig 8e",
    "Pleckstrin-homology-domain mutation driving constitutive membrane recruitment and PI3K"
    "\u2013AKT\u2013mTOR activation [[carpten2007]]; the one acquired lysine engaged selectively "
    "by a covalent campaign [[craven2025]].", "Established driver", "[[carpten2007]], [[craven2025]]"),
   ("TP53 R280K", "Fig 8d structure",
    "Reported as a gain-of-function allele promoting NOX4-dependent reactive-oxygen production and "
    "focal-adhesion-kinase activation, and so migration and invasion [[boudreau2014]], rather than "
    "acting only through the loss of DNA binding that its crystal structure explains "
    "[[gomes2018]].", "Functionally validated", "[[boudreau2014]], [[gomes2018]]"),
   ("PTPN11 E76K", "Fig 8d structure",
    "Charge reversal at the N-SH2\u2013PTP autoinhibitory interface destabilises the closed "
    "conformation of SHP2. Gain-of-function Ptpn11 alleles including E76K cause cytokine "
    "hypersensitivity in myeloid progenitors and myeloproliferative disease in mice [[xu2013]].",
    "Established driver", "[[xu2013]]"),
   ("BRAF V600K", "3 models (fig S10d); clears all three gates in two lines (fig S10c)",
    "Constitutive BRAF\u2013MEK\u2013ERK activation. Present in 10\u201330% of BRAF-mutant "
    "melanoma [[rubinstein2010]]; clinically distinct from V600E, with shorter progression-free "
    "survival on BRAF-directed therapy [[menzies2012]] and lower ERK but higher PI3K\u2013AKT "
    "pathway expression [[piresdasilva2019]].", "Established driver",
    "[[rubinstein2010]], [[menzies2012]], [[piresdasilva2019]]"),
   ("TP53 E285K", "6 models (fig S10d)",
    "Disrupts a salt-bridge network stabilising the DNA-binding domain, abolishing transcriptional "
    "activity and mislocalising the protein to the cytosol, with cisplatin resistance reversible "
    "by ferroptosis induction [[dunsche2023]]; reported as a late secondary event in lobular "
    "breast-cancer evolution, where the mutant is temperature-sensitive and its reconstitution "
    "restores CDKN1A, MDM2 and PHLDA3 [[christgen2012]].", "Functionally validated",
    "[[dunsche2023]], [[christgen2012]]"),
   ("TP53 E286K", "7 models (fig S10d)",
    "Transactivation-defective at every response element tested. In the arrayed yeast assay "
    "that scored the 2,314 p53 missense substitutions reachable by a single point mutation "
    "against eight response elements, E286K retained 18.3% of wild-type transactivation at the "
    "p21WAF1 promoter and 21.8% at p53R2, and between 6.0% and 25.9% across the eight, so the "
    "substitution is strongly hypomorphic rather than silent [[kato2003]].",
    "Loss of function", "[[kato2003]]"),
   ("MTOR E1799K", "Clears all three gates (fig S10c)",
    "Hyperactivating. E1799K is one of 33 pathway-activating MTOR mutations curated from "
    "tumour sequencing, and one of ten expressed and assayed directly; all ten raised "
    "phosphorylation of the mTORC1 and mTORC2 substrates S6K1, 4EBP1 or AKT1 [[grabiner2014]]. "
    "Two cancer lines carry it natively, HEC59 (endometrium) and SNU349 (kidney), with "
    "rapamycin IC50 values of 0.49 and 0.19 nmol/L against 500 nmol/L for a line carrying a "
    "non-activating MTOR mutation, and rapamycin halted HEC59 xenograft growth outright "
    "[[grabiner2014]]. HEC59 is the model in which this site clears all three gates here.",
    "Gain of function", "[[grabiner2014]]"),
   ("XRCC5 E621K", "Clears all three gates (fig S10c)",
    "No mutation-specific functional study identified. XRCC5 (Ku80) partners Ku70, whose K351 is "
    "one of the liganded lysines in the fragment screen [[abbasov2021]].", "Not assessed", "\u2014"),
   ("PIK3CD E1021K", "Clears all three gates (fig S10c)",
    "The one acquired lysine here with a proven gain of function, established in the germline "
    "rather than in cancer: E1021K is the canonical allele of activated PI3K\u03b4 syndrome, and "
    "the substitution increases p110\u03b4 membrane association and kinase activity "
    "[[angulo2013]], with germline gain of function at this subunit driving T-cell senescence and "
    "immunodeficiency [[lucas2014]]. Reported somatically in Burkitt lymphoma [[penther2019]], but "
    "not characterised as a somatic driver. PIK3CD already carries approved drugs directed "
    "elsewhere on the protein.", "Gain of function shown in the germline",
    "[[angulo2013]], [[lucas2014]], [[penther2019]]"),
  ],
  foot=[
   "Patient counts are distinct TCGA participants carrying the substitution in the MC3 call set, "
   "against the 10,224 participants with at least one MC3 call; percentages are of that cohort. "
   "MC3 is already one of the three variant sources behind Figure 1, so this adds a patient "
   "denominator rather than a new dataset. No PASS-only filter is applied, for the reason given in "
   "the Methods: it would delete the ovarian and myeloid cohorts entirely; the PASS-only "
   "counts are in tcga_site_frequency.json and are within a few patients of these throughout. "
   "Column 3 is what this paper measures; columns 4 to 6 are what others report, and nothing in "
   "them is a result of this work. The table exists so that the sites named in the Discussion "
   "carry their independent evidence with them, and so that the sites with none are visible "
   "rather than absent.",
   "__CITATION_STATUS__",
   "**BRAF V600K returns zero patients, and the zero is a property of the call set rather than of "
   "the disease.** MC3 contains no V600K call at all; at codon 600 it carries 588 V600E, 40 "
   "V600M, five V600V, four V600G and two indels. V600K is the one site in this table that cannot "
   "be reached by a single nucleotide substitution: valine\u2019s codons are two changes from "
   "lysine\u2019s, while the other twelve sites are one. MC3 calls single-nucleotide variants and "
   "does not merge adjacent ones into a dinucleotide event, so a substitution requiring two "
   "adjacent changes is invisible to it by construction, which is why DepMap reports V600K "
   "in three cell lines and this cohort reports it in none. The same reasoning applies to every "
   "acquired lysine created from A, C, D, F, G, H, L, P, S, V, W or Y, all of which need at least "
   "two nucleotide changes; only E, I, M, N, Q, R and T reach lysine in one. Patient counts in "
   "this column are therefore a floor for the class as a whole, not only for this row.",
  ])]

# The citation-status sentence is measured from the rows, not typed, so it cannot drift from the
# table it describes, the same discipline the Highlights character counts use.
def _status(t):
    cited   = [r for r in t["rows"] if "[[" in r[5]]
    pending = [r for r in t["rows"] if "[ref pending]" in r[5]]
    none_   = [r for r in t["rows"] if r[5].strip() == "\u2014"]
    miss = len(pending) + len(none_)
    def _s(n, one, many): return one if n == 1 else many
    parts = ["**%d of the %d rows carry a verified primary citation; %d %s not.**"
             % (len(cited), len(t["rows"]), miss, _s(miss, "does", "do"))]
    if pending:
        parts.append("%d %s marked [ref pending] %s a mechanism summarised from the literature "
                     "whose primary source has not been checked against PubMed."
                     % (len(pending), _s(len(pending), "row is", "rows are"),
                        _s(len(pending), "carries", "carry")))
    if none_:
        parts.append("%d %s no mutation-specific study we could identify at all, and %s the "
                     "specific thing to resolve before this table goes out."
                     % (len(none_), _s(len(none_), "site has", "sites have"),
                        _s(len(none_), "that one uncited mechanism is", "those uncited mechanisms are")))
    parts.append("Every reference that does appear was verified against PubMed, with its PMID "
                 "recorded in the bibliography file.")
    return " ".join(parts)

# Patient counts are inserted here from step53_tcga_patient_frequency.py rather than typed into
# the rows, so the table cannot drift from the pass that produced them. TCGA MC3 is already one of
# the three variant sources behind Figure 1, so this adds a patient denominator without adding a
# data source. BRAF V600K carries a footnote rather than a count: see the table's second note.
_TCGA = {'NRAS Q61K': (68,), 'PIK3CA E545K': (263,), 'PIK3CA E542K': (169,), 'XPO1 E571K': (11,), 'AKT1 E17K': (54,), 'TP53 R280K': (16,), 'PTPN11 E76K': (2,), 'BRAF V600K': (0,), 'TP53 E285K': (24,), 'TP53 E286K': (16,), 'MTOR E1799K': (11,), 'XRCC5 E621K': (7,), 'PIK3CD E1021K': (1,)}
_COHORT = 10224
for _t in SI_TABLES:
    _rows = []
    for _r in _t["rows"]:
        _k = _r[0]
        _f = _TCGA.get(_k)
        if _f is None:
            _c = "not looked up"
        elif _f[0] == 0:
            _c = "0, see note"
        else:
            _c = "%d  (%.2f%%)" % (_f[0], 100.0 * _f[0] / _COHORT)
        _rows.append((_r[0], _c) + tuple(_r[1:]))
    _t["rows"] = _rows
    _t["foot"] = [_status(_t) if f == "__CITATION_STATUS__" else f for f in _t["foot"]]
