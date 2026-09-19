# -*- coding: utf-8 -*-
"""Supporting Information: the nine supplementary figures with publishable captions.

Captions are the publishable form of the entries in FIGURE_LEGENDS_all.md — internal [VERIFY]
notes, decision records and drafting commentary are stripped and live in the section log instead.
"""

SI_FIGS = [
 dict(no="S1", fig="figS1.png", parent="Figure 1", title="The substitution landscape in full.", cap=[
  "**(a) The complete 20 × 20 substitution matrix.** Direct counts from the union, normalised over the whole matrix rather than per row, so cells are comparable across rows.",
  "**(b) Against a trinucleotide-context expectation, E→K is barely enriched.** Observed over expected for the twelve commonest substitutions: E→K = 1.12, against three CpG hotspots that dominate the panel (R→H 3.59, R→C 3.38, R→Q 3.16).",
  "**(c) gnomAD frequency bands as shares.** Absent 52.2%, below 10⁻⁵ 32.8%, 10⁻⁵–10⁻⁴ 6.5%, 10⁻⁴–10⁻³ 1.1%, above 10⁻³ 0.7%, not assessed 6.7%.",
  "**(d) Which residues are gained.** Sites gained per residue across the union. Lysine is fourth at 295,411; cysteine eleventh at 211,386.",
  "**(e) Net flux — gains minus losses.** Lysine gains 295,411 and loses 164,237, netting 131,174 and ranking fifth. Cysteine nets 144,683 and ranks second, because cysteine is rarely lost (66,703).",
  "**(f) Gains are one substitution class; losses are every other one.** Lysine gain is 85.5% C>T and 0% T>C; lysine loss is 0% C>T and spread across the other five. All missense sits between, at 51.1% C>T.",
  "**(g) The asymmetry is not inherited from the background spectrum.** The signature control, per lineage.",

  "**(h) Every residue, not only lysine \u2014 and cysteine on the identical calculation.** Somatic "
  "share of the missense burden divided by germline share, for all twenty residues, ordered by the "
  "ratio. Computed by step54_germline_by_residue.py on the same pass over dbNSFP 5.4a that produced "
  "Figure 1g: same file, same MANE restriction, same gnomAD v4.1 joint frequency field, same "
  "expected-diploid-count definition. It returns 471.1 acquired lysines per individual, which is "
  "Figure 1g\u2019s value, and that agreement is the check that generalising the calculation did "
  "not change it. Rust, lysine at 2.05, third of twenty; teal, cysteine at 1.76, fourth. The "
  "distribution runs from 0.43 (alanine) to 2.63 (phenylalanine) with a median of 1.07, and eight "
  "residues fall below 1. The panel exists because panel f shows all 32 lineages for the same "
  "reason: a single ratio cannot carry its own control. **It is not a measure of selection in "
  "cancer.** Germline shares are lowered by purifying selection, so residues whose introduction is "
  "most disruptive are depleted from inherited variation and score high here independently of any "
  "cancer-specific process \u2014 the likeliest explanation for phenylalanine, tyrosine and "
  "tryptophan at the top. Selection is Figure 1d\u2019s odds ratio of 1.68; the context-corrected "
  "enrichment for E\u2192K is panel b\u2019s 1.12-fold.",
 ]),
 dict(no="S2", fig="figS2.png", parent="Figure 2", title="The non-lysine control and what limits detection.", cap=[
  "**(a) The two gates, both classes.** 1,902 variant sites that do not create a lysine, through the identical pipeline. Acquired lysines (n = 8,470): gate 1 7.3% covered, gate 2 4.4% recovered. Other variants (n = 1,902): 4.6% and 6.9%. Acquired lysines are covered more and recovered less, because a non-lysine variant shares peptide boundaries with the wild-type peptide that established coverage while an acquired lysine produces a different peptide.",
  "**(b) Peptide length, both classes.** Medians 13 against 14 residues, so the difference in panel a is not a length artefact.",
  "**(c) Cleavage after arginine fails about twice as often as cleavage after lysine.** Peptides carrying a missed cleavage at the site, split by whether the first basic residue is K or R, PSM-weighted and by distinct sequence: K 15.0% / 23.6%; R 29.4% / 37.0%.",
  "**(d) Acquired-lysine peptides are digested like any other.** Missed cleavage across all identified peptides: 23.8% for peptides carrying an acquired lysine against 21.6% background.",
  "*Panels c and d report missed / (missed + cleaved) among identified peptides, as the field reports it. A cleavage site that produced no observable peptide appears in neither term, so each is a lower bound on how often the protease fails rather than a digestion efficiency in the absolute sense.*",
  "**(e) Score margin.** Hyperscore minus nextscore for acquired-lysine PSMs: 53 of the 89 PSMs clear margin ≥ 10, and they fall on 30 of the 53 sites plotted.",
  "**(f) Spectrum quality.** Hyperscore densities, variant against background: medians 33.6 against 31.2.",
  "**(g) An annotated spectrum.** GAPDH E106K in WM2664, scan 72116: 635 peaks, base peak 8378, 41 of 76 theoretical ions matched within 0.6 Da, and 18 of 19 backbone bonds observed, with y1–y3 directly evidencing the C-terminal lysine.",
 ]),
 dict(no="S3", fig="figS3.png", parent="Figure 2", title="Four explanations that are not the answer.", cap=[
  "**(a) The detections are not carried by subclonal sites.** Of the 53 sites where the mutant peptide was detected, the 27 whose wild-type region was also covered are 11% clonal (allele fraction ≥ 0.7) and the 26 seen only as the mutant are 27%, against an 11% background.",
  "**(b) Protein abundance does not explain which regions are covered.** Abundance from the CCLE quantification rather than a proxy; 1,769 of 8,470 observations have a value and the panel is drawn on those. Cliff's δ = 0.31 between found and never-covered — a real shift, far too small to account for roughly 93% of regions never being observed.",
  "**(c) Allele fraction does not either.** Detected sites have a median allele fraction of 0.48 against 0.44 across all acquired-lysine sites.",
  "**(d) What predicts identification is incomplete digestion, not charge.** Standardised effects: internal K/R (missed cleavages) −0.621; basic centres, all −0.584; histidine — basic, but never a cleavage site — −0.220; C-terminal K/R −0.004. Within peptides having no missed cleavage the apparent charge effect falls from −0.58 to −0.09.",
  "**(e) The evidence behind the 53 sites, at full depth.** 53 sites rest on at least one spectrum, 12 on two or more, 4 on two or more distinct peptides, and none on more than one experiment. Bars are drawn at the search’s reported 1% false-discovery rate, which is the population this panel bounds; of these, 30 sites clear the score margin ≥ 10 of Figure 3b, with 24 on a single spectrum, six on two or more and two on two or more distinct peptides.",

  "**(f) The mutant peptide is short, and the missed-cleavage allowance absorbs it.** Peptide "
  "length for the four species that can carry an acquired-lysine site, over the 8,020 sites that "
  "map onto a reference sequence, at perfect digestion. Coloured, shorter than seven residues; "
  "pale, seven to fifty and therefore searchable; dark, longer than fifty. The window is not a "
  "convention chosen here \u2014 both searches enforced it: digest_min_length = 7 and "
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
  "population\u2019s 4.4% recovery rate one event is expected, so the comparison \u2014 1 of 27 "
  "against 26 of 590 \u2014 has no power rather than a null result. What the panel does "
  "establish is that for half this class the only searchable species carries a missed cleavage, "
  "and panel d shows internal missed cleavages are the strongest negative predictor of "
  "identification measured here. Computed by step52_mutant_peptide_length.py; counts from the "
  "same per-observation table as panels a\u2013e.",
 ]),
 dict(no="S4", fig="figS4.png", parent="Figure 3", title="The entrapment design in full.", cap=[
  "**(a) The four-way design.** 589,588 variant-pair rows read from 8 plexes; entrapment sets of 147,397 each. 100% of entrapment sequences share residues with their target.",
  "**(b) Length matching.** Peptide-length densities for the two sets, identical by construction.",
  "**(c) What the shuffle preserves.** C-terminal residue composition. Each entrapment peptide is its target's own residues permuted with the C-terminal K or R kept, so the pair shares an exact precursor mass and only fragment ions distinguish them.",
  "**(d) Boundary control.** Identifications across 42 experiments, all against those surviving margin ≥ 10: downstream target 58 → 25; boundary control — the same region shifted one residue, tryptic in neither protein — 2 → 0; entrapment for the downstream route 0 → 0.",
  "**(e) Per-plex FDP.** 42 plexes ordered by identifications, 0.00% to 9.52%, pooled 0.58%.",
  "**(f) Nominal against observed FDR.** On plex 2 (79,365 identifications) the nominal q-value from decoys reaches 0.60% where the observed proportion from entrapment reaches 0.99% — decoys underestimate by 1.7-fold.",
  "**(g) All four sets.** Score-margin densities. Entrapment peptides are identified, and score like a population that cannot be present should.",
  "**(h) A published estimator, on the whole search.** The Wen et al. estimator over the entire plex-2 search — 79,365 identifications, not the variant subset — so it is not comparable with the 0.58% above.",
 ]),
 dict(no="S5", fig="figS5.png", parent="Figure 4", title="What the confidence gate does.", cap=[
  "**(a) The tiering null.** Five series. **(b) Every feature, both gates.** 15 features. **(c) Leave-one-out robustness.** 15 features \u2014 the same set as Figure 4, each dropped in turn.",
  "**(d) Dimensionality.** Principal components explain 46.3, 28.8, 18.5, 6.4 and 0.0% of variance for reactivity and 48.9, 14.5, 7.8, 7.3, 4.9 and 4.6% for ligandability — reactivity is close to three-dimensional, ligandability is not.",
  "**(e) Mutual adjustment.**",
  "**(f) What the gate removes, shown rather than asserted.** Two real engaged sites, each rendered whole-protein and as the lysine with everything within 8 Å, from the same AlphaFold models the features are measured on. PFKL K677 — pLDDT 95, 43 compounds at 100% occupancy — kept. LYSMD4 K72 (UniProt Q6MZR2, which the deposited site table carries under its TrEMBL clone name DKFZp779J2370) — pLDDT 28, 3 compounds at 95% occupancy — discarded. The same 8 Å sphere: in one of them there is no structure to measure.",
  "**(g) The gate is not a strength filter** — the control that excludes the objection that gating simply keeps the strongest sites.",
 ]),
 dict(no="S6", fig="figS6.png", parent="Figure 5", title="The ceiling, and the limit.", cap=[
  "**(a) The ceiling per dataset**, against how deeply each screen is split. Intervals are the 2.5th to 97.5th percentile across 50 random half-splits, the same construction as Figure 5a.",
  "**(b) Single features against the composite.** Each of the 14 features alone, against the composite, held out by protein. The composite reaches 0.617 against 0.598 for the best single feature, distance to the protein centroid — a gain of +0.019.",
  "**(c) Every coefficient.** The 14 composite coefficients, averaged over the 20 protein-held-out splits, with each feature’s own AUC on the twin axis. An AUC below 0.5 is an inverted feature rather than a weak one: distance to the centroid reads 0.402, which is 0.598 with the sign flipped.",
  "**(d) Train against held-out.** Train 0.619, held-out 0.617, optimism +0.002. Intervals are the 2.5th to 97.5th percentile across the 20 protein-held-out splits, the same construction as Figure 5a and 5c.",
  "**(e) The breadth gradient, adjusted for testing depth** — 21 breadth features, the control for the objection that broadly engaged sites are simply the sites tested most often.",
  "**(f) Two chemotypes agree.** 13 chemotype features.",
  "**(g) Across datasets — a limit, not a result.** This is the only cross-dataset evidence in the paper and the weakest evidence in it.",
 ]),
 dict(no="S7", fig="figS7.png", parent="Figure 6", title="The confound, and what no descriptor can see.", cap=[
  "**(a) The provenance confound named.** Pooled across the three screens, 21 of 84 descriptor × site-feature pairs clear q < 0.05 over 82 compounds — and every one of them falls on molecular weight (8), volume (7) or logP (6). The three pure-shape descriptors, NPR1, NPR2 and rod-likeness, produce none, so what survives pooling is the size and lipophilicity axis rather than shape. Within NatChem alone the count is 0 of 84 (62 compounds), and across the within-dataset slopes 0 of 72 (173 compounds). This is the same analysis as Figure 6d, at the same inclusion rule of at least fifteen engaged sites per compound; the panel names which descriptors carry that figure’s pooled column, and Figure 6f shows those are the descriptors that separate the three screens.",

  "**(b) Per-compound slopes across occupancy.**",
  "**(c) Three threshold-free continuous measures** — mean occupancy 14 features, maximum occupancy 12, normalised sum 15 — so the Figure 6 result does not depend on a 60% cut existing at all.",
  "**(d) The decile profile.** **(e) Where group sizes collapse**, which is why two of the five cuts in Figure 6e are not strictly monotonic.",
  "**(f) What no tested descriptor can see.** Charge: maximum |ρ| across tested descriptors is 0.193, against a median |ρ| of 0.552 for molecular weight. Chirality: NPR1 and NPR2 are ratios of principal moments of inertia and achiral by construction, so reflecting a conformer changes NPR1 by 0.000 while re-embedding an enantiomer changes it by up to 0.281 — against a between-compound standard deviation of 0.123, so that figure is conformer sampling rather than chirality. This panel bounds Figure 6's claim rather than supporting it.",
 ]),
 dict(no="S8", fig="figS8.png", parent="Figure 7", title="Each probe alone, and the pooling reversal.", cap=[
  "**(a) Each probe against unlabelled lysines.** 22 features per probe; 8 significant for the anionic probe and 8 for the neutral one, so the two-of-twenty-two difference in Figure 7b is a statement about the contrast rather than about power.",
  "**(b) Enclosure is not a restatement of solvent accessibility.** Enclosure against relative solvent accessibility across 167,057 sites: Spearman ρ = −0.862, so the two share about 74% of their variance and roughly a quarter is enclosure's own.",
  "**(c) Pooling reverses the sign of the enclosure effect.** Pooled across sites, labelled lysines look less enclosed (d = −0.032); within protein they are more (d = +0.073). The within-protein value is the one Figure 7e reports.",
  "**(d) What the neutral-probe assumption bounds.** If the neutral probe is in fact partly ionised, the true reagent effect on cationic environment is larger, not smaller: +0.233 under the neutral assumption, running to +0.46 if the probe is substantially ionised. The assumption makes the reported effect a lower bound. **This is a sensitivity range across assumed ionisation states, not a confidence interval**: it carries no sampling statement, and the two endpoints are the neutral and substantially-ionised limits rather than percentiles.",
 ]),
 dict(no="S9", fig="figS9.png", parent="Figure 8", title="The anchor tested rather than defended.", cap=[
  "**(a) Distances, four ways**, and **(b) each arm against its own control.** Three legitimate anchorings — both sides at Cβ, both sides modelled at Nζ, and both at Cβ under a symmetric convention — give acquired-versus-engaged distances of 0.571, 0.535 and 0.583 against controls of 0.410, 0.331 and 0.407. The mismatched anchoring, union positions at Cβ against measured lysines at Nζ, gives 2.256 against 0.330, a 6.84-fold ratio. The anchor does not change the answer; the mismatch does. The mismatched arm is retained as a negative control because it is the analysis this project originally ran.",
  "**(c) Overlap, four ways.**",
  "**(d) One site per protein.** The correction for proteins contributing many lysines; the ratio shrinks but the direction holds.",
  "**(e) The built side chain, checked against real lysines.** All-trans rotamers were constructed at every position and compared with the amine positions in the measured set: median error 4.26 Å across 30,805 measured lysines, against a Cβ-to-Nζ separation of 6.4 Å. All-trans is the weakest possible construction, so this is an upper bound on the construction error rather than an estimate of it — and the reference amine positions are themselves predicted, since pLDDT is a backbone measure that does not certify side-chain placement.",
  "**(f) Overlap at every cut.** Four quantile cuts (50, 75, 90, 95), both legitimate arms, both comparison populations. Acquired lysines sit at 42 / 69 / 82 / 88% inside the 50 / 75 / 90 / 95% engaged region under one anchoring and 45 / 67 / 83 / 89% under the other; labelled-but-unengaged lysines sit at 47 / 73 / 89 / 95% and 50 / 78 / 92 / 96%. The reading does not depend on where the cut is drawn.",
 ]),
]

SI_NOTES_HEADER = (
 "Audited against the finalised figures on 6 September 2026. Four topics are now carried by a "
 "panel and are struck rather than written; two rest on an argument that was withdrawn; one needs "
 "rewriting. The remainder stand."
)

SI_NOTES = [
 ("N1",  "keep",     "Site validation — 97.82% two-route agreement, 44,736 of 45,552 usable sites, and the marker rescue. Referenced from Methods and from every figure that names a site set."),
 ("N2",  "keep",     "Occupancy scale construction across the three datasets."),
 ("N3",  "keep",     "Screening depth by dataset — the basis for the statement that no screen has saturated, so every ligandability estimate is a lower bound."),
 ("N4",  "keep",     "Cross-screen concordance of site ligandability, and the compound-treatment "
                     "condition.", [
  "Compounds were applied to intact cells in the CyMPSA and meroterpenoid screens and to cell "
  "lysate in the NatChem screen. In all three, the probe was applied to lysate after compound "
  "treatment, so the readout is identical and the screens differ only in where the compound "
  "encountered the proteome.",

  "No compound in these data was assayed in both conditions, so no controlled comparison of the "
  "condition is available. What can be measured is whether ligandability agrees between screens. "
  "Among sites assayed in both members of a pair, ligandability in one screen predicts "
  "ligandability in the other in every pair: NatChem and CyMPSA, 6,811 shared sites, odds ratio "
  "1.60, Fisher exact p = 5.5 \u00d7 10\u207b\u2075; NatChem and the meroterpenoids, 3,712 shared "
  "sites, 3.12, p = 2.4 \u00d7 10\u207b\u2076; CyMPSA and the meroterpenoids, 8,211 shared sites, "
  "6.05, p = 2.0 \u00d7 10\u207b\u00b9\u00b9.",

  "The largest odds ratio belongs to the pair that shares a condition. That ordering is consistent "
  "with a condition effect and equally consistent with library similarity \u2014 CyMPSA and the "
  "meroterpenoids are 15 and 11 elaborated natural products, NatChem is 186 fragments \u2014 and "
  "across three screens the two explanations are collinear. No condition effect is claimed from the "
  "ordering, and none is excluded. The ligandable base rates carry the same warning: NatChem calls "
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
 ("N5",  "keep",     "Pocket membership survives the confidence gate, and the pocket axis is independent of charge."),
 ("N6",  "withdraw", "Cavity volumes and the covalent-proximity argument. WITHDRAWN — the figures it rested on (41.6% of lysines within 5 Å of a cavity, 10.3% of cavities large enough) appear in no step report and do not reproduce, and step 19's in_pocket IS dist_nearest_pocket ≤ 5, so the two halves of the sentence were one variable. Figure 7e replaces the argument with enclosure, measured independently of both fpocket and solvent accessibility."),
 ("N7",  "keep",     "Promiscuity control — the breadth gradient sharpens when promiscuous compounds are excluded."),
 ("N8",  "keep",     "Compound descriptors and 3D shape, the nulls — the full descriptor list behind figure S7, which is what makes figure S7f's bound checkable."),
 ("N9",  "withdraw", "Volume and shape against cavity. WITHDRAWN with N6 — same cavity argument."),
 ("N10", "cut",      "Compound-count tiers. CUT — carried by figure S7e, which shows where the group sizes collapse."),
 ("N11", "cut",      "Occupancy threshold, 40–80%. CUT — carried by figure S7c, three threshold-free continuous measures, and by Figure 6c."),
 ("N12", "cut",      "Every feature against continuous occupancy. CUT — carried by figure S7c."),
 ("N13", "cut",      "Mutant-only detections and clonality. CUT — carried by figure S3a, which reports the 27/26 split and the 11% background."),
 ("N14", "cut",      "Spectrum quality. CUT — carried by Figure 3d and figure S2f."),
 ("N15", "rewrite",  "Crystallographic validation candidates. REWRITTEN: step31_crystal_validation.py defines six sites — PTPN11 76, TP53 280, AKT1 17, PIK3CA 542, PIK3CA 545 and NRAS 61 — and all six were measured; none was excluded. The note must also record that the NRAS wild-type set is restricted to a single nucleotide state so every pair is state-matched, which is why its floor rests on one pair, and that the PI3Kα mutant structures 8GUA and 8GUD are cryo-EM. The earlier version of this note described two validated sites from ten candidates with eight excluded; that is the state before 4 September 2026 and no longer holds."),
 ("N16", "keep",     "PgxSAVy non-applicability."),
 ("N17", "keep",     "Entrapment prior art — the survey finding that 14 of 17 studies used foreign entrapment and only 3 of 17 a valid design. Cited, not claimed."),
 ("N18", "keep",     "Class-restricted entrapment — the false-discovery proportion measured on the acquired-lysine class itself rather than borrowed from the whole variant class.",
  ["Figure 3b measures the realised false-discovery proportion over every variant peptide in the "
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
   "3b exactly — 1,419 variant peptides against 70 entrapment at margin 0 (4.70%) and 832 "
   "against 7 at margin ≥ 10 (0.83%) — which is the check that it is slicing the "
   "published population and not a different one. The unit throughout is one best spectrum per "
   "distinct peptide sequence per plex, and the estimator is e / (t + e), both as in Figure 3b.",

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

  "**(a) Burden by lineage recovers the mutational processes of Figure 1e in a second dataset.** "
  "Acquired lysines per cell line, log scale, boxes at the quartiles with medians marked, lineages "
  "with at least 15 models, ordered by median. Skin leads at 35, then uterus 22, bowel 20 and lung "
  "19; eye sits second from last at 4. At subtype resolution the same contrast is sharper still: "
  "cutaneous melanoma 37 (n = 130) against ocular melanoma 3 (n = 15), which is the "
  "cutaneous-versus-uveal comparison of Figure 1e in cell lines rather than tumours. Metastatic "
  "models carry more than primary ones (median 13 against 9, n = 576 and 1,089, p = 7.7 "
  "\u00d7 10\u207b\u2078, Cliff\u2019s \u03b4 = \u22120.160). Rust, the three lineages the "
  "reanalysed screens sampled: cervix (HeLa), breast (MDA-MB-231, MDA-MB-436), lymphoid (Ramos).",

  "**(b) Every cell line in DepMap.** All 1,954 models carrying at least one missense lysine "
  "change, ranked by net gain; grey above the line, lysines gained, teal below, lysines lost. "
  "Symlog scale, linear within \u00b1 10. Across all models 39,390 lysines are gained and 25,531 "
  "lost, a ratio of 1.54, and 494 models (25.3%) are net-negative. The ratio is strongly "
  "lineage-dependent: 6.4 gained per lost in skin, 2.8 in bladder and soft tissue, 1.3 in "
  "oesophagus and stomach \u2014 the C>T asymmetry of figure S1f seen one cell line at a time, "
  "since lysine is gained 85.5% by C>T and never lost by it. Labelled: the three highest net "
  "gainers and the four cell lines the reanalysed screens were run in (rust). The extreme loss "
  "tail is a hypermutator losing 999 lysines.",

  "**(c) Four nested gates, and the lines that carry the class are the least characterised.** Of "
  "the 1,954 models, 1,630 have expression, 1,204 CRISPR and 372 CCLE proteomic data; only 308 "
  "have all three, and none of the ten highest net-gain lines does \u2014 MeWo has no CRISPR "
  "entry, SW684 none of the three. Bars are the twelve of those 308 with the most sites clearing "
  "every gate. Grey, unique acquired lysines; teal, gene at log\u2082(TPM+1) \u2265 1 in that "
  "line; charcoal, protein also quantified in that line\u2019s CCLE proteome; rust, gene also a "
  "CRISPR dependency at gene effect \u2264 \u22120.5. Aggregate over the 28 models examined: "
  "2,236 \u2192 1,035 expressed (46%) \u2192 686 also detected (31%) \u2192 132 also a "
  "dependency (6%).",

  "**(d) The recurrent driver-flagged acquired lysines are the sites Figure 8d already measured.** "
  "248 of 39,390 unique acquired lysines carry a hotspot or HessDriver flag (0.63%), in 210 of "
  "1,954 models. Ranked by the number of cell lines carrying each: NRAS Q61K 36, PIK3CA E545K 33, "
  "PIK3CA E542K 15, XPO1 E571K 9, TP53 E286K 7, TP53 E285K 6, AKT1 E17K 5. Rust marks the four "
  "with a structure in Figure 8d, and they are the four most recurrent. BRAF V600K, three models, "
  "is an acquired lysine on the most-drugged oncogene in melanoma. Driver-flagged sites clearing "
  "all three gates of panel b in the models drawn: MTOR E1799K (HEC59), PIK3CA E545K (MCF7), "
  "BRAF V600K (MDST8 and IGR1), PIK3CD E1021K (REH).",

  "**(e) Eleven acquired lysines in the four screen lines sit on a protein the matching screen "
  "probe-labels, and the measured lysine is never the acquired one.** A site is reachable when the "
  "screen run in that line labels at least one lysine on that protein, the precondition for "
  "measuring occupancy. Each bar is the highest occupancy at the best-engaged **probe-labelled** "
  "lysine on that protein in that screen, named on the axis with the acquired mutation in "
  "parentheses; the acquired lysine was outside the search space of all three screens and has no "
  "occupancy value. Labelled and acquired positions are a median of 271 residues apart (7 for "
  "FCHO1 to 2,391 for BIRC6), so a bar reports that the protein is reachable by a probe, not that "
  "the acquired site is. Filled bar, the one protein whose labelled lysine is engaged at "
  "\u2265 60% (SF3A3 K264, 72.8%). Cell lines: SF3A3, NARS2, SFXN4 HeLa; CAST Ramos; BIRC6 "
  "MDA-MB-436; the rest MDA-MB-231. TP53 R280K is one of the six sites measured in Figure 8d. "
  "Dashed line, the 60% engagement convention.",

  "Computed by figS10_v1.py from DepMap 26Q1 somatic mutations, expression, CRISPR gene effect and "
  "model annotation, the CCLE proteomes, and the per-site occupancy table behind Figures 4 to 8; "
  "every number above is printed by that script. Sites are matched to screens by cell line: HeLa "
  "to CyMPSA, MDA-MB-231 to the meroterpenoid and NatChem screens, MDA-MB-436 to the meroterpenoid "
  "screen, Ramos to NatChem.",
  ]))

# ================================================================ Supplementary tables
# Table S1 added 7 September 2026. It is provenance for sites this paper measures, not a survey:
# every row is an acquired lysine that appears in Figure 8d, in figure S10d's recurrence ranking,
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
    "Impairs intrinsic and GAP-mediated GTP hydrolysis, locking NRAS in the GTP-bound state; the "
    "commonest NRAS variant in melanoma.", "Established driver", "[ref pending]"),
   ("PIK3CA E545K", "Fig 8d structure; 33 models (fig S10d)",
    "Helical-domain hotspot. Oncogenic in vivo in the chick chorioallantoic membrane, with "
    "increased angiogenesis and Akt activation [[bader2006]]; separately, releases p85\u03b2 from "
    "p110\u03b1 to translocate to the nucleus, where it stabilises EZH1/EZH2 [[hao2022]].",
    "Established driver", "[[bader2006]], [[hao2022]]"),
   ("PIK3CA E542K", "Fig 8d structure; 15 models (fig S10d)",
    "As E545K: oncogenic in vivo [[bader2006]] and acting through the same p85\u03b2 "
    "nuclear-translocation route rather than canonical PI3K\u2013AKT alone [[hao2022]].",
    "Established driver", "[[bader2006]], [[hao2022]]"),
   ("XPO1 E571K", "9 models (fig S10d); clears all three gates (fig S10c)",
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
    "focal-adhesion-kinase activation, and so migration and invasion, rather than acting only "
    "through loss of DNA binding.", "Functionally validated", "[ref pending]"),
   ("PTPN11 E76K", "Fig 8d structure",
    "Charge reversal at the N-SH2\u2013PTP autoinhibitory interface destabilises the closed "
    "conformation of SHP2. Gain-of-function Ptpn11 alleles including E76K cause cytokine "
    "hypersensitivity in myeloid progenitors and myeloproliferative disease in mice [[xu2013]].",
    "Established driver", "[[xu2013]]"),
   ("BRAF V600K", "3 models (fig S10d); clears all three gates in two lines (fig S10c)",
    "Constitutive BRAF\u2013MEK\u2013ERK activation; a recurrent melanoma allele responsive to "
    "BRAF-targeted therapy, though with shorter responses and greater PI3K\u2013AKT use than "
    "V600E.", "Established driver", "[ref pending]"),
   ("TP53 E285K", "6 models (fig S10d)",
    "Disrupts a salt-bridge network stabilising the DNA-binding domain, with loss of CDKN1A, MDM2 "
    "and PHLDA3 transcription and cisplatin resistance; reported as a late secondary event in "
    "lobular breast-cancer evolution.", "Functionally validated", "[ref pending]"),
   ("TP53 E286K", "7 models (fig S10d)",
    "No mutation-specific functional study identified.", "Not assessed", "\u2014"),
   ("MTOR E1799K", "Clears all three gates (fig S10c)",
    "Not characterised at this residue. Other mTOR lysine-gain substitutions (T1977K, E2014K, "
    "Q2223K, E2288K, E2419K) are reported to hyperactivate mTORC1/mTORC2; E1799K is not among them.",
    "Not assessed at this site", "\u2014"),
   ("XRCC5 E621K", "Clears all three gates (fig S10c)",
    "No mutation-specific functional study identified. XRCC5 (Ku80) partners Ku70, whose K351 is "
    "one of the liganded lysines in the NatChem screen [[abbasov2021]].", "Not assessed", "\u2014"),
   ("PIK3CD E1021K", "Clears all three gates (fig S10c)",
    "No mutation-specific functional study identified in cancer. PIK3CD already carries approved "
    "drugs directed elsewhere on the protein.", "Not assessed", "\u2014"),
  ],
  foot=[
   "Patient counts are distinct TCGA participants carrying the substitution in the MC3 call set, "
   "against the 10,224 participants with at least one MC3 call; percentages are of that cohort. "
   "MC3 is already one of the three variant sources behind Figure 1, so this adds a patient "
   "denominator rather than a new dataset. No PASS-only filter is applied, for the reason given in "
   "the Methods \u2014 it would delete the ovarian and myeloid cohorts entirely; the PASS-only "
   "counts are in tcga_site_frequency.json and are within a few patients of these throughout. "
   "Column 3 is what this paper measures; columns 4 to 6 are what others report, and nothing in "
   "them is a result of this work. The table exists so that the sites named in the Discussion "
   "carry their independent evidence with them, and so that the sites with none are visible "
   "rather than absent.",
   "__CITATION_STATUS__",
   "**BRAF V600K returns zero patients, and the zero is a property of the call set rather than of "
   "the disease.** MC3 contains no V600K call at all \u2014 at codon 600 it carries 588 V600E, 40 "
   "V600M, five V600V, four V600G and two indels. V600K is the one site in this table that cannot "
   "be reached by a single nucleotide substitution: valine\u2019s codons are two changes from "
   "lysine\u2019s, while the other twelve sites are one. MC3 calls single-nucleotide variants and "
   "does not merge adjacent ones into a dinucleotide event, so a substitution requiring two "
   "adjacent changes is invisible to it by construction \u2014 which is why DepMap reports V600K "
   "in three cell lines and this cohort reports it in none. The same reasoning applies to every "
   "acquired lysine created from A, C, D, F, G, H, L, P, S, V, W or Y, all of which need at least "
   "two nucleotide changes; only E, I, M, N, Q, R and T reach lysine in one. Patient counts in "
   "this column are therefore a floor for the class as a whole, not only for this row.",
  ])]

# The citation-status sentence is measured from the rows, not typed, so it cannot drift from the
# table it describes — the same discipline the Highlights character counts use.
def _status(t):
    cited   = [r for r in t["rows"] if "[[" in r[5]]
    pending = [r for r in t["rows"] if "[ref pending]" in r[5]]
    none_   = [r for r in t["rows"] if r[5].strip() == "\u2014"]
    return ("**%d of the %d rows carry a verified primary citation; %d do not.** %d rows marked "
            "[ref pending] carry a mechanism summarised from the literature whose primary source "
            "has not been checked against PubMed, and %d sites have no mutation-specific study we "
            "could identify at all. Neither state is submittable as it stands, and the four "
            "uncited mechanisms are the specific thing to verify before this table goes out. Every "
            "reference that does appear was verified against PubMed on 7 September 2026, with its "
            "PMID recorded in the bibliography file."
            % (len(cited), len(t["rows"]), len(pending) + len(none_), len(pending), len(none_)))

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
            _c = "0 \u2014 see note"
        else:
            _c = "%d  (%.2f%%)" % (_f[0], 100.0 * _f[0] / _COHORT)
        _rows.append((_r[0], _c) + tuple(_r[1:]))
    _t["rows"] = _rows
    _t["foot"] = [_status(_t) if f == "__CITATION_STATUS__" else f for f in _t["foot"]]
