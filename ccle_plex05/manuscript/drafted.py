# -*- coding: utf-8 -*-
"""Drafted section content, ported from the per-section build scripts on 5 September 2026.

References are cited as [[key]] against manuscript/refs.py; the builder numbers them globally.
"""

# ---------------------------------------------------------------- 2.7

P_2_7 = [
"Figure 6 shows that no tested compound descriptor predicts which lysines a compound prefers, and figure S7f shows that charge was never among the descriptors tested. That leaves the question open rather than settled. Here it is put directly, and to the reagent rather than to the compound — because in a competitive experiment it is the probe, not the compound, that decides which lysines are visible to be competed for at all.",
"Two probes were compared that share a desthiobiotin-hexanoyl core and differ only in the acyl leaving group (Figure 7a). The first, dTB-STP, carries a 4-sulfo-2,3,5,6-tetrafluorophenyl ester — an aromatic sulfonate, anionic at physiological pH [[hacker2017]]. The second carries an N-(cyanomethyl)-N-(4-nitrobenzenesulfonyl)amide and is neutral overall [[ryan2024]]. Both were applied to the same HeLa proteome, and twenty-two site features were measured on 167,057 lysines with an AlphaFold confidence of at least 70 [[abbasov2021,jumper2021,varadi2022]], including relative solvent accessibility on the Tien scale [[tien2013]], predicted pKa from PROPKA3 [[olsson2011]], and pocket descriptors from fpocket [[leguilloux2009]]. Of those sites, 6,294 were labelled only by the anionic probe and 5,851 only by the neutral one, with 12,298 labelled by both; the discriminating comparison is drawn within the 383 proteins that carry sites of both kinds, so that a probe which happened to sample a different set of proteins cannot create a difference.",
"**Two features of twenty-two differ, and both are charge.** Lysines labelled only by the anionic probe sit among more cationic residues than those labelled only by the neutral probe (cationic_8A, d = +0.233, q = 1.6 × 10⁻⁴) and in more positive local environments overall (net charge within 8 Å, d = +0.180, q = 5.2 × 10⁻³). **Zero of five pocket features and zero of nine shape features differ** (Figure 7b). The next-largest difference, aromatic contacts at d = 0.118, does not clear the threshold.",
"**The specificity is the argument.** A reagent that shifted everything would be a confound rather than a finding, and a reagent that shifted nothing would leave the panel uninformative. This one moves charge and leaves shape and pocket geometry alone, which is what a steering effect should look like and what a general difference in proteome sampling should not. Neither probe is individually quiet: measured against unlabelled lysines, eight of the twenty-two features are significant for each — 4,612 anionic-only sites across 1,026 proteins and 4,574 neutral-only sites across 1,007 (Figure S8a). Two of twenty-two is therefore a statement about the contrast between the probes, not about the power to detect one.",
"**Where the charge effect belongs, feature by feature.** Bootstrap intervals over proteins separate the reagent's contribution from the site's (Figure 7c, d). Cationic attraction is the **reagent**: the neutral probe's interval brackets zero and excludes the anionic probe's estimate. Acidic avoidance, predicted pKa and salt bridges are the **site**: the two probes are indistinguishable on all three. Net local charge is **mixed**, and it is the one feature for which a partial attribution is the honest reading. We do not infer reagent attribution from one probe reaching significance while the other does not — a difference between significant and non-significant is not itself a significant difference — which is why the panels show intervals rather than asterisks.",
"This changes how an earlier reading of the same data has to be worded. Labelled lysines sit in cationic, acid-poor environments — fewer acidic neighbours (d = −0.46), more cationic ones (+0.40), fewer salt bridges (−0.32) and lower predicted pKa (−0.30) at pLDDT ≥ 70 — and the natural reading is nucleophilicity: fewer carboxylates stabilising the ammonium, more free amine available to react. The decomposition says that reading is right for the pKa and salt-bridge terms and wrong for cationic attraction, which tracks the reagent. Both effects are real; they are not the same effect.",
"**Concavity, by contrast, belongs to the site.** Labelled lysines sit in more enclosed environments than unlabelled ones (d = +0.073, q = 1.1 × 10⁻⁵), and the two probes do not differ on it (d = 0.019, q = 0.71) — so this is a property of where reactive lysines are, not of what reaches them (Figure 7e). Two bounds belong with that result. Enclosure and relative solvent accessibility share about 74% of their variance across the 167,057 sites (Spearman ρ = −0.862), so the concavity effect is largely, though not wholly, a solvent-accessibility effect (Figure S8b); and the effect is reported within protein because pooling across sites reverses its sign, from +0.073 to −0.032 (Figure S8c).",
"One assumption bounds the charge result rather than threatening it. The neutral probe carries hydrogens α to an electron-withdrawing acylsulfonamide and may enolise at physiological pH, acquiring partial anionic character; its pKa is not known. If it does, the two probes are less different than assumed and the measured reagent effect is an **underestimate** — the interval runs from +0.233 under the neutral assumption to +0.46 if the probe is substantially ionised (Figure S8d). The assumption therefore cannot explain the result away.",
"The consequence for the field is a design principle rather than a caveat. A chemoproteomic probe is not a neutral window onto lysine reactivity: part of what a lysine-directed experiment measures is where its reagent is drawn. Coverage and reactivity are entangled at the point of capture, and a probe panel is a sampling variable that should be varied deliberately rather than held fixed by convention.",
]

CAP_2_7 = [
("Figure 7. The probe is not a neutral reporter: reagent charge steers labelling.", True),
("Two probes sharing a desthiobiotin-hexanoyl core and differing only in the acyl leaving group, applied to the same HeLa proteome. Twenty-two site features on 167,057 lysines at pLDDT ≥ 70.", False),
("**(a) The design.** dTB-STP (4-sulfo-2,3,5,6-tetrafluorophenyl ester, anionic) against a neutral N-(cyanomethyl)-N-(4-nitrobenzenesulfonyl)amide probe. 6,294 sites labelled only by the anionic probe, 5,851 only by the neutral one, 12,298 by both; the discriminating contrast is drawn within the 383 proteins carrying both.", False),
("**(b) Two features of twenty-two differ, and both are charge.** cationic_8A d = +0.233 (q = 1.6 × 10⁻⁴); net charge within 8 Å d = +0.180 (q = 5.2 × 10⁻³). Zero of five pocket features and zero of nine shape features differ. Figure S8a shows each probe is individually informative — eight of twenty-two features significant for each against unlabelled sites — so this is a statement about the contrast rather than about power.", False),
("**(c, d) What each probe prefers, and which part of the charge effect belongs to which.** Effect sizes are standardised (mean divided by standard deviation across proteins); intervals are 95% percentile bootstrap confidence intervals, resampling proteins as the unit over 4,000 resamples. Cationic attraction is the reagent: the neutral probe's interval brackets zero and excludes the anionic estimate. Acidic avoidance, predicted pKa and salt bridges are the site: no probe difference. Net local charge is mixed, and only there does a partial attribution survive. Intervals rather than significance markers are shown, because a difference between significant and non-significant is not a significant difference.", False),
("**(e) Concavity belongs to the site, not the reagent.** Labelled lysines sit in concavities (enclosure d = +0.073, q = 1.1 × 10⁻⁵) and the two probes do not differ on it (d = 0.019, q = 0.71). Figure S8b bounds this: enclosure shares about 74% of its variance with relative solvent accessibility. Figure S8c is why it is reported within protein — pooled across sites the sign reverses.", False),
("Site sets are chemoproteomic assignments, concurring at 97.82% between two independent routes. Both probes were applied to HeLa cells; labelling is expressed as % occupancy, and nothing here is a statement about a tumour or a patient.", False),
]

# ---------------------------------------------------------------- 2.8

P_2_8 = [
"Figures 4 to 7 describe what governs engagement among lysines that compounds have actually been tested against: a property of the site rather than of the molecule, poorly predicted by the structural features available to us, and read out through a probe that is not itself neutral. The acquired lysines of Figure 1 have been tested against nothing. Before any of those rules can be expected to apply to them, one question has to be settled — do they sit inside the same population, or does mutation create lysines in environments unlike the ones the atlas has sampled?",
"We compared nine environment features at 156,781 union positions with an AlphaFold confidence of at least 70 [[jumper2021,varadi2022]], against the engaged and unengaged lysines measured in the three screens [[abbasov2021,bracken2024,ryan2024]]. The features describe a site's surroundings within 8 to 10 Å: neighbour count, atom density, hydrophobic, aromatic and polar contacts, acidic and basic residues, net local charge, and other lysines.",
"**The comparison had to be anchored before it could be made.** Our first version of it measured union positions at Cβ and measured lysines at Nζ, and returned a separation of 2.256 against a control of 0.330 — a 6.84-fold ratio that would have been the strongest structural result in the paper. It is an artefact. Cβ and Nζ lie about 6.4 Å apart in an extended lysine side chain, and the two populations were being measured from different atoms. Figure 8a shows the artefact beside its correction. Three legitimate anchorings — both sides at Cβ, both sides modelled at Nζ, and both at Cβ under a symmetric convention — give 0.571, 0.535 and 0.583 against controls of 0.410, 0.331 and 0.407 (Figure S9a, b). The mismatched arm is retained in that panel as a negative control, because it is the analysis we first ran. Where side chains had to be built, all-trans rotamers reproduce the amine positions of 30,805 measured lysines to a median of 4.26 Å (Figure S9e); all-trans is the weakest possible construction, so that figure is an upper bound on the error rather than an estimate of it.",
"**On a common anchor, the three populations are difficult to tell apart.** Of the nine features, five have identical medians across acquired, engaged and unengaged lysines — aromatic contacts, polar contacts, acidic residues and other lysines, with neighbour count identical between acquired and engaged sites at twelve (Figure 8b). Where the populations do differ, they differ by less than one standard deviation of the acquired distribution: basic residues by −0.73, net local charge by −0.56, hydrophobic contacts by +0.33. Reported without a threshold, the probability that a randomly chosen acquired lysine lies closer to the centre of the engaged population than a randomly chosen engaged lysine does is 0.443, against 0.5 for populations that cannot be distinguished at all (Figure 8c). The reading does not depend on where a boundary is drawn: acquired lysines fall inside the 50, 75, 90 and 95% engaged regions at 42, 69, 82 and 88% under one anchoring and 45, 67, 83 and 89% under another \u2014 marginally short of the labelled-but-unengaged control at every cut, which runs 47, 73, 89 and 95%, and 50, 78, 92 and 96% \u2014 so the two comparators bracket the same region rather than separating (Figure S9f).",
"**Six of these sites have crystal structures, and they behave the same way.** For each we measured local Cα displacement in a co-local shell and normalised it to that protein's own wild-type-versus-wild-type floor, so the comparison asks whether the acquired lysine moves its environment more than two structures of the same wild-type protein move relative to each other (Figure 8d). Five of the six sit at that floor: NRAS Q61K at 1.07, AKT1 E17K at 1.12, PI3Kα E542K at 1.02, PI3Kα E545K at 1.22 and p53 R280K at 0.89. SHP2 E76K is the exception at 3.53 — and it is the one site with a substitution series, where other substitutions at the same position give 3.21. Where an acquired lysine does perturb its surroundings, then, it does so no more than a different substitution at the same position would. Two of the six comparisons carry a caveat that belongs in the open: the PI3Kα mutant structures 8GUA and 8GUD are cryo-EM, so those pairs cross a method boundary the crystal-versus-crystal pairs do not, and the NRAS floor rests on a single pair because the wild-type set was deliberately restricted to one nucleotide state — residue 61 sits in switch II, and an earlier run that mixed GDP and GTP structures inflated that floor to 4.25 Å and flattered the mutant. These are mutation sites taken from variant calling rather than chemoproteomic assignments, so the provisional-assignment caveat that applies elsewhere in this work does not apply here.",
"**What follows from this, and what does not.** Acquired lysines are not structurally unusual, which removes the cheapest reason to expect the engagement rules measured in Figures 4 to 7 not to reach them. It does not establish that they can be engaged. Eight of the nine features compared here are among the fourteen from which the composite model of Figure 5 was built — every one except the count of neighbouring lysines — and that model reaches an AUC of 0.617 against a measured reproducibility ceiling of at least 0.72, recovering 11% of engaged sites at a 5% false-positive rate. Features that predict engagement that poorly cannot be turned around and used to infer engageability from similarity. The engaged and unengaged comparators come from the three chemoproteomic screens — CyMPSA in HeLa, meroterpenoids in MDA-MB-231 and MDA-MB-436 — and describe those cells, not a tumour and not a patient. What this comparison supports is narrower and, for the purposes of this paper, more useful: nothing about the chemistry of these sites explains why they are absent from the record. The barrier is detection.",
"One acquired lysine makes that concrete. AKT1 E17K, the transforming pleckstrin-homology-domain mutation identified in breast, colorectal and ovarian cancer [[carpten2007]], has since been engaged selectively in cells through lysine targeting and neo-zinc chelation, and its geometry is available from PDB 8UVY, 54 atoms in the local shell (Figure 8e) [[craven2025]]. The search databases carried it as three peptides — two spanning the new lysine and one immediately downstream of it — and it returns zero spectra. That negative rests on a single experiment: five DepMap models carry AKT1 E17K, one of them, KU-19-19 (RRID:CVCL_1344), is in the 375-line proteomic set, and its plex is therefore the only one of the 42 in which the mutation was present to be found. It is one informative observation, not a survey. It is worth stating anyway, because the one acquired lysine that has been drugged is the one this workflow did not see.",
]

CAP_2_8 = [
("Figure 8. Acquired lysines are structurally ordinary.", True),
("Environment features at 156,781 union positions (pLDDT ≥ 70) against measured engaged and unengaged lysines.", False),
("**(a) One reference atom, on both sides.** The first version of this comparison measured union positions at Cβ and measured lysines at Nζ and returned a separation of 2.256 against a control of 0.330 — a 6.84-fold ratio that is an artefact of the mismatch, not a result. Cβ and Nζ are about 6.4 Å apart in an extended lysine. The panel shows the artefact and the correction.", False),
("**(b) The three populations differ far less than the anchor artefact implied.** Feature medians for acquired, engaged and unengaged lysines side by side. Five of the nine features have identical medians across all three populations; the largest difference is 0.73 standard deviations of the acquired distribution.", False),
("**(c) The distance between acquired and engaged lysines, against permutation nulls.** Reported threshold-free: the probability that a randomly chosen acquired lysine sits closer to the engaged centre than a randomly chosen engaged lysine does is P = 0.443, against 0.5 for indistinguishable populations. The populations overlap heavily.", False),
("**(d) Crystallographic validation at six sites.** Local Cα RMSD in a co-local shell, each mutant normalised to its own protein's wild-type floor. Five of six sit at that floor (0.89–1.22); SHP2 E76K is the exception at 3.53, and is the only site carrying a substitution series, where other substitutions at the same position give 3.21. Six sites were defined and all six measured. The NRAS wild-type set is restricted to a single nucleotide state so that every pair is state-matched, which leaves its floor resting on one pair — annotated “floor n = 1” on the panel. The PI3Kα mutant structures 8GUA and 8GUD are cryo-EM. These are mutation sites from variant calling, not chemoproteomic assignments, so the provisional-assignment caveat does not apply to them. Structures used, all previously deposited and none determined here \u2014 PTPN11: 2SHP, 4DGP, 4DGX, 4H1O, 5EHR (wild type), 6CMS, 6CRF, 6CRG (E76K), 7EMN (E76A), 6CMR (E76D), 5IBS (E76Q); TP53: 2OCJ, 1TSR, 3KMD, 2AC0 (wild type), 6FF9 (R280K); AKT1: 1UNQ, 1H10, 1UNR, 1UNP (wild type), 2UZS, 8UVY (E17K); PIK3CA: 4OVU, 2RD0, 4L23, 5DXT (wild type), 8GUA (E542K), 8GUD (E545K); NRAS: 5UHV, 8TBI (wild type), 8VM2 (Q61K).", False),
("**(e) AKT1 E17K — engaged selectively in cells by others, and undetected here.** Site geometry from PDB 8UVY (54 atoms in the local shell). Craven et al. [[craven2025]] engaged this lysine selectively through lysine targeting and neo-zinc chelation; searched as the three peptides the databases carried — two spanning the new lysine, one immediately downstream — it returns zero PSMs. Five DepMap models carry this mutation and one of them, KU-19-19, is in the 375-line proteomic set, so the search had one informative chance rather than 42, and the zero is a single observation rather than a survey. The point is the contrast between what is chemically reachable and what this workflow can see.", False),
("The engaged and unengaged comparators come from the three chemoproteomic screens (CyMPSA in HeLa; meroterpenoids in MDA-MB-231 and MDA-MB-436); engagement is expressed as % occupancy. Nothing in this figure is a statement about a tumour or a patient.", False),
]

# ---------------------------------------------------------------- 2.6
P_2_6 = [
"If engagement is a property of the site rather than of the molecule, then the same lysines should be preferred by compounds that have little else in common. Figures 4 and 5 measured site properties; this section asks the complementary question of the compounds, using a design in which each compound is its own control.",
"Occupancy is the share of probe labelling a compound competes away, 100 \u00d7 (1 \u2212 compound / DMSO) \u2014 for a covalent compound, the share of that lysine's copies it has taken. Across the three screens [[abbasov2021,bracken2024,ryan2024]] there are 1,642,330 site \u00d7 compound measurements, of which 8,019 (0.49%) reach the 60% occupancy that defines engagement here; the median compound quantifies 6,332 lysines and engages 14 of them (Figure 6a). The distribution is smooth through 60%, so that threshold is a line drawn through a continuum rather than a boundary in the data. For each compound separately, the sites it engaged are compared with the sites it was quantified at and left below the cut, and each site feature is standardised within that compound. Only then are compounds compared, and what is compared is their preferences rather than their site lists \u2014 so a compound tested on a different set of lysines cannot manufacture a difference.",
"**Thirteen of fourteen site features point the same way in every compound.** Across n = 84 compounds with a usable contrast, engaged lysines are consistently more buried, more crowded, more hydrophobically packed and lower in predicted pKa (Figure 6b). The single exception is the count of nearby lysines within 8 \u00c5 (+0.038, q = 0.097), which is a property of the neighbourhood rather than of the lysine in question. The three charge features clear significance but are the three smallest of the thirteen that do, and figure S7b shows that two of them do not strengthen with occupancy.",
"**The preference is a gradient, not a threshold effect.** Across 958,108 site \u00d7 compound observations binned by occupancy, every one of six site features moves monotonically and accelerates through the 60% line (Figure 6c). The tiers are very unequal \u2014 734,394 observations in the 0\u201310% bin against 648 at \u2265 90% \u2014 so the direction is secure while the acceleration at the right-hand end rests on the two smallest bins. Three threshold-free continuous measures agree (figure S7c), so nothing here depends on a 60% cut existing at all.",
"**No tested compound descriptor predicts which lysines a compound prefers.** Correlating every compound descriptor against every site-feature preference gives 21 of 84 pairs clearing q < 0.05 when the three screens are pooled (82 compounds, against 4.2 expected), but **0 of 84 within NatChem alone** (62 compounds) and **0 of 72** in the larger within-dataset slope analysis (173 compounds, against 3.6 expected) (Figure 6d). Descriptors were computed with RDKit [[rdkit]].",
"**The pooled positives are provenance, not chemistry.** Molecular weight all but names the screen a compound came from \u2014 medians 306 Da for NatChem, 386 for the meroterpenoids and 480 for CyMPSA, Kruskal\u2013Wallis p = 1.3 \u00d7 10\u207b\u2079 \u2014 and each screen is a different cell line in a different sample state (Figure 6f). More directly, the descriptors that separate the three screens are the descriptors that produced the pooled false positives: across ten descriptors the two rank orders correlate at Spearman \u03c1 = 0.69, p = 0.027. Pooling therefore lets provenance masquerade as chemistry, which is why panel d has three columns rather than one.",
"**The breadth gradient of Figure 5e does not depend on where the occupancy cut is drawn.** Repeating it here on the site features, the most broadly engaged sites differ from a never-engaged site by 0.48 SD on the average feature, against 0.09 for shuffled labels, and the gradient holds at five cuts from 40% to 80% (Figure 6e). The two cuts that are not strictly monotonic, 40% and 80%, are the two where group sizes collapse (figure S7e).",
"**The boundary of this claim is measured rather than asserted, and it is narrower than it first appears.** The descriptors tested contain nothing that encodes charge: the maximum |\u03c1| any of them reaches with a charge-related preference is 0.193, against a median |\u03c1| of 0.552 for molecular weight. And every shape term is achiral by construction \u2014 NPR1 and NPR2 are ratios of principal moments of inertia [[sauer2003]], and reflecting a conformer changes them by exactly zero, while re-embedding an enantiomer changes them by up to 0.281 — larger than the 0.123 SD that separates one compound from another, so that number is conformer sampling and not chirality (figure S7f). The result is therefore that engagement is a property of the lysine **over the descriptors tested**, and it is silent on exactly the two axes this paper goes on to test \u2014 charge in Figure 7, and chirality in Figure 5f, both of which show molecular properties that do matter. It must not be read as a claim that no property of the molecule matters.",
]

CAP_2_6 = [
"Site sets are chemoproteomic assignments from three screens \u2014 NatChem in Ramos and MDA-MB-231, meroterpenoids in MDA-MB-231 and MDA-MB-436, CyMPSA in HeLa \u2014 concurring at 97.82% between two independent routes. Engagement is \u2265 60% occupancy.",
"**(a) What is measured, and what is done with it.** A probe labels lysines; a compound competes with the probe. Occupancy is the share of probe labelling the compound competes away, 100 \u00d7 (1 \u2212 compound / DMSO). The distribution shows every site \u00d7 compound measurement in the three screens (1,642,330); 8,019 of them (0.49%) reach 60%, and the three screens agree, so the cut is a line drawn through a continuum rather than a boundary in the data. The median compound quantifies 6,332 lysines and engages 14. The spike in the last bin is the clip at 100%. Right: the design \u2014 for one compound, the sites it blocked to \u2265 60% against the sites it was quantified at and left below 60%, standardised within that compound, so that compounds are compared on their preferences and not on their site lists.",
"**(b) Thirteen of fourteen site features point the same way in every compound.** Mean standardised shift between engaged and below-cut sites with 95% intervals, across n = 84 compounds. Engaged lysines are more buried, more crowded, more hydrophobically packed and lower in predicted pKa. Rust marks the single exception, nearby lysines within 8 \u00c5 (+0.038, q = 0.097): local lysine density is a property of the neighbourhood, not of this lysine. The three charge features clear significance but are the three smallest of the thirteen that do, and figure S7b shows two of them do not strengthen with occupancy.",
"**(c) The occupancy profile has no step at 60%.** 958,108 site \u00d7 compound observations binned by occupancy, six site features as standardised shifts from the 0\u201310% tier. Every curve moves monotonically and accelerates through the dashed 60% line. The tiers are very unequal \u2014 734,394 observations in 0\u201310% against 648 at \u2265 90% \u2014 so the direction is secure and the acceleration at the right-hand end rests on the two smallest bins.",
"**(d) No compound property predicts which lysines a compound prefers.** Every descriptor \u00d7 site-preference correlation, one point per pair. Pooled across the three screens (82 compounds), 21 of 84 pairs clear q < 0.05 against 4.2 expected, reaching |\u03c1| = 0.45. Within NatChem alone (62 compounds), 0 of 84, reaching 0.24. In the within-dataset slope analysis (173 compounds), 0 of 72 against 3.6 expected, reaching 0.18. Rust marks pairs clearing q < 0.05. This is a statement about the descriptors that were tested; figure S7f says which those are, and none encodes charge or chirality.",
"**(e) The gradient does not depend on where the cut is drawn.** Sites grouped by how many compounds engage them, against how far the group's average environment sits from a never-engaged site \u2014 the average shift per site feature, in SD, over fifteen features — the fourteen of panel b together with pocket occupancy. Repeated at five occupancy cuts. At 60% the most broadly engaged sites differ by 0.48 SD on the average feature, against 0.09 for shuffled labels (grey; the null shown is the one for the 60% cut). Rust marks 40% and 80%, the two cuts that are not strictly monotonic; figure S7e shows both are where group sizes collapse.",
"**(f) Why the analysis is run within a dataset.** Left: molecular weight all but names the screen a compound came from (medians NatChem 306, meroterpenoids 386, CyMPSA 480 Da; Kruskal\u2013Wallis p = 1.3 \u00d7 10\u207b\u2079), and each screen is a different cell line in a different sample state. Right: a wider sweep than panel d \u2014 ten RDKit descriptors against fifteen site features, so its counts are not panel d's \u2014 in which the descriptors that separate the three screens are the descriptors that produced the pooled false positives (Spearman \u03c1 = 0.69, p = 0.027 across the ten). TPSA is the labelled outlier: fourth in \u03b7\u00b2, it produced none. Pooling therefore lets provenance masquerade as chemistry.",
]


# ================================================================ 2.5  (Figure 5)
P_2_5 = [
"Figure 4 established that engaged lysines sit in measurably different environments from lysines "
"that are merely labelled. Whether that difference is large enough to act on is a separate "
"question: given a structure and no assay, can the ligandable lysines be picked out in advance? "
"Lysine ligandability has been profiled globally before [[hacker2017]], but the accuracy of a "
"structural prediction has never been measured against what the assay itself reproduces. Both "
"quantities come from one analysis, so this section reports them together.",

"A logistic model was fitted on fourteen features \u2014 thirteen describing the lysine's "
"environment and one describing how confidently the residue was modelled \u2014 across 30,805 "
"lysines on 6,201 proteins, of which 3,047 (9.9%) are engaged by at least one compound. It was "
"cross-validated by protein rather than by site, because lysines on one protein are correlated and "
"a site-level split leaks information between training and test. Held out, the composite reaches "
"**AUC 0.617 [0.594, 0.640]** against a train AUC of 0.619 \u2014 an optimism of +0.002, so the "
"number is not depressed by overfitting (Figure 5a; figure S6d). Refitting on the thirteen "
"environment features alone gives 0.611 [0.590, 0.632], which forecloses the objection that the "
"model is partly predicting how well a site was modelled rather than what the site is. This set is "
"*the composite* throughout; the fourteen site features of Figure 6 are a different set, and the "
"two are never pooled.",

"Whether 0.617 is disappointing depends on what a perfect model could reach, and that limit is set "
"by the assay rather than by chemistry. Engagement is the outcome of a measurement, not a property "
"a site simply has: if the same lysine is assayed twice and the two assays disagree, no predictor "
"can agree with both. Splitting each site's compound measurements into two independent halves, "
"taking the label from one half and the score from the other, gives the AUC a perfect structural "
"model would reach if structure captured everything the assay captures. At the deepest split the "
"data support \u2014 twenty compounds in each half \u2014 that ceiling is "
"**0.723 [0.706, 0.744]**, and the composite falls 0.106 short of it (Figure 5a).",

"**Two things make 0.106 the smallest shortfall the data permit.** Each half carries fewer "
"compounds than the full screen and is therefore noisier than the real label, so the estimate is a "
"lower bound on the true ceiling \u2014 and it is still rising at the deepest split. And 0.723 is "
"the lowest of the three screens: the CyMPSA and meroterpenoid ceilings sit between 0.90 and 0.95 "
"in every band where they are populated. The comparison is drawn against the weakest ceiling in "
"the data, and is quoted as \u2265 0.72.",

"**The ceiling is not an artefact of disordered sites.** Computed within five bands of AlphaFold confidence at the lysine "
"[[jumper2021,varadi2022]], the NatChem ceiling runs from 0.672 to 0.734 with no trend across "
"them, and reads 0.734 on the 8,436 sites at pLDDT \u2265 90, where the structure is most "
"trustworthy and the shortfall is therefore largest (Figure 5b). The other two screens are higher "
"in every populated band.",

"**Nor is the shortfall a feature-count problem.** The fourteen features do not carry fourteen "
"independent properties: the first principal component alone accounts for 48.9% of their variance, "
"and five components clear 80% (Figure 5d). Burial, crowding, hydrophobic packing and predicted "
"pKa [[tien2013,olsson2011]] are facets of one situation \u2014 a lysine in a pocket \u2014 "
"measured several ways, which is also why the composite improves on its best single feature, "
"distance to the protein centroid, by only +0.019 (figure S6b). More features of this kind would "
"not close the gap; a different kind of information would be needed.",

"**No operating point on the curve would change a decision.** An AUC of 0.617 reads as somewhat "
"predictive, which is why the receiver-operating curve matters more here than the summary "
"statistic. At a 5% false-positive rate the composite recovers 11% of engaged sites, and at 10% it "
"recovers 19% (Figure 5c). Nothing on this curve would alter which lysines a screen went after, "
"which is the practical content of the result: ligandability at these sites, in these cell lines, "
"has to be measured.",

"**What separates ligandable sites is breadth of engagement rather than a structural signature.** "
"Grouping sites by how many of the 34 tested chemotypes engage them at \u2265 60% occupancy, the "
"distance from a never-engaged site rises monotonically \u2014 0.49 for sites engaged by one "
"compound (n = 1,529), 0.91 for two or three (n = 624), 1.22 for four or five (n = 142) and 2.65 "
"for six or more (n = 82) \u2014 against nulls of 0.12, 0.18, 0.37 and 0.49 from 2,000 label "
"shuffles (Figure 5e). The direction is secure across the gradient; the magnitude at the far end "
"rests on 82 sites, and adjusting for how often each site was tested does not remove it "
"(figure S6e).",

"Two enantiomer pairs have enough shared coverage to compare. Guadial B engages nineteen sites "
"with both hands against 0.1 expected by chance and eight with a single hand; euglobal G2 engages "
"six with both against 0.0 expected and one with a single hand (Figure 5f). The dominant signal is "
"agreement: two molecules of identical constitution and opposite three-dimensional shape "
"overwhelmingly engage the same lysines. The minority engaged by a single hand cannot be intrinsic "
"reactivity, since the two hands present the same chemistry \u2014 but with n = 2 pairs that is a "
"bounded observation, not a selectivity claim.",

"Both halves of the heading are bounded. *Reproducible* is a measurement rather than an "
"assumption \u2014 the same split-half analysis supplies the ceiling. *Poorly predicted* is "
"bounded to these features and these cell lines, and Figure 6 asks the complementary question of "
"whether any tested property of the compound predicts which lysines it prefers. Carried across "
"screens rather than within one, this reasoning rests on 25 sites in 21 dimensions with two of six "
"comparisons short of their null (figure S6g), reported there as the weakest evidence in the paper "
"rather than as reinforcement.",
]

CAP_2_5 = [
"Ligandability measured in three chemoproteomic screens \u2014 NatChem in Ramos and MDA-MB-231, "
"meroterpenoids in MDA-MB-231 and MDA-MB-436, CyMPSA in HeLa \u2014 concurring at 97.82% between "
"two independent site-assignment routes. Engagement is \u2265 60% occupancy. *The composite* is a "
"logistic model on thirteen environment features plus AlphaFold confidence at the lysine, fitted "
"on 30,805 lysines across 6,201 proteins and cross-validated by protein. Intervals on model AUCs "
"are the 2.5th to 97.5th percentile across 20 protein-held-out splits; intervals on ceilings are "
"the same percentiles across 50 random half-splits.",

"**(a) The measured model against the reproducibility ceiling.** The composite reaches AUC 0.617 "
"[0.594, 0.640] against a ceiling of 0.723 [0.706, 0.744] derived from agreement between "
"independent halves of the same screen \u2014 a gap of 0.106. Train AUC 0.619, optimism +0.002, so "
"the gap is not overfitting. The ceiling drawn is NatChem's at the deepest split the data support, "
"twenty compounds per half; it is the lowest of the three screens and is still rising at that "
"depth, so \u2265 0.72 is a floor rather than an estimate. The conclusion is that ligandability "
"must be measured rather than predicted \u2014 from these features, in these cell lines.",

"**(b) The ceiling does not depend on structural order.** The ceiling for each screen across five "
"bands of AlphaFold confidence at the lysine \u2014 disordered (< 50), low (50\u201370), "
"confident (70\u201390), very high (\u2265 90) and no structure \u2014 with the composite's "
"0.617 drawn across as a reference line. Bands populated: NatChem 5, CyMPSA 5, meroterpenoids 2. "
"NatChem's ceiling runs 0.672 to 0.734 with no trend across the bands and reads 0.734 on the 8,436 "
"sites at pLDDT \u2265 90; CyMPSA and the meroterpenoids sit between 0.90 and 0.95 in every "
"populated band. Panel a therefore cannot be dismissed as an artefact of disordered sites.",

"**(c) There is no operating point at which the model is useful.** Receiver-operating curve for "
"the composite, every site predicted once, with the band spanning 95% of 20 held-out splits. 11% "
"of engaged sites are recovered at 5% false positives and 19% at 10%. An AUC of 0.617 reads as "
"somewhat predictive; this panel says what that buys in practice.",

"**(d) The shortfall is not a feature-count problem.** Variance captured against the number of "
"independent properties kept. The first component alone carries 48.9% and five clear the 80% line, "
"so the fourteen features carry about five independent properties. Adding more features of the "
"same kind would not close the gap in panel a \u2014 a different kind of information would be "
"needed.",

"**(e) Breadth across chemotypes.** Sites grouped by how many of the 34 tested chemotypes engage "
"them at \u2265 60% occupancy, against multivariate distance from a never-engaged site: one "
"compound, n = 1,529, 0.491; two or three, n = 624, 0.913; four or five, n = 142, 1.215; six or "
"more, n = 82, 2.650. Grey is the null over 2,000 label shuffles (0.120, 0.183, 0.369, 0.487); all "
"four bins clear it at p < 0.0005, the resolution limit of 2,000 shuffles. The gradient is "
"monotonic and the last bin is the smallest, so "
"the direction is secure while the magnitude at the right-hand end rests on 82 sites. Figure S6e "
"repeats the test adjusted for how often each site was tested.",

"**(f) Enantiomers.** The two pairs with enough shared coverage to compare, in the order plotted. "
"Guadial B, measured at 7,043 sites: nineteen sites engaged by both hands against 0.1 expected by "
"chance, eight by a single hand. Euglobal G2, measured at 7,560 sites: six by both against 0.0 "
"expected, one by a single hand. Each pair overwhelmingly engages the same sites; the minority "
"engaged by one hand only "
"cannot be intrinsic reactivity, since the two hands present the same chemistry. **n = 2 pairs**, "
"reported as a bounded observation rather than a selectivity claim.",

"*Engagement measured in cell lines. Nothing here is a statement about a tumour.*",
]


# ================================================================ 2.4  (Figure 4)
P_2_4 = [
"Two questions are routinely run together in work on covalent lysine chemistry: what makes a "
"lysine reactive enough to be labelled by a probe at all, and what makes a labelled lysine "
"blockable by a competing molecule [[hacker2017]]. One contrast cannot answer both, so two were "
"run against two different control groups (Figure 4a). Reactivity compares probe-labelled lysines "
"against every other lysine on the same protein; ligandability compares engaged lysines against "
"all labelled ones at \u2265 60% occupancy. The setting is 245,684 lysines across 6,201 proteins, "
"of which 30,805 (12.5%) are labelled. Fifteen site features were measured, six describing charge "
"and nine describing shape [[tien2013,olsson2011]]; the AlphaFold terms are excluded throughout, "
"because confidence is the gating variable here and not a property of the site.",

"Both contrasts use the same estimator \u2014 one difference per protein, tested against zero \u2014 "
"so the two sets of numbers are directly comparable, and both are run through the same ladder of "
"structural-confidence gates. The gate is always applied to both groups at once; gating one side "
"would manufacture a difference rather than reveal one.",

"**Ungated, the two contrasts are indistinguishable.** The mean distance between a feature's two "
"effects is 0.074 for charge features and 0.075 for shape features (Figure 4b), and within each "
"class the two contrasts track each other almost perfectly \u2014 0.98 for charge, 0.97 for shape. "
"An exact permutation over all 5,005 ways of splitting fifteen features six against nine gives an "
"observed separation of +0.009, p = 0.412. Read on its own, this panel says that reactivity and "
"ligandability are one property measured twice.",

"They are not, and the gate is what shows it. Every feature is measured on a predicted structure "
"[[jumper2021,varadi2022]], and where the model is uncertain there is no structure to measure: "
"figure S5f shows two genuinely engaged sites from the same models, PFKL K677 at pLDDT 95 kept "
"and LYSMD4 K72 at pLDDT 28 discarded, with the same 8 \u00c5 sphere drawn around each. Gating "
"costs sites \u2014 30,805 labelled lysines fall to roughly 17,000 at pLDDT \u2265 90 (Figure 4d) "
"\u2014 and some of those discarded sites are engaged, which is a real price and not a free "
"cleanup.",

"**At pLDDT \u2265 90 the two contrasts come apart, and only across the shape features.** The mean "
"pair distance for shape rises from 0.075 to 0.126 while charge is essentially unchanged at 0.068 "
"(Figure 4c). The class correlations move accordingly: charge from 0.98 to 0.99, shape from 0.97 "
"to 0.35. The permutation observed value rises to +0.639, p = 0.003.",

"**The packing features are where it happens, and the spread shows it is not range restriction.** "
"Crowding runs +0.197 to \u22120.042 in reactivity while holding at +0.143 to +0.100 in "
"ligandability; atom density +0.244 to +0.003 against +0.157 to +0.109; polar contacts +0.183 to "
"\u22120.002 against +0.185 to +0.197. In each case the feature's own standard deviation holds or "
"widens under the gate \u2014 crowding 4.56 to 4.70, polar 1.66 to 1.81, aromatic 1.08 to 1.25 "
"\u2014 so the reactivity effect is not shrinking because the spread shrank. Two shape features "
"are exceptions and are reported as such: exposure and distance to the centroid lose 12% and 29% "
"of their spread under the gate, so for those two the collapse is partly range restriction rather "
"than artefact removal.",

"**What Figure 4e measures is whether the two contrasts move together, not whether either is "
"stable.** Several charge effects change substantially under the gate \u2014 anionic contacts "
"strengthen from \u22120.307 to \u22120.454, salt bridges from \u22120.221 to \u22120.327 "
"\u2014 but they change in step on both sides, the ligandability arms running \u22120.210 to "
"\u22120.336 and \u22120.183 to \u22120.286. The shape features do not: their reactivity arms "
"collapse toward zero or reverse sign while their ligandability arms hold direction. That "
"asymmetry, not the size of any single effect, is the result.",

"**Both arms of the test have to be reported together, because either alone misleads.** Across the "
"ladder the shape coupling falls monotonically \u2014 0.97, 0.90, 0.77, 0.35 at gates of none, "
"\u2265 50, \u2265 70 and \u2265 90, with exact p of 0.412, 0.034, 0.009 and 0.003 \u2014 while "
"charge holds at 0.98 to 0.99 throughout (Figure 4e). The separation therefore does not exist "
"before the gate is applied, and it is not manufactured by the gate keeping only the strongest "
"sites, which is the control in figure S5g. Nor does it depend on the 60% cut existing: the shift "
"from the lowest occupancy tier rises monotonically across all ten tiers with no discontinuity "
"where the cut falls (Figure 4f).",

"The claim is conditional. Where predicted structure is "
"trustworthy, shape distinguishes ligandability from reactivity and charge does not; ungated, "
"neither does, because the apparent shape contribution to reactivity is carried by sites whose "
"geometry cannot be measured. This does not say that shape is irrelevant to reactivity in "
"general \u2014 it says that the evidence for it in these data comes from the sites where the "
"evidence is weakest. It also inherits the gate's own bias: AlphaFold confidence is not random "
"with respect to biology, and a claim that holds only above pLDDT 90 is a claim about ordered "
"protein.",
]

CAP_2_4 = [
"Reactivity and ligandability are measured on the same lysines and with the same estimator \u2014 "
"one difference per protein, tested against zero \u2014 so the two contrasts are directly "
"comparable. Site sets here and in Figures 5\u20137 are chemoproteomic assignments concurring at "
"97.82% between two independent routes (44,736 of 45,552 usable sites); the figure compares site "
"*sets*, and individual residue numbers within them are provisional. Every structural gate is "
"applied to both groups of a contrast at once.",

"**(a) The design.** Two comparisons against two control groups, run before and after a "
"structural-confidence gate so the gate's effect is visible rather than assumed. 245,684 lysines "
"across 6,201 proteins, 30,805 of them labelled (12.5%). Fifteen features, six charge and nine "
"shape, the identical set on both ladders; AlphaFold terms are excluded because confidence is the "
"gating variable rather than a site property.",

"**(b, c) The paired display, ungated and at pLDDT \u2265 90.** Each feature's two effects are "
"drawn on one line \u2014 filled circle for the effect on being labelled, open diamond for the "
"effect on being blocked \u2014 so that the movement under the gate reads directly instead of "
"having to be reconstructed from two separate panels. Average distance between the pair: ungated "
"charge 0.074, shape 0.075; gated charge 0.068, shape 0.126.",

"**(d) The gate, step by step.** Four representative features' reactivity effect at each rung, "
"with the lysines retained printed beneath: none 30k/214k, \u2265 50 28k/161k, \u2265 70 "
"25k/141k, \u2265 90 17k/80k (labelled / unlabelled). Charge features hold or strengthen across "
"the ladder; shape features fall.",

"**(e) The permutation test, over the gate ladder and against the null.** Two stacked axes under "
"one letter. The upper axis is how closely the two contrasts track as the gate tightens: charge "
"stays at 0.98\u20130.99 throughout while shape falls 0.97, 0.90, 0.77, 0.35 at gates of none, "
"\u2265 50, \u2265 70 and \u2265 90, with exact permutation p of 0.412, 0.034, 0.009 and 0.003. "
"The lower axis is the null itself \u2014 all 5,005 ways of splitting these fifteen features six "
"against nine \u2014 with both observed values marked: +0.009 ungated (p = 0.412) and +0.639 at "
"pLDDT \u2265 90 (p = 0.003). Both are stated because either alone misleads: the separation does "
"not exist until the gate is applied, and it is not manufactured by the gate keeping the "
"strongest sites (figure S5g).",

"**(f) Occupancy dose-response.** Shift from the 0\u201310% occupancy tier as occupancy rises, so "
"the result does not depend on the 60% engagement cut existing at all. Monotonic across all ten "
"tiers, with no discontinuity where the cut falls. The fully threshold-free version is figure S7c.",

"*Effects are within-protein standardised differences. Measured in cell lines; nothing here is a "
"statement about a tumour.*",
]


# ================================================================ 2.3  (Figure 3)
P_2_3 = [
"Figure 2 recovered a small number of mutant peptides, and figure S3e bounds what they rest on: "
"53 sites on at least one spectrum, 12 on two or more, four on two or more distinct peptides, and "
"none on more than one experiment; at margin \u2265 10 those become 30, six, two and none. A "
"result carried by single-spectrum identifications has to "
"survive an error-rate audit before anything is built on it, and the usual answer \u2014 that the "
"search was filtered at 1% false discovery \u2014 does not settle it. Target-decoy control is "
"calibrated over a whole search [[elias2007]], and variant peptides are a tiny, atypical subset of "
"one: a rate that holds globally need not hold within them. Chemoproteomic work on "
"mutation-created residues has answered this in two ways. DrugMap searches canonical and "
"cell-line-specific mutant sequences in one pass, takes consensus across four search engines, "
"and filters to 1% peptide false discovery at run level against reversed canonical decoys "
"[[takahashi2024]]. The chemoproteogenomic platform searches in two stages, removing spectra "
"confidently assigned to canonical sequences before the variant database is searched, so that "
"fewer variant peptides can be matched to wild-type spectra [[desai2024]]. Both control the "
"error. Neither reports a realised false discovery proportion measured within the variant "
"class, which is what follows \u2014 and what follows is a property of this search and this "
"database, not of theirs.",

"The test is burden-matched entrapment, and it is an application of a published method rather than "
"a new one [[wen2025]]. Each variant peptide is given a partner built from its own residues, "
"permuted, with the C-terminal lysine or arginine kept in place (Figure 3a). The pair therefore "
"shares an exact precursor mass and an identical length, and only fragment ions can separate them "
"\u2014 which is what makes this an upper-bound test rather than a decoy in a different costume. "
"Across eight plexes of the CCLE proteomes [[nusinow2020]], 589,588 variant-pair rows yielded "
"entrapment sets of 147,397 each, every one of them sharing residues with its target "
"(figure S4a\u2013c). An entrapment sequence cannot be in the sample, so every entrapment "
"identification is a false one, counted directly rather than estimated.",

"**At the reported filter, the realised error rate is roughly fivefold higher than the nominal "
"one.** Over 1,489 rows the false discovery proportion is **4.70%** with no score margin imposed "
"\u2014 1,419 variant peptides against 70 entrapment (Figure 3b). A search reporting 1% false "
"discovery is delivering about 4.7% within this class of peptide.",

"**A score margin recovers the intended rate.** The margin is the gap between a spectrum's best "
"interpretation and its next best, and requiring it to be large removes exactly the "
"identifications that a mass-matched impostor could have produced. The proportion falls to 3.91% "
"at a margin of 2, 2.36% at 5, and **0.83% at 10**, where 832 variant peptides survive against "
"seven entrapment. Beyond that the estimate stops being informative rather than continuing to "
"improve: only one entrapment peptide remains at a margin of 15, so the 0.18% there and the 0.31% "
"at 20 are the behaviour of a single count, not evidence of a cleaner set. Every error rate "
"quoted in this paper, and every identification quoted from the entrapment analysis, is taken "
"at margin \u2265 10; Figure 2\u2019s detection census is reported at the search\u2019s own 1% "
"filter and is given at both thresholds where it appears.",

"**That rate is measurable on the acquired-lysine class alone, and measured there it is not lower.** "
"A rate estimated on every variant peptide and then applied to one kind of substitution is an "
"inference across populations, and it does not have to hold: a substitution to lysine creates a "
"tryptic cleavage site, so the peptide the search must find is not the peptide it would have "
"looked for. The restriction is legitimate here because the entrapment design is pairwise \u2014 "
"each variant peptide has one paired entrapment peptide built from the same burden, and the pair "
"carries the same class label \u2014 so keeping only lysine-creating substitutions preserves the "
"r = 1 match instead of slicing an unmatched comparison (note N18). Restricted that way, 66 "
"acquired-lysine peptides compete with five entrapment at margin 0: **7.04%, against 4.58% for the "
"rest of the variant class**. The counts cannot resolve that difference (Fisher exact, odds ratio "
"1.58, p = 0.38), but the direction is the one that matters, because it means borrowing the "
"whole-class figure is not the conservative choice. At margin \u2265 10 the class contributes zero "
"entrapment peptides against the 38 variant peptides this paper quotes, which **bounds** its "
"false-discovery proportion at 7.9% by the rule of three rather than measuring it at zero. Seven "
"entrapment peptides across 42 plexes is the entire supply of counted errors, so no single class "
"of substitution can be resolved to better than a few percent on its own. The whole-class 0.83% "
"therefore remains the tighter estimate and is the number quoted; what the class-restricted run "
"adds is that the acquired-lysine set is not carrying a hidden excess of false discoveries past "
"the margin filter.",

"Two denominators are in use and they are not interchangeable. Figure 3b counts **distinct "
"peptides**, one best spectrum each, which is the unit the entrapment pairing is defined over, and "
"gives 0.83%. Figure 3c and figure S4e count **identifications**, and give six entrapment among "
"1,033 identifications \u2014 **0.58% pooled**, with entrapment appearing in four of 42 plexes "
"(Figure 3c). Six in 1,039 is 0.58% and seven in 839 is 0.83%; both are correct, and neither "
"number should be quoted with the other's counts.",

"**Decoys do not see this, which is why the margin is applied rather than the reported q-value "
"trusted.** On plex 2, where 79,365 identifications are available to compare, the nominal q-value "
"from decoys reaches 0.60% where the observed proportion from entrapment reaches 0.99% \u2014 an "
"underestimate of 1.7-fold (figure S4f).",

"**The surviving identifications are ordinary spectra, not marginal ones.** Judged on properties of "
"the spectrum alone, with no reference to the search engine or the database, the 38 acquired-lysine "
"hit spectra sit at 0.91 times the background peak count and 1.19 times its signal-to-noise across "
"802,067 spectra from the same searches (Figure 3d). Precursor mass errors are likewise ordinary: "
"91% fall within 0.05 Da, with medians of 0.0013 Da for variant identifications against 0.0009 Da "
"for background (Figure 3e). The one feature that invites suspicion is a tail at +1 Da, the "
"signature of the first isotope being selected instead of the monoisotopic peak; it accounts for "
"7.7% of variant PSMs and 7.3% of background, the same rate in both, so it is instrument behaviour "
"and not a property of variant peptides.",

"**The margin does work rather than simply discarding the weakest identifications of every kind.** "
"Entrapment peptides are identified, and they score like a population that cannot be present should "
"(Figure 3f); it is that separation the margin exploits. A boundary control makes the same point "
"from the other direction: the same genomic region shifted by one residue, tryptic in neither the "
"wild-type nor the mutant protein, produces two identifications and none after the margin, against "
"58 falling to 25 for the real downstream target (figure S4d).",

"**The error rates quoted in this paper are measured at margin \u2265 10 rather than assumed, "
"and Figure 2\u2019s census is reported at both thresholds so the two are never conflated.** "
"Neither makes 0.58% "
"a per-site confidence, and a single spectrum at margin \u2265 10 remains a single spectrum \u2014 "
"which is why Figure 2's detections are treated as a lower bound on coverage rather than a catalogue "
"of sites. Figure S4h applies the published estimator to the entire plex-2 search rather than to the "
"variant subset, so it is not comparable with 0.58% and is labelled in the supplement as such.",
]

CAP_2_3 = [
"Error rates measured on the CCLE proteomes by burden-matched entrapment at r = 1 (FDRBench "
"1.1.1). Entrapment sequences are each target's own residues permuted with the C-terminal lysine "
"or arginine kept, so a pair shares an exact precursor mass and an identical length and only "
"fragment ions separate them. Every entrapment identification is by construction false. "
"**Two denominators appear and are not interchangeable: panel b counts distinct peptides, one "
"best spectrum each; panel c counts identifications.**",

"**(a) The design.** Burden-matched entrapment at r = 1. Entrapment sequences are matched to the "
"variant set on the properties that drive identification, so the comparison is like for like; "
"figure S4a\u2013c shows the matching \u2014 589,588 variant-pair rows across eight plexes, "
"entrapment sets of 147,397 each, all of them sharing residues with their target.",

"**(b) The reported rate holds only once a score margin is imposed.** Margin dose-response over "
"1,489 rows, counted as distinct peptides. The false discovery proportion is **4.70% at margin 0** "
"(1,419 variant peptides against 70 entrapment) and falls to **0.83% at margin \u2265 10** (832 "
"against 7). Past a margin of 15 only one entrapment peptide remains, so the rate there reflects a "
"single count rather than a cleaner set. **A 1% reported FDR is not a 1% realised FDR for this "
"class of peptide unless a margin is applied.**",

"**(c) Per-plex consistency.** Identifications at margin \u2265 10, plex by plex, on a log axis; "
"rust marks the plexes in which an entrapment peptide was identified. Six entrapment among 1,033 "
"identifications \u2014 0.58% pooled \u2014 appearing in four of 42 plexes, so the estimate is "
"not carried by one experiment. This is the identification denominator, not panel b's peptide one.",

"**(d) Spectrum quality.** Median spectral-quality features for the 38 acquired-lysine "
"identifications expressed relative to all 802,067 spectra in the same searches, so 1.0 is the "
"background and the panel reads as a multiple of it: peak count 0.91\u00d7, signal-to-noise "
"1.19\u00d7. These features score a spectrum on its own properties, with no reference to the "
"search engine, the database or the FDR. SQS$_{GM}$ is omitted because the MGF export zeroed its "
"intensity term, which would collapse a geometric mean.",

"**(e) Precursor mass error.** Observed minus calculated precursor mass, variant against "
"background: 91% within 0.05 Da, medians 0.0013 against 0.0009 Da. The tail is stated because a "
"reader will ask about it: 7.7% of variant PSMs and 7.3% of background sit at +1 Da \u2014 the "
"first isotope selected instead of the monoisotopic peak, at the same rate in both, so it is an "
"instrument behaviour and not a property of variant peptides.",

"**(f) Paired score separation.** Score-margin densities for variant and entrapment "
"identifications. Entrapment peptides *are* identified, and they score like a population that "
"cannot be present should \u2014 which is what makes the margin in panel b do work rather than "
"simply discard the weakest identifications of every kind.",

"*Pooled false discovery proportion across all plexes is 0.58% as identifications (figure S4); "
"panel b's 0.83% is the same events counted as distinct peptides.*",
]


# ================================================================ 2.2  (Figure 2)
P_2_2 = [
"If mutation-created lysines are numerous, the next question is whether anyone has been looking at "
"them. The 42 TMT plexes of the CCLE proteomes [[nusinow2020]] cover 375 cell lines and are among "
"the deepest shotgun datasets available, so they are the fair test: re-searched with a "
"mutation-aware database, how many of these sites does conventional proteomics actually see? The "
"population here is not the 295,411-site union of Figure 1 but what these particular lines carry: "
"**8,470 site-in-cell-line observations** \u2014 8,055 distinct genomic sites across the 372 of "
"those lines that carry one, a few sites assayed in more than one plex \u2014 together with 1,902 "
"observations of variant "
"sites that do not create a lysine, carried through the identical pipeline as a control "
"(figure S2). Every rate below is per observation, because each line is an independent chance to "
"detect.",

"**An acquired lysine is not one sequence with a letter changed.** Because lysine is a tryptic "
"cleavage site, creating one redraws the peptide boundaries: four species can carry the site and "
"only two of them contain it \u2014 the mutant peptide and its missed-cleavage form \u2014 while "
"the wild-type peptide and the wild-type read-through do not (Figure 2a). For a glutamate-to-lysine "
"site the mutant peptide is 303.26 Da heavier once the new lysine carries TMTpro, with different "
"hydrophobicity and a different charge state. Detection therefore fails at two separate gates, and "
"they have to be counted separately.",

"**Gate 1: the region is almost never covered.** A conventional 20,431-protein search observes the "
"wild-type region at 640 of the 8,470 observations (7.56%), and the mutation-aware search over "
"10.17 "
"million peptide entries observes it at 617 (7.28%) (Figure 2b). The two architectures agree to "
"within 0.3 percentage points, so coverage is a property of the sample rather than of the "
"database. **The remaining 7,853 sites \u2014 92.7% \u2014 are never seen in any form, in any of "
"the 42 plexes.**",

"**Gate 2: where the region is covered, the mutant peptide is still usually missed.** The mutant "
"peptide was detected at 27 of the 617 covered observations (4.38%), against 26 of the 7,853 "
"uncovered ones (0.33%) by the mutant-only route \u2014 53 detections in total, 0.63%, falling on "
"52 distinct genomic sites (Figure 2c). Two thresholds are in use here and they are not "
"interchangeable: those counts are taken at the search\u2019s reported 1% false-discovery rate, "
"and 30 of the 53 also clear the score margin of \u2265 10 at which section 2.3 measures the "
"realised error rate \u2014 16 of 617 covered (2.59%) and 14 of 7,853 uncovered (0.18%), 0.35% "
"in total on 29 distinct genomic sites. Both are correct, and neither should be quoted with the "
"other\u2019s counts. The larger number is the census of what the searches returned; the smaller "
"is the subset over which the false-discovery proportion is 0.83% rather than the 4.70% that "
"holds at margin 0. "
"Under the conventional search the mutant count is zero, which is a passed control rather than a "
"failure: that database contains no variant sequences, so it can speak to coverage and to nothing "
"else.",

"**The misses include the driver hotspots.** Arraying the 617 covered sites by allele fraction and "
"splitting them into the 590 whose mutant peptide was missed and the 27 where it was found, all ten "
"driver hotspots in that population fall on the missed row \u2014 eight distinct gene\u2013mutation "
"labels across five genes: BRAF V600K, NRAS Q61K, PIK3CA E542K, PTEN I67K, and TP53 M133K, R280K, "
"E258K and E285K (Figure 2d). Driver status is the DepMap hotspot flag rather than a gene list "
"chosen here, which matters because recurrence alone does not establish that a site is a driver "
"[[hess2019]]. Allele fraction is on the axis so a reader can see the misses are not confined to "
"subclonal sites. This panel counts hotspots **among the covered sites**, which is the population "
"it plots: NRAS Q61K was in fact recovered in one line, HCC15, by the mutant-only route \u2014 the "
"peptide QVVIDGETCLLDILDTAGK ending at the acquired K61, score margin 16.08 \u2014 in a line where "
"the wild-type region was never observed. The same lysine is covered and missed in HT1080 and "
"SKNAS.",

"**The two gates fail for different reasons, and only gate 1 has a clean explanation.** Where the "
"region was never covered the shortest tryptic peptide spanning it is longer \u2014 mean 16.7 "
"residues against 13.3 where it was covered, Mann\u2013Whitney p = 3.5 \u00d7 10\u207b\u2078 "
"(Figure 2e). The two distributions diverge in the tail beyond 20 residues rather than shifting as "
"a whole, which points at a subset of regions needing peptides too long to observe rather than a "
"general depth problem. Digestion is not the cause: acquired-lysine peptides carry missed cleavages "
"at 23.8% against a 21.6% background, so the new lysine does not create an unusually hard substrate "
"(figure S2d).",

"**The mutant peptide has the opposite geometry, and the search absorbs it rather than being "
"stopped by it.** A new lysine is a new cleavage site, so the fully cleaved mutant peptide runs "
"from the previous lysine or arginine to the new one. Its median length is six residues against "
"thirteen for the wild-type peptide over the same site, and 52.3% of the class falls below the "
"seven-residue minimum that this search and the original one both enforced (figure S3f). That is "
"not a second hidden exclusion. Both searches allowed two missed cleavages, and an upstream "
"missed cleavage returns a too-short peptide to the window, so 8,066 of the 8,470 observations "
"still carry a searchable species ending at the acquired lysine. What the geometry costs is "
"quality rather than existence: for half this class the only searchable species is one carrying "
"a missed cleavage, and internal missed cleavages are the strongest negative predictor of "
"identification in these data (figure S3d). Whether that costs detections cannot be settled "
"here \u2014 only 27 of the 617 covered observations lack a fully cleaved species, which at a "
"4% recovery rate expects one event.",

"**The non-lysine control demonstrates the mechanism rather than controlling it away.** The 1,902 "
"variant sites that do not create a lysine are covered *less* (4.6% against 7.3%) and recovered "
"*more* (6.9% against 4.4%) than the acquired lysines (figure S2a). That is the expected direction: "
"a non-lysine variant leaves the peptide boundaries intact, so its mutant peptide is the same "
"species that established coverage in the first place, while an acquired lysine produces a peptide "
"that has to be found on its own terms. Median peptide lengths are 13 and 14 residues, so the "
"difference is not a length artefact (figure S2b).",

"**A second, independent route reaches the same class of site.** The tryptic peptide released "
"downstream of the new cleavage site contains no mutated residue but exists only if the cleavage "
"happened, and it identifies 15 sites, seven of them shared with the mutant-peptide route and eight "
"new \u2014 61 sites in the union (Figure 2f). The two routes are reported separately rather than "
"pooled because their error properties differ: the downstream route returned zero entrapment hits "
"across 42 experiments, against 0.58% for variant-containing peptides in the same searches, which "
"is the subject of Figure 3.",

"**What the evidence behind these counts will bear is limited, and Figure 3 exists because of it.** "
"Of the 53 mutant-peptide sites, 41 rest on a single spectrum and 12 on two or more; four rest on "
"two or more distinct peptides, and none on more than one experiment (figure S3e). Of the 30 that "
"clear margin \u2265 10, 24 rest on a single spectrum, six on two or more and two on two or more "
"distinct peptides. Three "
"alternative explanations were tested and none accounts for the gate-1 failure: the detections are "
"not carried by subclonal sites, protein abundance separates found from never-covered sites by only "
"Cliff's \u03b4 = 0.31, and allele fraction barely differs between detected sites and the class as "
"a whole, 0.48 against 0.44 (figure S3a\u2013c). What predicts identification is incomplete "
"digestion rather than charge: internal missed cleavages carry a standardised effect of "
"\u22120.621 while the C-terminal lysine or arginine \u2014 universal to every tryptic peptide "
"\u2014 carries \u22120.004, and histidine, basic but never a cleavage site, carries \u22120.220 "
"(figure S3d).",

"This figure states what has been looked for and found, in these cell lines, with this workflow. It "
"is not evidence that the peptides are absent from the sample, and the distinction is not "
"rhetorical: gate 1 is a statement about tryptic geometry and gate 2 about search sensitivity, and "
"neither licenses a claim about what the proteome contains.",
]

CAP_2_2 = [
"42 TMT plexes across 375 cell lines. **The unit is a site-in-cell-line observation, not a "
"distinct genomic site**: 8,470 acquired-lysine observations covering 8,055 distinct sites across "
"372 lines, with 1,902 non-lysine observations carried through the identical pipeline as the "
"control (figure S2). Every rate is per observation, since each line is an independent chance to "
"detect. Detection is counted at "
"two gates that fail for different reasons and are never combined: gate 1 is whether the wild-type "
"region was observed at all, gate 2 whether the mutant peptide was recovered where it was.",

"**(a) Four peptide species can carry an acquired-lysine site, and only two contain it.** The "
"mutant peptide and its missed-cleavage form carry the site; the wild-type peptide and the "
"wild-type read-through do not. An acquired lysine creates a tryptic cleavage site, so the mutant "
"peptide has different boundaries from the wild-type peptide \u2014 the two are not one sequence "
"with a letter changed. For an E\u2192K site the mutant is heavier by 303.26 Da once the new "
"lysine carries TMTpro, with different hydrophobicity and charge state.",

"**(b) Gate 1 \u2014 the region is almost never covered.** Both bars are the same 8,470 "
"acquired-lysine observations: a conventional 20,431-protein search covers 640 (7.56%), the "
"mutation-aware peptide-level search over 10.17 million entries covers 617 (7.28%). They agree to "
"within 0.3 percentage points, so coverage is a property of the sample rather than of the "
"database. 7,853 sites \u2014 92.7% \u2014 are never seen in any form in any of the 42 plexes.",

"**(c) Gate 2 \u2014 where the region is covered, the mutant is still usually missed.** The mutant "
"peptide was detected at 27 of the 617 covered sites (4.38%) and at 26 of the 7,853 uncovered ones "
"(0.33%), by the mutant-only route \u2014 53 sites in all, 0.63% of the class, at the "
"search\u2019s reported 1% false-discovery rate. **Thirty of those 53 also clear the score "
"margin \u2265 10 of Figure 3b** \u2014 16 of 617 covered and 14 of 7,853 uncovered, 0.35% of "
"the class. The panel is drawn at the reported filter; both counts are correct and neither "
"should be quoted with the other\u2019s denominators. Under the "
"conventional search the mutant count is zero: a passed control, not a failure, since that database "
"contains no variant sequences and can speak only to coverage. Detection is by a cleaved species, a "
"missed-cleavage species, or both.",

"**(d) Driver hotspots are missed too.** The 617 covered sites arrayed by allele fraction and split "
"into the 590 whose mutant peptide was missed and the 27 where it was found. All ten driver "
"hotspots in that set fall on the missed row \u2014 eight distinct gene\u2013mutation labels "
"across five genes: BRAF V600K (IGR1, MDST8), NRAS Q61K (HT1080, SKNAS), PIK3CA E542K (CAL51), "
"PTEN I67K (LS180), and TP53 M133K (HEL), R280K (MDAMB231), E258K (MEWO) and E285K (RPMI8226). "
"Driver status is the DepMap hotspot flag, not a gene list chosen here. Allele fraction is on the "
"axis so that misses can be seen not to be confined to subclonal sites. **The panel counts hotspots "
"among the covered sites, which is the population it plots**: NRAS Q61K was recovered in HCC15 by "
"the mutant-only route, where the wild-type region was never observed, so that site is not in this "
"population \u2014 it is one of the 26 mutant-only detections in figure S3a.",

"**(e) Where the region is never covered, the wild-type tryptic peptide is the reason.** The "
"shortest tryptic peptide covering the site is longer where the region was never observed: mean "
"16.7 residues (n = 7,853) against 13.3 (n = 617) where it was covered, Mann\u2013Whitney "
"p = 3.5 \u00d7 10\u207b\u2078. The two distributions diverge in the tail beyond 20 residues "
"rather than shifting as a whole, pointing at a subset of regions requiring peptides too long to "
"observe. A gate-1 result, kept separate from gate 2 because the two failures have different "
"causes.",

"**(f) Two independent routes to the same class of site**, reported separately because their error "
"properties differ. The mutant peptide itself gives 53 sites at the reported 1% filter, 30 of "
"them at margin \u2265 10. The downstream peptide \u2014 "
"released by the new cleavage site and containing no mutated residue \u2014 gives 15, of which 7 "
"are also found by the first route and 8 are new: 61 sites in the union. The downstream route "
"returned zero entrapment hits across 42 experiments (figure S4d) against 0.58% for "
"variant-containing peptides in the same searches, which is why the two are kept apart rather than "
"pooled. The panel also carries the depth behind the counts: of the 53 mutant-peptide sites, 41 "
"rest on a single spectrum and 12 on two or more, and of the 30 at margin \u2265 10, 24 and six "
"(figure S3e). That depth is what Figure 3 "
"addresses.",

"*This figure states what has been looked for and found, in these cell lines, with this workflow. "
"It is not evidence that the peptides are absent from the sample.*",
]


# ================================================================ 2.1  (Figure 1)
P_2_1 = [
"Covalent drugs are written around the residue they modify, and for lysine the usual question is "
"which of the roughly six hundred thousand lysines already in a proteome can be reached. A "
"different question is available in cancer: how many lysines does the disease itself create? "
"Missense mutation writes new residues into proteins, and a lysine written by mutation is present "
"in tumour cells and absent from every normal cell in the same patient. Whether that is a class "
"worth building chemistry around depends first on how large it is.",

"**It is large, and it is not the product of one dataset.** Pooling missense sites from DepMap "
"[[tsherniak2017]], the TCGA MC3 call set [[ellrott2018]] and COSMIC [[tate2019]] gives 4,395,403 "
"missense sites, of which **295,411 create a lysine, across 20,365 genes** (Figure 1a). Of those, "
"93,769 (31.7%) are supported by more than one source. The three sets are nested rather than "
"independent \u2014 DepMap sits almost entirely inside COSMIC, and the region in DepMap and TCGA "
"but not COSMIC holds 147 sites \u2014 so the union is a catalogue of what has been reported, not "
"three replicate surveys.",

"**Glutamate to lysine is the commonest substitution in cancer, and that is a statement about "
"glutamate and about C>T mutation rather than about lysine.** E\u2192K ranks first of the 363 "
"substitution pairs observed, at 181,252 sites or 4.12% of all missense (Figure 1b; figure S1a). "
"But against a trinucleotide-context expectation the enrichment is only **1.12-fold**, against "
"3.59, 3.38 and 3.16 for the CpG hotspots R\u2192H, R\u2192C and R\u2192Q (figure S1b). Lysine "
"is fourth of twenty residues by sites gained, behind leucine, serine and valine, and fifth by net "
"flux \u2014 gaining 295,411 and losing 164,237 \u2014 where cysteine ranks second, netting "
"144,683 because it is so rarely lost (figure S1d, e). What does distinguish lysine is the "
"asymmetry of the routes: lysine is gained 85.5% by C>T and never by T>C, and lost 0% by C>T "
"(figure S1f). A residue that a C>T-dominated process creates and does not destroy accumulates. "
"The class is abundant for a mechanistic and unremarkable reason.",

"**The sites are somatic rather than common polymorphism.** Across all 295,411 sites, gnomAD v4.1 "
"[[karczewski2020]] records 154,188 as absent and 96,892 below a frequency of 10\u207b\u2075; "
"19,194 fall between 10\u207b\u2075 and 10\u207b\u2074, 3,297 between 10\u207b\u2074 and "
"10\u207b\u00b3, and 1,971 above 10\u207b\u00b3, with 19,869 never assessed (Figure 1c). "
"**91.12% of assessed sites are absent or below 10\u207b\u2075.** The 6.73% never assessed is "
"counted separately rather than folded into a denominator that would flatter the result. Frequency "
"cannot exclude a private germline variant; only matched normals can, and these data have none.",

"**Acquired lysines fall on driver hotspots more often than other substitutions do.** Among DepMap "
"variants, 1.024% of 19,239 acquired lysines carry the Hess driver-hotspot flag against 0.6138% of "
"311,166 other substitutions \u2014 odds ratio **1.68**, Fisher exact p = 1.3 \u00d7 "
"10\u207b\u00b9\u2070 (Figure 1d). The flag is the one carried in the variant table rather than "
"a gene list chosen here, and it is conservative by construction: the method exists to strip "
"recurrent passenger sites from hotspot lists [[hess2019]], so it should if anything understate "
"the enrichment. Hotspot recurrence is a statement about where mutation is selected, not about "
"whether the resulting protein is expressed or reachable.",

"**The excess tracks mutational processes that make lysine, not mutation in general.** Splitting "
"driver hotspots by their assigned signature, acquired lysines (n = 315) against other "
"substitutions (n = 3,465), the four processes that generate the relevant base changes are all "
"enriched: ultraviolet 23.5% against 9.6%, microsatellite instability 20.6% against 6.3%, POLE "
"15.2% against 6.0%, APOBEC 8.3% against 4.3% (Figure 1e). The two that do not make lysine run the "
"other way \u2014 smoking 18.7% against 26.2%, CpG deamination 10.5% against 28.2%. Both "
"directions are labelled.",

"**Burden is lineage-specific, and the ranking carries its own control.** Cutaneous melanoma "
"carries a median of 39 acquired lysines per tumour, thirteen times the pan-TCGA median, and it "
"leads all 32 lineages (Figure 1f). **Uveal melanoma carries zero.** Both are melanocytic and only "
"the first is ultraviolet-driven, so if the excess were a property of the lineage rather than of "
"the mutational process both would be high. All 32 lineages are shown because a truncated ranking "
"cannot carry that control.",

"**The excess survives a germline reference chosen to remove the obvious objection.** Acquired lysines are 6.72% of "
"all missense in the somatic union (295,411 of 4,395,403) against **3.29%** in the germline, a "
"**2.04-fold** difference (Figure 1g). The germline figure is an expected diploid count per "
"individual \u2014 471.1 acquired lysines among 14,342 missense variants, summed as twice the "
"allele frequency over dbNSFP rows carrying a MANE annotation, using gnomAD v4.1 joint frequencies "
"from about 807,000 individuals \u2014 and not a count of distinct variants. The obvious objection "
"is that gnomAD contains cancer-affected individuals, so the same analysis was run over the gnomAD "
"v2.1.1 **non-cancer** exome subset (about 118,000 individuals): 508.0 acquired lysines per "
"individual, **3.41%**. The two agree to within 0.12 percentage points. They are not independent "
"cohorts \u2014 both are gnomAD, with heavy overlap \u2014 so what the second controls is "
"specifically the presence of cancer-affected individuals, not reference bias in general.",

"**2.04-fold means little without the other nineteen residues, so they were measured too.** The "
"same pass over dbNSFP was repeated for every residue rather than for lysine alone, on the "
"identical file, filter and frequency field; the lysine value it returns is 471.1 per individual, "
"which is the figure above. Ranked by the ratio of somatic to germline share, lysine is third of "
"twenty at 2.05 and **cysteine fourth at 1.76**, in a distribution running from 0.43 to 2.63 with "
"a median of 1.07 (figure S1h). Two things follow. The comparison the field draws with acquired "
"cysteine can now be made on one calculation instead of two, and on that calculation lysine is "
"not the weaker case. And the ratio is not a measure of selection: germline shares are depressed "
"by purifying selection, so a residue that is disruptive to introduce scores high here for a "
"reason that has nothing to do with cancer \u2014 which is the likeliest reading of phenylalanine, "
"tyrosine and tryptophan above lysine. The selection evidence is the driver-hotspot odds ratio "
"below, and the context-corrected enrichment for E\u2192K is 1.12-fold (figure S1b).",

"Together these exclude the cheapest alternative to the premise of this paper \u2014 that acquired "
"lysines are simply what mutation does everywhere. They are created at twice the germline rate, "
"they are enriched on driver hotspots, and the processes that make them are identifiable. None of "
"that says the resulting proteins are expressed, detectable or reachable by a molecule, and the "
"rest of the paper is about how far each of those holds.",
]

CAP_2_1 = [
"Sites are missense positions in the union of three variant sources; a site is an acquired lysine "
"if the substituted residue is lysine. **Panels a\u2013c are union-wide (DepMap, TCGA and "
"COSMIC), panels d and e are DepMap only** \u2014 because the driver and signature flags are "
"carried in that table \u2014 **and panel f is TCGA.** Panel g compares the union against two "
"germline references.",

"**(a) The class, assembled from three variant sources.** 295,411 acquired-lysine sites across "
"20,365 genes, from 4,395,403 missense sites. 93,769 sites (31.7%) are supported by more than one "
"source. Circle areas and overlaps are both computed from the seven region counts rather than "
"positioned by hand, so the structure shown is the data's: DepMap sits almost entirely inside "
"COSMIC, and the DepMap-with-TCGA-but-not-COSMIC region is 147 sites. Region fills carry no "
"meaning \u2014 set identity is in the outlines.",

"**(b) Glutamate-to-lysine is the single commonest substitution, and lysine is the fourth "
"most-gained residue.** Codon routes into and out of lysine. E\u2192K ranks first of 363 observed "
"substitution pairs; 380 are possible (20 \u00d7 19), and the 17 never seen all require more than "
"one nucleotide change. Three things the panel cannot show: the full 20 \u00d7 20 ranking is "
"figure S1a; lysine is **not** uniquely protected from C>T loss \u2014 seven residues (I, F, K, N, "
"Q, W, Y) have no C>T exit route; and figure S1b is the conservative reading of the "
"85.5%-versus-51.1% contrast, because against a trinucleotide-context expectation E\u2192K is "
"enriched only 1.12-fold. The abundance of this class follows from how common glutamate is and "
"from a C>T-dominated mutational process, not from an unusual per-site enrichment.",

"**(c) The sites are somatic, not common polymorphism.** gnomAD v4.1 frequency bands over all "
"295,411 sites: absent 154,188; below 10\u207b\u2075 96,892; 10\u207b\u2075\u201310\u207b\u2074 "
"19,194; 10\u207b\u2074\u201310\u207b\u00b3 3,297; above 10\u207b\u00b3 1,971; not assessed "
"19,869. 91.12% of assessed sites are absent from gnomAD or below 10\u207b\u2075. The 6.73% never "
"assessed is counted separately rather than folded into a denominator. Frequency cannot exclude a "
"private germline variant; only matched normals can, and these data have none.",

"**(d) Acquired lysines fall on driver hotspots more often than other substitutions do.** 1.024% "
"of 19,239 acquired lysines are Hess driver hotspots, against 0.6138% of 311,166 other "
"substitutions \u2014 odds ratio 1.68, Fisher exact p = 1.3 \u00d7 10\u207b\u00b9\u2070. "
"Driver status is the flag carried in the DepMap variant table, not a gene list chosen for this "
"analysis, and the method it comes from exists to remove recurrent passenger sites, so it is "
"conservative. Hotspot recurrence is a statement about where mutation is selected, not about "
"whether the resulting protein is expressed or reachable \u2014 Figure 2 addresses the first of "
"those.",

"**(e) The excess tracks mutational processes that make lysine, not mutation in general.** Share "
"of driver hotspots by signature, acquired lysines (n = 315) against other substitutions "
"(n = 3,465): UV 23.5% vs 9.6%, MSI 20.6% vs 6.3%, POLE 15.2% vs 6.0%, APOBEC 8.3% vs 4.3%. The "
"two signatures that do not make lysine run the other way \u2014 smoking 18.7% vs 26.2%, CpG "
"10.5% vs 28.2%. Enrichments and depletions are both labelled.",

"**(f) Burden across every TCGA lineage, with the control the ranking exists for.** Median "
"acquired lysines per tumour across 32 lineages. Cutaneous melanoma carries 39, thirteen times the "
"pan-TCGA median. **Uveal melanoma carries zero.** Both are melanocytic; only the first is "
"ultraviolet-driven. If the excess were intrinsic to the lineage rather than to the mutational "
"process, both would be high. All 32 lineages are shown because a truncated list cannot carry that "
"control.",

"**(g) Somatic mutation creates acquired lysines at twice the germline rate, on two gnomAD "
"references.** Acquired lysines as a share of all missense: somatic union 6.72% (295,411 "
"of 4,395,403) against germline 3.29% \u2014 2.04-fold. The germline figure is an expected diploid "
"count per individual, the sum of 2 \u00d7 allele frequency over dbNSFP 5.4a rows carrying a MANE "
"annotation and an unambiguous single-letter amino-acid change, using gnomAD v4.1 joint allele "
"frequencies (~807,000 individuals): 471.1 acquired lysines against 14,342 missense variants per "
"individual. It is not a count of distinct variants. The control is a second, independent germline "
"reference: the same analysis over the gnomAD v2.1.1 non-cancer exome subset (~118,000 "
"individuals) gives 508.0 per individual, 3.41%. The two agree to within 0.12 percentage points, "
"so the germline rate is not an artefact of cancer-affected individuals being present in the "
"reference cohort. The two cohorts overlap heavily and are not independent samples of human "
"variation; that specific objection is what this control removes.",

"*Counts are of sites, not of tumours or of individuals, except in panels f and g where the unit "
"is stated.*",
]


# ================================================================ 1. Introduction
P_INTRO = [
"Covalent chemistry directed at lysine has reached the point where the limiting question is no "
"longer whether a lysine can be engaged but which one to engage. Proteome-wide profiling "
"established that lysine reactivity and ligandability vary enormously between sites "
"[[hacker2017]], and three subsequent screens have mapped thousands of ligandable lysines across "
"several cell lines with chemically distinct probes [[abbasov2021,bracken2024,ryan2024]]. What "
"none of that work supplies is a principled way to choose a site. Selectivity in a covalent "
"campaign is usually argued from the chemistry of the molecule; the alternative is to choose a "
"residue that exists only in the diseased cell.",

"Cancer supplies such residues. Missense mutation writes new amino acids into proteins, and a "
"residue created that way is present in tumour cells and absent from every normal cell in the "
"same patient \u2014 which is the logic that made KRAS G12C tractable, where the drug depends on "
"a cysteine the mutation created and cannot touch the wild-type protein [[ostrem2013]]. That "
"argument has since been generalised for cysteine: acquired cysteines have been catalogued across "
"cancer genomes, shown to be widespread, and pursued chemoproteomically as a target class "
"[[desai2024]]. Lysine has not been counted the same way, although mutation creates more of them: "
"the union of three variant catalogues holds 295,411 acquired lysines against 211,386 acquired "
"cysteines, and glutamate-to-lysine is the single commonest substitution in cancer. On net flux "
"the order reverses, because cysteine is rarely lost and lysine is not; neither residue is "
"exceptional, and the point is only that one of the two has been catalogued as a target class "
"and the other has not.",

"That E\u2192K is common is not a discovery of this work, and the reason it is common is "
"unremarkable \u2014 glutamate is abundant and the dominant mutational process in cancer is C>T. "
"Against a trinucleotide-context expectation the substitution is enriched only 1.12-fold. The "
"question that matters is not whether these lysines are numerous but whether they are worth "
"targeting, and the obvious counter-position is that they are passengers: residues that mutation "
"produces incidentally, on proteins under no selection, carrying no more therapeutic promise than "
"any other missense change. Figure 1d addresses that directly by asking whether acquired lysines "
"fall on driver hotspots more often than other substitutions do.",

"A second question is why nobody has counted them. The answer is a measurement problem with two "
"parts, and both are specific to this class rather than general to proteomics. First, lysine is a "
"tryptic cleavage site, so creating one redraws the peptide boundaries: the mutant peptide is not "
"the wild-type peptide with a letter changed but a different species with a different mass, "
"length and charge, which a reference database cannot contain and a conventional search cannot "
"find. Second, the error control that would license the identifications it does make has mostly "
"been tested in the wrong place. Target-decoy false-discovery rates are calibrated over a whole "
"search [[elias2007]], and variant peptides are a small, atypical subset of one. Two-stage search "
"strategies that control error separately for the variant space have been developed "
"[[desai2024]]; atlases of mutation-created residues have instead filtered a combined "
"canonical-and-mutant search at a global 1% [[takahashi2024]]; and entrapment methods that "
"*measure* the realised rate within a subset now "
"exist [[wen2025]]. What is rarely done is to turn the entrapment measurement on a variant "
"search and report what comes back. A class that is hard to detect and whose detections are hard to trust will not be "
"catalogued by accident.",

"This paper is a reanalysis. It generates no new data: it re-searches 42 published TMT plexes of "
"the CCLE proteomes [[nusinow2020]] with a mutation-aware database, re-examines three published "
"chemoproteomic screens, and asks of both what can be said with the error rates measured rather "
"than assumed. The claims it makes are correspondingly bounded. It counts the class and tests "
"whether it is somatic and driver-enriched; it measures how much of it conventional proteomics "
"can see; it measures the realised error rate for variant peptides instead of quoting the "
"nominal one; and it asks what distinguishes an engaged lysine from a labelled one, and how well "
"that distinction can be predicted from structure. It does not report a new probe, a new "
"selectivity claim, or a disease model. Every script and intermediate is published, so each "
"number can be regenerated rather than taken on trust.",
]


# ================================================================ Abstract
# Seven sentences, each bound to a named panel. The binding is a drafting discipline, not text:
# BMC abstracts carry no figure citations. In order:
#   1 Fig 1a, g   2 Fig 1c, d   3 Fig 1e, f   4 Fig 2b, c   5 Fig 3b
#   6 Fig 5a, c and Fig 6       7 Fig 7b, e       8 Fig 8c, e
# Sentence 7 (the reagent) was added 7 September. The order of 6 and 7 is deliberate: 6 says no
# compound descriptor tested predicts preference, 7 names the one variable that does. Paid for by
# compression elsewhere, not by cutting a control; the build measures the word count.
P_ABSTRACT = [
"Missense mutation in cancer creates lysines present in tumour cells and absent from every "
"normal cell in the same patient: 295,411 sites across 20,365 genes, 6.72% of all missense against "
"3.29% in the germline, on a reference and on its cancer-free subset. The sites are somatic \u2014 91.12% of "
"assessed sites are absent from gnomAD or below 10\u207b\u2075 in frequency \u2014 and they fall on "
"driver hotspots more often than other substitutions (odds ratio 1.68, "
"p = 1.3 \u00d7 10\u207b\u00b9\u2070). The excess tracks mutational processes that make lysine "
"rather than mutation in general: cutaneous melanoma carries a median of 39 per tumour while uveal "
"melanoma, melanocytic but not ultraviolet-driven, carries none. Re-searching 42 published "
"proteomic experiments with a mutation-aware database recovers the wild-type region in 7.28% of "
"cases and the mutant peptide in 0.63%. Burden-matched entrapment puts the realised false-discovery "
"proportion for variant peptides at 4.70% where 1% was reported, falling to 0.83% once a score "
"margin is imposed, at which margin the mutant-peptide recovery is 0.35%. Across three "
"chemoproteomic screens engagement is a property of the lysine "
"rather than of any compound descriptor tested, and structure predicts it poorly: a site-feature "
"composite reaches AUC 0.617 against a reproducibility ceiling of at least 0.72. The reagent is "
"the exception: two probes differing only in the acyl leaving group separate on two of twenty-two "
"site features, and both are charge. Acquired lysines "
"are structurally ordinary, so the barrier is detection rather than chemistry \u2014 and one of "
"them, AKT1 E17K, has been drugged elsewhere while yielding no spectra in the single line here "
"that carries it.",
]


# ================================================================ 3. Discussion
P_DISCUSSION = [
"Mutation-created lysines are numerous, largely absent from common germline variation, enriched on "
"driver hotspots, written by identifiable mutational processes, and almost entirely missing from "
"the proteomic record. Each of those is measured: 295,411 distinct genomic sites across 20,365 "
"genes, at twice the germline rate (Figure 1a, g); 91.12% of assessed sites absent from gnomAD or "
"below a frequency of 10\u207b\u2075 (Figure 1c); an odds ratio of 1.68 on driver hotspots "
"(Figure 1d); ultraviolet, microsatellite-instability, POLE and APOBEC signatures enriched among "
"those hotspots while smoking and CpG deamination run the other way (Figure 1e); and, across 42 "
"TMT plexes of the published CCLE proteomes re-searched with a mutation-aware database, the "
"wild-type region observed in 7.28% of 8,470 site-in-cell-line observations and the mutant peptide "
"in 0.63%, or 0.35% at margin \u2265 10 (Figure 2b, c) \u2014 a different unit from Figure 1's "
"distinct-site count. The error "
"rate that would license even those few identifications does not hold at the reported filter, and "
"a score margin is what restores it (Figure 3b).",

"**The closest prior work approached this from two directions, and neither reaches the question "
"asked here.** The three chemoproteomic screens reanalysed in this paper mapped ligandable lysines "
"across several proteomes [[abbasov2021,bracken2024,ryan2024]], building on the reactivity and "
"ligandability landscape established earlier [[hacker2017]] \u2014 but they searched reference "
"databases, so a mutation-created lysine was outside their search space by construction. From the "
"other direction, AKT1 E17K has been engaged selectively in cells through lysine targeting and "
"neo-zinc chelation [[craven2025]]: an acquired lysine drugged by target-based design, one site at "
"a time, without any claim about the class it belongs to. The measurability question has been "
"asked and addressed for cysteine, where a purpose-built chemoproteogenomic pipeline with "
"variant-aware databases and staged error control catalogued acquired cysteines directly "
"[[desai2024]]. Lysine has had no equivalent, and it is the harder case for a reason that is "
"specific rather than incidental: an acquired cysteine leaves the tryptic peptide it sits in "
"intact, while an acquired lysine creates a cleavage site and produces a different peptide "
"altogether. The variant-aware database that solves the cysteine problem does not, on its own, "
"solve this one.",

"**The obvious next experiment is a re-search of those three screens, and it is bounded before it "
"starts.** The four cell lines they ran in carry 56 unique acquired lysines between them, against "
"a median of 10 and a maximum of 623 across 1,954 DepMap models: MDA-MB-231 24, MDA-MB-436 19, "
"HeLa seven and Ramos six (figure S10b). Eleven of the 56 fall on a protein the matching screen "
"probe-labels at one or more lysines, and one of them, SF3A3, already yields a lysine engaged "
"above 60% (figure S10e). The lysine measured is not the acquired one and usually sits far from "
"it, a median of 271 residues away, so what the eleven establish is that the protein is reachable "
"by a probe rather than that the acquired site is. One of the eleven is TP53 R280K in "
"MDA-MB-231, whose structure is measured in Figure 8d. The re-search is worth running. It is not "
"this paper: it would test eleven sites, not a class.",

"**Which cells carry the class is not arbitrary, and it can be chosen.** Burden by lineage runs "
"skin, uterus, bowel and lung at the top and eye second from last; at subtype resolution "
"cutaneous melanoma carries a median of 37 against ocular melanoma\u2019s three, which is the "
"cutaneous-versus-uveal contrast of Figure 1e recovered in cell lines rather than tumours "
"(figure S10a). The asymmetry that makes the class is visible one cell line at a time: across "
"1,954 models 39,390 lysines are gained against 25,531 lost, and the ratio runs from 6.4 in skin "
"to 1.3 in oesophagus and stomach, which is figure S1f\u2019s C>T asymmetry at cell-line "
"resolution (figure S10b). A quarter of models are net-negative. Burden is not the same as "
"opportunity either. Passing an acquired lysine through three "
"gates \u2014 the gene expressed in that line, the protein quantified in that line\u2019s CCLE "
"proteome, and the gene a CRISPR dependency \u2014 leaves 46%, then 31%, then 6% of 2,236 sites "
"(figure S10c). Six percent is not a small number in absolute terms: 132 sites in 27 cell lines "
"are on a protein that is present, detectable and required. The gates also expose a gap in the "
"reference data: only 308 of 1,954 models carry expression, proteomic and CRISPR data "
"together, and none of the ten highest net-gain lines does. The cells that carry most of the "
"class are the least characterised.",

"**The sites that recur are the sites this paper already measured.** Of the 39,390 unique "
"acquired lysines, 248 carry a hotspot or driver flag, and ranked by how many cell lines carry "
"them the leaders are NRAS Q61K in 36, PIK3CA E545K in 33, PIK3CA E542K in 15, XPO1 E571K in "
"nine and AKT1 E17K in five (figure S10d). Four of the six sites whose structures Figure 8d "
"measures are on that list, and they are the four most recurrent \u2014 which was not how they "
"were chosen. BRAF V600K, in three models, is an acquired lysine on the most-drugged oncogene in "
"melanoma. Five driver-flagged sites clear all three gates in the models examined: MTOR E1799K, "
"PIK3CA E545K, BRAF V600K in two lines, and PIK3CD E1021K. PIK3CD already "
"carries an approved drug against another site. What is independently known about these sites is "
"collected in table S1, which is provenance rather than a survey. Twelve of the thirteen are "
"carried by TCGA patients as well as by cell lines \u2014 PIK3CA E545K by 263 of the 10,224 "
"participants in the MC3 call set, E542K by 169, NRAS Q61K by 68, AKT1 E17K by 54 \u2014 and "
"nine carry a functional consequence measured by someone else, among them impaired GTP "
"hydrolysis at NRAS Q61K and accelerated leukemogenesis at XPO1 E571K. The thirteenth is "
"instructive rather than absent: BRAF V600K returns no MC3 patients because valine reaches "
"lysine only through two adjacent nucleotide changes and that call set does not merge adjacent "
"variants, so the site is invisible to it by construction while DepMap reports it in three "
"lines. These are the prospective screen\u2019s first targets, and naming them, with the gaps "
"marked, is what a catalogue is for.",

"Two things follow, and the first corrects the account just given of why lysine is the harder "
"case. An acquired lysine creates a tryptic cleavage site in a shotgun experiment but not in a "
"lysine-chemoproteomic one, because the probe acylates the \u03b5-amine and trypsin does not "
"cleave after an acylated lysine \u2014 so a labelled acquired lysine leaves the peptide "
"boundaries near wild-type while the peptide is streptavidin-enriched. The obstacle that makes "
"this class harder than acquired cysteine is partly dissolved by the experiment that would be "
"used to look for it. The second is a choice of cells. Lysine chemoproteomics has been run in "
"lines that barely carry the class, while the chemoproteogenomic work on acquired cysteine ran in "
"high-burden lines [[desai2024]]. Ranked on the sites a screen could actually reach \u2014 "
"acquired lysines on a gene the line expresses \u2014 the leaders are MAMEL46 with 248, "
"MCC13 with 178 and HEC251 with 176, against three in HeLa. A prospective screen in any of "
"them would put the same chemistry against roughly sixty times the reachable sites.",

"**The passenger reading of this class deserves a direct answer.** "
"Glutamate-to-lysine is the commonest substitution in cancer for reasons that carry no biological "
"weight \u2014 glutamate is abundant, the dominant process is C>T, and against a "
"trinucleotide-context expectation the substitution is enriched only 1.12-fold (figure S1b). It "
"would be reasonable to conclude from that alone that these are incidental residues. One "
"observation argues otherwise: acquired lysines carry the driver-hotspot flag at an odds ratio of "
"1.68 against other substitutions, and that flag is conservative by construction, since the method "
"behind it exists to strip recurrent passenger sites from hotspot lists [[hess2019]]. The "
"signature composition of those hotspots (Figure 1e) is **not** independent support and should not "
"be offered as such: it compares lysine-creating hotspots against other hotspots, and since lysine "
"is made 85.5% by C>T and never by T>C (figure S1f), any lysine-creating set would show the same "
"enrichment whether or not it were under selection. The test that would separate the two "
"\u2014 acquired-lysine hotspots against acquired-lysine non-hotspots \u2014 is not in this "
"paper. What the odds ratio establishes is that these sites are selected more often than other "
"substitutions; it does not establish that the affected proteins are expressed, detectable, or "
"reachable by a molecule, and the rest of the paper is about how far those three hold.",

"**For probe design, Figure 7 yields a principle rather than a recommendation.** Two probes sharing "
"a desthiobiotin-hexanoyl core and differing only in the acyl leaving group label overlapping but "
"distinct sets of lysines in the same proteome, and of twenty-two site features only two separate "
"them \u2014 both charge. A chemoproteomic probe is therefore not a neutral window onto lysine "
"reactivity: part of what any lysine-directed experiment measures is where its reagent is drawn. "
"Coverage and reactivity are entangled at the point of capture, which makes the probe a sampling "
"variable that should be varied deliberately rather than held fixed by convention. The paper stops "
"there, and deliberately: establishing which panel of probes would maximise coverage requires "
"prospective experiments a reanalysis cannot contain, and a two-probe comparison supports the "
"principle without licensing a design. One result does travel with the site rather than the "
"reagent \u2014 labelled lysines sit in more enclosed environments than unlabelled ones, and the "
"two probes do not differ on that \u2014 though enclosure and solvent accessibility share about "
"three-quarters of their variance, so it is largely, though not wholly, a burial effect under "
"another name, and it holds only within protein: pooling across sites reverses its sign "
"(figure S8b, c).",

"**For measurement, three things follow, and the most obvious of them is wrong.** The first is that "
"a score margin is not optional for this peptide class: at the reported 1% filter the realised "
"false-discovery proportion is 4.70%, and only a margin brings it to 0.83% (Figure 3b). The second "
"is that entrapment must be burden-matched to be informative \u2014 an entrapment sequence built "
"from its target's own residues, sharing an exact precursor mass, tests what a mass-matched "
"impostor could have produced, which is the question a decoy does not answer; on plex 2 decoys "
"underestimate the realised rate by 1.7-fold (figure S4f). The third concerns what would actually "
"recover this class, and here the intuitive answer does not survive the data. Deeper or "
"enrichment-based proteomics addresses abundance, but abundance does not look like the limiting "
"gate: it separates sites whose mutant peptide was found from those whose region was never covered "
"by Cliff's \u03b4 = 0.31, on the 1,769 of 8,470 observations that carry an abundance value and "
"only 40 found sites (figure S3b). A subset of the failures is instead associated with long "
"tryptic peptides, the two distributions diverging in the tail beyond twenty residues rather than "
"shifting as a whole (Figure 2e); no effect size is available for that comparison, so neither test "
"is decisive. What can be said is that changing which proteins are sampled leaves the peptide "
"boundaries untouched, and that three measured results point at the same lever. Half this class "
"has no fully cleaved mutant peptide above the seven-residue floor the searches enforce, so its "
"only searchable species is one carrying a missed cleavage (figure S3f). An internal missed "
"cleavage is the strongest negative predictor of identification measured here, at a standardised "
"effect of \u22120.621 against \u22120.004 for the C-terminal lysine or arginine every tryptic "
"peptide carries (figure S3d). And the cleavage that fails is specifically at arginine: 29.4% of "
"identified peptides whose first basic residue is arginine carry a missed cleavage against 15.0% "
"for lysine (figure S2c). What those three identify is completeness at the lysine and arginine "
"sites trypsin already targets, not a different specificity \u2014 a protease cutting at "
"additional residues would shorten peptides further, and short is where this class already "
"fails. Digesting the same proteome more completely is a prediction this paper does not test.",

"**The limitations are specific, and several of them bound results this paper otherwise leans on.** "
"Everything measured here is measured in cell lines \u2014 375 in the proteomic re-search, four "
"cell lines across the three chemistry screens \u2014 and no result in this paper is a statement "
"about a tumour or a patient. The screens are not saturated, so every ligandability figure is a "
"lower bound: a site recorded as unengaged was unengaged by the compounds tested. The structural "
"separation of reactivity from ligandability emerges progressively as the confidence gate tightens "
"and is strongest above pLDDT 90 (Figure 4e), which makes it a claim about ordered protein, and "
"AlphaFold confidence is not random with respect to biology. The compound descriptors tested "
"encode neither charge nor chirality, so the finding that none of them predicts which lysines a "
"compound prefers is silent on exactly the two axes Figures 7 and 5f show to matter. Figure 8's "
"nine features describe a site's surroundings within 8 to 10 \u00c5; the lysine's own N\u03b6 "
"chemistry is not compared, and at most acquired positions that N\u03b6 is a built all-trans "
"rotamer carrying a median error of 4.26 \u00c5 against a reference that is itself predicted "
"(figure S9e) \u2014 which is why the six-structure check matters, and that check crosses a "
"method boundary at two PI3K\u03b1 pairs and rests on a single wild-type pair at NRAS. The "
"cross-dataset version of the breadth result is the weakest evidence in the paper \u2014 25 sites "
"in 21 dimensions, with two of six comparisons short of their null (figure S6g) \u2014 and is "
"reported as a limit rather than as support. And the detections themselves are thin: of 53, none "
"rests on more than one experiment (figure S3e), which is why Figure 3 exists and why Figure 2's "
"counts are a lower bound on coverage rather than a catalogue of sites.",

"Taken together the results describe a class that is large, selected more often than other "
"substitutions, chemically unremarkable and effectively unmeasured, and they locate the obstacle "
"precisely. It is not that acquired lysines are structurally unusual: nine environment features "
"place them within the ordinary range of the lysines chemoproteomics already engages (Figure 8b, "
"c). Nor is the chemistry unimaginable \u2014 one such lysine has been engaged selectively, by a "
"bespoke campaign against a single site [[craven2025]], which establishes possibility there rather "
"than tractability for the class. It is that the same site returns zero spectra in the one plex "
"of the 42 whose cell line carries it, and that a class which cannot be seen cannot be prioritised or argued about "
"with data. Screening does not require proteomic visibility, as the three reanalysed screens show; "
"knowing which sites exist, in which cells, and with what confidence does. The measurement problem "
"is the tractable one.",
]


# ================================================================ 4. Conclusion
P_CONCLUSION = [
"Mutation writes 295,411 new lysines into cancer proteomes, and the class carries the "
"driver-hotspot flag more often than other substitutions do. It is nonetheless close to "
"unmeasured: a mutation-aware re-search of 42 published TMT plexes recovers the mutant peptide in "
"0.63% of site-in-cell-line observations, and the realised error rate on those recoveries does not "
"match the reported one until a score margin is imposed. The obstacle is not the chemistry "
"\u2014 acquired lysines sit within the ordinary range of the lysines chemoproteomics already "
"engages, and one has been engaged selectively by a campaign built for that site. What is missing "
"is the ability to see them. One of the two things that would change that is measured here: error "
"control within the variant subset rather than assumed from the whole search. The other, "
"digesting the same proteome more completely at the lysine and arginine sites trypsin already "
"targets rather than cutting at different residues, is a prediction this paper does not test.",
]


# ================================================================ 5. Methods
# Subheadings are carried as bold lead-ins inside the paragraph, matching the house style used
# for the Results. No new analysis is described here; everything is provenance and parameters.
P_METHODS = [
"**Variant catalogues and the union.** Missense sites were pooled from three sources: DepMap "
"release 26Q1 [[tsherniak2017]], the TCGA MC3 controlled call set [[ellrott2018]], and COSMIC "
"Cancer Mutation Census v103 [[tate2019]]. The union holds 4,395,403 distinct missense sites, "
"keyed on gene, position, reference residue and alternate residue, of which 295,411 substitute a "
"lysine. A site is counted once however many sources report it, and 93,769 are reported by more "
"than one. Two TCGA cohorts, ovarian serous cystadenocarcinoma and acute myeloid leukaemia, were "
"generated from whole-genome-amplified material; a PASS-only filter removes them entirely, so "
"filtering was applied without that restriction and the two cohorts are retained. Driver-hotspot "
"and mutational-signature flags are those carried in the DepMap variant table, derived by the "
"method of Hess et al. [[hess2019]]; no gene list was chosen for this work.",

"**Germline comparison.** Germline rates were computed from dbNSFP 5.4a [[liu2020]], restricted to "
"rows carrying a MANE annotation and an unambiguous single-letter amino-acid change. The reported "
"figure is an expected diploid count per individual \u2014 the sum of twice the allele frequency "
"over qualifying rows \u2014 and is not a count of distinct variants. Allele frequencies are "
"gnomAD v4.1 joint frequencies, approximately 807,000 individuals [[karczewski2020]]; 91,571,221 "
"rows were scanned and 14,084,288 carried a frequency above zero. The control repeats the same "
"calculation over the gnomAD v2.1.1 non-cancer exome subset, approximately 118,000 individuals. "
"The two cohorts overlap heavily and are not independent samples of human variation; the control "
"exists to remove one specific objection, that cancer-affected individuals inflate the germline "
"baseline.",

"**Proteomic reanalysis.** The 42 sample plexes of the CCLE proteomes [[nusinow2020]] \u2014 nine "
"cell lines per plex and a pooled bridge in the tenth channel, 378 samples in all, three lines run "
"twice, so 375 distinct cell lines \u2014 were re-searched with MSFragger 4.4.1 within FragPipe "
"v24.0, with Percolator and "
"Philosopher. Two architectures were run. The conventional search used a 20,431-protein MANE v1.2 "
"database with contaminants and reversed decoys. The mutation-aware search used a peptide-level "
"database of 10.17 million entries built from the variant union. Detection was scored at two "
"gates that are never combined: gate 1 asks whether the wild-type region was observed at all, "
"gate 2 whether the mutant peptide was recovered where it was. **The unit throughout Figures 2 and "
"3 is a site-in-cell-line observation, not a distinct genomic site**: 8,470 acquired-lysine "
"observations cover 8,055 distinct sites across the 372 lines carrying one, in 10,372 rows, with "
"1,902 non-lysine variant observations carried through the identical pipeline as a control. "
"Figure 1's 295,411 is a distinct-site count and is a different denominator.",

"**Entrapment and error control.** Realised error rates for variant peptides were measured by "
"burden-matched entrapment at r = 1 using FDRBench v1.1.1 [[wen2025]]. Each variant peptide was "
"paired with an entrapment sequence built by permuting its own residues with the C-terminal lysine "
"or arginine held in place, so the pair shares an exact precursor mass and an identical length and "
"only fragment ions distinguish them; 589,588 variant-pair rows across eight plexes yielded "
"entrapment sets of 147,397 each, all sharing residues with their target. Because every entrapment "
"identification is false by construction, the proportion of entrapment hits estimates the realised "
"false-discovery proportion directly, without the distributional assumptions target-decoy control "
"requires [[elias2007]]. The score margin is the difference between a spectrum's best and "
"second-best interpretation; margin \u2265 10 is used throughout. **Two denominators are in use "
"and are not interchangeable:** Figure 3b counts distinct peptides, one best spectrum each, and "
"gives 0.83% at that margin; Figure 3c and figure S4e count identifications and give 0.58%.",

"**Site assignment.** Chemoproteomic site assignments in Figures 4 to 7 were reconciled between "
"two independent routes, which concur at 97.82% (44,736 of 45,552 usable sites). Assignments are "
"described as provisional throughout: the figures compare site *sets*, and an individual residue "
"number within a set may be wrong even where the set is right. Figure 8 is the exception, since "
"its sites come from variant calling rather than from chemoproteomic assignment.",

"**Occupancy and engagement.** Occupancy is the share of probe labelling a compound competes away, "
"100 \u00d7 (1 \u2212 compound / DMSO). A site is engaged at \u2265 60% occupancy. Where source "
"screens report competition ratios rather than occupancies, a ratio of 4 corresponds to 75% "
"engagement; the 60% cut is a reporting convention shared with those screens rather than a "
"boundary in the data, and Figures 4f, 6a and 6c together with figure S7c show the results do not "
"depend on it.",

"**Structural features.** Site features were measured on AlphaFold Monomer v2.0 models retrieved "
"from the AlphaFold Protein Structure Database [[jumper2021,varadi2022]], files tried at v6, v5 "
"and v4. Positions are anchored at C\u03b2 unless stated; where an N\u03b6 is required at a "
"position that carries no lysine in the model it is built as an all-trans rotamer, which "
"reproduces measured amine positions to a median 4.26 \u00c5 and is an upper bound rather than a "
"placement. Relative solvent accessibility follows Tien et al. [[tien2013]], predicted pKa is from "
"PROPKA3 [[olsson2011]] and pocket descriptors from fpocket [[leguilloux2009]]. **Two feature sets "
"are used and must not be conflated.** Figure 4 pairs 15 features across both contrasts \u2014 6 "
"charge and 9 shape \u2014 and excludes the AlphaFold terms, because confidence is the gating "
"variable there rather than a site property. Figure 5's composite is 14 features: 13 environment "
"features plus model confidence, and the composite is reported both with it (AUC 0.617) and "
"without it (0.611). **Two gate rungs are also in use:** Figures 4b and 4c display pLDDT \u2265 "
"90, while Figure 8 is gated at \u2265 70.",

"**Statistics.** Models were fitted with scikit-learn [[pedregosa2011]]. The composite classifier "
"is logistic regression cross-validated by protein rather than by site, since lysines on one "
"protein are correlated; intervals on model AUCs are the 2.5th to 97.5th percentile across 20 "
"protein-held-out splits, and intervals on reproducibility ceilings the same percentiles across 50 "
"random half-splits. Permutation tests in Figure 4e are exact, enumerating all C(15,6) = 5,005 "
"ways of splitting the features. Null distributions for the breadth analyses use 2,000 label "
"shuffles, which sets a resolution limit of p < 0.0005. Compound descriptors were computed with "
"RDKit [[rdkit]]; NPR1 and NPR2 are principal-moment ratios as defined by Sauer and Schwarz "
"[[sauer2003]]. Random seeds are fixed in the deposited scripts.",

"**Cell lines.** Figures 2 and 3 are the CCLE proteomics reanalysis across 375 cell lines. Figures "
"4 to 8 are the three published chemoproteomic screens, in four cell lines: CyMPSA in HeLa "
"(RRID:CVCL_0030) [[ryan2024]]; the meroterpenoid screen in MDA-MB-231 (RRID:CVCL_0062) and "
"MDA-MB-436 (RRID:CVCL_0623) [[bracken2024]]; and the NatChem screen in Ramos (RRID:CVCL_0597) and "
"MDA-MB-231 [[abbasov2021]]. Compounds were applied to intact cells in the CyMPSA and "
"meroterpenoid screens and to cell lysate in the NatChem screen; in all three, probe labelling "
"was performed on lysate after compound treatment, so the readout is identical across screens "
"and they differ only in where the compound encountered the proteome (note N4). No cells were "
"cultured for this work.",

"**Code and reproducibility.** No new data were generated. Every analysis and figure script is "
"deposited with its version and random seed, and every number in this paper is printed by the "
"script that produces the corresponding panel rather than transcribed from the artwork.",
]
