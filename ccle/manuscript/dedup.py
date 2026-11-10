# -*- coding: utf-8 -*-
"""Remove repeated text from the counted main text. Applied at build time, not in drafted.py,
so drafted.py stays the record of what was written and every removal is reviewable as a list.

Measured 19 October 2026 across the live tree: the 98.21% site-assignment sentence stood in four
figure legends, the "nothing here is a statement about a tumour" disclaimer in five places, the
margin >= 10 rule in ten, and Figure 2's section closer was its own legend footer verbatim.
Stating each once is not a cut to the argument. It removes a lap.

Two of these are free rather than cheap: Cell Press excludes STAR Methods from the word count,
so a method statement moved out of a legend and into Methods costs nothing at all.

Every bound is kept. The tumour disclaimer is an evidentiary bound, not boilerplate, so it stays
where a reader meets the acquired-lysine comparison (Figure 8) and in the Discussion's
limitations paragraph; only the three legends that repeat it lose it.
"""

Q  = "’"       # right single quote, as drafted.py uses it
GE = "≥"       # greater-than-or-equal
EN = "–"       # en dash

EDITS = [   # (block, old, new, label) -- located by content, never by index
 # ---- the 98.21% site-assignment statement, four legends to zero (it moves to Methods)
 ("CAP_2_4",
  "Site sets here and in Figures 5" + EN + "7 are chemoproteomic assignments concurring at 98.21% "
  "between two independent routes (44,736 of 45,552 usable sites); the figure compares site "
  "*sets*, and individual residue numbers within them are provisional.",
  "The figure compares site *sets*; individual residue numbers within them are provisional "
  "(Methods).", "98.21% out of the Figure 4 legend"),
 ("CAP_2_5",
  "concurring at 98.21% between two independent site-assignment routes. ", "",
  "98.21% out of the Figure 5 legend"),
 ("CAP_2_6",
  "Site sets are chemoproteomic assignments from three screens (the fragment screen in Ramos and "
  "MDA-MB-231, meroterpenoids in MDA-MB-231 and MDA-MB-436, CyMPSA in HeLa) concurring at 98.21% "
  "between two "
  "independent routes.",
  "Site sets are the three chemoproteomic screens of Figure 5 (Methods).",
  "98.21% out of the Figure 6 legend"),
 ("CAP_2_7",
  "Site sets are chemoproteomic assignments, concurring at 98.21% between two independent routes. "
  "Both probes were applied to HeLa cells; labelling",
  "Both probes were applied to HeLa cells; labelling",
  "98.21% out of the Figure 7 legend"),

 # ---- the tumour disclaimer, five places to two
 ("CAP_2_4",
  "*Effects are within-protein standardised differences. Measured in cell lines; nothing here is "
  "a statement about a tumour.*",
  "*Effects are within-protein standardised differences, measured in cell lines.*",
  "disclaimer out of the Figure 4 legend"),
 ("CAP_2_5",
  "*Engagement measured in cell lines. Nothing here is a statement about a tumour.*",
  "*Engagement measured in cell lines.*", "disclaimer out of the Figure 5 legend"),
 ("CAP_2_7",
  ", and nothing here is a statement about a tumour or a patient.", ".",
  "disclaimer out of the Figure 7 legend"),

 # ---- Figure 2's closer opened with the sentence its own legend footer carries verbatim.
 # The gate-1 against gate-2 bound is unique to it and stays.
 ("P_2_2",
  "This figure states what has been looked for and found, in these cell lines, with this "
  "workflow. It is not evidence that the peptides are absent from the sample, and the distinction "
  "is not rhetorical: gate 1 is a statement about tryptic geometry and gate 2 about search "
  "sensitivity, and neither licenses a claim about what the proteome contains.",
  "Gate 1 is a statement about tryptic geometry and gate 2 about search sensitivity, and neither "
  "licenses a claim about what the proteome contains.",
  "Figure 2 closer, first sentence duplicated its legend"),

 # ---- Figure 2's opener restated the unit definition its legend header gives in full
 ("P_2_2",
  "The population here is not the 295,411-site union of Figure 1 but what these particular lines "
  "carry: **8,470 site-in-cell-line observations** (8,055 distinct genomic sites across the 372 of "
  "those lines that carry one, a few sites assayed in more than one plex) together with 1,902 "
  "observations of variant sites that do not create a lysine, carried through the identical "
  "pipeline as a control (figure S2). Every rate below is per observation, because each line is "
  "an independent chance to detect.",
  "The population here is not the 295,411-site union of Figure 1 but what these particular lines "
  "carry, and the unit is the site-in-cell-line observation rather than the distinct genomic "
  "site.", "Figure 2 opener, units already in its legend"),

 # ---- the estimator stood in both the body and the legend of Figure 4; the legend keeps it,
 # because a legend has to stand alone
 ("P_2_4",
  "Both contrasts use the same estimator (one difference per protein, tested against zero), so "
  "the two sets of numbers are directly comparable, and both are run through the same ladder of "
  "structural-confidence gates. The gate is always applied to both groups at once; gating one "
  "side would manufacture a difference rather than reveal one.",
  "Both contrasts run through the same ladder of structural-confidence gates, always applied to "
  "both groups at once, since gating one side would manufacture a difference rather than reveal "
  "one.", "Figure 4, estimator stated in body and legend"),

 # ---- Figure 3's opener re-listed the spectrum counts Figure 2 had just given
 ("P_2_3",
  "Figure 2 recovered a small number of mutant peptides, and figure S3e bounds what they rest on: "
  "53 sites on at least one spectrum, 12 on two or more, four on two or more distinct peptides, "
  "and none on more than one experiment; at margin " + GE + " 10 those become 30, six, two and "
  "none. A result carried by single-spectrum identifications",
  "Figure 2 recovered a small number of mutant peptides, and figure S3e bounds what they rest on. "
  "A result carried by single-spectrum identifications",
  "Figure 3 opener, spectrum counts restated"),

 # ---- the margin rule was restated in full in Figure 3's closer, three paragraphs after it
 # was established
 ("P_2_3",
  "**The error rates quoted in this paper are measured at margin " + GE + " 10 rather than "
  "assumed, and Figure 2" + Q + "s census is reported at both thresholds so the two are never "
  "conflated.** Neither makes 0.58% a per-site confidence, and a single spectrum at margin "
  + GE + " 10 remains a single spectrum,",
  "Neither 0.83% nor 0.58% is a per-site confidence, and a single spectrum at margin " + GE +
  " 10 remains a single spectrum,", "Figure 3 closer, margin rule restated"),

 # ---- Figure 8's opener recapped Figures 4 to 7 before asking its own question
 ("P_2_8",
  "Figures 4 to 7 describe what governs engagement among lysines that compounds have actually "
  "been tested against: a property of the site rather than of the molecule, poorly predicted by "
  "the structural features available to us, and read out through a probe that is not itself "
  "neutral. The acquired lysines of Figure 1 have been tested against nothing. Before any of "
  "those rules can be expected to apply to them, one question has to be settled: do they sit "
  "inside the same population, or does mutation create lysines in environments unlike the ones "
  "the atlas has sampled?",
  "The acquired lysines of Figure 1 have been tested against nothing, so before the engagement "
  "rules measured above can be expected to reach them, one question has to be settled. Do they "
  "sit inside the same population, or does mutation create lysines in environments unlike the "
  "ones the atlas has sampled?", "Figure 8 opener, recap of Figures 4 to 7"),

 # ---- Figure 1's closer recapped the section before stating its bound. The bound stays.
 ("P_2_1",
  "Together these exclude the cheapest alternative to the premise of this paper, that acquired "
  "lysines are simply what mutation does everywhere. They are created at twice the germline rate, "
  "they are enriched on driver hotspots, and the processes that make them are identifiable. None "
  "of that says",
  "Together these exclude the cheapest alternative to the premise of this paper, that acquired "
  "lysines are simply what mutation does everywhere. None of that says",
  "Figure 1 closer, recap before its bound"),

 # ---- Figure 7's opener announced its question before asking it
 ("P_2_7",
  "Figure 6 shows that no tested compound descriptor predicts which lysines a compound prefers, "
  "and figure S7f shows that charge was never among the descriptors tested. That leaves the "
  "question open rather than settled. Here it is put directly, and to the reagent rather than to "
  "the compound, because in a competitive experiment it is the probe, not the compound, that "
  "decides which lysines are visible to be competed for at all.",
  "Charge was never among the descriptors tested (figure S7f), so the question is open rather "
  "than settled. It is put here to the reagent rather than to the compound, because in a "
  "competitive experiment it is the probe that decides which lysines are visible to be competed "
  "for at all.", "Figure 7 opener, announcement"),
 # ---- a third pass, 19 October: the remaining margin-rule appearances are mostly USES of the
 # rule rather than restatements of it, so the 180 words first estimated here were not there.
 # These three are genuine duplicates between a body paragraph and its own legend.
 ("P_2_2",
  " Of the 30 that clear margin ≥ 10, 24 rest on a single spectrum, six on two or more and "
  "two on two or more distinct peptides.",
  "", "Figure 2 body, spectrum depth already in its legend"),
 ("CAP_2_2",
  "The mutant peptide itself gives 53 sites at the reported 1% filter, 30 of them at margin "
  + GE + " 10. ",
  "", "Figure 2 legend, the 53-and-30 count given twice in one legend"),
 ("P_2_3",
  " Beyond that the estimate stops being informative rather than continuing to improve: only one "
  "entrapment peptide remains at a margin of 15, so the 0.18% there and the 0.31% at 20 are the "
  "behaviour of a single count, not evidence of a cleaner set.",
  "", "Figure 3 body, the margin-15 point is in its legend"),
 # ---- method text moved to STAR Methods, which Cell Press excludes from the word count. These
 # are free saves: nothing leaves the paper. Cell Press REQUIRES legends to define intervals,
 # error bars and replicate counts, so every interval definition stays in the legend; only the
 # construction of the reagent, the database and the model moves.
 ("CAP_2_3",
  "Error rates measured on the CCLE proteomes by burden-matched entrapment at r = 1 (FDRBench "
  "1.1.1). Entrapment sequences are each target's own residues permuted with the C-terminal "
  "lysine or arginine kept, so a pair shares an exact precursor mass and an identical length and "
  "only fragment ions separate them. Every entrapment identification is by construction false.",
  "Burden-matched entrapment at r = 1 (Methods). Every entrapment identification is by "
  "construction false.", "Figure 3 legend, entrapment construction to Methods"),
 ("CAP_2_5",
  "*The composite* is a logistic model on thirteen environment features plus AlphaFold confidence "
  "at the lysine, fitted on 30,805 lysines across 6,201 proteins and cross-validated by protein. "
  "Intervals",
  "*The composite* is a logistic model on fourteen site features, cross-validated by protein "
  "(Methods). Intervals", "Figure 5 legend, model fitting to Methods"),
 ("P_2_8",
  "We compared nine environment features at 156,781 union positions with an AlphaFold confidence "
  "of at least 70 [[jumper2021,varadi2022]], against the engaged and unengaged lysines measured "
  "in the three screens [[abbasov2021,bracken2024,ryan2024]]. The features describe a site's "
  "surroundings within 8 to 10 \u00c5: neighbour count, atom density, hydrophobic, aromatic and "
  "polar contacts, acidic and basic residues, net local charge, and other lysines.",
  "We compared nine features of a site's surroundings at 156,781 union positions with an "
  "AlphaFold confidence of at least 70 [[jumper2021,varadi2022]], against the engaged and "
  "unengaged lysines measured in the three screens [[abbasov2021,bracken2024,ryan2024]] "
  "(Methods).", "Figure 8 body, the feature list to Methods"),
]

# Stated once, in Methods, which Cell Press excludes from the word count.
METHODS_ADDITION = (
 "**Site assignment.** Site sets in Figures 4 to 7 are chemoproteomic assignments from three "
 "screens (the fragment screen in Ramos and MDA-MB-231 [[abbasov2021]], the meroterpenoid screen "
 "in MDA-MB-231 and MDA-MB-436 [[bracken2024]], the CyMPSA screen "
 "in HeLa), concurring at 98.21% between two independent assignment routes, 44,736 of 45,552 "
 "usable sites. Engagement is 60% occupancy or above throughout. Engagement was measured in cell "
 "lines, and no result in this paper is a statement about a tumour or a patient. "
 "**Entrapment.** Error rates were measured on the CCLE proteomes by burden-matched entrapment at "
 "r = 1 using FDRBench 1.1.1. Entrapment sequences are each target's own residues permuted with "
 "the C-terminal lysine or arginine kept, so a pair shares an exact precursor mass and an "
 "identical length and only fragment ions separate them. "
 "**Composite model.** The composite is a logistic model on thirteen environment features plus "
 "AlphaFold confidence at the lysine, fitted on 30,805 lysines across 6,201 proteins and "
 "cross-validated by protein. "
 "**Environment features.** The nine features compared at union positions describe a site's "
 "surroundings within 8 to 10 \u00c5: neighbour count, atom density, hydrophobic, aromatic and "
 "polar contacts, acidic and basic residues, net local charge, and other lysines.")

def apply(block_name, items):
    """Return items with every edit for this block applied. An edit is located by its own text,
    so it cannot be aimed at the wrong paragraph if the block gains or loses an entry. Raises if
    an edit matches no paragraph or more than one, so a drafted.py revision can never silently
    drop or double-apply one of these."""
    out = list(items)
    def txt(x): return x[0] if isinstance(x, tuple) else x
    for name, old, new, label in EDITS:
        if name != block_name:
            continue
        hits = [i for i, x in enumerate(out) if old in txt(x)]
        assert len(hits) == 1, "dedup edit matched %d paragraphs: %s" % (len(hits), label)
        i = hits[0]
        t2 = txt(out[i]).replace(old, new, 1)
        out[i] = (t2,) + out[i][1:] if isinstance(out[i], tuple) else t2
    return out
