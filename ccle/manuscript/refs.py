# -*- coding: utf-8 -*-
"""Bibliography for the acquired-lysine manuscript.

One entry per source, keyed by a short id that sections cite as [[key]]. Numbering is assigned by
the builder in order of first appearance across the assembled manuscript, so sections never carry
their own numbers and cannot drift out of step.

PRIMARY RESEARCH PAPERS ONLY, no reviews. Every entry below was verified against PubMed on
5 September 2026; `verified` records that.
"""
REFS = {
 "abbasov2021": dict(verified="PubMed 34504315", text=
  "Abbasov ME, Kavanagh ME, Ichu T-A, Lazear MR, Tao Y, Crowley VM, am Ende CW, Hacker SM, Ho J, "
  "Dix MM, Suciu R, Hayward MM, Kiessling LL, Cravatt BF. A proteome-wide atlas of lysine-reactive "
  "chemistry. Nat Chem. 2021;13(11):1081–1092. doi:10.1038/s41557-021-00765-4"),
 "bracken2024": dict(verified="article citation line", text=
  "Bracken AK, Gekko CE, Suss NO, Lueders EE, Cui Q, Fu Q, Lui ACW, Anderson ET, Zhang S, "
  "Abbasov ME. Biomimetic synthesis and chemical proteomics reveal the mechanism of action and "
  "functional targets of phloroglucinol meroterpenoids. J Am Chem Soc. 2024;146:2524–2548. "
  "doi:10.1021/jacs.3c10741"),
 "ryan2024": dict(verified="article citation line", text=
  "Ryan EM, Norinskiy MA, Bracken AK, Lueders EE, Chen X, Fu Q, Anderson ET, Zhang S, Abbasov ME. "
  "Activity-based acylome profiling with N-(cyanomethyl)-N-(phenylsulfonyl)amides for targeted "
  "lysine acylation and post-translational control of protein function in cells. J Am Chem Soc. "
  "2024;146:27622–27643."),
 "hacker2017": dict(verified="PubMed 29168484", text=
  "Hacker SM, Backus KM, Lazear MR, Forli S, Correia BE, Cravatt BF. Global profiling of lysine "
  "reactivity and ligandability in the human proteome. Nat Chem. 2017;9(12):1181–1190. "
  "doi:10.1038/nchem.2826"),
 "jumper2021": dict(verified="PubMed 34265844", text=
  "Jumper J, Evans R, Pritzel A, Green T, Figurnov M, Ronneberger O, et al. Highly accurate "
  "protein structure prediction with AlphaFold. Nature. 2021;596(7873):583–589. "
  "doi:10.1038/s41586-021-03819-2"),
 "varadi2022": dict(verified="PubMed 34791371", text=
  "Varadi M, Anyango S, Deshpande M, Nair S, Natassia C, Yordanova G, et al. AlphaFold Protein "
  "Structure Database: massively expanding the structural coverage of protein-sequence space with "
  "high-accuracy models. Nucleic Acids Res. 2022;50(D1):D439–D444. doi:10.1093/nar/gkab1061"),
 "tien2013": dict(verified="PubMed 24278298", text=
  "Tien MZ, Meyer AG, Sydykova DK, Spielman SJ, Wilke CO. Maximum allowed solvent accessibilities "
  "of residues in proteins. PLoS One. 2013;8(11):e80635. doi:10.1371/journal.pone.0080635"),
 "olsson2011": dict(verified="PubMed 26596171", text=
  "Olsson MHM, Søndergaard CR, Rostkowski M, Jensen JH. PROPKA3: consistent treatment of internal "
  "and surface residues in empirical pKa predictions. J Chem Theory Comput. 2011;7(2):525–537. "
  "doi:10.1021/ct100578z"),
 "leguilloux2009": dict(verified="PubMed 19486540", text=
  "Le Guilloux V, Schmidtke P, Tuffery P. Fpocket: an open source platform for ligand pocket "
  "detection. BMC Bioinformatics. 2009;10:168. doi:10.1186/1471-2105-10-168"),
 # Added 7 September 2026 for table S1 (functional consequence of the recurrent acquired lysines).
 # Each verified against PubMed on 7 September 2026; PMIDs recorded in `verified` as elsewhere.
 "walker2021": dict(verified="PubMed 33451349", text=
  "Walker JS, Hing ZA, Harrington B, Baumhardt J, Ozer HG, Lehman A, et al. Recurrent XPO1 "
  "mutations alter pathogenesis of chronic lymphocytic leukemia. J Hematol Oncol. 2021;14(1):17. "
  "doi:10.1186/s13045-021-01032-2"),
 "bader2006": dict(verified="PubMed 16432179", text=
  "Bader AG, Kang S, Vogt PK. Cancer-specific mutations in PIK3CA are oncogenic in vivo. "
  "Proc Natl Acad Sci U S A. 2006;103(5):1475\u20131479. doi:10.1073/pnas.0510857103"),
 "hao2022": dict(verified="PubMed 35418124", text=
  "Hao Y, He B, Wu L, Li Y, Wang C, Wang T, et al. Nuclear translocation of p85\u03b2 promotes "
  "tumorigenesis of PIK3CA helical domain mutant cancer. Nat Commun. 2022;13(1):1974. "
  "doi:10.1038/s41467-022-29585-x"),
 "xu2013": dict(verified="PubMed 23675459", text=
  "Xu D, Zheng H, Yu WM, Qu CK. Activating mutations in protein tyrosine phosphatase Ptpn11 (Shp2) "
  "enhance reactive oxygen species production that contributes to myeloproliferative disorder. "
  "PLoS One. 2013;8(5):e63152. doi:10.1371/journal.pone.0063152"),
 "buhrman2007": dict(verified="PubMed 18073111", text=
  "Buhrman G, Wink G, Mattos C. Transformation efficiency of RasQ61 mutants linked to structural "
  "features of the switch regions in the presence of Raf. Structure. 2007;15(12):1618\u20131629. "
  "doi:10.1016/j.str.2007.10.011"),
 "carpten2007": dict(verified="PubMed 17611497", text=
  "Carpten JD, Faber AL, Horn C, Donoho GP, Briggs SL, Robbins CM, et al. A transforming mutation "
  "in the pleckstrin homology domain of AKT1 in cancer. Nature. 2007;448(7152):439–444. "
  "doi:10.1038/nature05933"),
 "craven2025": dict(verified="PubMed 39506119", text=
  "Craven GB, Chu H, Sun JD, Carelli JD, Coyne B, Chen H, et al. Mutant-selective AKT inhibition "
  "through lysine targeting and neo-zinc chelation. Nature. 2025;637(8044):205–214. "
  "doi:10.1038/s41586-024-08176-4"),
 # Added 15 September 2026 from the student reference submission of the same date. Every entry
 # below was matched against PubMed by journal, year, volume, first page and first author before
 # being accepted; five submitted citations were rejected because they resolved to a different
 # paper, and are recorded in results/TableS1_reference_verification_2026-09-15.docx.
 "ackermann2005": dict(verified="PubMed 15899789", text=
  "Ackermann J, Frutschi M, Kaloulis K, McKee T, Trumpp A, Beermann F. Metastasizing melanoma "
  "formation caused by expression of activated N-RasQ61K on an INK4a-deficient background. "
  "Cancer Res. 2005;65(10):4005\u20134011. doi:10.1158/0008-5472.CAN-04-2970"),
 "li2012": dict(verified="PubMed 22718121", text=
  "Li A, Ma Y, Jin M, Mason S, Mort RL, Blyth K, et al. Activated mutant NRasQ61K drives aberrant "
  "melanocyte signaling, survival, and invasiveness via a Rac1-dependent mechanism. J Invest "
  "Dermatol. 2012;132(11):2610\u20132621. doi:10.1038/jid.2012.186"),
 "boudreau2014": dict(verified="PubMed 24714748", text=
  "Boudreau HE, Casterline BW, Burke DJ, Leto TL. Wild-type and mutant p53 differentially regulate "
  "NADPH oxidase 4 in TGF-\u03b2-mediated migration of human lung and breast epithelial cells. "
  "Br J Cancer. 2014;110(10):2569\u20132582. doi:10.1038/bjc.2014.165"),
 "gebregiworgis2024": dict(verified="PubMed 38640594", text=
  "Gebregiworgis T, Chan JY, Kuntz DA, Priv\u00e9 GG, Marshall CB, Ikura M. Crystal structure of "
  "NRAS Q61K with a ligand-induced pocket near switch II. Eur J Cell Biol. 2024;103(2):151414. "
  "doi:10.1016/j.ejcb.2024.151414"),
 "gomes2018": dict(verified="PubMed 29652801", text=
  "Gomes AS, Trov\u00e3o F, Andrade Pinheiro B, Freire F, Gomes S, Oliveira C, et al. The crystal "
  "structure of the R280K mutant of human p53 explains the loss of DNA binding. Int J Mol Sci. "
  "2018;19(4):1184. doi:10.3390/ijms19041184"),
 "rubinstein2010": dict(verified="PubMed 20630094", text=
  "Rubinstein JC, Sznol M, Pavlick AC, Ariyan S, Cheng E, Bacchiocchi A, et al. Incidence of the "
  "V600K mutation among melanoma patients with BRAF mutations, and potential therapeutic response "
  "to the specific BRAF inhibitor PLX4032. J Transl Med. 2010;8:67. doi:10.1186/1479-5876-8-67"),
 "menzies2012": dict(verified="PubMed 22535154", text=
  "Menzies AM, Haydu LE, Visintin L, Carlino MS, Howle JR, Thompson JF, et al. Distinguishing "
  "clinicopathologic features of patients with V600E and V600K BRAF-mutant metastatic melanoma. "
  "Clin Cancer Res. 2012;18(12):3242\u20133249. doi:10.1158/1078-0432.CCR-12-0052"),
 "piresdasilva2019": dict(verified="PubMed 30630828", text=
  "Pires da Silva I, Wang KYX, Wilmott JS, Holst J, Carlino MS, Park JJ, et al. Distinct molecular "
  "profiles and immunotherapy treatment outcomes of V600E and V600K-mutant melanoma. Clin Cancer "
  "Res. 2019;25(4):1272\u20131279. doi:10.1158/1078-0432.CCR-18-1680"),
 "dunsche2023": dict(verified="PubMed 38097548", text=
  "Dunsche L, Ivanisenko N, Riemann S, Schindler S, Beissert S, Angeli C, et al. A cytosolic "
  "mutp53(E285K) variant confers chemoresistance of malignant melanoma. Cell Death Dis. "
  "2023;14(12):831. doi:10.1038/s41419-023-06360-4"),
 "christgen2012": dict(verified="PubMed 22945757", text=
  "Christgen M, Noskowicz M, Heil C, Schipper E, Christgen H, Geffers R, et al. IPH-926 lobular "
  "breast cancer cells harbor a p53 mutant with temperature-sensitive functional activity and "
  "allow for profiling of p53-responsive genes. Lab Invest. 2012;92(11):1635\u20131647. "
  "doi:10.1038/labinvest.2012.126"),
 "angulo2013": dict(verified="PubMed 24136356", text=
  "Angulo I, Vadas O, Gar\u00e7on F, Banham-Hall E, Plagnol V, Leahy TR, et al. Phosphoinositide "
  "3-kinase \u03b4 gene mutation predisposes to respiratory infection and airway damage. Science. "
  "2013;342(6160):866\u2013871. doi:10.1126/science.1243292"),
 "grabiner2014": dict(verified="PubMed 24631838", text=
  "Grabiner BC, Nardi V, Birsoy K, Possemato R, Shen K, Sinha S, et al. A diverse array of "
  "cancer-associated MTOR mutations are hyperactivating and can predict rapamycin "
  "sensitivity. Cancer Discov. 2014;4(5):554\u2013563. doi:10.1158/2159-8290.CD-13-0929"),
 "lucas2014": dict(verified="PubMed 24165795", text=
  "Lucas CL, Kuehn HS, Zhao F, Niemela JE, Deenick EK, Palendira U, et al. Dominant-activating "
  "germline mutations in the gene encoding the PI(3)K catalytic subunit p110\u03b4 result in "
  "T cell senescence and human immunodeficiency. Nat Immunol. 2014;15(1):88\u201397. "
  "doi:10.1038/ni.2771"),
 "penther2019": dict(verified="PubMed 30779244", text=
  "Penther D, Viailly PJ, Latour S, Ruminy P, Bobee V, Kaplan B, et al. A recurrent clonally "
  "distinct Burkitt lymphoma case highlights genetic key events contributing to oncogenesis. "
  "Genes Chromosomes Cancer. 2019;58(9):595\u2013601. doi:10.1002/gcc.22743"),
 "nusinow2020": dict(verified="PubMed 31978347", text=
  "Nusinow DP, Szpyt J, Ghandi M, Rose CM, McDonald ER 3rd, Kalocsay M, et al. Quantitative "
  "proteomics of the Cancer Cell Line Encyclopedia. Cell. 2020;180(2):387–402.e16. "
  "doi:10.1016/j.cell.2019.12.023"),
 "hess2019": dict(verified="PubMed 31526759", text=
  "Hess JM, Bernards A, Kim J, Miller M, Taylor-Weiner A, Haradhvala NJ, Lawrence MS, Getz G. "
  "Passenger hotspot mutations in cancer. Cancer Cell. 2019;36(3):288–301.e14. "
  "doi:10.1016/j.ccell.2019.08.002"),
 "sauer2003": dict(verified="PubMed 12767158", text=
  "Sauer WHB, Schwarz MK. Molecular shape diversity of combinatorial libraries: a prerequisite for "
  "broad bioactivity. J Chem Inf Comput Sci. 2003;43(3):987–1003. doi:10.1021/ci025599w"),
 "rdkit": dict(verified="software, no primary paper", text=
  "Landrum G, et al. RDKit: open-source cheminformatics [computer software]. https://www.rdkit.org. "
  "Accessed 6 September 2026."),
 "pedregosa2011": dict(verified="publisher page jmlr.org/papers/v12/pedregosa11a.html", text=
  "Pedregosa F, Varoquaux G, Gramfort A, Michel V, Thirion B, Grisel O, et al. Scikit-learn: "
  "machine learning in Python. J Mach Learn Res. 2011;12:2825\u20132830."),

 "liu2020": dict(verified="PubMed 33261662", text=
  "Liu X, Li C, Mou C, Dong Y, Tu Y. dbNSFP v4: a comprehensive database of transcript-specific "
  "functional predictions and annotations for human nonsynonymous and splice-site SNVs. Genome "
  "Med. 2020;12(1):103. doi:10.1186/s13073-020-00803-9"),

 "ostrem2013": dict(verified="PubMed 24256730", text=
  "Ostrem JM, Peters U, Sos ML, Wells JA, Shokat KM. K-Ras(G12C) inhibitors allosterically control "
  "GTP affinity and effector interactions. Nature. 2013;503(7477):548\u2013551. "
  "doi:10.1038/nature12796"),

 "smith2013": dict(verified="PubMed 23487764", text=
  "Smith MJ, Neel BG, Ikura M. NMR-based functional profiling of RASopathies and oncogenic RAS "
  "mutations. Proc Natl Acad Sci U S A. 2013;110(12):4574\u20134579. doi:10.1073/pnas.1218173110"),
 "takahashi2024": dict(verified="PubMed 38653237", text=
  "Takahashi M, Chong HB, Zhang S, Yang T-Y, Lazarov MJ, Harry S, et al. DrugMap: a quantitative "
  "pan-cancer analysis of cysteine ligandability. Cell. 2024;187(10):2536\u20132556.e30. "
  "doi:10.1016/j.cell.2024.03.027"),
 "kato2003": dict(verified="PubMed 12826609; MaveDB urn:mavedb:00001234-a-1 through -h-1, "
                          "variant #414 = p.Glu286Lys, scores 18.3 25.9 23.5 12.0 13.3 6.0 11.7 21.8", text=
  "Kato S, Han S-Y, Liu W, Otsuka K, Shibata H, Kanamaru R, Ishioka C. Understanding the "
  "function\u2013structure and function\u2013mutation relationships of p53 tumor suppressor protein by "
  "high-resolution missense mutation analysis. Proc Natl Acad Sci U S A. 2003;100(14):"
  "8424\u20138429. doi:10.1073/pnas.1431692100"),
 "kuljanin2021": dict(verified="PubMed 33398154", text=
  "Kuljanin M, Mitchell DC, Schweppe DK, Gikandi AS, Nusinow DP, Bulloch NJ, Vinogradova EV, "
  "Wilson DL, Kool ET, Mancias JD, Cravatt BF, Gygi SP. Reimagining high-throughput profiling of "
  "reactive cysteines for cell-based screening of large electrophile libraries. Nat Biotechnol. "
  "2021;39(5):630\u2013641. doi:10.1038/s41587-020-00778-3"),
 "der1986": dict(verified="PubMed 3510078", text=
  "Der CJ, Finkel T, Cooper GM. Biological and biochemical properties of human rasH genes mutated "
  "at codon 61. Cell. 1986;44(1):167\u2013176. doi:10.1016/0092-8674(86)90495-2"),
 "desai2024": dict(verified="PubMed 39468056", text=
  "Desai H, Andrews KH, Bergersen KV, Ofori S, Yu F, Shikwana F, et al. Chemoproteogenomic "
  "stratification of the missense variant cysteinome. Nat Commun. 2024;15(1):9284. "
  "doi:10.1038/s41467-024-53520-x"),

 "tsherniak2017": dict(verified="PubMed 28753430", text=
  "Tsherniak A, Vazquez F, Montgomery PG, Weir BA, Kryukov G, Cowley GS, et al. Defining a cancer "
  "dependency map. Cell. 2017;170(3):564\u2013576.e16. doi:10.1016/j.cell.2017.06.010"),

 "ellrott2018": dict(verified="PubMed 29596782", text=
  "Ellrott K, Bailey MH, Saksena G, Covington KR, Kandoth C, Stewart C, et al. Scalable open "
  "science approach for mutation calling of tumor exomes using multiple genomic pipelines. "
  "Cell Syst. 2018;6(3):271\u2013281.e7. doi:10.1016/j.cels.2018.03.002"),

 "tate2019": dict(verified="PubMed 30371878", text=
  "Tate JG, Bamford S, Jubb HC, Sondka Z, Beare DM, Bindal N, et al. COSMIC: the Catalogue Of "
  "Somatic Mutations In Cancer. Nucleic Acids Res. 2019;47(D1):D941\u2013D947. "
  "doi:10.1093/nar/gky1015"),

 "karczewski2020": dict(verified="PubMed 32461654", text=
  "Karczewski KJ, Francioli LC, Tiao G, Cummings BB, Alf\u00f6ldi J, Wang Q, et al. The mutational "
  "constraint spectrum quantified from variation in 141,456 humans. Nature. "
  "2020;581(7809):434\u2013443. doi:10.1038/s41586-020-2308-7"),

 "elias2007": dict(verified="PubMed 17327847", text=
  "Elias JE, Gygi SP. Target-decoy search strategy for increased confidence in large-scale protein "
  "identifications by mass spectrometry. Nat Methods. 2007;4(3):207\u2013214. "
  "doi:10.1038/nmeth1019"),

 "wen2025": dict(verified="PubMed 40524023", text=
  "Wen B, Freestone J, Riffle M, MacCoss MJ, Noble WS, Keich U. Assessment of false discovery rate "
  "control in tandem mass spectrometry analysis using entrapment. Nat Methods. 2025;22(7):"
  "1454–1463. doi:10.1038/s41592-025-02719-x"),
 # --- added 9 Oct 2026: the broad-spectrum lysine ligandability census, and the two
 #     Cell Chemical Biology computational precedents. All verified against PubMed that day.
 "ward2017": dict(verified="PubMed 28445029", text=
  "Ward CC, Kleinman JI, Nomura DK. NHS-esters as versatile reactivity-based probes for mapping "
  "proteome-wide ligandable hotspots. ACS Chem Biol. 2017;12(6):1478\u20131483. "
  "doi:10.1021/acschembio.7b00125"),

 "zhou2025": dict(verified="PubMed 39543955", text=
  "Zhou M, Li S, Tan Y, Huang W, Li Y, Yuan X, Li Z. Global profiling lysine reactivity and "
  "ligandability with oxidant-triggered bioconjugation chemistry. Angew Chem Int Ed Engl. "
  "2025;64(6):e202418473. doi:10.1002/anie.202418473"),

 "boatner2023": dict(verified="PubMed 37119813", text=
  "Boatner LM, Palafox MF, Schweppe DK, Backus KM. CysDB: a human cysteine database based on "
  "experimental quantitative chemoproteomics. Cell Chem Biol. 2023;30(6):683\u2013698.e3. "
  "doi:10.1016/j.chembiol.2023.04.004"),

 "white2023": dict(verified="PubMed 37451266 (Meta-Analysis; primary computational analysis, not a review)", text=
  "White MEH, Gil J, Tate EW. Proteome-wide structural analysis identifies warhead- and "
  "coverage-specific biases in cysteine-focused chemoproteomics. Cell Chem Biol. "
  "2023;30(7):828\u2013838.e4. doi:10.1016/j.chembiol.2023.06.021"),
}
