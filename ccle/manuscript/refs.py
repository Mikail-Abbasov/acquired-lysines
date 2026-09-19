# -*- coding: utf-8 -*-
"""Bibliography for the acquired-lysine manuscript.

One entry per source, keyed by a short id that sections cite as [[key]]. Numbering is assigned by
the builder in order of first appearance across the assembled manuscript, so sections never carry
their own numbers and cannot drift out of step.

PRIMARY RESEARCH PAPERS ONLY — no reviews. Every entry below was verified against PubMed on
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
 "carpten2007": dict(verified="PubMed 17611497", text=
  "Carpten JD, Faber AL, Horn C, Donoho GP, Briggs SL, Robbins CM, et al. A transforming mutation "
  "in the pleckstrin homology domain of AKT1 in cancer. Nature. 2007;448(7152):439–444. "
  "doi:10.1038/nature05933"),
 "craven2025": dict(verified="PubMed 39506119", text=
  "Craven GB, Chu H, Sun JD, Carelli JD, Coyne B, Chen H, et al. Mutant-selective AKT inhibition "
  "through lysine targeting and neo-zinc chelation. Nature. 2025;637(8044):205–214. "
  "doi:10.1038/s41586-024-08176-4"),
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

 "takahashi2024": dict(verified="PubMed 38653237", text=
  "Takahashi M, Chong HB, Zhang S, Yang T-Y, Lazarov MJ, Harry S, et al. DrugMap: a quantitative "
  "pan-cancer analysis of cysteine ligandability. Cell. 2024;187(10):2536\u20132556.e30. "
  "doi:10.1016/j.cell.2024.03.027"),
 "kuljanin2021": dict(verified="PubMed 33398154", text=
  "Kuljanin M, Mitchell DC, Schweppe DK, Gikandi AS, Nusinow DP, Bulloch NJ, Vinogradova EV, "
  "Wilson DL, Kool ET, Mancias JD, Cravatt BF, Gygi SP. Reimagining high-throughput profiling of "
  "reactive cysteines for cell-based screening of large electrophile libraries. Nat Biotechnol. "
  "2021;39(5):630\u2013641. doi:10.1038/s41587-020-00778-3"),
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
}
