#!/usr/bin/env python3
"""
build_release.py -- assemble the shipping tree for the code+data deposit.

WHY A STAGING COPY AND NOT A CLEANUP
The working tree holds ~270 scripts: eighteen build the published figures, the rest are earlier
versions of those eighteen, R61 grant figures, and one-off exploration. Deleting or moving them
would rewrite a working analysis directory to serve packaging, which is the wrong trade -- and
the earlier versions are the record of how the figures got to where they are. Instead this
computes what the eighteen actually need and copies exactly that into release/. Nothing in the
working tree moves.

HOW THE SHIP LIST IS DERIVED
Transitive closure, not a hand list:

    the 18 published figure scripts
      + every script that writes a file they read
      + every script that writes a file THOSE read, and so on
      + the shared modules and the manuscript build

Any file at the edge of that closure with no producer in the tree is an EXTERNAL input. Those
are not allowed to be anonymous: each must appear in EXTERNAL below with its source and version,
and the build FAILS if one does not. That is the whole point -- an undeclared input is how a
"reproducible" deposit turns out not to be.

RUN IT FROM repo/. Importing build_manifest regenerates the manifest as a side effect, so the
release can never be assembled against a stale one.
"""
import os, re, sys, json, shutil, hashlib
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_manifest as bm          # regenerates manifest.csv / manifest.md on import

TOP = bm.TOP
OUT = os.path.join(TOP, "release")

# ---------------------------------------------------------------------------------------------
# Every input at the edge of the closure, with where it comes from. Version strings are the ones
# the manuscript's Methods state; if a version changes there it must change here.
EXTERNAL = {
    # --- DepMap 26Q1
    "Model.csv":                        ("DepMap Public 26Q1", "depmap.org/portal/download"),
    "OmicsSomaticMutations.csv":        ("DepMap Public 26Q1", "depmap.org/portal/download"),
    "CRISPRGeneEffect.csv":             ("DepMap Public 26Q1", "depmap.org/portal/download"),
    "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv":
                                        ("DepMap Public 26Q1", "depmap.org/portal/download"),
    # --- somatic variant catalogues
    "mc3.v0.2.8.PUBLIC.maf":            ("TCGA MC3 v0.2.8-PUBLIC", "gdc.cancer.gov/about-data/publications/mc3-2017"),
    "cmc_export.tsv":                   ("COSMIC Cancer Mutation Census v103", "cancer.sanger.ac.uk/cmc"),
    # --- population frequency
    "dbNSFP5.4a":                       ("dbNSFP 5.4a", "database.liulab.science/dbNSFP"),
    "gnomad.joint.v4.1":                ("gnomAD v4.1 joint", "gnomad.broadinstitute.org"),
    "gnomad.exomes.r2.1.1.non_cancer":  ("gnomAD v2.1.1 non-cancer exomes", "gnomad.broadinstitute.org"),
    # --- reference proteome and annotation
    "MANE.GRCh38.v1.2":                 ("MANE Select v1.2", "ncbi.nlm.nih.gov/refseq/MANE"),
    "uniprot_gene_map.tsv":             ("UniProt release used at analysis time", "uniprot.org"),
    "uniprot_annotation.tsv":           ("UniProt release used at analysis time", "uniprot.org"),
    "uniprot_keyword_categories.tsv":   ("UniProt release used at analysis time", "uniprot.org"),
    "uniprot_hugo_entrez_id_mapping.csv": ("UniProt / HGNC / Entrez cross-reference", "uniprot.org"),
    # --- CCLE proteomes (Nusinow et al. 2020, Cell 180:387-402)
    "harmonized_MS_CCLE_Gygi.csv":      ("Nusinow et al. 2020 Table S2", "doi.org/10.1016/j.cell.2019.12.023"),
    "Table_S1_Sample_Information.xlsx": ("Nusinow et al. 2020 Table S1", "doi.org/10.1016/j.cell.2019.12.023"),
    # --- the three chemoproteomic screens, compound structures from their SI
    # Redistributed here (20 KB) so Figure 6 runs from the deposit. Compound structures are
    # factual data taken from the three papers' supporting information and are cited as such.
    "SMILES.xlsx":                      ("Abbasov 2021 Nat Chem; Bracken 2024 JACS; Ryan 2024 JACS -- SI compound tables, redistributed here with citation",
                                         "see manuscript references"),
    # --- entrapment
    "uniprot_entrap_decoy.txt":         ("FDRBench v1.1.1", "github.com/Noble-Lab/FDRBench"),
    "uniprot_entrapment_pep.fasta":     ("FDRBench v1.1.1", "github.com/Noble-Lab/FDRBench"),
}

EXTERNAL.update({
    "CRISPRInferredCommonEssentials.csv": ("DepMap Public 26Q1", "depmap.org/portal/download"),
    "gdc_barcode_project.csv":            ("GDC sample-to-project barcode map", "portal.gdc.cancer.gov"),
    "CyMPSA":         ("Ryan et al. 2024 JACS, SI site tables", "see manuscript references"),
    "Meroterpenoid":  ("Bracken et al. 2024 JACS, SI site tables", "see manuscript references"),
    "NatChem":        ("Abbasov et al. 2021 Nat Chem, SI site tables", "see manuscript references"),
})

# ---------------------------------------------------------------------------------------------
# ANALYSIS INTERMEDIATES THAT SHIP AS DATA.
# Each of these is produced inside this project, but its producing script either assembles the
# filename at runtime or writes it through a config constant, so the manifest's pattern scan
# cannot tie it to a writer. Rather than let them read as "external" -- which would be false --
# they are deposited alongside the code with the note below, and the figures therefore run from
# the deposit even where the producing step cannot be re-executed. This is the whole reason the
# deposit has two tiers.
DEPOSITED = {
    "hits62.tsv":                     "per-PSM detail for the acquired-lysine detections",
    "fig3_margins.csv":               "per-peptide score margins, cached from the PSM tables",
    "fig3_deltamass.csv":             "precursor mass errors, cached from the PSM tables (97 MB)",
    "region_compare_pep.tsv":         "per-site search outcome, mutation-aware peptide database",
    "region_compare_conv.tsv":        "per-site search outcome, conventional MANE database",
    "fdp_p2.csv":                     "plex-2 q-value vs entrapment curve (9.5 MB)",
    "step10_roc_pooled.csv":          "pooled ROC of the ligandability composite",
    "step10_roc_splits.csv":          "per-split ROC band of the ligandability composite",
    "step37_within_NatChem.csv":      "within-NatChem descriptor x site-preference correlations",
    "step43_enclosure.csv":           "per-site enclosure, from the AlphaFold models (4.3 MB)",
    "step27_modnz_measured.csv":      "measured lysine NZ positions (2.3 MB)",
    "step27_modnz_union.csv":         "modelled NZ positions for union sites (9.0 MB)",
    "proteomics_subset.csv":          "the 375-line CCLE proteomic subset (61 MB)",
    "acquired_K_searched_vs_detected.csv": "per-site searched/detected table for acquired lysines",
    "depmap_variants_annotated.parquet":   "annotated DepMap variant table (46 MB)",
    "sites_covalent_residues.csv":    "covalent-residue site table across DepMap (28 MB)",
    "residue_flux_by_line.csv":       "per-line residue gain/loss counts",
    "features_final_full.csv":        "per-site structural features, ungated (11 MB)",
    "mc3_nt.tsv":                     "MC3 nucleotide-level table (31 MB)",
    "mc3_bysample.tsv":               "MC3 per-sample table (34 MB)",
}

# Files that are outputs of the search pipeline rather than of a Python script in this tree:
# MSFragger/FragPipe/Philosopher write them, and the pipeline configuration ships instead.
PIPELINE = {
    "psm.tsv", "peptide.tsv", "protein.tsv", "ion.tsv",
    "hit_spectra_mgf_features.tsv",       # SpecQuality, run over the mzML exports
}

# Shared modules and build machinery that no figure "reads" as data but every figure needs.
ALWAYS = [
    "ccle/abbasov_style.py",
    "ccle/figure_data.py",
    "ccle/manuscript/build.py",
    "ccle/manuscript/sections.py",
    "ccle/manuscript/drafted.py",
    "ccle/manuscript/si.py",
    "ccle/manuscript/refs.py",
    "ccle/manuscript/seams_check.py",
    "repo/build_manifest.py",
    "repo/build_release.py",
    "repo/verify.py",
    "repo/reverse_check.py",
]


# Files that belong at the root of the deposit rather than inside repo/.
# Files reached through a helper rather than a literal open()/read_csv(), so no static scan can
# see them. Figure 1g calls _read_germline("germline_lysine.out"); the acceptance test caught
# their absence as five NaNs in panel g. Add to this list when a figure reads through indirection.
EXTRA_DATA = [
    "ccle/germline_lysine.out",
    "ccle/germline_nc.out",
    "SpecQuality/hit_spectra_mgf_features.tsv",
    "ccle/SMILES.xlsx",
]

# Globs, for inputs a script discovers rather than names. Figure 8e counts the AKT1 E17K peptides
# by walking the search databases; without them it reports zero, which reads as "never searched"
# and is the opposite of the claim. 120 MB to make the paper's closing move checkable.
EXTRA_GLOBS = [
    "ccle/plex*_varpep.fasta",
    "ccle/plex*_downpep.fasta",
    # 319 MB of entrapment pairings, so figure S4 -- the entrapment figure -- re-runs from the
    # deposit. Entrapment is what the paper's false-discovery claim rests on, so it is the last
    # thing that should be un-checkable.
    "ccle/plex*_varpair.txt",
]

ROOT_FILES = ["README.md", "LICENSE", "LICENSE-DATA", "CITATION.cff", ".zenodo.json",
              ".gitignore"]


def _copy(src, dst):
    """Copy unless an identical file is already there. Identical means same size AND same mtime:
    size alone would keep a truncated file that happens to match, and always-copy would re-write
    ~800 MB over the mount on every build. A 1-byte fig3_margins.csv survived an earlier
    skip-if-exists and killed every run of that figure until it was noticed."""
    if os.path.exists(dst):
        a, b = os.stat(src), os.stat(dst)
        if a.st_size == b.st_size and int(a.st_mtime) == int(b.st_mtime):
            return False
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copy2(src, dst)
    return True


def _ver(path):
    """Sort key that ranks foo_v3.py > foo_v2.py > foo.py, then prefers the shorter name."""
    b = os.path.basename(path)
    m = re.search(r"_v(\d+)\.py$", b)
    return (int(m.group(1)) if m else 0, -len(b))


def closure():
    """Walk from the 18 figure scripts back through producers. Returns (scripts, external)."""
    want = {s for _, s in bm.FIGURES}
    seen, ext = set(), defaultdict(set)
    while want:
        name = want.pop()
        if name in seen:
            continue
        seen.add(name)
        p = bm.find_script(name)
        if p is None:
            print("  !! figure script not found: %s" % name)
            continue
        reads, _, _, _ = bm.scan(p)
        for r in reads:
            b = os.path.basename(r)
            src = bm.producers.get(b, [])
            # A script that writes its own report and then re-reads it is NOT an external input.
            # Excluding self here made 30 in-tree files look external on the first run.
            if src:
                # Several files are written by both a script and its own later version
                # (step6_features.py and step6_features_v2.py both write features.csv). Pulling
                # in both would ship a superseded script AND report a false collision, so take
                # the highest version -- the convention this tree uses everywhere.
                best = max((x for x in src if os.path.basename(x) != name),
                           key=_ver, default=None)
                if best is not None:
                    want.add(os.path.basename(best))
            else:
                ext[b].add(name)
    return seen, ext


def declared(b):
    """Is this edge input accounted for -- external, pipeline output, or deposited intermediate?"""
    if b in EXTERNAL or b in PIPELINE or b in DEPOSITED:
        return True
    return any(b.startswith(k) for k in EXTERNAL)


def source_of(b):
    """Where an edge input comes from, for the report and for external_inputs.json."""
    if b in EXTERNAL:
        return EXTERNAL[b][0]
    if b in DEPOSITED:
        return "deposited intermediate -- " + DEPOSITED[b]
    if b in PIPELINE:
        return "search-pipeline output (MSFragger / FragPipe / Philosopher)"
    for k in EXTERNAL:
        if b.startswith(k):
            return EXTERNAL[k][0]
    return "UNDECLARED"


def main():
    print("\n" + "=" * 78)
    print("BUILD RELEASE")
    print("=" * 78)
    scripts, ext = closure()

    # ---- pre-flight: two scripts writing one filename
    # step31_crystal_validation.py silently overwrote step31_compound_site_matching.py's report,
    # and the within-dataset arm of that analysis was lost until it was re-run on 2026-09-07.
    # Nothing warned. This does.
    clash = {out: sorted(set(w)) for out, w in bm.producers.items()
             if len({os.path.basename(x) for x in w}) > 1
             and not out.endswith((".png", ".pdf"))}
    shipping = {os.path.basename(x) for x in scripts}
    hot = {out: w for out, w in clash.items()
           if sum(1 for x in w if os.path.basename(x) in shipping) > 1}
    if hot:
        print("\nCOLLISIONS AMONG SHIPPING SCRIPTS -- resolve before depositing:")
        for out, w in sorted(hot.items()):
            print("   %-30s %s" % (out, ", ".join(w)))
    else:
        print("output collisions among shipping scripts: none")
    cold = {k: v for k, v in clash.items() if k not in hot}
    if cold:
        print("  (%d further collisions involve superseded scripts that do not ship: %s)"
              % (len(cold), ", ".join(sorted(cold)[:4]) + (" ..." if len(cold) > 4 else "")))
    print("scripts in the closure of the 18 published figures: %d" % len(scripts))

    undeclared = sorted(b for b in ext if not declared(b))
    if undeclared:
        print("\nUNDECLARED EXTERNAL INPUTS -- add each to EXTERNAL with its source and version:")
        for b in undeclared:
            print("   %-46s  read by %s" % (b, ", ".join(sorted(ext[b]))))
        print("\n%d undeclared. Release NOT built." % len(undeclared))
        return 1

    print("edge inputs, all accounted for: %d" % len(ext))
    for b in sorted(ext):
        print("   %-46s  %s" % (b, source_of(b)))

    # ---- assemble
    if os.path.isdir(OUT):
        shutil.rmtree(OUT, ignore_errors=True)
        if os.path.isdir(OUT):
            # Some environments (this one included) refuse deletes. The build then OVERLAYS the
            # previous release, so a file that should have disappeared silently persists and the
            # size report undercounts. Say so rather than pretending the tree is clean.
            print("  NOTE: release/ could not be cleared; this build overlays the previous one.")
            print("        File counts below count only what THIS run copied.")
    n_scripts = n_root = n_panel = 0
    inter_n, inter_bytes = 0, 0
    b_scripts = b_root = b_panel = 0
    ship = sorted({os.path.join(d, s) for s in scripts
                   for d in bm.SEARCH_DIRS
                   if os.path.exists(os.path.join(TOP, d, s))} | set(ALWAYS))
    for rel in ship:
        src = os.path.join(TOP, rel)
        if not os.path.exists(src):
            print("   missing, skipped: %s" % rel); continue
        dst = os.path.join(OUT, rel)
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copy2(src, dst); n_scripts += 1; b_scripts += os.path.getsize(src)

    for rel in EXTRA_DATA:
        src = os.path.join(TOP, rel)
        dst = os.path.join(OUT, rel)
        if os.path.exists(src) and not os.path.exists(dst):
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            shutil.copy2(src, dst); n_scripts += 1; b_scripts += os.path.getsize(src)
        elif not os.path.exists(src):
            print("   EXTRA_DATA missing: %s" % rel)

    import glob as _glob
    for pat in EXTRA_GLOBS:
        hits = sorted(_glob.glob(os.path.join(TOP, pat)))
        if not hits:
            print("   EXTRA_GLOBS matched nothing: %s" % pat)
        for src in hits:
            rel = os.path.relpath(src, TOP)
            if _copy(src, os.path.join(OUT, rel)):
                inter_n += 1; inter_bytes += os.path.getsize(src)

    for f in ROOT_FILES:
        src = os.path.join(TOP, "repo", f)
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(OUT, f)); n_root += 1
            b_root += os.path.getsize(src)
        else:
            print("   root file missing: %s" % f)

    # the deposited panel data and the ledgers travel with the code
    pd_src = os.path.join(TOP, "results", "panel_data")
    if os.path.isdir(pd_src):
        for f in sorted(os.listdir(pd_src)):
            dst = os.path.join(OUT, "panel_data", f)
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            shutil.copy2(os.path.join(pd_src, f), dst); n_panel += 1
            b_panel += os.path.getsize(os.path.join(pd_src, f))

    # ---- every in-tree file the closure reads, at the path its reader expects.
    # Shipping the producing SCRIPTS alone is not enough: a fresh clone has none of their
    # outputs, so a figure either dies on the first missing input or -- worse -- skips a panel
    # in silence. Deriving this set rather than listing it means it cannot fall behind the code.
    SUB = bm.SEARCH_DIRS + ["Datasets/pipeline_out", "Datasets/pandepmap_out",
                            "ligand_rebuild/superseded", "SpecQuality"]
    inter_missing = []
    for name in sorted(scripts):
        sp = bm.find_script(name)
        if not sp:
            continue
        reads, _, _, _ = bm.scan(sp)
        for r in reads:
            b = os.path.basename(r)
            if b in EXTERNAL or b in DEPOSITED or b in PIPELINE:
                continue
            hit = None
            for d in SUB:
                cand = os.path.join(TOP, d, b)
                if os.path.exists(cand):
                    hit = cand; break
            if hit is None:
                inter_missing.append(b); continue
            rel = os.path.relpath(hit, TOP)
            dst = os.path.join(OUT, rel)
            if _copy(hit, dst):
                inter_n += 1; inter_bytes += os.path.getsize(hit)
    if inter_missing:
        print("\n  intermediates not on disk (regenerate before depositing): %s"
              % ", ".join(sorted(set(inter_missing))))

    # the deposited intermediates, so the figures run from the deposit even where the
    # producing step cannot be re-executed
    dep_n, dep_bytes, dep_missing, dep_paths = 0, 0, [], {}
    for b in sorted(DEPOSITED):
        hit = None
        for d in bm.SEARCH_DIRS + ["Datasets/pipeline_out", "Datasets/pandepmap_out",
                                   "ligand_rebuild/superseded"]:
            cand = os.path.join(TOP, d, b)
            if os.path.exists(cand):
                hit = cand; break
        if hit is None:
            dep_missing.append(b); continue
        # Keep each intermediate at the relative path its readers expect. A flat data/ directory
        # was tidier but broke tier 2: the scripts locate inputs with a find() that searches "."
        # and "..", so a file moved out of ccle/ or ligand_rebuild/ is simply not found.
        rel = os.path.relpath(hit, TOP)
        dep_paths[b] = rel
        dst = os.path.join(OUT, rel)
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        if _copy(hit, dst):
            dep_n += 1; dep_bytes += os.path.getsize(hit)

    # a machine-readable record of where every edge input comes from
    rec = {b: {"kind": "external", "source": EXTERNAL[b][0], "where": EXTERNAL[b][1]}
           for b in sorted(EXTERNAL)}
    rec.update({b: {"kind": "deposited intermediate", "note": DEPOSITED[b],
                    "path": dep_paths.get(b, "NOT FOUND")} for b in sorted(DEPOSITED)})
    json.dump(rec, open(os.path.join(OUT, "external_inputs.json"), "w"), indent=1)
    if dep_missing:
        print("\n  deposited intermediates NOT FOUND: %s" % ", ".join(dep_missing))

    print("\nrelease/  %d files, %.0f MB"
          % (n_scripts + n_root + n_panel + dep_n + inter_n,
             (b_scripts + b_root + b_panel + dep_bytes + inter_bytes) / 1e6))
    print("  scripts             %3d  (%.2f MB)" % (n_scripts, b_scripts / 1e6))
    print("  root files          %3d  (%.2f MB)" % (n_root, b_root / 1e6))
    print("  panel data+ledgers  %3d  (%.2f MB)   <- tier 1, goes to git" % (n_panel, b_panel / 1e6))
    print("  deposited data      %3d  (%.1f MB)   <- Zenodo only" % (dep_n, dep_bytes / 1e6))
    print("  intermediates       %3d  (%.1f MB)   <- Zenodo only" % (inter_n, inter_bytes / 1e6))
    print("\n  git tier   : %.2f MB   (code + panel data; every figure redrawable)" %
          ((b_scripts + b_root + b_panel) / 1e6))
    print("  Zenodo tier: %.0f MB   (+ all inputs; every figure re-runnable)" %
          ((b_scripts + b_root + b_panel + dep_bytes + inter_bytes) / 1e6))
    # The counters above are what THIS run copied. When the tree could not be cleared they
    # undercount, so report what is actually on disk as well -- that is what gets uploaded.
    tot_n = tot_b = 0
    for root, _, files in os.walk(OUT):
        for f in files:
            tot_n += 1; tot_b += os.path.getsize(os.path.join(root, f))
    print("\n  release/ on disk: %d files, %.0f MB" % (tot_n, tot_b / 1e6))
    return 0


if __name__ == "__main__":
    sys.exit(main())
