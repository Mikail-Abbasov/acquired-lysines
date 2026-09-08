#!/usr/bin/env python3
"""
build_atlas_v2.py — lysine ligandability atlas, with protein class and pathway.

CHANGES FROM v1
  ADDED    protein family, EC number, functional class, Reactome pathways, GO process
  REMOVED  the cancer cross-reference columns. They rested on a gene-symbol to UniProt
           mapping that disagrees with the variant callers' transcripts for 8.2% of checkable
           sites, concentrated in large multi-isoform proteins (UNC79 91% mismatch, ARID1B 98%,
           MYH11 73%). They were also uninformative: `protein_acquires_lysine` was true for 97%
           of sites, and a spatial test found ligandable lysines no closer to acquired ones than
           chance (15.9 A observed vs 15.0 A permuted null, 45.5% closer than their own null).
           Restoring that cross-reference needs a transcript-aware mapping, not a patch.

WHAT IS UNAFFECTED BY THAT PROBLEM
Occupancy, compounds tested, chemotype breadth and every structural feature are keyed on UniProt
accession throughout and never pass through gene symbols. Those columns were always sound.

PROTEIN CLASS is derived from UniProt keywords and EC number, in a fixed precedence order, so a
kinase that is also ATP-binding is classed as a kinase. The rule is crude but explicit, and the
underlying keywords are kept so anyone can reclassify.

USAGE
    python3 build_atlas_v2.py
"""
import os, re, sys, time, urllib.request
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")

TIERS = (60, 80, 90)
CACHE = "uniprot_annotation.tsv"


def find(n):
    for d in (".", "..", os.path.expanduser("~/Desktop/DepMap/ligand_rebuild"),
              os.path.expanduser("~/Desktop/DepMap")):
        p = os.path.join(d, n)
        if os.path.exists(p):
            return p
    return None


occ = pd.read_csv(find("occ_long.csv"), low_memory=False)
print("measurements: %d" % len(occ))

pat = re.compile(r"^(\d+)([a-z])")
def chemotype(ds, c):
    if ds == "NatChem":
        m = pat.match(str(c))
        if not m:
            return None
        n = m.group(1)
        return "cympsa_family" if n == "7" else "mero_family" if n == "28" else "nat_%s" % n
    return "cympsa_family" if ds == "CyMPSA" else "mero_family" if ds == "Meroterpenoid" else None
occ["ct"] = [chemotype(d, c) for d, c in zip(occ.dataset, occ.compound)]

g = occ.groupby(["uniprot", "position"])
A = g.agg(n_compounds_tested=("compound", "nunique"),
          n_datasets=("dataset", "nunique"),
          max_occupancy=("occupancy_pct", "max"),
          median_occupancy=("occupancy_pct", "median")).reset_index()
idx = A.set_index(["uniprot", "position"]).index
for t in TIERS:
    e = occ[occ.occupancy_pct >= t].groupby(["uniprot", "position"]).compound.nunique()
    A["n_engaged_%d" % t] = idx.map(e).fillna(0).astype(int)
ct = occ[occ.occupancy_pct >= 60].groupby(["uniprot", "position"]).ct.nunique()
A["n_chemotypes_engaged"] = idx.map(ct).fillna(0).astype(int)
best = occ.loc[occ.groupby(["uniprot", "position"]).occupancy_pct.idxmax(),
               ["uniprot", "position", "compound", "dataset"]]
A = A.merge(best.rename(columns={"compound": "best_compound", "dataset": "best_dataset"}),
            on=["uniprot", "position"], how="left")
A["datasets"] = idx.map(occ.groupby(["uniprot", "position"]).dataset
                        .apply(lambda s: "+".join(sorted(s.unique()))))

feat = find("features_expanded.csv")
if feat:
    F = pd.read_csv(feat, low_memory=False).drop_duplicates(["uniprot", "position"])
    if "in_pocket" in F:
        F["in_pocket"] = F.in_pocket.astype(str).str.lower().isin(["true", "1", "1.0"]).astype(int)
    cols = [c for c in ("rsa_tien", "nz_sasa", "sasa_sidechain", "crowding_8A", "atoms_10A",
                        "hydrophobic_8A", "aromatic_8A", "polar_8A", "anionic_8A", "cationic_8A",
                        "net_charge_8A", "n_lys_8A", "dist_centroid", "propka_pKa", "in_pocket",
                        "pocket_volume", "plddt_site") if c in F.columns]
    A = A.merge(F[["uniprot", "position"] + cols], on=["uniprot", "position"], how="left")

# ---- annotation from UniProt
if os.path.exists(CACHE):
    ann = pd.read_csv(CACHE, sep="\t")
    print("annotation loaded from cache: %d accessions" % len(ann))
else:
    accs = sorted(A.uniprot.dropna().astype(str).unique())
    print("annotating %d accessions from UniProt ..." % len(accs))
    rows = []
    fields = ("accession,gene_primary,protein_name,protein_families,ec,keyword,go_f,go_p")
    for i in range(0, len(accs), 100):
        url = ("https://rest.uniprot.org/uniprotkb/accessions?accessions=" +
               ",".join(accs[i:i + 100]) + "&fields=" + fields + "&format=tsv")
        for attempt in range(3):
            try:
                txt = urllib.request.urlopen(url, timeout=60).read().decode()
                for ln in txt.rstrip("\n").split("\n")[1:]:
                    f = ln.split("\t")
                    if len(f) >= 8:
                        rows.append(f[:8])
                break
            except Exception as e:
                if attempt == 2:
                    print("  chunk %d failed: %s" % (i, str(e)[:60]))
                time.sleep(2)
        if (i // 100) % 10 == 0:
            print("  %d/%d" % (i, len(accs)))
        time.sleep(0.15)
    ann = pd.DataFrame(rows, columns=["uniprot", "gene", "protein_name", "protein_family",
                                      "ec_number", "keywords", "go_function", "go_process"])
    ann = ann.drop_duplicates("uniprot")
    ann.to_csv(CACHE, sep="\t", index=False)
    print("annotated %d of %d" % (len(ann), len(accs)))

# ---- classification comes from the curated sources, not from rules written here.
# Three vocabularies, each reported as its authors assigned it:
#   UniProt keywords, split by the CATEGORY UniProt itself gives them (Molecular function,
#     Biological process, Ligand) — 197, 495 and 69 controlled terms respectively
#   UniProt protein family, the curated family/superfamily assignment
#   Reactome pathways, plus the top-level pathway each belongs to, obtained by walking
#     Reactome's own hierarchy
# Nothing here is a judgement of ours about what a protein "is".

def fetch(url):
    return urllib.request.urlopen(
        urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (research script)"}),
        timeout=180).read().decode("utf8", errors="replace")

KWCACHE = "uniprot_keyword_categories.tsv"
if os.path.exists(KWCACHE):
    kw = pd.read_csv(KWCACHE, sep="\t")
else:
    print("fetching the UniProt keyword vocabulary ...")
    d = fetch("https://rest.uniprot.org/keywords/stream?query=*&format=tsv&fields=id,name,category")
    rows = [l.split("\t") for l in d.strip().split("\n")[1:] if len(l.split("\t")) >= 3]
    kw = pd.DataFrame(rows, columns=["kw_id", "keyword", "category"])
    kw.to_csv(KWCACHE, sep="\t", index=False)
    print("  %d keywords in %d categories" % (len(kw), kw.category.nunique()))
KWCAT = dict(zip(kw.keyword, kw.category))

def split_keywords(kwstring, category):
    """Return the protein's keywords belonging to one UniProt category, as UniProt assigned them."""
    if not isinstance(kwstring, str) or not kwstring.strip():
        return ""
    hits = [k.strip() for k in kwstring.split(";") if KWCAT.get(k.strip()) == category]
    return "; ".join(hits)

ann["molecular_function"] = [split_keywords(k, "Molecular function") for k in ann.keywords]
ann["biological_process"] = [split_keywords(k, "Biological process") for k in ann.keywords]
ann["ligand"] = [split_keywords(k, "Ligand") for k in ann.keywords]
for c in ("molecular_function", "biological_process", "ligand"):
    ann[c] = ann[c].replace("", "not annotated")
ann["protein_family"] = ann.protein_family.fillna("no family assigned")

# ---- Reactome: every pathway, and the top-level pathway each belongs to
RCACHE, HCACHE = "reactome_pathways.tsv", "reactome_toplevel.tsv"
# the cache is only reused if it has the columns this version expects — an earlier version
# wrote the same filename with a different schema, which read back cleanly and then failed
if os.path.exists(RCACHE):
    rx = pd.read_csv(RCACHE, sep="\t")
    if "pw_id" not in rx.columns:
        print("  cached %s predates this version; refetching" % RCACHE)
        rx = None
else:
    rx = None
if rx is None:
    print("fetching Reactome pathway names ...")
    txt = fetch("https://reactome.org/download/current/UniProt2Reactome.txt")
    rows = []
    for ln in txt.split("\n"):
        f = ln.split("\t")
        if len(f) >= 6 and f[5].strip() == "Homo sapiens":
            rows.append((f[0], f[1], f[3]))
    rx = pd.DataFrame(rows, columns=["uniprot", "pw_id", "pathway"]).drop_duplicates()
    rx.to_csv(RCACHE, sep="\t", index=False)
    print("  %d human protein-pathway pairs" % len(rx))

top = None
if os.path.exists(HCACHE):
    top = pd.read_csv(HCACHE, sep="\t")
    if "top_level_pathway" not in top.columns:
        top = None
if top is None:
    print("walking the Reactome hierarchy to top-level pathways ...")
    rel = fetch("https://reactome.org/download/current/ReactomePathwaysRelation.txt")
    parent = {}
    for ln in rel.split("\n"):
        f = ln.split("\t")
        if len(f) == 2 and f[0].startswith("R-HSA"):
            parent.setdefault(f[1].strip(), []).append(f[0].strip())
    names = fetch("https://reactome.org/download/current/ReactomePathways.txt")
    pname = {}
    for ln in names.split("\n"):
        f = ln.split("\t")
        if len(f) >= 3 and f[2].strip() == "Homo sapiens":
            pname[f[0].strip()] = f[1].strip()
    def root_of(pid, seen=None):
        seen = seen or set()
        if pid in seen:
            return pid
        seen.add(pid)
        ps = parent.get(pid)
        return pid if not ps else root_of(ps[0], seen)
    ids = sorted(set(rx.pw_id.dropna()))
    top = pd.DataFrame([(p, pname.get(root_of(p), "")) for p in ids],
                       columns=["pw_id", "top_level_pathway"])
    top.to_csv(HCACHE, sep="\t", index=False)
    print("  %d pathways mapped to %d top-level categories"
          % (len(top), top.top_level_pathway.nunique()))

rx = rx.merge(top, on="pw_id", how="left")
agg = rx.groupby("uniprot").agg(
    pathways=("pathway", lambda s: "; ".join(sorted(set(s)))),
    n_pathways=("pathway", "nunique"),
    pathway_categories=("top_level_pathway", lambda s: "; ".join(sorted({x for x in s if x}))))
ann = ann.merge(agg, on="uniprot", how="left")
ann["n_pathways"] = ann.n_pathways.fillna(0).astype(int)
for c in ("pathways", "pathway_categories"):
    ann[c] = ann[c].fillna("not in Reactome")

A = A.merge(ann, on="uniprot", how="left")
# ---- SMILES for the best compound at each site, and a companion compound table.
# The compound names in occ_long carry a concentration for the NatChem series ("17g (500uM)");
# the SMILES sheet lists the compound itself ("17g"), so the concentration is stripped for the
# join and kept in the atlas.
SM = None
for p in ("SMILES.xlsx", os.path.expanduser("~/Desktop/SMILES.xlsx"),
          os.path.expanduser("~/Desktop/DepMap/SMILES.xlsx")):
    if os.path.exists(p):
        SM = p; break
if SM:
    sheets = {"NatChem": "NatChem", "Meroterpenoids": "Meroterpenoid", "CyMPSA": "CyMPSA"}
    rows = []
    xl = pd.ExcelFile(SM)
    for sh, ds in sheets.items():
        if sh not in xl.sheet_names:
            continue
        d = xl.parse(sh).dropna(subset=["SMILES"])
        for r in d.itertuples():
            rows.append((ds, str(r.Compounds).strip(), str(r.SMILES).strip()))
    smi = pd.DataFrame(rows, columns=["dataset", "compound_base", "SMILES"]).drop_duplicates(
        ["dataset", "compound_base"])
    print("SMILES loaded: %d compounds from %s" % (len(smi), os.path.basename(SM)))

    conc = re.compile(r"\s*\(([^)]*)\)\s*$")
    def base_and_conc(c):
        c = str(c)
        m = conc.search(c)
        return (conc.sub("", c).strip(), m.group(1) if m else "")
    A["best_compound_base"], A["best_compound_conc"] = zip(*A.best_compound.map(base_and_conc))
    A = A.merge(smi.rename(columns={"dataset": "best_dataset", "compound_base": "best_compound_base",
                                    "SMILES": "best_compound_SMILES"}),
                on=["best_dataset", "best_compound_base"], how="left")
    print("  SMILES matched for %d of %d sites" % (int(A.best_compound_SMILES.notna().sum()), len(A)))
    miss = sorted(set(A.loc[A.best_compound_SMILES.isna(), "best_compound_base"].dropna()))[:8]
    if miss:
        print("  unmatched compound names:", ", ".join(miss))

    # companion table: every compound, its structure, and what it engaged
    occ["compound_base"] = [base_and_conc(c)[0] for c in occ.compound]
    cs = occ.groupby(["dataset", "compound_base"]).agg(
        n_sites_tested=("occupancy_pct", "size"),
        n_sites_engaged_60=("occupancy_pct", lambda s: int((s >= 60).sum())),
        max_occupancy=("occupancy_pct", "max"),
        median_occupancy=("occupancy_pct", "median")).reset_index()
    cs = cs.merge(smi, on=["dataset", "compound_base"], how="left")
    cs = cs.sort_values("n_sites_engaged_60", ascending=False)
    with open("compound_table.csv", "w") as fh:
        fh.write("# Compounds in the lysine ligandability atlas\n")
        fh.write("# One row per compound. n_sites_engaged_60 counts site-observations at or\n")
        fh.write("# above 60% occupancy, so a compound tested at several concentrations is\n")
        fh.write("# pooled here; the atlas keeps the concentration of the best observation.\n")
        cs.to_csv(fh, index=False)
    print("  wrote compound_table.csv — %d compounds" % len(cs))
else:
    print("SMILES.xlsx not found — best_compound_SMILES omitted")

A = A.sort_values(["max_occupancy", "n_chemotypes_engaged"], ascending=False)

HEADER = """# Lysine ligandability atlas v2
# One row per lysine, merged from three chemoproteomic datasets on a common occupancy scale.
#
# READ BEFORE FILTERING
#   occupancy is COMPETITION, not affinity. It does not establish a binding constant, and does
#     not separate a molecule binding at the site from one reacting there more readily.
#   a zero in n_engaged_60 means nothing unless n_compounds_tested is large. Check that column.
#   structural columns are AlphaFold-derived. Below plddt_site 70 the geometric ones are
#     unreliable; the charge ones hold up better.
#   classification is taken from curated sources as their authors assigned it, not from rules
#     written for this table:
#       molecular_function / biological_process / ligand  UniProt keywords, split by the
#         category UniProt itself gives each keyword
#       protein_family                                    UniProt curated family assignment
#       pathways / pathway_categories                     ALL Reactome pathways for the protein,
#         and the top-level Reactome pathway each belongs to, from Reactome's own hierarchy
#     GO molecular-function and process strings are kept alongside for anything not covered.
#   best_compound_SMILES gives the structure of the single strongest observation at each site.
#     compound_table.csv is the companion file: every compound, its SMILES, and how many sites
#     it engaged.
#
# A CAVEAT WORTH KNOWING BEFORE USING THIS TO PICK TARGETS
#   Sites in ribosomal proteins engage at 19.0% against 9.3% elsewhere, and this survives
#   matching on how deeply each site was screened. But their structural profile is the OPPOSITE
#   of what predicts engagement elsewhere — more exposed, less crowded, fewer hydrophobic
#   neighbours, less often in a pocket — while carrying the cationic, acid-poor environment that
#   governs probe REACTIVITY. Competition measurements cannot separate a molecule binding at a
#   site from one reacting there more readily, so high occupancy at these sites may reflect
#   hyperreactivity rather than ligandability.
#
# NOT INCLUDED, deliberately: any cancer-mutation cross-reference. The gene-symbol mapping it
#   would need disagrees with the variant callers' transcripts for 8% of checkable sites, and a
#   spatial test found ligandable lysines no closer to cancer-acquired lysines than chance.
"""
with open("lysine_ligandability_atlas_v2.csv", "w") as fh:
    fh.write(HEADER)
    A.to_csv(fh, index=False)

print("\nwrote lysine_ligandability_atlas_v2.csv — %d sites, %d columns" % (len(A), A.shape[1]))
for col, label in (("molecular_function", "UniProt molecular function"),
                   ("pathway_categories", "Reactome top-level pathway")):
    print("\nengagement by %s (a site counts once per term it carries):" % label)
    ex = A.assign(eng=A.n_engaged_60 > 0)
    ex = ex.assign(term=ex[col].astype(str).str.split(";")).explode("term")
    ex["term"] = ex.term.str.strip()
    ex = ex[~ex.term.isin(["not annotated", "not in Reactome", ""])]
    t = ex.groupby("term").agg(sites=("eng", "size"), engaged=("eng", "sum"))
    t = t[t.sites >= 150]
    t["pct"] = 100 * t.engaged / t.sites
    print(t.sort_values("engaged", ascending=False).head(12).to_string())
print("\npathways per protein: median %d, max %d"
      % (A.n_pathways.median(), A.n_pathways.max()))
print("\ntop sites:")
show = ["gene", "position", "molecular_function", "max_occupancy", "n_chemotypes_engaged",
        "n_compounds_tested"]
print(A.head(10)[[c for c in show if c in A.columns]].to_string(index=False))
