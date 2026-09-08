#!/usr/bin/env python3
"""
union_all_missense.py — count EVERY missense substitution in the union, not only those creating
lysine.

WHY THIS EXISTS
`union_sites_v2.py` filters to `alt_aa == 'K'` at the point of reading each source (line 34 for
DepMap, line 58 for TCGA), so the total missense count was never computed by it. The figure of
4,395,403 total union missense sites — and therefore the claim that E→K accounts for 4.12% of
all missense substitutions and ranks first of 380 — has no derivation in the pipeline.

That claim appears in the second paragraph of the manuscript. It needs a verified denominator
rather than a remembered one.

WHAT THIS DOES
The same union logic over the same three sources, but retaining every ref→alt pair. Output is a
380-cell substitution matrix with counts, from which the rank and percentage follow directly.

DEDUPLICATION must match the lysine pipeline exactly or the denominator will not be comparable:
sites are keyed on GENE|REF|POSITION|ALT with the gene symbol upper-cased and any suffix after
an underscore removed, and a site seen in several sources counts once. Synonymous entries
(ref == alt) are excluded, as are anything failing the amino-acid alphabet.

CHECK BUILT IN
The lysine-creating subset of this matrix must reproduce 295,411 sites. If it does not, the two
pipelines have diverged and neither number should be used until that is resolved.

USAGE
    python3 union_all_missense.py

OUTPUT
    union_subst_matrix.csv     ref_aa, alt_aa, n_sites
    union_all_missense.tsv     gene, ref_aa, position, alt_aa, source flags
    union_all_missense_report.txt
"""
import os, re, sys, csv, gzip, time
import numpy as np, pandas as pd

AA = set("ACDEFGHIKLMNPQRSTVWY")
MC3 = os.environ.get("MC3", "TCGA_mc3-v028-PUBLIC.maf.gz")
PARQ = os.environ.get("PARQ", "depmap_variants_annotated.parquet")
COSMIC = os.environ.get("COSMIC", "COSMIC-Cancer-Mutation-Census-All-Data-v103-GRCh37.tsv")
pAA = re.compile(r"^p\.([A-Z])(\d+)([A-Z])$")
csv.field_size_limit(10 ** 7)

log = []
def say(m):
    print(m, flush=True); log.append(m)


def key(gene, ref, pos, alt):
    return "%s|%s|%d|%s" % (str(gene).upper().split("_")[0], ref, int(pos), alt)


say("=" * 74)
say("UNION SUBSTITUTION MATRIX — all missense, not only lysine-creating")
say("=" * 74)

sites = {}          # key -> [ref, alt, in_depmap, in_tcga, in_cosmic, gene, pos]
t0 = time.time()

# ------------------------------------------------------------------ DepMap
if os.path.exists(PARQ):
    say("")
    say("DepMap ...")
    v = pd.read_parquet(PARQ, columns=["gene", "ref_aa", "position", "alt_aa", "vclass"])
    d = v[v.vclass == "missense"]
    say("  %d missense rows" % len(d))
    n = 0
    for g, r, p, a in zip(d.gene, d.ref_aa, d.position, d.alt_aa):
        if r in AA and a in AA and r != a and p == p:
            k = key(g, r, p, a)
            e = sites.get(k)
            if e is None:
                sites[k] = [r, a, 1, 0, 0, str(g).upper().split("_")[0], int(p)]
            else:
                e[2] = 1
            n += 1
    say("  %d usable -> %d distinct union sites so far" % (n, len(sites)))
else:
    say("!! DepMap parquet not found at %s — set PARQ=" % PARQ)

# ------------------------------------------------------------------ TCGA MC3
if os.path.exists(MC3):
    say("")
    say("TCGA MC3 ... (streaming, this is the slow step)")
    n = kept = 0
    with gzip.open(MC3, "rt", encoding="utf8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            n += 1
            if n % 1_000_000 == 0:
                say("  ...%d rows | %d sites | %.0f min" % (n, len(sites), (time.time() - t0) / 60))
            if row.get("Variant_Classification") != "Missense_Mutation":
                continue
            if (row.get("FILTER") or "").strip() not in ("PASS", "", "."):
                continue
            m = pAA.match((row.get("HGVSp_Short") or "").strip())
            if not m:
                continue
            r, p, a = m.group(1), m.group(2), m.group(3)
            if r not in AA or a not in AA or r == a:
                continue
            g = row.get("Hugo_Symbol") or ""
            k = key(g, r, p, a)
            e = sites.get(k)
            if e is None:
                sites[k] = [r, a, 0, 1, 0, str(g).upper().split("_")[0], int(p)]
            else:
                e[3] = 1
            kept += 1
    say("  %d rows read, %d usable -> %d distinct union sites" % (n, kept, len(sites)))
else:
    say("!! MC3 not found at %s — set MC3=" % MC3)

# ------------------------------------------------------------------ COSMIC
# Column names and filtering copied EXACTLY from union_sites_v2.py. COSMIC does not use HGVS
# p.X###Y notation here; it stores the wild-type and mutant residues and the position in
# separate fields. The length and alphabet checks are what corrected the earlier 295,598 count
# to 295,411 — they drop '*', '', multi-residue entries like 'TQ', indels, and X>X
# non-substitutions. Reproducing them is required for the totals to be comparable.
if COSMIC and os.path.exists(COSMIC):
    say("")
    say("COSMIC ...")
    n = kept = dropped = 0
    with open(COSMIC, newline="", encoding="utf8", errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            n += 1
            if n % 1_000_000 == 0:
                say("  ...%d rows | %d sites | %d dropped | %.0f min"
                    % (n, len(sites), dropped, (time.time() - t0) / 60))
            a = (row.get("AA_MUT_ALLELE_SEQ") or "").strip()
            r = (row.get("AA_WT_ALLELE_SEQ") or "").strip()
            p = (row.get("AA_MUT_START") or "").strip()
            if (r not in AA or len(r) != 1 or a not in AA or len(a) != 1
                    or not p.isdigit() or r == a):
                dropped += 1
                continue
            g = row.get("GENE_NAME") or ""
            k = key(g, r, p, a)
            e = sites.get(k)
            if e is None:
                sites[k] = [r, a, 0, 0, 1, str(g).upper().split("_")[0], int(p)]
            else:
                e[4] = 1
            kept += 1
    say("  %d rows, %d usable, %d dropped -> %d distinct union sites"
        % (n, kept, dropped, len(sites)))
else:
    say("")
    say("!! COSMIC not supplied — set COSMIC=. Totals will NOT match the published figures.")

# ------------------------------------------------------------------ matrix
say("")
say("=" * 74)
rows = [(v[0], v[1]) for v in sites.values()]
mat = pd.DataFrame(rows, columns=["ref_aa", "alt_aa"]).groupby(["ref_aa", "alt_aa"]).size()
mat = mat.rename("n_sites").reset_index().sort_values("n_sites", ascending=False)
total = int(mat.n_sites.sum())
mat["pct_of_all_missense"] = 100 * mat.n_sites / total
mat["rank"] = range(1, len(mat) + 1)
mat.to_csv("union_subst_matrix.csv", index=False)

say("TOTAL DISTINCT UNION MISSENSE SITES: %s" % f"{total:,}")
say("distinct substitution pairs observed: %d (of 380 possible)" % len(mat))
say("")
say("TOP 12 SUBSTITUTIONS")
say("  %-6s %-8s %14s %10s" % ("rank", "subst", "sites", "% of all"))
for _, r in mat.head(12).iterrows():
    mark = "  <-- lysine-creating" if r.alt_aa == "K" else ""
    say("  %-6d %-8s %14s %9.2f%%%s" % (r["rank"], "%s>%s" % (r.ref_aa, r.alt_aa),
                                        f"{int(r.n_sites):,}", r.pct_of_all_missense, mark))

lys = mat[mat.alt_aa == "K"]
say("")
say("LYSINE-CREATING SUBSTITUTIONS: %s sites across %d pairs"
    % (f"{int(lys.n_sites.sum()):,}", len(lys)))
say("  E>K: rank %d of %d, %s sites, %.2f%% of all missense"
    % (int(mat.loc[(mat.ref_aa == 'E') & (mat.alt_aa == 'K'), 'rank'].iloc[0]),
       len(mat),
       f"{int(mat.loc[(mat.ref_aa=='E')&(mat.alt_aa=='K'),'n_sites'].iloc[0]):,}",
       float(mat.loc[(mat.ref_aa=='E')&(mat.alt_aa=='K'),'pct_of_all_missense'].iloc[0])))

say("")
say("CONSISTENCY CHECK against the lysine pipeline")
say("  this run, lysine-creating sites : %s" % f"{int(lys.n_sites.sum()):,}")
say("  union_sites_v2.py reported      : 295,411")
if abs(int(lys.n_sites.sum()) - 295411) > 0:
    say("  !! MISMATCH. The two pipelines have diverged — most likely because COSMIC was not")
    say("     supplied here. Do NOT use either number until this reconciles.")
else:
    say("  MATCH — the denominator is verified and the percentage can be quoted.")

# write the site table for downstream use
out = pd.DataFrame([[v[5], v[0], v[6], v[1], v[2], v[3], v[4]] for v in sites.values()],
                   columns=["gene", "ref_aa", "position", "alt_aa",
                            "in_depmap", "in_tcga", "in_cosmic"])
out.to_csv("union_all_missense.tsv", sep="\t", index=False)
open("union_all_missense_report.txt", "w").write("\n".join(log) + "\n")
say("")
say("wrote union_subst_matrix.csv, union_all_missense.tsv, union_all_missense_report.txt")
say("elapsed %.0f min" % ((time.time() - t0) / 60))
