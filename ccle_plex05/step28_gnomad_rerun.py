#!/usr/bin/env python3
"""
step28_gnomad_rerun.py — settle the two remaining Section 1 numbers in one pass over dbNSFP.

TASK 1 — a single, defensible somatic percentage
The previous annotation wrote dbNSFP's missing-data marker (".") for two different situations:
a site that was looked up and found absent from gnomAD, and a site that was never matched to a
dbNSFP record at all. Both landed in the "absent" band, so "91.0% confidently somatic" is an
upper bound that silently includes 21,582 sites (7.3%) which were never assessed. The lower
bound is 83.7%. This pass distinguishes the two by recording match status explicitly:

    not_in_dbnsfp        no record for this gene/position/substitution   -> NOT ASSESSED
    in_dbnsfp_no_gnomad  record exists, no gnomAD allele frequency       -> confidently absent
    in_dbnsfp_with_AF    record exists, gnomAD frequency present         -> banded as before

TASK 2 — the germline baseline, with a denominator that matches the somatic one
The draft compares a somatic rate of 6.72% against a germline rate of 3.41% and states ~500
lysine-creating variants per healthy individual. Those numbers are only comparable if the
denominators are. They are computed here as an explicit parallel:

    somatic  = acquired-lysine sites / all union missense sites
             = 295,411 / 4,395,403 = 6.72%
    germline = expected lysine-creating variants per individual
               / expected missense variants per individual
             where expected count = sum over variants of 2 x allele frequency
               (2 because each individual carries two alleles)

Both are then "what fraction of the missense burden creates a lysine", one somatic and one
inherited, and the twofold comparison is legitimate. If the germline figure cannot be derived
cleanly it should be CUT rather than defended — the somatic claim rests on the per-site
annotation (Task 1), not on the population rate.

METHOD
One streaming pass over dbNSFP. For every missense row: accumulate 2 x AF into the all-missense
total, and if the variant creates a lysine, into the lysine total and into a lookup keyed by
(gene, protein position) for the union match. dbNSFP stores one row per genomic variant, so
several rows can encode the same protein-level substitution; the union is protein-level, so
rows are collapsed to (gene, aapos, aaalt) taking the maximum AF, and the germline sums are
computed over GENOMIC rows, which is correct for a per-individual expectation.

USAGE
    DBNSFP=/path/to/dbNSFP5.4a_grch38.gz python3 step28_gnomad_rerun.py

OUTPUT
    union_sites_gnomad_v2.tsv   union sites with match_status and AF
    step28_report.txt
"""
import os, sys, gzip, time
import numpy as np, pandas as pd
import warnings
warnings.filterwarnings("ignore")

DB = os.environ.get("DBNSFP", "dbNSFP5.4a_grch38.gz")
UNION = os.environ.get("UNION", "union_sites.tsv")
# TWO frequency fields, deliberately, because they trade off against each other:
#   gnomAD4.1_joint_AF            ~800k individuals — best power to detect rare germline
#                                 variants, but INCLUDES cancer cohorts, so a somatic mutation
#                                 present in a tumour sample appears germline and would be
#                                 wrongly excluded from the somatic count
#   gnomAD2.1.1_exomes_non_cancer_AF  ~118k individuals — no cancer cohorts, so no tumour
#                                 contamination, but misses rarer germline variants and will
#                                 therefore OVERSTATE how somatic a site is
# A site absent from BOTH is confidently somatic. Sites where they disagree are reported
# separately rather than resolved by picking one.
AF_PRIMARY = ["gnomAD4.1_joint_AF", "gnomAD_joint_AF", "gnomAD4.1_exomes_AF"]
AF_NONCANCER = ["gnomAD2.1.1_exomes_non_cancer_AF", "gnomAD_exomes_non_cancer_AF"]

log = []
def say(m):
    print(m, flush=True); log.append(m)

say("=" * 78)
say("STEP 28 — SINGLE SOMATIC PERCENTAGE + GERMLINE BASELINE")
say("=" * 78)

if not os.path.exists(DB):
    say("dbNSFP not found at %s — set DBNSFP=/path/to/dbNSFP5.4a_grch38.gz" % DB)
    sys.exit(1)

u = pd.read_csv(UNION, sep="\t")
say("union: %d sites, %d genes" % (len(u), u.gene.nunique()))
want = set(zip(u.gene.astype(str), u.position.astype(int)))

# ---- one streaming pass
say("")
say("streaming dbNSFP (this is a single pass over the whole file)...")
t0 = time.time()
with gzip.open(DB, "rt", errors="replace") as f:
    header = f.readline().rstrip("\n").split("\t")
    idx = {c: i for i, c in enumerate(header)}
    need = ["genename", "aaref", "aaalt", "aapos"]
    missing = [c for c in need if c not in idx]
    if missing:
        say("!! required columns missing from dbNSFP header: %s" % missing)
        say("   header sample: %s" % header[:25]); sys.exit(1)
    af_col = next((c for c in AF_PRIMARY if c in idx), None)
    nc_col = next((c for c in AF_NONCANCER if c in idx), None)
    if af_col is None:
        cands = [c for c in header if "AF" in c and "gnomAD" in c]
        say("!! no primary AF field found. gnomAD AF-like columns present:")
        say("   %s" % cands[:12]); sys.exit(1)
    say("  primary AF (all cohorts) : %s" % af_col)
    say("  non-cancer AF            : %s" % (nc_col or "NOT FOUND — proceeding with primary only"))
    gi, ri, ai, pi = idx["genename"], idx["aaref"], idx["aaalt"], idx["aapos"]
    fi = idx[af_col]
    ni = idx[nc_col] if nc_col else None

    sum2af_all = 0.0
    sum2af_lys = 0.0
    n_missense_rows = 0
    n_lys_rows = 0
    best = {}                       # (gene, aapos) -> max AF among K-creating rows
    seen_prot_all = set()
    seen_prot_lys = set()

    for n, line in enumerate(f, 1):
        if n % 5_000_000 == 0:
            say("    %d M rows | %.0f min" % (n // 1_000_000, (time.time() - t0) / 60))
        p = line.rstrip("\n").split("\t")
        if len(p) <= fi:
            continue
        aaref, aaalt = p[ri], p[ai]
        if aaref == "." or aaalt == "." or aaref == aaalt:
            continue
        n_missense_rows += 1
        def parse_af(field_i):
            if field_i is None or len(p) <= field_i:
                return None
            r = p[field_i]
            if r in (".", "", "NA"):
                return None
            try:
                return max(float(x) for x in r.split(";") if x not in (".", ""))
            except ValueError:
                return None
        af_v = parse_af(fi)
        nc_v = parse_af(ni)
        raw = p[fi] if len(p) > fi else "."
        af = af_v if af_v is not None else 0.0
        sum2af_all += 2.0 * af
        genes = [g for g in p[gi].split(";") if g and g != "."]
        pos_s = [x for x in p[pi].split(";") if x not in (".", "")]
        for g in genes:
            for ps in pos_s:
                try:
                    pos = int(ps)
                except ValueError:
                    continue
                seen_prot_all.add((g, pos, aaalt))
                if aaalt == "K":
                    seen_prot_lys.add((g, pos))
                    if (g, pos) in want:
                        k = (g, pos)
                        prev = best.get(k)
                        cand = ("in_dbnsfp_with_AF" if af_v is not None else "in_dbnsfp_no_gnomad",
                                af_v if af_v is not None else 0.0,
                                nc_v if nc_v is not None else 0.0)
                        if prev is None or cand[1] > prev[1]:
                            best[k] = cand
        if aaalt == "K":
            n_lys_rows += 1
            sum2af_lys += 2.0 * af

say("  pass complete in %.0f min | %d missense rows | %d lysine-creating rows"
    % ((time.time() - t0) / 60, n_missense_rows, n_lys_rows))

# ---- Task 1
keys = list(zip(u.gene.astype(str), u.position.astype(int)))
u["match_status"] = [best.get(k, ("not_in_dbnsfp", np.nan, np.nan))[0] for k in keys]
u["gnomad_AF"] = [best.get(k, ("", np.nan, np.nan))[1] for k in keys]
u["gnomad_noncancer_AF"] = [best.get(k, ("", np.nan, np.nan))[2] for k in keys]

def band(row):
    if row.match_status == "not_in_dbnsfp":
        return "not_assessed"
    if row.match_status == "in_dbnsfp_no_gnomad":
        return "absent"
    a = row.gnomad_AF
    if a is None or a != a or a == 0:
        return "absent"
    if a < 1e-5:  return "lt_1e-5"
    if a < 1e-4:  return "1e-5_1e-4"
    if a < 1e-3:  return "1e-4_1e-3"
    return "gt_1e-3"

u["band"] = u.apply(band, axis=1)
u.to_csv("union_sites_gnomad_v2.tsv", sep="\t", index=False)

say("")
say("TASK 1 — MATCH STATUS")
for k, v in u.match_status.value_counts().items():
    say("  %-22s %8d  %5.2f%%" % (k, v, 100 * v / len(u)))
say("")
say("BANDS")
for k in ("absent", "lt_1e-5", "1e-5_1e-4", "1e-4_1e-3", "gt_1e-3", "not_assessed"):
    v = int((u.band == k).sum())
    say("  %-14s %8d  %5.2f%%" % (k, v, 100 * v / len(u)))

assessed = u[u.band != "not_assessed"]
som = ((assessed.band.isin(["absent", "lt_1e-5"])).sum())
say("")
say("  confidently somatic (absent or AF<1e-5), OF ASSESSED SITES:")
say("    %d of %d = %.2f%%" % (som, len(assessed), 100 * som / max(len(assessed), 1)))
say("  as a fraction of ALL union sites: %.2f%%" % (100 * som / len(u)))
say("  not assessed: %d (%.2f%%)" % (int((u.band == 'not_assessed').sum()),
                                     100 * (u.band == 'not_assessed').mean()))
say("")
say("  -> report the percentage OF ASSESSED SITES, with the unassessed fraction stated.")
say("     That is one number with a stated scope, not a range.")

# agreement between the two frequency sources
a = pd.to_numeric(assessed.gnomad_AF, errors="coerce").fillna(0)
b = pd.to_numeric(assessed.gnomad_noncancer_AF, errors="coerce").fillna(0)
say("")
say("  AGREEMENT BETWEEN FREQUENCY SOURCES (assessed sites only)")
say("    absent from BOTH (confidently somatic)      : %d (%.2f%%)"
    % (int(((a == 0) & (b == 0)).sum()), 100 * ((a == 0) & (b == 0)).mean()))
say("    present in all-cohort only (cancer cohorts?): %d" % int(((a > 0) & (b == 0)).sum()))
say("    present in non-cancer only                  : %d" % int(((a == 0) & (b > 0)).sum()))
say("    present in both                             : %d" % int(((a > 0) & (b > 0)).sum()))
say("    A site seen ONLY in the all-cohort set may be a somatic variant captured from a")
say("    tumour sample, not a germline polymorphism. Report that count.")

# ---- Task 2
say("")
say("TASK 2 — GERMLINE BASELINE")
say("  expected variants per individual = sum over variants of 2 x allele frequency")
say("    all missense        : %10.1f" % sum2af_all)
say("    lysine-creating     : %10.1f" % sum2af_lys)
if sum2af_all > 0:
    frac = 100 * sum2af_lys / sum2af_all
    say("    lysine-creating share of germline missense burden: %.2f%%" % frac)
    say("")
    say("  COMPARISON, matched denominators:")
    say("    somatic  : 295,411 / 4,395,403 union missense sites = 6.72%")
    say("    germline : %.2f%% of a healthy individual's missense burden" % frac)
    say("    ratio    : %.2fx" % (6.72 / frac) if frac > 0 else "")
say("")
say("  CAVEAT: the somatic figure counts DISTINCT SITES observed across many tumours; the")
say("  germline figure is an EXPECTED COUNT within one individual. They answer 'what fraction")
say("  of missense change creates a lysine' in each setting, which is the intended comparison,")
say("  but they are not the same kind of count. If this distinction cannot be made cleanly in")
say("  the text, CUT the germline comparison — the somatic claim rests on Task 1, not on it.")

open("step28_report.txt", "w").write("\n".join(log) + "\n")
say("")
say("wrote union_sites_gnomad_v2.tsv and step28_report.txt")
