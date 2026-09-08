#!/usr/bin/env python3
"""
step54_germline_by_residue.py -- the germline baseline for EVERY residue, not just lysine.

WHY
Figure 1g compares one somatic share against one germline share: acquired lysines are 6.72% of
union missense sites against 3.29% of a healthy individual's missense burden, a 2.04-fold
difference. The same comparison has never been run for the other nineteen residues, and it needs
to be, for two reasons.

  1. CYSTEINE IS THE COMPARATOR THE FIELD ALREADY USES. Acquired cysteines are the precedent this
     paper's premise rests on, and the cysteine literature quotes a somatic enrichment of about
     2.6-fold over baseline frequency. That figure is not computed the way Figure 1g computes
     lysine's, so the two cannot be set side by side as they stand. Running the identical
     calculation for cysteine makes them comparable for the first time.
  2. A 2.04-FOLD ENRICHMENT MEANS NOTHING WITHOUT ITS DISTRIBUTION. If most residues sit near
     2-fold, lysine is unremarkable and the paper should say so. If lysine and cysteine stand out,
     that is a result. Either way the answer belongs in the paper, and at present neither is known.

METHOD -- identical to germline_lysine.py, generalised
One streaming pass over dbNSFP 5.4a. For every missense row on a MANE transcript with a gnomAD
allele frequency above zero, accumulate 2 x AF (diploid expectation) into a running total and into
a bucket keyed by the residue CREATED. Nothing else changes: same file, same filters, same column
indices, same frequency field. germline_lysine.py's lysine total is reproduced as an internal
check, and the script fails loudly if it does not come back at 471.1 +- 1.

RUNTIME -- this is the reason the script exists rather than the analysis
dbNSFP 5.4a GRCh38 is 55 GB compressed, 91,571,221 rows. A full pass takes hours and cannot be run
inside a short-lived shell. Run it detached and come back to it:

    cd ~/Desktop/DepMap/ccle_plex05
    nohup python3 step54_germline_by_residue.py > step54_report.txt 2>&1 &
    tail -f step54_report.txt        # progress, with an ETA, every 10M rows

It prints progress with an ETA every 10M rows and writes germline_by_residue.csv at the end. It
also checkpoints every 10M rows to step54_partial.json, so a run killed part-way still yields an
interpretable partial answer -- the residue SHARES converge long before the absolute totals do.
The germline results are written to disk BEFORE the somatic join, and pandas is imported
defensively, so neither a missing library nor a missing union_subst_matrix.csv can lose the pass.

OUTPUT
  germline_by_residue.csv   res, expected_per_individual, germline_pct, somatic_sites,
                            somatic_pct, ratio
  step54_partial.json       checkpoint, overwritten every 10M rows
"""
import subprocess, collections, sys, os, json, shutil, time
# pandas is used only to join the somatic side at the very end. It is imported defensively so that
# a missing library cannot destroy a multi-hour pass: without it the germline numbers are still
# written, as JSON and as a hand-rolled CSV, and the join can be done later anywhere.
try:
    import pandas as pd
except ImportError:
    pd = None
    print("NOTE: pandas not available. The germline pass will still run and write "
          "germline_by_residue_raw.json; the somatic join is skipped.", file=sys.stderr, flush=True)

DB = "dbNSFP5.4a_grch38.gz"
CAT = "gzcat" if shutil.which("gzcat") else "zcat"
AF_COL, MANE_COL, REF_COL, ALT_COL = 427, 28, 4, 5   # 0-based, as in germline_lysine.py
CHECK_K, CHECK_TOL = 471.1, 1.0

if not os.path.exists(DB):
    sys.exit("dbNSFP not found: %s" % DB)

TOTAL_ROWS = 91571221          # from germline_lysine.py's completed run, used only for the ETA
T0 = time.time()
exp = collections.defaultdict(float)
expAll = 0.0
n = obs = 0
p = subprocess.Popen([CAT, DB], stdout=subprocess.PIPE, bufsize=1 << 22)
p.stdout.readline()
for line in p.stdout:
    n += 1
    f = line.rstrip(b"\n").split(b"\t")
    if len(f) <= AF_COL:
        continue
    ref, alt = f[REF_COL], f[ALT_COL]
    if len(ref) != 1 or len(alt) != 1 or ref == alt:
        continue
    if not ref.isalpha() or not alt.isalpha():
        continue
    if f[MANE_COL] in (b".", b""):
        continue
    try:
        a = float(f[AF_COL])
    except ValueError:
        continue
    if a <= 0:
        continue
    obs += 1
    expAll += 2 * a
    exp[alt.decode()] += 2 * a
    if n % 10000000 == 0:
        json.dump(dict(rows=n, observed=obs, expAll=expAll,
                       exp={k: v for k, v in exp.items()}),
                  open("step54_partial.json", "w"), indent=1)
        el = time.time() - T0
        eta = el * (TOTAL_ROWS - n) / n / 60.0
        print("  ...%d of ~%d rows (%.0f%%), %.0f min elapsed, ~%.0f min left | "
              "expAll=%.0f K=%.1f C=%.1f  [K/C share %.2f%% / %.2f%%]"
              % (n, TOTAL_ROWS, 100.0 * n / TOTAL_ROWS, el / 60.0, eta, expAll,
                 exp.get("K", 0), exp.get("C", 0),
                 100 * exp.get("K", 0) / expAll if expAll else 0,
                 100 * exp.get("C", 0) / expAll if expAll else 0),
              file=sys.stderr, flush=True)

print("dbNSFP 5.4a, MANE transcripts, gnomAD v4.1 joint (~807k individuals)")
print("  rows scanned             : %d" % n)
print("  observed missense (AF>0) : %d" % obs)
print("  expected missense variants per individual: %.0f" % expAll)
gotK = exp.get("K", 0.0)
print("\n  CHECK against germline_lysine.py: lysine %.1f (expected %.1f)" % (gotK, CHECK_K))
if abs(gotK - CHECK_K) > CHECK_TOL:
    print("  !! MISMATCH -- the generalisation has changed the lysine number. Do not use this "
          "output until the difference is explained.")

# Germline numbers are written before anything that could fail, so the pass is never wasted.
json.dump(dict(rows=n, observed=obs, expAll=expAll, exp=dict(exp)),
          open("germline_by_residue_raw.json", "w"), indent=1)
with open("germline_by_residue_germline_only.csv", "w") as fh:
    fh.write("res,expected_per_individual,germline_pct\n")
    for a in sorted(exp, key=lambda k: -exp[k]):
        fh.write("%s,%.4f,%.5f\n" % (a, exp[a], 100 * exp[a] / expAll if expAll else 0))
print("\nwrote germline_by_residue_raw.json and germline_by_residue_germline_only.csv")
if pd is None or not os.path.exists("union_subst_matrix.csv"):
    print("somatic join skipped (pandas missing or union_subst_matrix.csv not found); "
          "the germline half is safe on disk.")
    sys.exit(0)

S = pd.read_csv("union_subst_matrix.csv")
som = S.groupby("alt_aa").n_sites.sum()
somtot = int(S.n_sites.sum())
rows = []
for a in sorted(exp, key=lambda k: -exp[k]):
    if a not in som.index:
        continue
    gpct = 100 * exp[a] / expAll
    spct = 100 * som[a] / somtot
    rows.append(dict(res=a, expected_per_individual=round(exp[a], 2), germline_pct=round(gpct, 3),
                     somatic_sites=int(som[a]), somatic_pct=round(spct, 3),
                     ratio=round(spct / gpct, 3) if gpct else None))
R = pd.DataFrame(rows).sort_values("ratio", ascending=False)
R.to_csv("germline_by_residue.csv", index=False)
print("\nSOMATIC SHARE AGAINST GERMLINE SHARE, EVERY RESIDUE, RANKED BY RATIO")
print(R.to_string(index=False))
print("\n  lysine and cysteine, side by side on the identical calculation:")
for a in ("K", "C"):
    r = R[R.res == a]
    if len(r):
        r = r.iloc[0]
        print("    %s: somatic %.2f%%, germline %.2f%%, ratio %.2f"
              % (a, r.somatic_pct, r.germline_pct, r.ratio))
print("\nwrote germline_by_residue.csv")
