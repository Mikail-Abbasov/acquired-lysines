#!/usr/bin/env python3
"""
tcga_burden_by_lineage.py — acquired-lysine burden per tumour, split by TCGA cancer type.

WHY THIS EXISTS
The manuscript claims melanoma carries a median of 39 acquired lysines per tumour against 3
across TCGA, a 13-fold excess. **That figure has no derivation anywhere in the pipeline.** It
appears only in a session summary. The pan-TCGA half checks out — mc3_perK.pkl gives a median
of 3 across all 9,092 tumours — but the melanoma half cannot be reproduced, because nothing in
the working directory maps a tumour barcode to a cancer type.

Two other numbers this week had the same problem. One (the 4.12% substitution denominator) was
recomputed and proved correct; another (the AUC ceiling) was recomputed and proved wrong. This
one is being recomputed rather than trusted.

BARCODE -> CANCER TYPE
Rather than rely on a memorised tissue-source-site table, the mapping is retrieved from the GDC
API, which returns submitter_id and project_id for every TCGA case. That makes the mapping
verifiable and reproducible: anyone can re-run this and get the same assignment.

WHAT IS COMPUTED
  - median and mean acquired lysines per tumour, per cancer type
  - the pan-TCGA figure over all tumours, as the comparator
  - SKCM (cutaneous melanoma) against UVM (uveal melanoma)

THE UVM CONTROL IS THE POINT
Cutaneous melanoma is ultraviolet-driven; uveal melanoma is not. If the acquired-lysine excess
is caused by UV mutagenesis rather than being an intrinsic property of melanocytes, SKCM should
be high and UVM should not. That contrast is stronger evidence than the SKCM number alone, and
it is a control the original claim did not include.

NOTE ON DENOMINATORS
mc3_perK.pkl holds two dictionaries: `perK` (7,403 tumours carrying at least one acquired
lysine) and `per` (9,092 tumours in total). Medians computed over `perK` alone would exclude
tumours with zero and inflate every figure. All statistics here use the full 9,092, with zeros
included, and the script reports both so the difference is visible.

USAGE
    python3 tcga_burden_by_lineage.py

OUTPUT
    tcga_burden_by_lineage.csv
    tcga_burden_report.txt
"""
import os, json, time, urllib.parse, urllib.request
import numpy as np, pandas as pd
import warnings
warnings.filterwarnings("ignore")

PKL = os.environ.get("PKL", "mc3_perK.pkl")
CACHE = "gdc_barcode_project.csv"
log = []
def say(m):
    print(m, flush=True); log.append(m)

say("=" * 74)
say("ACQUIRED-LYSINE BURDEN PER TUMOUR, BY TCGA CANCER TYPE")
say("=" * 74)

# ---- barcode -> project, from the GDC API (cached after first run)
if os.path.exists(CACHE):
    m = pd.read_csv(CACHE)
    say("barcode->project mapping loaded from cache: %d cases" % len(m))
else:
    say("querying the GDC API for TCGA case -> project mapping ...")
    rows, size, frm = [], 1000, 0
    flt = json.dumps({"op": "in",
                      "content": {"field": "project.program.name", "value": ["TCGA"]}})
    while True:
        url = ("https://api.gdc.cancer.gov/cases?filters=" + urllib.parse.quote(flt) +
               "&fields=submitter_id,project.project_id&format=json"
               "&size=%d&from=%d" % (size, frm))
        r = json.loads(urllib.request.urlopen(url, timeout=120).read().decode())
        hits = r["data"]["hits"]
        if not hits:
            break
        for h in hits:
            rows.append((h["submitter_id"], h["project"]["project_id"]))
        frm += size
        say("  ...%d" % len(rows))
        if frm >= r["data"]["pagination"]["total"]:
            break
        time.sleep(0.2)
    m = pd.DataFrame(rows, columns=["barcode", "project"]).drop_duplicates("barcode")
    m.to_csv(CACHE, index=False)
    say("retrieved %d cases -> %s" % (len(m), CACHE))

m["cancer"] = m.project.str.replace("TCGA-", "", regex=False)

# ---- burden
d = pd.read_pickle(PKL)
perK, per = d["perK"], d["per"]
say("")
say("mc3_perK.pkl: %d tumours with >=1 acquired lysine, %d tumours total"
    % (len(perK), len(per)))

b = pd.DataFrame({"barcode": list(per.keys())})
b["n_acquired_K"] = b.barcode.map(lambda x: perK.get(x, 0))
b = b.merge(m[["barcode", "cancer"]], on="barcode", how="left")
say("mapped to a cancer type: %d of %d (%.1f%%)"
    % (b.cancer.notna().sum(), len(b), 100 * b.cancer.notna().mean()))
unm = b[b.cancer.isna()]
if len(unm):
    say("  unmapped examples: %s" % ", ".join(unm.barcode.head(4)))

say("")
say("PAN-TCGA, all %d tumours: median %.0f, mean %.1f"
    % (len(b), b.n_acquired_K.median(), b.n_acquired_K.mean()))
say("  restricted to tumours with >=1: median %.0f, mean %.1f  <-- do NOT quote this"
    % (b.loc[b.n_acquired_K > 0, "n_acquired_K"].median(),
       b.loc[b.n_acquired_K > 0, "n_acquired_K"].mean()))

g = (b.dropna(subset=["cancer"]).groupby("cancer").n_acquired_K
     .agg(n="size", median="median", mean="mean",
          q25=lambda s: s.quantile(.25), q75=lambda s: s.quantile(.75))
     .reset_index().sort_values("median", ascending=False))
g["fold_vs_pan"] = g["median"] / max(b.n_acquired_K.median(), 1e-9)
g.to_csv("tcga_burden_by_lineage.csv", index=False)

say("")
say("BY CANCER TYPE — top 12 by median")
say("  %-8s %6s %8s %8s %10s %8s" % ("type", "n", "median", "mean", "IQR", "fold"))
for _, r in g.head(12).iterrows():
    say("  %-8s %6d %8.0f %8.1f %4.0f-%-5.0f %7.1fx"
        % (r.cancer, r.n, r["median"], r["mean"], r.q25, r.q75, r.fold_vs_pan))

say("")
say("THE ULTRAVIOLET CONTRAST")
for code, label in (("SKCM", "cutaneous melanoma (UV-driven)"),
                    ("UVM", "uveal melanoma (NOT UV-driven)")):
    s = g[g.cancer == code]
    if len(s):
        r = s.iloc[0]
        say("  %-5s %-32s n=%4d  median %5.0f  mean %7.1f  %.1fx pan-TCGA"
            % (code, label, r.n, r["median"], r["mean"], r.fold_vs_pan))
    else:
        say("  %-5s not present in the mapped set" % code)
say("")
say("  If SKCM is high and UVM is not, the excess follows the mutagen rather than the cell type.")
say("  That contrast is the control the original claim lacked.")

open("tcga_burden_report.txt", "w").write("\n".join(log) + "\n")
say("")
say("wrote tcga_burden_by_lineage.csv and tcga_burden_report.txt")
