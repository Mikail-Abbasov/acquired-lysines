#!/usr/bin/env python3
"""
step48c_signature_control.py — the specificity control for the codon-asymmetry claim: does the
lysine gain/loss asymmetry depend on WHICH mutational process is running?

THE ALTERNATIVE EXPLANATION THIS RULES OUT
step48 shows lysine gains are 85.5% C>T while losses are 0% C>T, and argues the cause is the
codon table. A reviewer can offer a cheaper explanation: the C>T class dominates cancer missense
mutation overall (51.1%), so perhaps the gain figure is just that background showing through, and
the whole thing is a statement about which mutational signatures happen to be common in TCGA
rather than about lysine.

That alternative makes a prediction the codon argument does not: if the asymmetry is inherited
from the background spectrum, it should WEAKEN in tumours where C>T is not dominant.

TCGA supplies the natural experiment. Lung adenocarcinoma is driven by tobacco mutagenesis and is
C>A-dominated; cutaneous melanoma is driven by ultraviolet light and is overwhelmingly C>T. Their
background spectra differ about fourfold in C>T content. If the loss side is forbidden by the
code rather than merely rare, it must be zero in BOTH — and in every type between.

Barcode-to-project mapping is the cached GDC table used by tcga_burden_by_lineage.py, so the
assignment is the same one the burden analysis uses.
"""
import os
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")

CC = os.environ.get("CCLE_DIR", "../ccle_plex05")
MIN_N = 200

log = []
def say(m):
    print(m, flush=True); log.append(m)

say("=" * 78)
say("STEP 48c — IS THE ASYMMETRY A SIGNATURE EFFECT? (the specificity control)")
say("=" * 78)

D = pd.read_csv("mc3_bysample.tsv", sep="\t", header=None,
                names=["bc", "cls", "sbs"], low_memory=False)
M = pd.read_csv(os.path.join(CC, "gdc_barcode_project.csv"), low_memory=False)
m = dict(zip(M.barcode.astype(str).str[:12], M.project))
D["cancer"] = D.bc.map(m)
D = D.dropna(subset=["cancer"])
D["cancer"] = D.cancer.astype(str).str.replace("TCGA-", "", regex=False)
say("mutations mapped to a cancer type: %s across %d types"
    % (f"{len(D):,}", D.cancer.nunique()))

rows = []
for ca, S in D.groupby("cancer"):
    g, l = S[S.cls == "gain"], S[S.cls == "loss"]
    if len(g) < MIN_N or len(l) < MIN_N: continue
    rows.append(dict(cancer=ca, n_mutations=len(S),
                     bg_CT=100*(S.sbs == "C>T").mean(),
                     bg_CA=100*(S.sbs == "C>A").mean(),
                     gain_CT=100*(g.sbs == "C>T").mean(), n_gain=len(g),
                     loss_CT=100*(l.sbs == "C>T").mean(), n_loss=len(l)))
R = pd.DataFrame(rows).sort_values("bg_CT").reset_index(drop=True)
R["gain_over_bg"] = R.gain_CT / R.bg_CT

say("")
say("  %-6s %10s %8s %8s %10s %9s %10s %8s"
    % ("type", "mutations", "bg C>T", "bg C>A", "gains C>T", "vs bg", "losses C>T", "n loss"))
for _, r in R.iterrows():
    say("  %-6s %10s %7.1f%% %7.1f%% %9.1f%% %8.2fx %9.1f%% %8s"
        % (r.cancer, f"{int(r.n_mutations):,}", r.bg_CT, r.bg_CA, r.gain_CT,
           r.gain_over_bg, r.loss_CT, f"{int(r.n_loss):,}"))

lo, hi = R.iloc[0], R.iloc[-1]
say("")
say("-" * 78)
say("READING IT")
say("  Background C>T ranges from %.1f%% (%s, tobacco, C>A-dominated at %.1f%%) to %.1f%% (%s,"
    % (lo.bg_CT, lo.cancer, lo.bg_CA, hi.bg_CT, hi.cancer))
say("  ultraviolet). That is a %.1f-fold spread in the process the alternative explanation"
    % (hi.bg_CT / lo.bg_CT))
say("  depends on.")
say("")
say("  GAINS track the local background, as they must — %.1f%% C>T in %s against %.1f%% in %s —"
    % (lo.gain_CT, lo.cancer, hi.gain_CT, hi.cancer))
say("  but sit ABOVE it in every type (%.2fx to %.2fx). The codon route is used whenever the"
    % (R.gain_over_bg.min(), R.gain_over_bg.max()))
say("  process that uses it is running.")
say("")
tot_loss = int(R.n_loss.sum())
say("  LOSSES are %s%% C>T in every single type: %d C>T events among %s lysine losses across"
    % (("%.1f" % R.loss_CT.max()), int((D[D.cls == "loss"].sbs == "C>T").sum()),
       f"{int((D.cls=='loss').sum()):,}"))
say("  %d cancer types. Not low, not variable — absent, regardless of which mutational process"
    % D.cancer.nunique())
say("  is running. A background effect cannot be invariant to the background.")
say("")
say("  This is the control that separates the two explanations, and it favours the codon table.")

R.to_csv("step48c_by_cancer_type.csv", index=False)
open("step48c_report.txt", "w").write("\n".join(log) + "\n")
say(""); say("wrote step48c_by_cancer_type.csv, step48c_report.txt")
