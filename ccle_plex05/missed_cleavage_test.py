#!/usr/bin/env python3
"""
missed_cleavage_test.py — is "basic centres" measuring charge, or measuring missed cleavages?

THE CONFLICT WITH THE LITERATURE
Our analysis found that peptides with FEWER basic centres are identified more often. The
proteomics literature says the opposite: reagents that add basic groups are designed to RAISE
charge state, and higher charge is generally considered better for fragmentation. Tertiary-amine
tags, DiPyrO and SOT are all sold on exactly that basis.

THE LIKELY EXPLANATION
A correctly cleaved tryptic peptide has about two basic centres: the free N-terminus and the
C-terminal lysine or arginine. Five or six means several INTERNAL K/R residues that trypsin did
not cut. So the variable may be counting missed cleavages rather than charge — and missed
cleavages are independently associated with poorer identification for reasons unrelated to
charge state.

THE TEST
Split the basic centres into two parts that mean different things:
  terminal      the N-terminus plus a C-terminal K/R — present in every properly cleaved peptide
  internal      K/R/H that are not at the C-terminus — each one is a missed cleavage or a His
If the effect lives in the internal count, the variable is digestion completeness, and calling
it charge would be wrong. If it lives in the terminal count, or in histidine specifically, then
it really is about charge.
"""
import os, re
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.metrics import roc_auc_score

log = []
def say(m):
    print(m, flush=True); log.append(m)

D = pd.read_csv("peptide_chem.csv")
D = D[D.workflow == "tmt"].copy()

def decompose(p):
    """Split the basic residues into terminal and internal."""
    c_term_kr = 1 if p[-1] in "KR" else 0
    internal_kr = sum(1 for a in p[:-1] if a in "KR")
    his = p.count("H")
    return pd.Series(dict(n_term=1, c_term_kr=c_term_kr, internal_kr=internal_kr, his=his,
                          missed=internal_kr))

D = pd.concat([D, D.peptide.apply(decompose)], axis=1)
say("=" * 78)
say("IS IT CHARGE, OR IS IT MISSED CLEAVAGES?")
say("=" * 78)
say("%d peptides, %d identified" % (len(D), int(D.found.sum())))

say("")
say("MEDIAN COUNTS")
say("  %-22s %12s %12s" % ("", "identified", "not identified"))
for c, nm in (("basic", "basic centres (all)"), ("c_term_kr", "C-terminal K/R"),
              ("internal_kr", "internal K/R (missed)"), ("his", "histidine"),
              ("L", "length")):
    say("  %-22s %12.1f %12.1f" % (nm, D.loc[D.found, c].median(), D.loc[~D.found, c].median()))

def delta(a, b):
    a, b = np.asarray(a), np.asarray(b)
    if len(a) < 3 or len(b) < 10:
        return np.nan
    return (a[:, None] > b[None, :]).mean() - (a[:, None] < b[None, :]).mean()

say("")
say("EFFECT ON IDENTIFICATION, unadjusted")
for c, nm in (("basic", "basic centres"), ("internal_kr", "internal K/R (missed cleavages)"),
              ("his", "histidine"), ("c_term_kr", "C-terminal K/R")):
    say("  %-34s %+.3f" % (nm, delta(D.loc[D.found, c], D.loc[~D.found, c])))

say("")
say("EFFECT OF BASIC CENTRES WITHIN STRATA OF MISSED CLEAVAGES")
say("  if the effect disappears here, the variable was counting missed cleavages")
for k in sorted(D.internal_kr.unique())[:5]:
    s = D[D.internal_kr == k]
    if s.found.sum() < 4 or (~s.found).sum() < 20:
        continue
    d = delta(s.loc[s.found, "basic"], s.loc[~s.found, "basic"])
    say("    %d missed cleavage%-2s n=%4d, %2d identified   effect %+.3f"
        % (k, "s" if k != 1 else " ", len(s), int(s.found.sum()), d))

say("")
say("AND THE REVERSE: missed cleavages within strata of basic centres")
for k in sorted(D.basic.unique())[:6]:
    s = D[D.basic == k]
    if s.found.sum() < 4 or (~s.found).sum() < 20:
        continue
    d = delta(s.loc[s.found, "internal_kr"], s.loc[~s.found, "internal_kr"])
    say("    %d basic centres   n=%4d, %2d identified   effect %+.3f"
        % (k, len(s), int(s.found.sum()), d))

def auc(cols):
    X = D[cols].values; y = D.found.astype(int).values
    out = []
    for r in range(5):
        cv = StratifiedKFold(5, shuffle=True, random_state=r)
        pred = np.zeros(len(y), float)
        for tr, te in cv.split(X, y):
            m = make_pipeline(StandardScaler(),
                              LogisticRegression(max_iter=2000, class_weight="balanced"))
            m.fit(X[tr], y[tr]); pred[te] = m.predict_proba(X[te])[:, 1]
        out.append(roc_auc_score(y, pred))
    return float(np.mean(out))

say("")
say("DISCRIMINATION")
for cols, nm in ((["basic"], "basic centres alone"),
                 (["internal_kr"], "missed cleavages alone"),
                 (["basic", "internal_kr"], "both"),
                 (["mw"], "mass alone"),
                 (["internal_kr", "mw"], "missed cleavages + mass"),
                 (["basic", "mw"], "basic centres + mass")):
    say("  %-30s AUC %.3f" % (nm, auc(cols)))

open("missed_cleavage_report.txt", "w").write("\n".join(log) + "\n")
say("")
say("wrote missed_cleavage_report.txt")
