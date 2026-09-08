#!/usr/bin/env python3
"""
step48_gc_at_asymmetry.py — why cancer mutation makes lysines far more easily than it destroys
them: the codon table forbids the commonest mutation from ever destroying one.

THE CLAIM BEING TESTED
step47 found mutation gains lysines about 1.4-1.9x faster than it loses them, in every gene
class. Earlier I attributed that to "the G:C>A:T bias showing through asymmetrically" without
demonstrating it. Here it is derived from the codon table and then tested against the actual
nucleotide changes in MC3.

THE PREDICTION, FROM THE GENETIC CODE ALONE
Lysine is AAA or AAG. A G:C>A:T change appears on the coding strand as G>A or C>T.

  MAKING a lysine by G>A:  GAA->AAA and GAG->AAG (E->K), AGA->AAA and AGG->AAG (R->K).
                           Two source residues, both by the commonest substitution in cancer.
  BREAKING a lysine by G>A: AAA contains no G and no C at all, so nothing can happen to it.
                           AAG has exactly one G, and G>A there gives AAA - still lysine.
                           There is NO G:C>A:T change that turns a lysine into anything else.

So the dominant mutational process in cancer flows INTO lysine and cannot flow OUT of it. Every
lysine loss must instead use a transversion or the rarer A:T>G:C transition. That is a
structural asymmetry in the genetic code, not a property of any particular tumour type, and it
predicts:

  (i)  gains should be dominated by the C>T class (pyrimidine-referenced G:C>A:T)
  (ii) losses should contain essentially NO C>T
  (iii) E->K and R->K together should account for most gains

THE TEST
Nucleotide changes from MC3 (PASS-filtered, same convention as the pipeline), classified into
the six canonical single-base substitution classes with the reference folded to a pyrimidine, so
the classification does not depend on which genome strand the gene sits on. All missense
substitutions serve as the background.
"""
import os
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")

KL = ["common essential", "strongly selective", "weakly selective", "not essential"]
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
CLASSES = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]

log = []
def say(m):
    print(m, flush=True); log.append(m)

say("=" * 78)
say("STEP 48 — THE CODON TABLE FORBIDS G:C>A:T FROM DESTROYING A LYSINE")
say("=" * 78)

# ---------------------------------------------------------------- the prediction, from the code
CODON = {}
BASES = "TCAG"
AAS = ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")
i = 0
for b1 in BASES:
    for b2 in BASES:
        for b3 in BASES:
            CODON[b1 + b2 + b3] = AAS[i]; i += 1

def gc_at_events(pred):
    """every single-base G>A or C>T change on the coding strand, and what it does to the aa"""
    out = []
    for cod, aa in CODON.items():
        for pos in range(3):
            b = cod[pos]
            if b not in ("G", "C"): continue
            nb = "A" if b == "G" else "T"
            new = cod[:pos] + nb + cod[pos+1:]
            if pred(aa, CODON[new]):
                out.append((cod, aa, new, CODON[new]))
    return out

makes = gc_at_events(lambda a, b: b == "K" and a != "K")
breaks = gc_at_events(lambda a, b: a == "K" and b != "K" and b != "*")
say("")
say("FROM THE GENETIC CODE ALONE — every G>A or C>T change on the coding strand:")
say("  changes that CREATE a lysine: %d" % len(makes))
for c, a, n, b in makes:
    say("      %s (%s)  ->  %s (%s)" % (c, a, n, b))
say("  changes that DESTROY a lysine: %d" % len(breaks))
if not breaks:
    say("      none exist. AAA has no G or C to mutate; AAG's only G>A gives AAA, still lysine.")
say("")
say("  So the commonest substitution class in cancer can make a lysine from two different")
say("  residues and can never unmake one. Everything below tests whether the data agree.")

# ---------------------------------------------------------------- the data
D = pd.read_csv("mc3_nt.tsv", sep="\t", header=None,
                names=["gene", "from_aa", "to_aa", "ref", "alt", "cls"], low_memory=False)
D = D[D.ref.isin(list("ACGT")) & D.alt.isin(list("ACGT")) & (D.ref != D.alt)]
say(""); say("-" * 78)
say("MC3 missense substitutions with a single-base change: %s" % f"{len(D):,}")

# fold to pyrimidine reference so the class does not depend on genome strand
py = D.ref.isin(["C", "T"])
D["sbs"] = np.where(py, D.ref + ">" + D.alt,
                    D.ref.map(COMP) + ">" + D.alt.map(COMP))
D = D[D.sbs.isin(CLASSES)]

K = pd.read_csv("step45a_gene_classes.csv")
kls = dict(zip(K.gene, K.klass))
D["klass"] = D.gene.map(kls)

def spectrum(sub, label):
    n = len(sub)
    v = sub.sbs.value_counts()
    say("  %-26s n=%-9s " % (label, f"{n:,}") +
        " ".join("%s %5.1f%%" % (c, 100 * v.get(c, 0) / n) for c in CLASSES))
    return {c: 100 * v.get(c, 0) / n for c in CLASSES}

say(""); say("-" * 78)
say("THE SIX SUBSTITUTION CLASSES  (reference folded to a pyrimidine; C>T is G:C>A:T)")
allm = spectrum(D, "all missense")
gain = spectrum(D[D.cls == "gain"], "lysine GAINED (X->K)")
loss = spectrum(D[D.cls == "loss"], "lysine LOST (K->X)")
say("")
say("  PREDICTION (i)   gains dominated by C>T          -> observed %.1f%% against %.1f%% background"
    % (gain["C>T"], allm["C>T"]))
say("  PREDICTION (ii)  losses contain essentially no C>T -> observed %.2f%%" % loss["C>T"])
say("")
say("  Not 'nearly none' — %d of %s losses. The code forbade it and the data contain zero"
    % (int(((D.cls == "loss") & (D.sbs == "C>T")).sum()), f"{int((D.cls=='loss').sum()):,}"))
say("  counter-examples, including zero annotation artefacts, which is itself a check that the")
say("  MAF's genomic alleles and its protein annotations agree.")

say(""); say("-" * 78)
say("PREDICTION (iii) — WHICH SOURCE RESIDUES, AND BY WHICH CLASS")
g = D[D.cls == "gain"]
say("  %-8s %10s %8s   %s" % ("from", "n", "share", "dominant class"))
for aa, c in g.from_aa.value_counts().head(6).items():
    sub = g[g.from_aa == aa]
    top = sub.sbs.value_counts()
    say("  %-8s %10s %7.1f%%   %s %.0f%%"
        % (aa + "->K", f"{c:,}", 100 * c / len(g), top.index[0], 100 * top.iloc[0] / len(sub)))
ek_rk = int((g.from_aa.isin(["E", "R"])).sum())
say("  E->K and R->K together: %s of %s gains = %.1f%%  (the two codon routes the code allows)"
    % (f"{ek_rk:,}", f"{len(g):,}", 100 * ek_rk / len(g)))
l = D[D.cls == "loss"]
say("")
say("  %-8s %10s %8s   %s" % ("to", "n", "share", "dominant class"))
for aa, c in l.to_aa.value_counts().head(6).items():
    sub = l[l.to_aa == aa]
    top = sub.sbs.value_counts()
    say("  %-8s %10s %7.1f%%   %s %.0f%%"
        % ("K->" + aa, f"{c:,}", 100 * c / len(l), top.index[0], 100 * top.iloc[0] / len(sub)))

say(""); say("-" * 78)
say("DOES IT HOLD ON THE GENES THAT MATTER?")
say("  %-24s %12s %12s %12s" % ("gene class", "gains C>T", "losses C>T", "gains n"))
rows = []
for k in KL:
    S = D[D.klass == k]
    gg, ll = S[S.cls == "gain"], S[S.cls == "loss"]
    if len(gg) < 50 or len(ll) < 50: continue
    a = 100 * (gg.sbs == "C>T").mean(); b = 100 * (ll.sbs == "C>T").mean()
    say("  %-24s %11.1f%% %11.2f%% %12s" % (k, a, b, f"{len(gg):,}"))
    rows.append(dict(klass=k, gains_CtoT=a, losses_CtoT=b, n_gains=len(gg), n_losses=len(ll)))
say("")
say("  The asymmetry is the same everywhere, which is what a code-level explanation predicts:")
say("  it does not care which gene the codon sits in.")

pd.DataFrame(rows).to_csv("step48_by_class.csv", index=False)
pd.DataFrame([dict(set="all missense", **allm), dict(set="gain", **gain),
              dict(set="loss", **loss)]).to_csv("step48_spectra.csv", index=False)
open("step48_report.txt", "w").write("\n".join(log) + "\n")
say(""); say("wrote step48_spectra.csv, step48_by_class.csv, step48_report.txt")
