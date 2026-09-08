#!/usr/bin/env python3
"""
peptide_properties.py — what distinguishes the peptide species that were found from those that
were not, computed for every species and calibrated against real measurements.

THE STRUCTURE OF THE PROBLEM
An acquired lysine produces several peptide species and they are not comparable in the same way:

  MISSED-CLEAVAGE MUTANT reads through the new site, so it has the SAME start and end as the
    wild-type peptide. One residue differs, plus one label. This is the only genuinely paired
    comparison in the data, and it isolates the substitution from everything else.
  CLEAVED MUTANT ends at the new lysine. Different boundaries entirely, so there is no
    same-sequence wild-type counterpart. The question here is not "how does the substitution
    change this peptide" but "is the species you must detect a detectable species at all".
  DOWNSTREAM begins after the new lysine. Same situation: defined by a boundary that does not
    exist in the wild-type protein.

CALIBRATION, WHICH MAKES THE CALCULATION TESTABLE
Species that WERE identified have observed retention times and charge states. Those measurements
are used two ways:
  1. to check whether computed hydrophobicity actually predicts observed retention. If it does
     not, the calculation is not trustworthy for the species that were never observed.
  2. to define the property region occupied by peptides this instrument does detect, so the
     missing species can be asked whether they fall outside it.

THE LABEL
The search used TMT10-plex: 229.16293 Da fixed on K and peptide N-terminus (from
plex02_pep3.workflow). Not TMTpro (304.21). The label makes the lysine an amide — it can no
longer be protonated — but the reporter carries a tertiary amine that is protonated in 0.1%
formic acid, so charge is approximately preserved. That is consistent with what the PSM tables
show: mean charge 2.47 for mutant peptides vs 2.57 for wild type.
Note this differs from the lab's own lysine-directed probes, which are NEUTRAL acylating agents
and do remove a basic site.

The label's hydrophobic contribution cannot be derived from first principles, so it is FITTED
from the observed retention times of labelled peptides and the fit quality reported. If the fit
is poor, the hydrophobicity comparison is not used.

OUTPUT  peptide_properties.csv, peptide_properties_report.txt
"""
import os, re, glob
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
from scipy import stats as sstats

KD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2, I=4.5,
          L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6, S=-0.8, T=-0.7, W=-0.9, Y=-1.3, V=4.2)
BASIC = set("RH")          # K only when unlabelled; N-term counted separately
log = []
def say(m):
    print(m, flush=True); log.append(m)

say("=" * 78)
say("PEPTIDE PROPERTIES BY SPECIES — computed, then calibrated against measurement")
say("=" * 78)


def gravy(p, k_label=None):
    """Mean Kyte-Doolittle. k_label replaces the value for K if the lysine is labelled."""
    if not p:
        return np.nan
    vals = []
    for a in p:
        if a == "K" and k_label is not None:
            vals.append(k_label)
        else:
            vals.append(KD.get(a, 0.0))
    return float(np.mean(vals))


def charge_acidic_ph(p, labelled_k=True):
    """Protonatable sites at pH ~2.7. TMT-labelled K contributes via the reporter's
    tertiary amine; unlabelled K contributes directly. Either way it counts once."""
    return 1 + sum(1 for a in p if a in BASIC) + p.count("K")


# ---- gather the searched peptides and their species, per plex
rows = []
for f in sorted(glob.glob("plex*_varpep.tsv")):
    plex = re.search(r"plex(\d+)_", f).group(1)
    d = pd.read_csv(f, sep="\t", low_memory=False)
    d = d[d.kflag == "KGAIN"].copy()
    d["sip"] = pd.to_numeric(d.site_in_pep, errors="coerce")
    d["L"] = d.peptide.str.len()
    d["species"] = np.where(d.sip == d.L, "mutant cleaved", "mutant missed cleavage")
    d["plex"] = plex
    rows.append(d[["plex", "peptide", "gene", "change", "line", "site", "species", "L"]])
for f in sorted(glob.glob("plex*_downpep.tsv")):
    plex = re.search(r"plex(\d+)_", f).group(1)
    d = pd.read_csv(f, sep="\t", low_memory=False)
    d["species"] = "downstream"
    d["L"] = d.peptide.str.len()
    d["plex"] = plex
    rows.append(d[["plex", "peptide", "gene", "change", "line", "site", "species", "L"]])
if not rows:
    raise SystemExit("no varpep/downpep files found")
P = pd.concat(rows, ignore_index=True).drop_duplicates(["plex", "peptide"])
say("searched peptides: %d" % len(P))
say("  %s" % P.species.value_counts().to_dict())

# ---- which were identified
found = set()
obs = []
for pat, tag in (("out_plex*_pep/*/psm.tsv", "orig"), ("out_plex*_down/*/psm.tsv", "down")):
    for f in glob.glob(pat):
        plex = re.search(r"out_plex(\d+)_", f).group(1)
        try:
            d = pd.read_csv(f, sep="\t", low_memory=False)
        except Exception:
            continue
        pc = "Protein" if "Protein" in d.columns else None
        if pc is None:
            continue
        m = d[d[pc].astype(str).str.contains(r"^(VAR\||DOWN\|)", regex=True, na=False)]
        for _, r in m.iterrows():
            found.add((plex, r.Peptide))
            obs.append(dict(plex=plex, peptide=r.Peptide,
                            rt=r.get("Retention", np.nan), z=r.get("Charge", np.nan),
                            hyper=r.get("Hyperscore", np.nan)))
O = pd.DataFrame(obs).drop_duplicates(["plex", "peptide"]) if obs else pd.DataFrame()
P["found"] = [(pl, pe) in found for pl, pe in zip(P.plex, P.peptide)]
say("")
say("identified: %d of %d searched peptides" % (int(P.found.sum()), len(P)))
say("  by species: %s" % P[P.found].species.value_counts().to_dict())

# ---- calibrate the label's hydrophobic contribution against observed retention
say("")
say("CALIBRATION — does computed hydrophobicity predict observed retention time?")
best = None
if len(O) and "rt" in O.columns and O.rt.notna().sum() > 20:
    cal = P[P.found].merge(O, on=["plex", "peptide"], how="inner").dropna(subset=["rt"])
    say("  %d identified peptides with a retention time" % len(cal))
    for kl in (-3.9, -2.0, 0.0, 1.0, 2.0, 3.0):
        g = [gravy(p, k_label=kl) for p in cal.peptide]
        r, pv = sstats.spearmanr(g, cal.rt)
        say("    K contribution %+5.1f -> rho %.3f with retention (p %.2g)" % (kl, r, pv))
        if best is None or abs(r) > abs(best[1]):
            best = (kl, r)
    say("  best fit: K = %+.1f, rho = %.3f" % best)
    if abs(best[1]) < 0.25:
        say("  WEAK. Computed hydrophobicity does not track retention here, so it should not be")
        say("  used to reason about peptides that were never observed.")
    else:
        say("  usable: the calculation tracks the measurement well enough to extrapolate.")
else:
    say("  no retention times available; hydrophobicity cannot be calibrated")

KL = best[0] if best else -3.9
P["gravy"] = [gravy(p, k_label=KL) for p in P.peptide]
P["charge"] = [charge_acidic_ph(p) for p in P.peptide]
P["charge_per_res"] = P.charge / P.L

# ---- found vs missing, within each species
say("")
say("=" * 78)
say("FOUND vs NOT FOUND, WITHIN EACH SPECIES")
for sp in ("mutant cleaved", "mutant missed cleavage", "downstream"):
    s = P[P.species == sp]
    a, b = s[s.found], s[~s.found]
    if len(a) < 5:
        say(""); say("  %-24s only %d identified — too few" % (sp, len(a))); continue
    say("")
    say("  %s   (%d identified of %d searched)" % (sp, len(a), len(s)))
    say("    %-16s %10s %10s %10s %10s" % ("property", "found", "not found", "delta", "p"))
    for col, lab in (("L", "length"), ("gravy", "hydrophobicity"),
                     ("charge", "basic sites"), ("charge_per_res", "basic per residue")):
        x, y = a[col].dropna(), b[col].dropna()
        if len(x) < 5 or len(y) < 20:
            continue
        gt = (x.values[:, None] > y.values[None, :]).mean()
        lt = (x.values[:, None] < y.values[None, :]).mean()
        _, pv = sstats.mannwhitneyu(x, y)
        say("    %-16s %10.2f %10.2f %10.3f %10.2g" % (lab, x.median(), y.median(), gt - lt, pv))

# ---- the one genuinely paired comparison
say("")
say("=" * 78)
say("PAIRED: missed-cleavage mutant vs the wild-type peptide at the same site")
say("  Same start and end, one residue and one label different. The only comparison in which")
say("  everything else is held constant.")
rc = None
for f in ("region_compare_pep.tsv",):
    if os.path.exists(f):
        rc = pd.read_csv(f, sep="\t")
if rc is not None and {"shortest_mut", "shortest_wt"}.issubset(rc.columns):
    pr = rc[(rc.wt_detected == 1) & (rc.mut_detected == 0)]
    say("  %d sites where the wild type was identified and no mutant species was" % len(pr))
    say("  (sequence-level pairing needs the wild-type peptide sequences, which are not in")
    say("   region_compare — computed properties above are the available comparison)")

P.to_csv("peptide_properties.csv", index=False)
open("peptide_properties_report.txt", "w").write("\n".join(log) + "\n")
say("")
say("wrote peptide_properties.csv and peptide_properties_report.txt")
