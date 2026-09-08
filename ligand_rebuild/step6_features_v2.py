#!/usr/bin/env python3
"""
step6_features_v2.py — locate every lysine site by TWO independent routes, then compute
geometry only where the site is established.

WHY THIS REPLACES THE OFFSET APPROACH
The first version inferred one numbering offset per (dataset, protein) and forced every site
through it. That is wrong whenever conventions are mixed inside a protein, and beta-actin
proved they are: P60709 positions 68, 113 and 215 are genuine lysines in the AlphaFold
sequence, but a -1 offset inferred from that protein's OTHER sites moved them to 67, 112 and
214 and marked all three failed. The offset was damaging sites that were already correct.

Inferring a convention is the mistake. Each site carries its own evidence, so each site is
resolved on its own evidence.

TWO INDEPENDENT ROUTES TO THE SAME COORDINATE
  route 1  the PEPTIDE. The identified sequence is the observation. Locating it in the
           AlphaFold sequence fixes a stretch of the protein without reference to any
           numbering convention.
  route 2  the REPORTED LYSINE NUMBER, which every dataset provides.

Neither is assumed correct. They are compared, and the comparison is the result.

VERDICTS
  AGREE            peptide matches once, contains a lysine, and that lysine is at the
                   reported position. Confirmed twice over — the strongest class.
  PEPTIDE_ONLY     peptide matches once and contains exactly one lysine, but not at the
                   reported position. The peptide is the observation and the number is
                   bookkeeping, so the peptide wins; the discrepancy is recorded.
  PEPTIDE_PICK     peptide matches once but contains several lysines, and the reported
                   number selects one of them. Both routes used jointly. Flagged.
  NUMBER_ONLY      peptide cannot be located (isoform difference, wrong accession, sequence
                   mismatch) but the reported position IS a lysine. Single-route evidence,
                   kept and flagged so it can be excluded in sensitivity analysis.
  MULTI_MATCH      peptide occurs more than once (internal repeats — the C2H2 zinc-finger
                   TGEKP linker is the common case). Position undetermined. Excluded.
  AMBIGUOUS        peptide matches once, several lysines, reported number selects none of
                   them. Excluded.
  NO_EVIDENCE      no peptide match and reported position is not a lysine. Excluded.

A third, free check is recorded but not used to exclude: whether the residue preceding the
peptide is K or R, i.e. whether the match is consistent with tryptic cleavage.

OUTPUT
  features.csv        one row per (dataset, uniprot, reported position) with verdict,
                      resolved position, and geometry where resolved
  step6_report.txt

Geometry is measured from the lysine NZ atom, the nitrogen that does the chemistry.
"""
import os, re, sys, time, urllib.request
import numpy as np, pandas as pd
import warnings
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

AF_DIR = os.environ.get("AF_DIR", "af")
AF_VERSIONS = (6, 5, 4)
AF_URL = "https://alphafold.ebi.ac.uk/files/AF-%s-F1-model_v%d.pdb"
HYDRO, AROM, POLAR = set("AVLIMFWPC"), set("FWYH"), set("STNQCY")
ANION, CATION = set("DE"), set("KRH")
T3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
      "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
      "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
os.makedirs(AF_DIR, exist_ok=True)

log = []
def say(m):
    print(m, flush=True); log.append(m)


def fetch(acc):
    p = os.path.join(AF_DIR, "%s.pdb" % acc)
    if os.path.exists(p) and os.path.getsize(p) > 1000:
        return p
    for v in AF_VERSIONS:
        try:
            urllib.request.urlretrieve(AF_URL % (acc, v), p)
            if os.path.getsize(p) > 1000:
                return p
        except Exception:
            pass
        if os.path.exists(p):
            os.remove(p)
    return None


def bare(p):
    if not isinstance(p, str):
        return None
    s = re.sub(r"^\[?[A-Z]+\]?\.", "", p.strip())
    s = re.sub(r"\.\[?[A-Z-]+\]?$", "", s)
    return re.sub(r"[^A-Z]", "", s) or None


sites = pd.read_csv("sites.csv", low_memory=False)
sites["pep"] = sites.peptide.map(bare)
sites["pep_marked"] = sites.get("peptide_marked", pd.Series(index=sites.index, dtype=object))
want = (sites.sort_values("pep", na_position="last")
             .drop_duplicates(["dataset", "uniprot", "position"])
             [["dataset", "uniprot", "position", "pep"]])
by_prot = {a: g for a, g in want.groupby("uniprot")}
accs = sorted(by_prot)
say("=" * 74)
say("STEP 6 (v2) — TWO-WAY SITE VALIDATION: peptide match AND reported lysine number")
say("=" * 74)
say("%d proteins, %d (dataset, protein, position) rows" % (len(accs), len(want)))

parser, sr = PDBParser(QUIET=True), ShrakeRupley()
rows, no_struct, t0 = [], [], time.time()

for i, acc in enumerate(accs, 1):
    if i % 250 == 0:
        say("  %d/%d proteins | %d rows | %.0f min" % (i, len(accs), len(rows),
                                                       (time.time() - t0) / 60))
    path = fetch(acc)
    if not path:
        no_struct.append(acc)
        continue
    try:
        st = parser.get_structure(acc, path)
        chain = next(st[0].get_chains())
        res = {r.id[1]: r for r in chain if r.id[0] == " "}
        order = sorted(res)
        seq = "".join(T3.get(res[i2].get_resname(), "X") for i2 in order)
        sr.compute(chain, level="A")
    except Exception:
        no_struct.append(acc)
        continue

    heavy = [(a.coord, a.get_parent()) for a in chain.get_atoms() if a.element != "H"]
    coords = np.array([h[0] for h in heavy]); parents = [h[1] for h in heavy]
    centroid = coords.mean(axis=0)
    cys = np.array([x["SG"].coord for x in res.values()
                    if x.get_resname() == "CYS" and "SG" in x]) if res else np.empty((0, 3))

    for t in by_prot[acc].itertuples():
        pos, pep = int(t.position), t.pep
        rec = {"dataset": t.dataset, "uniprot": acc, "position_reported": pos,
               "peptide": pep}
        num_is_K = (T3.get(res[pos].get_resname(), "X") == "K") if pos in res else False
        rec["reported_is_lysine"] = num_is_K

        n = seq.count(pep) if pep else 0
        rec["peptide_matches"] = n
        resolved, verdict = None, None
        if n == 1:
            j = seq.find(pep)
            rec["tryptic_nterm"] = (seq[j - 1] in "KR") if j > 0 else True
            ks = [order[j + k] for k, a in enumerate(pep) if a == "K"]
            if not ks:
                verdict = "NO_EVIDENCE" if not num_is_K else "NUMBER_ONLY"
                resolved = pos if num_is_K else None
            elif pos in ks:
                verdict, resolved = "AGREE", pos
            elif len(ks) == 1:
                verdict, resolved = "PEPTIDE_ONLY", ks[0]
            else:
                near = [k for k in ks if abs(k - pos) <= 3]
                if len(near) == 1:
                    verdict, resolved = "PEPTIDE_PICK", near[0]
                else:
                    verdict, resolved = "AMBIGUOUS", None
        elif n > 1:
            verdict, resolved = "MULTI_MATCH", None
        else:
            verdict = "NUMBER_ONLY" if num_is_K else "NO_EVIDENCE"
            resolved = pos if num_is_K else None
        rec["verdict"] = verdict
        rec["position"] = resolved
        rec["shift"] = (resolved - pos) if resolved is not None else np.nan

        r0 = res.get(resolved) if resolved is not None else None
        if r0 is None or "NZ" not in r0:
            rec["qc"] = "unresolved" if r0 is None else "no_NZ"
            rows.append(rec); continue

        nz = r0["NZ"].coord
        d = np.linalg.norm(coords - nz, axis=1)
        near8 = {parents[k] for k in np.where(d <= 8)[0]} - {r0}
        aas = [T3.get(x.get_resname(), "X") for x in near8]
        pl = [x["CA"].get_bfactor() for x in near8 if "CA" in x]
        rec.update(
            qc="ok",
            crowding_8A=len(near8), atoms_10A=int((d <= 10).sum()),
            hydrophobic_8A=sum(a in HYDRO for a in aas),
            aromatic_8A=sum(a in AROM for a in aas),
            polar_8A=sum(a in POLAR for a in aas),
            anionic_8A=sum(a in ANION for a in aas),
            cationic_8A=sum(a in CATION for a in aas),
            net_charge_8A=sum(a in CATION for a in aas) - sum(a in ANION for a in aas),
            n_lys_8A=sum(a == "K" for a in aas),
            sasa_total=float(sum(a.sasa for a in r0)),
            sasa_sidechain=float(sum(a.sasa for a in r0
                                     if a.get_name() not in ("N", "CA", "C", "O"))),
            nz_sasa=float(r0["NZ"].sasa),
            rsa=float(sum(a.sasa for a in r0)) / 200.0,
            dist_centroid=float(np.linalg.norm(nz - centroid)),
            dist_nearest_cys=float(np.min(np.linalg.norm(cys - nz, axis=1)))
                             if len(cys) else np.nan,
            plddt_site=float(r0["CA"].get_bfactor()) if "CA" in r0 else np.nan,
            plddt_8A=float(np.mean(pl)) if pl else np.nan)
        rows.append(rec)

feat = pd.DataFrame(rows)
feat.to_csv("features.csv", index=False)

say("")
say("TWO-WAY VALIDATION")
n = len(feat)
order_v = ["AGREE", "PEPTIDE_ONLY", "PEPTIDE_PICK", "NUMBER_ONLY",
           "AMBIGUOUS", "MULTI_MATCH", "NO_EVIDENCE"]
for v in order_v:
    c = int((feat.verdict == v).sum())
    if c:
        say("  %-14s %6d  %5.2f%%" % (v, c, 100 * c / n))
usable = feat[feat.qc == "ok"]
say("  %-14s %6d  %5.2f%%" % ("USABLE", len(usable), 100 * len(usable) / n))

say("")
say("  agreement between the two routes, where both are available:")
both = feat[feat.verdict.isin(["AGREE", "PEPTIDE_ONLY", "PEPTIDE_PICK"])]
if len(both):
    say("    both routes point to the same residue : %d (%.1f%%)"
        % ((both.verdict == "AGREE").sum(),
           100 * (both.verdict == "AGREE").sum() / len(both)))
    say("    peptide overrides the reported number : %d (%.1f%%)"
        % ((both.verdict != "AGREE").sum(),
           100 * (both.verdict != "AGREE").sum() / len(both)))
    s = both.loc[both.verdict != "AGREE", "shift"].dropna()
    if len(s):
        say("    shift where they disagree: median %+.0f, |shift|<=1 in %d of %d"
            % (s.median(), (s.abs() <= 1).sum(), len(s)))
say("")
for ds, g in feat.groupby("dataset"):
    say("  %-14s AGREE %5.1f%% | peptide-corrected %4.1f%% | unusable %4.1f%%"
        % (ds, 100 * (g.verdict == "AGREE").mean(),
           100 * g.verdict.isin(["PEPTIDE_ONLY", "PEPTIDE_PICK"]).mean(),
           100 * (g.qc != "ok").mean()))
if "tryptic_nterm" in feat:
    tn = feat[feat.verdict == "AGREE"].tryptic_nterm.dropna()
    if len(tn):
        say("")
        say("  third check (not used to exclude): peptide preceded by K/R in %.1f%% of"
            " AGREE sites" % (100 * tn.mean()))
say("")
say("STRUCTURES UNAVAILABLE: %d of %d proteins (%.1f%%)"
    % (len(no_struct), len(accs), 100 * len(no_struct) / max(len(accs), 1)))
open("step6_report.txt", "w").write("\n".join(log) + "\n")
say("wrote features.csv (%d rows)" % len(feat))
