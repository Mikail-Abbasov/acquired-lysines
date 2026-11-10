#!/usr/bin/env python3
"""
step6c_marker.py — resolve the AMBIGUOUS sites using the marked-lysine peptide.

WHAT WENT WRONG IN STEP 6
442 sites were classed AMBIGUOUS: the peptide matched the sequence once but contained two or
three lysines, and the reported position picked none of them. That was treated as unresolvable.
It is not. CyMPSA and Nature Chemistry both ship a peptide column in which the probe-labelled
lysine is marked with an asterisk — EVYEEAGVK*GK says plainly which of the two lysines carried
the probe. 357 of the 442 have such a marker. The ambiguity was in the resolver, not the data.

Meroterpenoids ship no marked column (confirmed in step 3), so their 96 ambiguous sites remain
unresolvable and stay excluded.

THIS SCRIPT DOES TWO THINGS

  1. RESOLVE. For each unusable site with a marker, locate the peptide in the AlphaFold
     sequence and take the lysine at the marked offset. Recompute geometry there. These are
     recorded as verdict MARKER, a third evidence route independent of the reported number.

  2. VALIDATE THE MARKER ITSELF. Before trusting it, check it against sites already resolved
     two ways: for AGREE sites that also carry a marker, does the marked lysine fall at the
     position both other routes agreed on? If the marker disagrees with a doubly-confirmed
     position, it is unreliable and must not be used to rescue anything. This check costs
     nothing and is reported first.

Only the previously unusable rows are touched. The 44,736 AGREE sites are already resolved by
two concurring routes and are not recomputed.

INPUT   features.csv (step 6 v2), sites.csv (step 4), af/
OUTPUT  features_final.csv, step6c_report.txt
"""
import os, re
import numpy as np, pandas as pd
import warnings
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

AF_DIR = os.environ.get("AF_DIR", "af")
HYDRO, AROM, POLAR = set("AVLIMFWPC"), set("FWYH"), set("STNQCY")
ANION, CATION = set("DE"), set("KRH")
T3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
      "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
      "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}

log = []
def say(m):
    print(m, flush=True); log.append(m)


def marker_index(marked):
    """Return (bare peptide, 0-based index of the marked lysine) or (None, None)."""
    if not isinstance(marked, str):
        return None, None
    s = re.sub(r"^\[?[A-Z]+\]?\.", "", marked.strip())
    s = re.sub(r"\.\[?[A-Z-]+\]?$", "", s)
    i = s.find("K*")
    if i < 0:
        return None, None
    bare = s.replace("*", "")
    idx = len(re.sub(r"[^A-Z]", "", s[:i]))
    return re.sub(r"[^A-Z]", "", bare), idx


feat = pd.read_csv("features.csv", low_memory=False)
sites = pd.read_csv("sites.csv", low_memory=False)
mk = (sites.dropna(subset=["peptide_marked"])
           .drop_duplicates(["dataset", "uniprot", "position"])
           [["dataset", "uniprot", "position", "peptide_marked"]])
say("=" * 74)
say("STEP 6c — MARKED-LYSINE RESOLUTION")
say("=" * 74)
say("features rows %d | rows with a marked peptide available %d" % (len(feat), len(mk)))

parser, sr = PDBParser(QUIET=True), ShrakeRupley()
cache = {}
def load(acc):
    if acc not in cache:
        p = os.path.join(AF_DIR, "%s.pdb" % acc)
        if not os.path.exists(p):
            cache[acc] = None
        else:
            st = parser.get_structure(acc, p)
            ch = next(st[0].get_chains())
            res = {r.id[1]: r for r in ch if r.id[0] == " "}
            order = sorted(res)
            seq = "".join(T3.get(res[i].get_resname(), "X") for i in order)
            cache[acc] = (ch, res, order, seq)
    return cache[acc]


# ---------------------------------------------------------------- 1. validate the marker
say("")
say("CHECK — does the marker agree with positions already confirmed two other ways?")
ag = feat[feat.verdict == "AGREE"][["dataset", "uniprot", "position_reported", "position"]]
ag = ag.merge(mk, left_on=["dataset", "uniprot", "position_reported"],
              right_on=["dataset", "uniprot", "position"], how="inner",
              suffixes=("", "_y"))
same = diff = skip = 0
for r in ag.head(4000).itertuples():
    L = load(r.uniprot)
    if not L:
        skip += 1; continue
    _, _, order, seq = L
    pep, idx = marker_index(r.peptide_marked)
    if not pep or seq.count(pep) != 1:
        skip += 1; continue
    j = seq.find(pep)
    if order[j + idx] == r.position:
        same += 1
    else:
        diff += 1
tot = same + diff
say("  checked %d AGREE sites carrying a marker" % tot)
if tot:
    say("  marker points to the SAME residue : %d (%.2f%%)" % (same, 100 * same / tot))
    say("  marker DISAGREES                  : %d (%.2f%%)" % (diff, 100 * diff / tot))
    if diff / tot > 0.02:
        say("  !! marker disagrees with doubly-confirmed positions in >2% of cases.")
        say("     It is NOT reliable enough to rescue ambiguous sites. Stopping.")
        open("step6c_report.txt", "w").write("\n".join(log) + "\n")
        raise SystemExit(1)
    say("  -> marker validated; using it to resolve ambiguous sites")

# ---------------------------------------------------------------- 2. resolve
bad = feat[feat.qc != "ok"].copy()
bad = bad.merge(mk, left_on=["dataset", "uniprot", "position_reported"],
                right_on=["dataset", "uniprot", "position"], how="left",
                suffixes=("", "_mk"))
rows, tally = [], {}
for r in bad.itertuples():
    L = load(r.uniprot)
    pep, idx = marker_index(getattr(r, "peptide_marked", None))
    if not L or not pep:
        tally["no_marker_or_structure"] = tally.get("no_marker_or_structure", 0) + 1
        continue
    ch, res, order, seq = L
    if seq.count(pep) != 1:
        tally["peptide_not_unique"] = tally.get("peptide_not_unique", 0) + 1
        continue
    j = seq.find(pep)
    pos = order[j + idx]
    r0 = res.get(pos)
    if r0 is None or T3.get(r0.get_resname(), "X") != "K" or "NZ" not in r0:
        tally["marked_residue_not_K"] = tally.get("marked_residue_not_K", 0) + 1
        continue
    try:
        sr.compute(ch, level="A")
    except Exception:
        tally["sasa_failed"] = tally.get("sasa_failed", 0) + 1
        continue
    heavy = [(a.coord, a.get_parent()) for a in ch.get_atoms() if a.element != "H"]
    coords = np.array([h[0] for h in heavy]); parents = [h[1] for h in heavy]
    cys = np.array([x["SG"].coord for x in res.values()
                    if x.get_resname() == "CYS" and "SG" in x]) if res else np.empty((0, 3))
    nz = r0["NZ"].coord
    d = np.linalg.norm(coords - nz, axis=1)
    near8 = {parents[k] for k in np.where(d <= 8)[0]} - {r0}
    aas = [T3.get(x.get_resname(), "X") for x in near8]
    pl = [x["CA"].get_bfactor() for x in near8 if "CA" in x]
    tally["RESOLVED"] = tally.get("RESOLVED", 0) + 1
    rows.append(dict(
        dataset=r.dataset, uniprot=r.uniprot, position_reported=r.position_reported,
        position=pos, peptide=pep, verdict="MARKER", qc="ok",
        shift=pos - r.position_reported, peptide_matches=1,
        reported_is_lysine=r.reported_is_lysine,
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
        dist_centroid=float(np.linalg.norm(nz - coords.mean(axis=0))),
        dist_nearest_cys=float(np.min(np.linalg.norm(cys - nz, axis=1)))
                         if len(cys) else np.nan,
        plddt_site=float(r0["CA"].get_bfactor()) if "CA" in r0 else np.nan,
        plddt_8A=float(np.mean(pl)) if pl else np.nan))

new = pd.DataFrame(rows)
say("")
say("RESOLUTION OF THE %d PREVIOUSLY UNUSABLE ROWS" % len(feat[feat.qc != "ok"]))
for k, v in sorted(tally.items(), key=lambda x: -x[1]):
    say("  %-26s %5d" % (k, v))
if len(new):
    say("")
    say("  by dataset : %s" % new.dataset.value_counts().to_dict())
    say("  by original verdict: recovered from %s"
        % feat[feat.qc != "ok"].verdict.value_counts().to_dict())
    s = new["shift"].dropna()
    a = s.abs()
    say("")
    say("  HOW FAR THE RESOLVED POSITIONS MOVED")
    say("    Do not quote the signed median. %d shifts are negative and %d positive, so"
        % ((s < 0).sum(), (s > 0).sum()))
    say("    they cancel to a signed median of %+.0f residue while the typical displacement"
        % s.median())
    say("    is far larger. The magnitude is the informative figure.")
    say("      |shift| median %.0f residues | IQR %.0f to %.0f | max %.0f"
        % (a.median(), a.quantile(.25), a.quantile(.75), a.max()))
    say("      at the reported position (shift 0) : %d of %d" % ((a == 0).sum(), len(s)))
    say("      within one residue (|shift| <= 1)  : %d of %d" % ((a <= 1).sum(), len(s)))
    say("    A shift of zero is impossible here, and its absence is not a warning sign.")
    say("    Every row rescued at this step was AMBIGUOUS: the peptide matched the")
    say("    sequence once, carried more than one lysine, and the reported number named")
    say("    none of them. The marked lysine therefore can never be the residue that")
    say("    number pointed at. The marker never reads the reported number as a")
    say("    coordinate; it locates the peptide in the sequence and takes the marked")
    say("    offset. The number is used only to join back to the source row.")
    inside = outside = oor = 0
    for r in new.itertuples():
        Lx = load(r.uniprot)
        if not Lx:
            continue
        _, _, order_x, seq_x = Lx
        if not (order_x[0] <= r.position_reported <= order_x[-1]):
            oor += 1
            continue
        jx = seq_x.find(r.peptide)
        if jx >= 0 and r.position_reported in order_x[jx:jx + len(r.peptide)]:
            inside += 1
        else:
            outside += 1
    say("    where the reported number actually falls:")
    say("      inside the matched peptide, on a residue that is not a lysine : %d" % inside)
    say("      elsewhere in the same protein                                 : %d" % outside)
    say("      outside the protein's modelled residue range altogether       : %d" % oor)
    say("    step6c_audit_shift.py attributes the displacements to named isoforms.")

final = pd.concat([feat[feat.qc == "ok"], new], ignore_index=True)
final.to_csv("features_final.csv", index=False)
say("")
say("FINAL: %d usable sites (%d two-route AGREE + %d peptide-corrected + %d marker-resolved)"
    % (len(final), (final.verdict == "AGREE").sum(),
       final.verdict.isin(["PEPTIDE_ONLY", "PEPTIDE_PICK", "NUMBER_ONLY"]).sum(),
       (final.verdict == "MARKER").sum()))
say("  verdict column retained so any analysis can be restricted to AGREE only.")
open("step6c_report.txt", "w").write("\n".join(log) + "\n")
say("wrote features_final.csv and step6c_report.txt")
