#!/usr/bin/env python3
"""
step6d_expand.py — expand the per-lysine parameter set, with reactivity and ligandability
kept as separate, labelled blocks.

THE POINT OF THE SPLIT
The lab's central claim is that lysine selectivity comes from BINDING, not from intrinsic
reactivity: a lysine is unreactive until a bound molecule holds an electrophile against it.
That claim is testable here, but only if the two properties are measured separately and never
pooled.

  BLOCK R — REACTIVITY. How nucleophilic is this lysine in itself?
      propka_pKa            predicted pKa of the NZ amine (PROPKA 3)
      pKa_shift             propka_pKa - 10.5, the model lysine value. A negative shift means
                            a more deprotonated, more nucleophilic amine.
      net_charge_8A, anionic_8A, cationic_8A   local electrostatics
      hbond_partners_NZ     donors/acceptors within 3.5 A of NZ that pin the amine

  BLOCK L — LIGANDABILITY. Can a molecule bind here?
      in_pocket             does fpocket place this lysine in a detected cavity
      pocket_volume         volume of that cavity (A^3)
      pocket_score          fpocket's own pocket score
      pocket_drugscore      fpocket druggability estimate
      dist_nearest_pocket   NZ to the nearest alpha-sphere centre, for lysines outside a pocket
      burial / crowding / aromatic contacts / secondary structure

  BLOCK S — STRUCTURAL CONTEXT (covariates, not hypotheses)
      dssp_ss               secondary structure (H / E / C)
      rsa_tien              RSA normalised by Tien et al. empirical maxima rather than a flat
                            200 A^2 as in step 6
      sasa_freesasa         independent SASA implementation, as a cross-check on Shrake-Rupley
      plddt_site, plddt_8A  AlphaFold confidence

THE TEST THIS ENABLES
If block L predicts engagement rate and block R does not, engagement is driven by binding, not
by how nucleophilic the lysine happens to be. If block R dominates, the opposite. The blocks
are modelled separately and corrected for multiple testing WITHIN block, so neither borrows
significance from the other.

CAVEAT THAT MUST TRAVEL WITH ANY pKa RESULT
These three datasets measure LIGANDABILITY (competition with a probe), not reactivity (probe
dose-response). A pKa prediction therefore cannot be validated against them. pKa is reported
as a descriptor and as a test of whether reactivity explains engagement — it is never treated
as a validated reactivity measurement.

RESTARTABLE. Per-protein results are cached as JSON in cache6d/, so an interrupted run resumes
without recomputing.

INPUT   features_final.csv, af/*.pdb
OUTPUT  features_expanded.csv, step6d_report.txt
"""
import os, sys, json, glob, time, shutil, subprocess, tempfile, re
import numpy as np, pandas as pd
import warnings
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

AF_DIR = os.environ.get("AF_DIR", "af")
CACHE = os.environ.get("CACHE_DIR", "cache6d")
os.makedirs(CACHE, exist_ok=True)
T3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
      "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
      "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
# Tien et al. 2013 empirical maximum accessibility, A^2
TIEN = {"A": 121, "R": 265, "N": 187, "D": 187, "C": 148, "Q": 214, "E": 214, "G": 97,
        "H": 216, "I": 195, "L": 191, "K": 230, "M": 203, "F": 228, "P": 154, "S": 143,
        "T": 163, "W": 264, "Y": 255, "V": 165}
HAVE = {}
for mod in ("freesasa", "pydssp", "propka"):
    try:
        __import__(mod); HAVE[mod] = True
    except Exception:
        HAVE[mod] = False
HAVE["fpocket"] = shutil.which("fpocket") is not None

log = []
def say(m):
    print(m, flush=True); log.append(m)


def run_fpocket(pdb, workdir):
    """Return list of pockets, each with volume/score/drugscore and its alpha-sphere centres."""
    dst = os.path.join(workdir, os.path.basename(pdb))
    shutil.copy(pdb, dst)
    try:
        subprocess.run(["fpocket", "-f", dst], capture_output=True, timeout=300)
    except Exception:
        return []
    out = dst[:-4] + "_out"
    info = os.path.join(out, os.path.basename(dst)[:-4] + "_info.txt")
    pockets, cur = [], None
    if os.path.exists(info):
        for line in open(info):
            m = re.match(r"Pocket\s+(\d+)", line)
            if m:
                cur = {"id": int(m.group(1))}
                pockets.append(cur)
            elif cur is not None and ":" in line:
                k, v = [x.strip() for x in line.split(":", 1)]
                try:
                    cur[k] = float(v)
                except ValueError:
                    pass
    # alpha-sphere centres per pocket, from the *_pockets.pqr / pocketN_atm.pdb files
    centres = {}
    for f in glob.glob(os.path.join(out, "pockets", "pocket*_vert.pqr")):
        pid = int(re.search(r"pocket(\d+)_", os.path.basename(f)).group(1)) + 1
        pts = []
        for line in open(f):
            if line.startswith(("ATOM", "HETATM")):
                try:
                    pts.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                except ValueError:
                    pass
        if pts:
            centres[pid] = np.array(pts)
    for p in pockets:
        p["centres"] = centres.get(p["id"])
    return pockets


def propka_pkas(pdb, workdir):
    """Return {resnum: pKa} for lysines, or {} on failure."""
    dst = os.path.join(workdir, os.path.basename(pdb))
    if not os.path.exists(dst):
        shutil.copy(pdb, dst)
    try:
        subprocess.run([sys.executable, "-m", "propka", dst],
                       capture_output=True, timeout=1800, cwd=workdir)
    except Exception:
        return {}
    pka_file = dst[:-4] + ".pka"
    out = {}
    if not os.path.exists(pka_file):
        return out
    # PROPKA writes fixed-width columns: "LYS  27 A  10.45". The residue-number field is
    # three characters wide, so at position >=1000 it overflows and abuts the residue name
    # ("LYS1006 A"). Splitting on whitespace then fails to match "LYS" and the site is
    # silently lost. This cost 2,214 of 2,378 missing pKa values, ALL of them at position
    # >=1000 - a systematic loss of high-position lysines, i.e. biased against large
    # proteins. Parse by column instead of by token.
    for line in open(pka_file):
        if not line.startswith("LYS"):
            continue
        try:
            num = int(line[3:8].strip())          # residue number, overflow-tolerant
            pka = float(line[10:16].strip())      # predicted pKa
        except (ValueError, IndexError):
            continue
        out.setdefault(num, pka)                  # detailed section first, SUMMARY repeats
    return out


def protein_expand(acc, positions):
    """Compute the expanded parameters for one protein. Cached."""
    cf = os.path.join(CACHE, "%s.json" % acc)
    if os.path.exists(cf):
        try:
            return json.load(open(cf))
        except Exception:
            pass
    pdb = os.path.join(AF_DIR, "%s.pdb" % acc)
    if not os.path.exists(pdb):
        return {}
    parser, sr = PDBParser(QUIET=True), ShrakeRupley()
    st = parser.get_structure(acc, pdb)
    chain = next(st[0].get_chains())
    res = {r.id[1]: r for r in chain if r.id[0] == " "}
    sr.compute(chain, level="A")

    work = tempfile.mkdtemp()
    pockets = run_fpocket(pdb, work) if HAVE["fpocket"] else []
    pkas = propka_pkas(pdb, work) if HAVE["propka"] else {}

    fs_area = {}
    if HAVE["freesasa"]:
        try:
            import freesasa
            s = freesasa.Structure(pdb)
            r = freesasa.calc(s)
            for i in range(s.nAtoms()):
                fs_area.setdefault(int(s.residueNumber(i)), 0.0)
                fs_area[int(s.residueNumber(i))] += r.atomArea(i)
        except Exception:
            pass

    ss_map = {}
    if HAVE["pydssp"]:
        try:
            import pydssp, torch
            order = sorted(res)
            need = ["N", "CA", "C", "O"]
            coords, keep = [], []
            for i in order:
                r0 = res[i]
                if all(a in r0 for a in need):
                    coords.append([r0[a].coord for a in need]); keep.append(i)
            if coords:
                arr = torch.tensor(np.array(coords, dtype=np.float32))[None]
                ss = pydssp.assign(arr, out_type="c3")[0]
                ss_map = {k: str(v) for k, v in zip(keep, ss)}
        except Exception:
            pass
    shutil.rmtree(work, ignore_errors=True)

    out = {}
    for pos in positions:
        r0 = res.get(int(pos))
        if r0 is None or "NZ" not in r0:
            continue
        nz = np.array(r0["NZ"].coord)
        rec = {}
        # --- block L: pockets
        best = None
        dmin = np.inf
        for p in pockets:
            c = p.get("centres")
            if c is None:
                continue
            dd = np.min(np.linalg.norm(c - nz, axis=1))
            if dd < dmin:
                dmin, best = dd, p
        rec["dist_nearest_pocket"] = float(dmin) if np.isfinite(dmin) else np.nan
        rec["in_pocket"] = bool(np.isfinite(dmin) and dmin <= 5.0)
        if best is not None:
            rec["pocket_volume"] = best.get("Volume", np.nan)
            rec["pocket_score"] = best.get("Score", np.nan)
            rec["pocket_drugscore"] = best.get("Druggability Score", np.nan)
        # --- block R: reactivity
        rec["propka_pKa"] = pkas.get(int(pos), np.nan)
        rec["pKa_shift"] = (rec["propka_pKa"] - 10.5) if rec["propka_pKa"] == rec["propka_pKa"] else np.nan
        nzc = nz
        hb = 0
        for r2 in res.values():
            for a in r2:
                if a.element in ("N", "O") and a.get_parent() is not r0:
                    if np.linalg.norm(a.coord - nzc) <= 3.5:
                        hb += 1
        rec["hbond_partners_NZ"] = hb
        # --- block S
        rec["dssp_ss"] = ss_map.get(int(pos), "")
        aa = T3.get(r0.get_resname(), "X")
        rec["rsa_tien"] = float(sum(a.sasa for a in r0)) / TIEN.get(aa, 200)
        rec["sasa_freesasa"] = fs_area.get(int(pos), np.nan)
        out[str(pos)] = rec
    json.dump(out, open(cf, "w"))
    return out


feat = pd.read_csv("features_final.csv", low_memory=False)
by_prot = feat.groupby("uniprot").position.apply(lambda s: sorted(set(int(x) for x in s))).to_dict()
accs = sorted(by_prot)
say("=" * 74)
say("STEP 6d — EXPANDED PARAMETERS (reactivity / ligandability / context)")
say("tools: %s" % ", ".join("%s=%s" % (k, "yes" if v else "NO") for k, v in HAVE.items()))
say("=" * 74)
say("%d proteins, %d sites" % (len(accs), len(feat)))
done = len(glob.glob(os.path.join(CACHE, "*.json")))
if done:
    say("resuming: %d proteins already cached" % done)

rows, t0 = [], time.time()
for i, acc in enumerate(accs, 1):
    if i % 100 == 0:
        say("  %d/%d proteins | %.0f min | eta %.0f min"
            % (i, len(accs), (time.time() - t0) / 60,
               (time.time() - t0) / 60 / i * (len(accs) - i)))
    try:
        got = protein_expand(acc, by_prot[acc])
    except Exception as e:
        continue
    for pos, rec in got.items():
        rows.append(dict(uniprot=acc, position=int(pos), **rec))

exp = pd.DataFrame(rows)
exp.to_csv("expanded_raw.csv", index=False)
merged = feat.merge(exp, on=["uniprot", "position"], how="left")
merged.to_csv("features_expanded.csv", index=False)

say("")
say("COVERAGE OF THE NEW PARAMETERS (%d sites)" % len(merged))
for c in ["in_pocket", "pocket_volume", "pocket_drugscore", "dist_nearest_pocket",
          "propka_pKa", "hbond_partners_NZ", "dssp_ss", "rsa_tien", "sasa_freesasa"]:
    if c in merged:
        n = merged[c].notna().sum()
        say("  %-22s %6d  %5.1f%%" % (c, n, 100 * n / len(merged)))
if "in_pocket" in merged:
    say("")
    say("  lysines placed in a detected pocket: %d (%.1f%%)"
        % (merged.in_pocket.sum(), 100 * merged.in_pocket.mean()))
if "propka_pKa" in merged and merged.propka_pKa.notna().any():
    p = merged.propka_pKa.dropna()
    say("  predicted pKa: median %.2f | 10th %.2f | 90th %.2f | below 9.0: %d"
        % (p.median(), p.quantile(.1), p.quantile(.9), (p < 9).sum()))
if "sasa_freesasa" in merged and merged.sasa_freesasa.notna().any():
    both = merged[["sasa_total", "sasa_freesasa"]].dropna()
    if len(both) > 10:
        say("  SASA cross-check, Shrake-Rupley vs FreeSASA: r = %.4f"
            % both.corr().iloc[0, 1])
open("step6d_report.txt", "w").write("\n".join(log) + "\n")
say("")
say("wrote features_expanded.csv and step6d_report.txt")
say("Blocks are labelled in the docstring; model them separately, correct within block.")
