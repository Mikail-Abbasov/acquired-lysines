#!/usr/bin/env python3
"""
step19_pocket_reactivity.py — are probe-labelled lysines found at cavities?

THE GAP THIS FILLS
Step 6d ran fpocket on all 6,201 proteins, but recorded pocket membership only for the 45,552
probe-labelled lysines — that was the site list it was given. Step 9 then computed geometry for
every lysine on those proteins (245,684 sites, labelled and unlabelled) using a separate
function that did not call fpocket. So the reactivity comparison in steps 11-12 could test
burial, crowding and charge, but never pocket membership. This closes that gap.

WHY IT IS NOT JUST ANOTHER SHAPE FEATURE
Crowding and burial count how much material is nearby. A pocket is different: a partially
enclosed cavity. A lysine can be crowded because it is packed into the core with nowhere for a
molecule to sit, or crowded because it lines the wall of a binding site — neighbour counts
cannot separate those, fpocket can.

And a pocket is a low-dielectric environment. Low dielectric destabilises charge, which would
LOWER a lysine's pKa and RAISE its reactivity. So pocket membership may act through
electrostatics rather than through accessibility — in which case it could survive the pLDDT
gate that collapsed every other shape feature (step 11: crowding +0.197 -> +0.014, atom density
+0.244 -> +0.003, polar +0.183 -> +0.043, all with flat SD).

WHAT EACH OUTCOME MEANS
  pocket membership SURVIVES the gate  -> reactive lysines sit preferentially at cavities.
      Probe labelling is then already enriching for sites where a molecule could bind, so
      reactivity and ligandability are partly coupled at the point of capture.
  pocket membership COLLAPSES like the other shape features -> the two properties are cleanly
      independent, which is the stronger form of the lab's standing distinction.

Both are publishable; they say different things about what the platform selects for.

DESIGN
  - fpocket re-run per protein (the step 6d cache kept only per-site values, not the pocket
    coordinates, so the cavities must be recomputed)
  - every lysine in the structure annotated, labelled and unlabelled alike
  - same within-protein paired test as step 11: one difference per protein, then tested against
    zero, so protein-level confounders cancel
  - same pLDDT gates (none / 50 / 70 / 90) with feature SD reported at each, so a decline in
    effect size can be checked against a decline in spread before it is called artefact removal
  - restartable: per-protein results cached as JSON

INPUT   step9_all_lysines.csv (which lysines are labelled), af/*.pdb
OUTPUT  step19_pockets.csv (per-site), step19_report.txt
"""
import os, re, sys, glob, json, time, shutil, subprocess, tempfile
import numpy as np, pandas as pd
import warnings
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser
from scipy import stats as sstats
from statsmodels.stats.multitest import multipletests

AF_DIR = os.environ.get("AF_DIR", "af")
CACHE = os.environ.get("CACHE19", "cache19")
GATES = [(0, "none"), (50, ">=50"), (70, ">=70"), (90, ">=90")]
POCKET_FEATS = ["in_pocket", "dist_nearest_pocket", "pocket_volume",
                "pocket_score", "pocket_drugscore"]
os.makedirs(CACHE, exist_ok=True)

log = []
def say(m):
    print(m, flush=True); log.append(m)


def pockets_for(acc):
    """Run fpocket and return, for every lysine NZ, its pocket annotation."""
    cf = os.path.join(CACHE, "%s.json" % acc)
    if os.path.exists(cf):
        try:
            return json.load(open(cf))
        except Exception:
            pass
    pdb = os.path.join(AF_DIR, "%s.pdb" % acc)
    if not os.path.exists(pdb):
        return {}
    work = tempfile.mkdtemp()
    dst = os.path.join(work, "%s.pdb" % acc)
    shutil.copy(pdb, dst)
    out = {}
    try:
        subprocess.run(["fpocket", "-f", dst], capture_output=True, timeout=600)
        odir = dst[:-4] + "_out"
        info = os.path.join(odir, "%s_info.txt" % acc)
        pockets, cur = [], None
        if os.path.exists(info):
            for line in open(info):
                m = re.match(r"Pocket\s+(\d+)", line)
                if m:
                    cur = {"id": int(m.group(1))}; pockets.append(cur)
                elif cur is not None and ":" in line:
                    k, v = [x.strip() for x in line.split(":", 1)]
                    try:
                        cur[k] = float(v)
                    except ValueError:
                        pass
        centres = {}
        for f in glob.glob(os.path.join(odir, "pockets", "pocket*_vert.pqr")):
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

        ch = next(PDBParser(QUIET=True).get_structure(acc, pdb)[0].get_chains())
        for r in ch:
            if r.id[0] != " " or r.get_resname() != "LYS" or "NZ" not in r:
                continue
            nz = np.array(r["NZ"].coord)
            best, dmin = None, np.inf
            for p in pockets:
                c = centres.get(p["id"])
                if c is None:
                    continue
                dd = float(np.min(np.linalg.norm(c - nz, axis=1)))
                if dd < dmin:
                    dmin, best = dd, p
            rec = {"dist_nearest_pocket": dmin if np.isfinite(dmin) else None,
                   "in_pocket": bool(np.isfinite(dmin) and dmin <= 5.0)}
            if best is not None:
                rec["pocket_volume"] = best.get("Volume")
                rec["pocket_score"] = best.get("Score")
                rec["pocket_drugscore"] = best.get("Druggability Score")
            out[str(r.id[1])] = rec
    except Exception:
        pass
    shutil.rmtree(work, ignore_errors=True)
    json.dump(out, open(cf, "w"))
    return out


say("=" * 78)
say("STEP 19 — POCKET MEMBERSHIP IN THE REACTIVITY COMPARISON")
say("=" * 78)
if shutil.which("fpocket") is None:
    say("fpocket not on PATH — install with: brew install brewsci/bio/fpocket")
    sys.exit(1)

allK = pd.read_csv("step9_all_lysines.csv", low_memory=False)
say("%d lysines on %d proteins | labelled %d (%.1f%%)"
    % (len(allK), allK.uniprot.nunique(), allK.labelled.sum(), 100 * allK.labelled.mean()))

accs = sorted(allK.uniprot.unique())
done = len(glob.glob(os.path.join(CACHE, "*.json")))
if done:
    say("resuming: %d proteins cached" % done)
rows, t0 = [], time.time()
for i, acc in enumerate(accs, 1):
    if i % 200 == 0:
        say("  %d/%d | %.0f min | eta %.0f min"
            % (i, len(accs), (time.time() - t0) / 60,
               (time.time() - t0) / 60 / i * (len(accs) - i)))
    for pos, rec in pockets_for(acc).items():
        rows.append(dict(uniprot=acc, position=int(pos), **rec))
pk = pd.DataFrame(rows)
say("pocket annotation computed for %d lysines" % len(pk))

d = allK.merge(pk, on=["uniprot", "position"], how="inner", suffixes=("", "_pk"))
for c in POCKET_FEATS:
    if c in d:
        d[c] = pd.to_numeric(d[c].replace({True: 1.0, False: 0.0}), errors="coerce")
say("merged: %d lysines with both geometry and pocket data" % len(d))
d.to_csv("step19_pockets.csv", index=False)

say("")
say("OVERALL POCKET OCCUPANCY")
for lab, sub in (("probe-labelled", d[d.labelled]), ("unlabelled", d[~d.labelled])):
    say("  %-16s %6d lysines | in a pocket %5.1f%% | median dist to pocket %.2f A"
        % (lab, len(sub), 100 * sub.in_pocket.mean(), sub.dist_nearest_pocket.median()))


def within(df, feats):
    out = []
    for f in feats:
        diffs = []
        for _, g in df.groupby("uniprot"):
            a = pd.to_numeric(g.loc[g.labelled, f], errors="coerce").dropna()
            b = pd.to_numeric(g.loc[~g.labelled, f], errors="coerce").dropna()
            if len(a) >= 2 and len(b) >= 2:
                diffs.append(a.mean() - b.mean())
        if len(diffs) < 30:
            continue
        v = np.array(diffs)
        t, p = sstats.ttest_1samp(v, 0)
        out.append(dict(feature=f, n_proteins=len(v), mean_diff=float(v.mean()),
                        d=float(v.mean() / v.std()) if v.std() else np.nan, p=float(p),
                        sd_all=float(pd.to_numeric(df[f], errors="coerce").std())))
    r = pd.DataFrame(out)
    if len(r):
        r["q"] = multipletests(r.p, method="fdr_bh")[1]
    return r


say("")
say("WITHIN-PROTEIN COMPARISON ACROSS pLDDT GATES")
say("  (effect / feature-SD at each gate; SD reported so a decline in effect can be checked")
say("   against a decline in spread before it is called artefact removal)")
res = []
for cut, lab in GATES:
    sub = d if cut == 0 else d[d.plddt_site >= cut]
    if sub.labelled.sum() < 200:
        continue
    r = within(sub, [f for f in POCKET_FEATS if f in sub.columns])
    if not len(r):
        continue
    r["gate"] = cut
    res.append(r)
    say("")
    say("  GATE %s — labelled %d, unlabelled %d"
        % (lab, int(sub.labelled.sum()), int((~sub.labelled).sum())))
    say("  %-22s %9s %10s %8s %10s %8s"
        % ("feature", "proteins", "mean diff", "d", "q", "SD"))
    for _, x in r.sort_values("q").iterrows():
        say("  %-22s %9d %10.4f %8.3f %10.3g %8.2f%s"
            % (x.feature[:22], x.n_proteins, x.mean_diff, x.d, x.q, x.sd_all,
               " *" if x.q < 0.05 else ""))

if res:
    R = pd.concat(res, ignore_index=True)
    R.to_csv("step19_gates.csv", index=False)
    say("")
    say("DOSE-RESPONSE")
    say("  %-22s %s" % ("feature", "  ".join("%13s" % g[1] for g in GATES)))
    for f in POCKET_FEATS:
        g = R[R.feature == f]
        if len(g) < 2:
            continue
        cells = []
        for cut, _ in GATES:
            row = g[g.gate == cut]
            cells.append("%6.3f/%-6.2f" % (row.d.iloc[0], row.sd_all.iloc[0])
                         if len(row) else "%13s" % "-")
        say("  %-22s %s" % (f[:22], "  ".join(cells)))
    say("")
    say("READ AGAINST STEP 11, where every shape feature collapsed at >=70 with flat SD")
    say("  (crowding +0.197 -> +0.014, atoms +0.244 -> +0.003, polar +0.183 -> +0.043)")
    say("  while electrostatics strengthened (anionic -0.307 -> -0.464).")
    say("")
    say("  If pocket membership HOLDS at >=70: reactive lysines sit at cavities, and probe")
    say("  capture is already enriching for sites a molecule could occupy — reactivity and")
    say("  ligandability are coupled at the point of capture.")
    say("  If it COLLAPSES: the two properties are cleanly independent.")

open("step19_report.txt", "w").write("\n".join(log) + "\n")
say("")
say("wrote step19_pockets.csv, step19_gates.csv, step19_report.txt")
