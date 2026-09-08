#!/usr/bin/env python3
"""
step31_crystal_validation.py — the experimental leg of Figure 8, in real crystal structures.

THE QUESTION
Everything else in Figure 8 is computed on AlphaFold models. This asks the same question of
experimental structures: when a cancer mutation puts a lysine where a glutamate was, does the
surrounding structure move?

THE DESIGN, from step 26b
PTPN11 (SHP2) position 76 is the strongest case in the whole PDB: three independent E76K
structures, plus E76A, E76D and E76Q at the SAME position, plus wild-type structures. That
substitution series separates two things that mutant-vs-wild-type alone cannot:
    "mutating this residue perturbs the site"   vs   "introducing a LYSINE perturbs the site"
TP53 R280K (6FF9) is carried as a second, independent site.

THE MEASUREMENT, and why it is built this way
For each pair of structures:
  1. superpose on CA atoms of a SURROUNDING FRAME, 12-25 A from the site -- never on the local
     shell itself, since superposing on the shell would minimise exactly the quantity measured
  2. report CA RMSD over the LOCAL shell (<= 8 A of the site)
  3. separately, superpose on EVERY residue > 12 A and record that RMSD too. When the whole-
     protein alignment fails while the local frame aligns, the two structures differ by a
     DOMAIN REARRANGEMENT rather than by anything local -- which for SHP2 E76K is a real and
     well-known effect (E76K is activating; it opens the autoinhibited conformation), not a
     defect in the data. The first version of this script used the global frame only and threw
     those pairs away, which would have hidden the most interesting thing in the set.
Comparisons: WT vs WT gives the noise floor -- how much two crystals of the same protein differ
for reasons that have nothing to do with the mutation. Every substitution is read against it.

INPUT   RCSB (downloaded to xtal/)
OUTPUT  step31_crystal.csv, step31_report.txt

NUMBERING NOTE. This was first written as step29 and its report overwrote step29_report.txt,
which fig5 parses for the label-reproducibility ceiling -- fig5 failed on the next full build
with KeyError: 'NatChem'. Renamed to step31 (step30 is the compound-promiscuity analysis) and
step29_report.txt regenerated from step29_ceiling_fast.py. Check for a free number before
adding a step.
"""
import os, itertools, urllib.request
import numpy as np, pandas as pd
import warnings
warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser, Superimposer

XT = "xtal"
os.makedirs(XT, exist_ok=True)
log = []
def say(s=""):
    print(s, flush=True); log.append(s)

SITES = {
    "PTPN11 76": dict(pos=76, wt_aa="GLU", groups={
        "wild type":  ["2SHP", "4DGP", "4DGX", "4H1O", "5EHR"],
        "E76K":       ["6CMS", "6CRF", "6CRG"],
        "E76A":       ["7EMN"],
        "E76D":       ["6CMR"],
        "E76Q":       ["5IBS"],
    }),
    "TP53 280": dict(pos=280, wt_aa="ARG", groups={
        "wild type":  ["2OCJ", "1TSR", "3KMD", "2AC0"],
        "R280K":      ["6FF9"],
    }),
    # added 2026-09-04: the two acquired lysines that matter most to the paper's argument.
    # AKT1 E17K is the site Craven et al. drugged (8UVY is that structure, "AKT1(E17K) with
    # compound 3"); 2UZS is Carpten's original E17K PH domain. PIK3CA E542K / E545K are the
    # helical-domain hotspots -- 8GUA and 8GUD are cryo-EM, so their comparison carries a
    # method difference that crystal-vs-crystal pairs do not, and that is stated rather than
    # hidden. Both were on step 26b's list and neither had been measured.
    "AKT1 17": dict(pos=17, wt_aa="GLU", groups={
        "wild type":  ["1UNQ", "1H10", "1UNR", "1UNP"],
        "E17K":       ["2UZS", "8UVY"],
    }),
    "PIK3CA 542": dict(pos=542, wt_aa="GLU", groups={
        "wild type":  ["4OVU", "2RD0", "4L23", "5DXT"],
        "E542K":      ["8GUA"],
    }),
    "PIK3CA 545": dict(pos=545, wt_aa="GLU", groups={
        "wild type":  ["4OVU", "2RD0", "4L23", "5DXT"],
        "E545K":      ["8GUD"],
    }),
    # NRAS Q61K, added 2026-09-04. Residue 61 sits in switch II, which is the single most
    # nucleotide-state-dependent loop in a RAS protein, so this site's wild-type floor is large
    # BY NATURE and not because the measurement is poor. That is precisely why every site is
    # read against its own floor rather than against zero. WT set deliberately spans states
    # (5UHV GppNHp, 8TBI GppNHp in a tricomplex, 6ZIO GDP) so the floor reflects the real
    # conformational range of the loop. Only one mutant structure exists (8VM2, Q61K.GTP).
    #
    # STATE-MATCHED ON PURPOSE. A first run included 6ZIO (WT.GDP) among the wild types. That
    # put a GDP-vs-GTP pair into the floor, inflating it to 4.25 A and making Q61K (1.81 A)
    # look far better than it is -- an error in the direction that flatters the paper. The set
    # is now GppNHp/GTP only, so every pair compares the same nucleotide state. The cost is
    # that the floor rests on ONE pair (5UHV vs 8TBI); the legend says so.
    "NRAS 61": dict(pos=61, wt_aa="GLN", groups={
        "wild type":  ["5UHV", "8TBI"],
        "Q61K":       ["8VM2"],
    }),
}

def fetch(pid):
    p = os.path.join(XT, "%s.pdb" % pid)
    if os.path.exists(p) and os.path.getsize(p) > 5000:
        return p
    try:
        urllib.request.urlretrieve("https://files.rcsb.org/download/%s.pdb" % pid, p)
        return p if os.path.getsize(p) > 5000 else None
    except Exception:
        return None

parser = PDBParser(QUIET=True)

def load(pid, pos, want=None):
    """Return EVERY chain that carries the site. Many of these structures hold several copies
    in the asymmetric unit, and matching the wrong pair of copies produces a superposition that
    fails silently -- the first run of this script reported 24-27 A 'local RMSD' for exactly
    that reason. Returning all candidate chains lets the caller pick the pair that actually
    aligns, and the far-atom RMSD guard below rejects any pair that does not."""
    p = fetch(pid)
    if p is None:
        return None, None, "download failed"
    try:
        model = parser.get_structure(pid, p)[0]
    except Exception:
        return None, None, "parse failed"
    out = []
    for ch in model:
        res = {r.id[1]: r for r in ch if r.id[0] == " " and "CA" in r}
        if pos in res and len(res) >= 60:
            if want is not None and res[pos].get_resname() != want:
                continue          # e.g. 8UVY chain B is the nanobody, which also has a 17
            out.append(res)
    if not out:
        return None, None, "position %d absent" % pos
    out.sort(key=len, reverse=True)
    return out, out[0][pos], None

say("=" * 78)
say("STEP 31 - DOES AN ACQUIRED LYSINE MOVE THE STRUCTURE AROUND IT?")
say("=" * 78)
say("Experimental crystal structures, not models. Superposition on residues > 12 A from the")
say("site; RMSD reported over CA atoms <= 8 A of it. WT-vs-WT is the noise floor.")

rows = []
for site, cfg in SITES.items():
    pos = cfg["pos"]
    say()
    say("-" * 78)
    say("%s  (wild type %s)" % (site, cfg["wt_aa"]))
    say("-" * 78)
    loaded = {}
    for grp, pids in cfg["groups"].items():
        for pid in pids:
            want_aa = cfg["wt_aa"] if grp == "wild type" else {
                "K": "LYS", "A": "ALA", "D": "ASP", "Q": "GLN"}[grp[-1]]
            res, r, err = load(pid, pos, want_aa)
            if err:
                say("  %-6s %-10s SKIPPED - %s" % (pid, grp, err)); continue
            loaded[pid] = (grp, res, r)
            say("  %-6s %-10s residue at %d is %s, %d chain(s) carry it, largest %d residues"
                % (pid, grp, pos, r.get_resname(), len(res), len(res[0])))
    # identity check: the group label must match what is in the file
    for pid, (grp, res, r) in list(loaded.items()):
        want = cfg["wt_aa"] if grp == "wild type" else {
            "K": "LYS", "A": "ALA", "D": "ASP", "Q": "GLN"}[grp[-1]]
        if r.get_resname() != want:
            say("  !! %s labelled %s but carries %s at %d - DROPPED"
                % (pid, grp, r.get_resname(), pos))
            del loaded[pid]
    if len(loaded) < 2:
        say("  too few usable structures"); continue

    say()
    say("  %-16s %-16s %7s %7s %8s %8s %8s %-7s" % ("A", "B", "frameCA", "shellCA", "frameRMS", "SHELL", "shared", "frame"))
    FAR_GUARD = 2.0          # A; above this the two chains are not the same thing
    # MIN_SHELL was 5 and silently excluded 8UVY -- the drugged AKT1 E17K structure -- because
    # its co-local shell holds 4 residues, not 5. AKT1 residue 17 sits on a small PH-domain
    # loop with few neighbours, so 4 is what the site has. Shell size is printed for every
    # pair; an RMSD over 4 residues is noisier than one over 10 and should be read that way.
    MIN_SHELL = 4
    n_reject = 0
    for a, b in itertools.combinations(sorted(loaded), 2):
        ga, chains_a = loaded[a][0], loaded[a][1]
        gb, chains_b = loaded[b][0], loaded[b][1]
        best = None
        for resa in chains_a:
            for resb in chains_b:
                common = sorted(set(resa) & set(resb))
                if len(common) < 60:
                    continue
                ca_a = np.array(resa[pos]["CA"].coord)
                d_a = {i: float(np.linalg.norm(np.array(resa[i]["CA"].coord) - ca_a))
                       for i in common}
                d_b = {i: float(np.linalg.norm(np.array(resb[i]["CA"].coord) -
                                               np.array(resb[pos]["CA"].coord)))
                       for i in common}
                # CO-LOCAL definition. The first attempt defined the shell by structure A's
                # geometry alone. For SHP2 that is wrong: residue 76 sits at the N-SH2 / PTP
                # interface, so in the CLOSED structures its 8 A neighbourhood contains PTP
                # residues that are nowhere near it once the protein OPENS (6CRF). Comparing
                # those positions across states measures the domain opening, not the site.
                # Shell and frame are therefore residues that are local in BOTH structures.
                shell = [i for i in common if d_a[i] <= 8 and d_b[i] <= 8 and i != pos]
                glob = [i for i in common if d_a[i] > 12 and d_b[i] > 12]
                # and the turnover itself is reported, because it is the more direct measure
                near_a = {i for i in common if d_a[i] <= 8 and i != pos}
                near_b = {i for i in common if d_b[i] <= 8 and i != pos}
                jac = len(near_a & near_b) / max(len(near_a | near_b), 1)
                frame, tier = None, None
                for lo, hi, tag in ((12, 25, "12-25A"), (8, 16, "8-16A")):
                    cand = [i for i in common
                            if lo < d_a[i] <= hi and lo < d_b[i] <= hi]
                    if len(cand) < 25 or len(shell) < MIN_SHELL:
                        continue
                    t = Superimposer()
                    t.set_atoms([resa[i]["CA"] for i in cand], [resb[i]["CA"] for i in cand])
                    if frame is None or float(t.rms) < 2.0:
                        frame, tier, sup, frame_rms = cand, tag, t, float(t.rms)
                    if float(t.rms) < 2.0:
                        break
                if frame is None:
                    continue
                if best is None or frame_rms < best[0]:
                    moved = [resb[i]["CA"].copy() for i in shell]
                    sup.apply(moved)
                    rmsd = float(np.sqrt(np.mean(
                        [np.sum((np.array(m.coord) - np.array(resa[i]["CA"].coord)) ** 2)
                         for m, i in zip(moved, shell)])))
                    gsup = Superimposer()
                    if len(glob) >= 40:
                        gsup.set_atoms([resa[i]["CA"] for i in glob], [resb[i]["CA"] for i in glob])
                        g = float(gsup.rms)
                    else:
                        g = float("nan")
                    best = (frame_rms, rmsd, len(frame), len(shell), g, tier, jac)
        if best is None:
            continue
        far_rms, rmsd, nf, ns, glob_rms, tier, jac = best
        if far_rms > FAR_GUARD:
            n_reject += 1
            say("  %-16s %-16s   REJECTED - even the local frame will not align (%.2f A)"
                % ("%s %s" % (a, ga), "%s %s" % (b, gb), far_rms))
            continue
        domain = (jac < 0.80) or (glob_rms == glob_rms and glob_rms > 3.0)
        kind = ("WT vs WT" if ga == gb == "wild type" else
                ("K vs WT" if "wild type" in (ga, gb) and (ga.endswith("K") or gb.endswith("K"))
                 else ("other sub vs WT" if "wild type" in (ga, gb) else "mutant vs mutant")))
        rows.append(dict(site=site, a=a, b=b, group_a=ga, group_b=gb, kind=kind,
                         n_frame=nf, n_shell=ns, frame_rms=far_rms, rmsd=rmsd,
                         global_rms=glob_rms, domain_shift=domain, frame_tier=tier,
                         shared_neighbours=jac))
        say("  %-16s %-16s %7d %7d %8.3f %8.3f %7.0f%% %-7s%s"
            % ("%s %s" % (a, ga), "%s %s" % (b, gb), nf, ns, far_rms, rmsd, 100 * jac, tier,
               "  <- DOMAIN SHIFT" if domain else ""))
    if n_reject:
        say("  (%d pair(s) rejected by the far-atom alignment guard)" % n_reject)

D = pd.DataFrame(rows)
D.to_csv("step31_crystal.csv", index=False)

say()
say("=" * 78)
say("SUMMARY - local CA RMSD (A) within 8 A of the site")
say("=" * 78)
say("  %-14s %-22s %5s %8s %8s %8s %10s"
    % ("site", "comparison", "n", "median", "min", "max", "neighbours"))
for site in D.site.unique():
    ss = D[D.site == site]
    for kind in ["WT vs WT", "K vs WT", "other sub vs WT", "mutant vs mutant"]:
        k = ss[ss.kind == kind]
        if len(k):
            say("  %-14s %-22s %5d %8.3f %8.3f %8.3f %9.0f%%"
                % (site, kind, len(k), k.rmsd.median(), k.rmsd.min(), k.rmsd.max(),
                   100 * k.shared_neighbours.median()))
say()
say("DOMAIN REARRANGEMENT, reported separately from the local question")
dm = D[D.domain_shift]
if len(dm):
    say("  %d pair(s) align locally but NOT globally (whole-protein CA RMSD > 3 A):" % len(dm))
    for _, r in dm.iterrows():
        say("     %-6s %-10s vs %-6s %-10s  local %.2f A, %.0f%% of neighbours shared"
            % (r.a, r.group_a, r.b, r.group_b, r.rmsd, 100 * r.shared_neighbours))
    say("  6CRF is deposited as 'SHP2 E76K GOF mutant in the OPEN conformation' -- its own")
    say("  title. SHP2 E76K is an ACTIVATING mutation: it releases the autoinhibited, closed")
    say("  conformation. A large global difference alongside an ordinary local one is that")
    say("  effect, and it is an argument for the figure rather than against it -- the site")
    say("  itself is unremarkable even when the protein as a whole has moved.")
else:
    say("  none - every pair that aligns locally also aligns globally")
say()
say("READING IT")
say("  If 'K vs WT' sits inside the range of 'WT vs WT', the acquired lysine moves the")
say("  surrounding structure no more than two crystals of the same protein differ from each")
say("  other -- i.e. not at all, at this resolution. If 'other sub vs WT' is similar to")
say("  'K vs WT', then whatever movement there is belongs to mutating the residue, not to")
say("  introducing a lysine specifically.")
say()
say("LIMITS, stated because they are real")
say("  Crystal structures differ for crystallographic reasons - space group, resolution,")
say("  ligands, temperature - and this measure cannot separate those from mutation effects.")
say("  That is exactly why WT-vs-WT is the comparator rather than zero.")

open("step31_report.txt", "w").write("\n".join(log) + "\n")
print("\nwrote step31_crystal.csv, step31_report.txt")
