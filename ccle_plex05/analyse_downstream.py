#!/usr/bin/env python3
"""
analyse_downstream.py — collect the 42-plex downstream search and answer three questions.

  1. Does the tryptic boundary carry information?
     Downstream identifications vs BOUNDARY controls (same regions, shifted one residue,
     tryptic in neither protein). If boundary controls match at similar rates, the approach
     fails and that must be reported.

  2. What is the false-discovery proportion for this peptide class?
     From the paired shuffled entrapment, computed separately from the variant class. They are
     not assumed to behave alike: a downstream peptide is wild-type in sequence and can be
     confused with wild-type protein in ways a variant peptide cannot.

  3. Does it add sites, or confirm existing ones?
     Downstream identifications are cross-referenced against the mutant-peptide detections.
     Confirmation is valuable in its own right — no acquired lysine is currently detected in
     more than one experiment — but it is a different claim from discovery and is reported
     separately.

The A30 margin filter (Hyperscore - Nextscore >= 10) is applied throughout, as everywhere else
in this work. Plex 23 showed why: an unfiltered 7-residue match had expectation 0.89.

OUTPUT  downstream_summary.csv, downstream_sites.csv, downstream_report.txt
"""
import os, re, glob
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")

MARGIN = float(os.environ.get("MARGIN", 10))
log = []
def say(m):
    print(m, flush=True); log.append(m)

say("=" * 78)
say("DOWNSTREAM PEPTIDE SEARCH — 42 PLEXES")
say("margin filter: Hyperscore - Nextscore >= %.0f" % MARGIN)
say("=" * 78)

dirs = sorted(glob.glob("out_plex*_down"))
say("search directories found: %d" % len(dirs))

rows, per_plex = [], []
for d in dirs:
    plex = re.search(r"out_plex(\d+)_down", d).group(1)
    f = glob.glob(os.path.join(d, "*", "psm.tsv"))
    if not f:
        say("  plex%s: no psm.tsv" % plex); continue
    p = pd.read_csv(f[0], sep="\t", low_memory=False)
    pc = "Protein" if "Protein" in p.columns else [c for c in p.columns if "rotein" in c][0]
    p["margin"] = p.Hyperscore - p.Nextscore
    p["prot"] = p[pc].astype(str)
    # each spectrum appears twice (uncalibrated and uncalibrated_uncalibrated); collapse
    p["scan"] = p.Spectrum.astype(str).str.replace("_uncalibrated", "", regex=False)
    p = p.drop_duplicates(["scan", "Peptide"])
    keep = p[p.margin >= MARGIN]
    cnt = {}
    for lab, pre in (("down", "DOWN|"), ("bound", "BOUNDARY|"),
                     ("down_entrap", "PTARGET|DOWN|"), ("var", "VAR|"),
                     ("var_entrap", "PTARGET|VAR|")):
        cnt[lab] = int(keep.prot.str.startswith(pre).sum())
        cnt[lab + "_all"] = int(p.prot.str.startswith(pre).sum())
    per_plex.append(dict(plex=plex, **cnt))
    for _, r in keep[keep.prot.str.startswith("DOWN|")].iterrows():
        m = re.match(r"DOWN\|([^|]+)\|([^|]+)\|([^|]+)\|([^|]+)\|(\d+)\|", r.prot)
        if m:
            rows.append(dict(plex=plex, gene=m.group(1), change=m.group(2), line=m.group(3),
                             site=int(m.group(5)), peptide=r.Peptide, margin=r.margin,
                             hyper=r.Hyperscore, expect=r.Expectation))

S = pd.DataFrame(per_plex)
D = pd.DataFrame(rows)
if len(S):
    S.to_csv("downstream_summary.csv", index=False)
if len(D):
    D.to_csv("downstream_sites.csv", index=False)

say("")
say("1. DOES THE BOUNDARY CARRY INFORMATION?")
if len(S):
    say("  %-26s %10s %10s" % ("class", "margin>=%d" % MARGIN, "unfiltered"))
    for lab, name in (("down", "downstream target"), ("bound", "BOUNDARY control"),
                      ("down_entrap", "entrapment, downstream"),
                      ("var", "variant target"), ("var_entrap", "entrapment, variant")):
        say("  %-26s %10d %10d" % (name, S[lab].sum(), S[lab + "_all"].sum()))
    nd, nb = S["down"].sum(), S["bound"].sum()
    say("")
    if nd and not nb:
        say("  %d downstream identifications, ZERO boundary controls." % nd)
        say("  The tryptic boundary carries the information.")
    elif nd and nb >= 0.5 * nd:
        say("  %d downstream vs %d boundary — the boundary carries little." % (nd, nb))
        say("  THE APPROACH FAILS and should be reported as such.")
    elif nd:
        say("  %d downstream vs %d boundary — partial support, state the ratio." % (nd, nb))
    else:
        say("  No downstream peptides identified at this margin.")

say("")
say("2. FALSE DISCOVERY FOR THIS PEPTIDE CLASS")
if len(S):
    for lab, name in (("down", "downstream"), ("var", "variant")):
        t = S[lab].sum(); e = S[lab + "_entrap"].sum()
        if t:
            say("  %-12s %4d targets, %3d entrapment -> lower-bound FDP %.2f%%"
                % (name, t, e, 100 * e / t))
        else:
            say("  %-12s no targets at this margin" % name)
    say("  (variant-class FDP is reported at 4.70%% in the published analysis; this search is")
    say("   a subset of plexes and the two are not directly comparable)")

say("")
say("3. NEW SITES, OR CONFIRMATION OF EXISTING ONES?")
if len(D):
    say("  downstream identifications: %d PSMs at %d distinct sites"
        % (len(D), D.groupby(["gene", "change", "line", "site"]).ngroups))
    if os.path.exists("hits62.tsv"):
        h = pd.read_csv("hits62.tsv", sep="\t")
        known = set(zip(h.gene, h.line, h.site))
        D["already_detected"] = [(g, l, s) in known for g, l, s in zip(D.gene, D.line, D.site)]
        sites = D.groupby(["gene", "change", "line", "site"]).already_detected.max()
        say("  confirming a site the mutant peptide already found : %d" % int(sites.sum()))
        say("  NEW sites, not detected by any mutant peptide      : %d" % int((~sites).sum()))
        say("")
        new = D[~D.already_detected]
        if len(new):
            say("  NEW SITES")
            say("  %-8s %-10s %-12s %-12s %8s %10s"
                % ("plex", "gene", "change", "line", "margin", "expect"))
            for _, r in new.sort_values("margin", ascending=False).iterrows():
                say("  %-8s %-10s %-12s %-12s %8.1f %10.2e"
                    % (r.plex, r.gene, str(r.change)[:12], str(r.line)[:12], r.margin, r.expect))
    say("")
    say("  peptide length: median %d, min %d, max %d"
        % (D.peptide.str.len().median(), D.peptide.str.len().min(), D.peptide.str.len().max()))
    short = D[D.peptide.str.len() <= 8]
    if len(short):
        say("  %d PSMs from peptides of 8 residues or fewer — check these individually,"
            % len(short))
        say("  short peptides match more readily by chance even above the margin threshold")
else:
    say("  no downstream identifications at this margin")

open("downstream_report.txt", "w").write("\n".join(log) + "\n")
say("")
say("wrote downstream_summary.csv, downstream_sites.csv, downstream_report.txt")
