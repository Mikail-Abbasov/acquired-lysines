#!/usr/bin/env python3
"""
build_manifest.py -- generate the figure -> script -> data manifest by reading the code.

WHY IT IS GENERATED AND NOT TYPED
The working tree holds 255 Python scripts and 95 .bak_* files; eighteen of them build the
published figures. A hand-written manifest would be wrong within a week and there would be no way
to tell. This reads the scripts themselves, so re-running it after any edit re-derives the truth.

WHAT IT RECORDS, per published figure
  script        the file that draws it
  reads         every data file the script names
  writes        every file the script writes
  produced_by   for each input, the script elsewhere in the tree that writes that filename
  numbers       the headline values the script prints into its HEAD dict, which is the paper's
                own mechanism for "every number is printed by the script that makes the panel"

WHAT IT DELIBERATELY DOES NOT DO
It does not execute anything. Inputs are collected by pattern, so a filename built at runtime is
missed; those are listed under `unresolved` rather than silently dropped. Verifying that the
numbers still match the manuscript is a separate job (verify.py) and this file is its input.
"""
import os, re, sys, csv, json
from collections import defaultdict

TOP = os.path.expanduser("~/mnt/DepMap")
SEARCH_DIRS = ["ccle", "ligand_rebuild", "Datasets"]

# The published figure set, in the order the manuscript presents it.
FIGURES = [
    ("Figure 1",  "fig1_v3.py"),   ("Figure 2",  "fig2_v3.py"),
    ("Figure 3",  "fig3_v3.py"),   ("Figure 4",  "fig4_v3.py"),
    ("Figure 5",  "fig5_v1.py"),   ("Figure 6",  "fig6_v1.py"),
    ("Figure 7",  "fig7_v1.py"),   ("Figure 8",  "fig8_v1.py"),
    ("Figure S1", "figS1_v3.py"),  ("Figure S2", "figS2_v3.py"),
    ("Figure S3", "figS3_v3.py"),  ("Figure S4", "figS4_v3.py"),
    ("Figure S5", "figS5_v3.py"),  ("Figure S6", "figS6_v1.py"),
    ("Figure S7", "figS7_v1.py"),  ("Figure S8", "figS8_v1.py"),
    ("Figure S9", "figS9_v1.py"),  ("Figure S10","figS10_v1.py"),
]

DATA_EXT = (".csv", ".tsv", ".txt", ".json", ".xlsx", ".pkl", ".parquet", ".fasta", ".faa",
            ".gz", ".maf", ".pdb", ".mzML",
            # .out is a data extension in this tree: germline_lysine.py and germline_lysine_nc.py
            # write their results as .out, and Figure 1g reads both. Without this they were
            # invisible to the manifest and were left out of the release, which the acceptance
            # test caught as five NaN values in panel g.
            ".out")
READ  = re.compile(r'(?:find|read_csv|read_excel|read_table|read_pickle|read_parquet|open|load)'
                   r'\s*\(\s*["\']([^"\']+)["\']')
# Several scripts build paths with os.path.join(DIR, "name.csv"); the filename is the last
# string literal inside the join. Without this, fig8, figS9 and figS10 report no inputs at all.
READ_J  = re.compile(r'(?:read_csv|read_excel|read_table|read_pickle|read_parquet|open)'
                     r'\s*\(\s*os\.path\.join\(([^)]*)\)', re.S)
WRITE_J = re.compile(r'(?:to_csv|to_json|to_pickle|to_parquet|savefig|dump)'
                     r'\s*\(\s*os\.path\.join\(([^)]*)\)', re.S)
LIT     = re.compile(r'["\']([^"\']+)["\']')
READ2 = re.compile(r'(?:read_csv|read_excel|read_table)\s*\(\s*find\(\s*["\']([^"\']+)["\']')
WRITE = re.compile(r'(?:to_csv|to_json|to_pickle|to_parquet|savefig|dump)\s*\(\s*["\']([^"\']+)["\']')
WRITE2= re.compile(r'open\s*\(\s*["\']([^"\']+)["\']\s*,\s*["\']w')
# open(os.path.join(DIR, "name.txt"), "w") -- the join form of WRITE2. Without this, a script
# that writes its report through os.path.join looks like it never writes it, and the file then
# reads as an external input. figS10_report.txt was exactly that.
WRITE_J2 = re.compile(r'open\s*\(\s*os\.path\.join\(([^)]*)\)\s*,\s*["\']w')
HEADK = re.compile(r'HEAD\[\s*["\']([^"\']+)["\']\s*\]')
SAVE  = re.compile(r'st\.save\s*\(\s*\w+\s*,\s*["\']([^"\']+)["\']')

def find_script(name):
    for d in SEARCH_DIRS:
        p = os.path.join(TOP, d, name)
        if os.path.exists(p):
            return p
    return None

def _joined(pattern, src):
    """Filenames assembled with os.path.join(...): take the last string literal in the call."""
    out = []
    for grp in pattern.findall(src):
        lits = [x for x in LIT.findall(grp) if x.endswith(DATA_EXT)]
        if lits:
            out.append(lits[-1])
    return out

def scan(path):
    s = open(path, encoding="utf-8", errors="replace").read()
    reads = {m for m in READ.findall(s) + READ2.findall(s) + _joined(READ_J, s)
             if m.endswith(DATA_EXT)}
    writes = {m for m in WRITE.findall(s) + WRITE2.findall(s) + _joined(WRITE_J, s)
              + _joined(WRITE_J2, s)
              if m.endswith(DATA_EXT)}
    for stem in SAVE.findall(s):
        writes.add(stem + ".png"); writes.add(stem + ".pdf")
    return sorted(reads), sorted(writes), sorted(set(HEADK.findall(s))), s.count("print(")

# ---- index every script in the tree by what it writes, so inputs can be traced to a producer
producers = defaultdict(list)
all_scripts = []
for d in SEARCH_DIRS:
    root = os.path.join(TOP, d)
    if not os.path.isdir(root):
        continue
    for f in sorted(os.listdir(root)):
        if not f.endswith(".py") or ".bak" in f:
            continue
        p = os.path.join(root, f)
        all_scripts.append(os.path.join(d, f))
        try:
            _, w, _, _ = scan(p)
        except Exception:
            continue
        for out in w:
            producers[os.path.basename(out)].append(os.path.join(d, f))

rows, unresolved = [], []
for figname, script in FIGURES:
    p = find_script(script)
    if p is None:
        rows.append(dict(figure=figname, script=script + "  [NOT FOUND]", reads="", writes="",
                         produced_by="", numbers="", prints=0))
        continue
    reads, writes, heads, nprint = scan(p)
    prod = []
    for r in reads:
        b = os.path.basename(r)
        src = [x for x in producers.get(b, []) if os.path.basename(x) != script]
        if src:
            prod.append("%s <- %s" % (b, src[0]))
        else:
            prod.append("%s <- EXTERNAL or unresolved" % b)
            unresolved.append((figname, b))
    rows.append(dict(figure=figname, script=os.path.relpath(p, TOP),
                     reads="; ".join(reads), writes="; ".join(writes),
                     produced_by="; ".join(prod), numbers="; ".join(heads), prints=nprint))

with open("manifest.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["figure", "script", "reads", "writes", "produced_by", "numbers", "prints"])
    w.writeheader(); w.writerows(rows)

# ---- readable version
with open("manifest.md", "w") as fh:
    fh.write("# Figure manifest\n\nGenerated by `build_manifest.py` from the scripts themselves. "
             "Re-run after any edit.\n\n")
    fh.write("Scripts scanned: %d across %s. Published figures: %d.\n\n"
             % (len(all_scripts), ", ".join(SEARCH_DIRS), len(FIGURES)))
    for r in rows:
        fh.write("## %s\n\n**Script:** `%s`\n\n" % (r["figure"], r["script"]))
        if r["reads"]:
            fh.write("**Reads**\n\n")
            for x in r["produced_by"].split("; "):
                fh.write("- %s\n" % x)
            fh.write("\n")
        if r["writes"]:
            fh.write("**Writes:** %s\n\n" % ", ".join(r["writes"].split("; ")))
        if r["numbers"]:
            fh.write("**Headline numbers printed:** %s\n\n" % r["numbers"])

print("=" * 78)
print("FIGURE MANIFEST")
print("=" * 78)
print("  scripts in tree (excluding .bak): %d" % len(all_scripts))
print("  published figures                : %d" % len(FIGURES))
missing = [r["figure"] for r in rows if "NOT FOUND" in r["script"]]
print("  figure scripts not found         : %s" % (", ".join(missing) if missing else "none"))
print()
print("  %-11s %-16s %5s %6s %8s %7s" % ("figure","script","reads","writes","HEAD nums","prints"))
for r in rows:
    print("  %-11s %-16s %5d %6d %8d %7d"
          % (r["figure"], os.path.basename(r["script"]),
             len(r["reads"].split("; ")) if r["reads"] else 0,
             len(r["writes"].split("; ")) if r["writes"] else 0,
             len(r["numbers"].split("; ")) if r["numbers"] else 0, r["prints"]))
# Fallback: an input with no detected writer may still be produced by a script that assembles the
# name at runtime. Grep the tree for the literal and report candidates rather than calling it
# external, so a human resolves it instead of the manifest quietly under-reporting.
mentions = defaultdict(list)
for rel in all_scripts:
    try:
        txt = open(os.path.join(TOP, rel), encoding="utf-8", errors="replace").read()
    except Exception:
        continue
    for b in {b for _, b in unresolved}:
        if b in txt and os.path.basename(rel) not in [s2 for _, s2 in FIGURES]:
            mentions[b].append(rel)
json.dump({k: v for k, v in mentions.items()}, open("unresolved_inputs.json", "w"), indent=1)
nin = len({b for _, b in unresolved})
print()
print("  distinct inputs with no producer in the tree (external or built at runtime): %d" % nin)
for b in sorted({b for _, b in unresolved})[:25]:
    cand = mentions.get(b, [])
    print("    %-42s %s" % (b, ("mentioned in " + cand[0]) if cand else "EXTERNAL INPUT"))
print("\nwrote manifest.csv and manifest.md")
