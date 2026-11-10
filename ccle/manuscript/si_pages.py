#!/usr/bin/env python3
"""Fill the cached page numbers in the Supporting Information contents.

Each contents line ends in a PAGEREF field with a cached result, so it reads correctly the
moment the file is opened rather than showing "Error! Bookmark not defined." until someone
presses F9. Word recomputes the true page on an update-fields; this only has to get the cache
right. The pagination comes from converting the built file and reading which page each legend
lands on, so the numbers are the ones the PDF will carry.

    python3 si_pages.py      # convert, read the pages, write results/panel_data/si_pages.json

Run build.py again afterwards to bake them in.
"""
import os, re, json, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
SUBMIT = os.path.abspath(os.path.join(HERE, "..", "..", "CCB-submission"))
OUT = os.path.join(HERE, "..", "results", "panel_data", "si_pages.json")
DOCX = os.path.join(SUBMIT, "Supporting_Information.docx")

import si as SI

def pages():
    tmp = os.path.join(os.path.expanduser("~"), "_sipages")
    os.makedirs(tmp, exist_ok=True)
    subprocess.run(["soffice", "--headless", "-env:UserInstallation=file://%s/lo" % tmp,
                    "--convert-to", "pdf", "--outdir", tmp, DOCX],
                   check=True, capture_output=True, timeout=600)
    pdf = os.path.join(tmp, "Supporting_Information.pdf")
    import pymupdf
    d = pymupdf.open(pdf)
    text = [p.get_text() for p in d]
    found = {}
    # The legend's own opening is the anchor. The contents line is a different wording, so a
    # search for it would land on the contents page every time.
    for f in SI.SI_FIGS:
        key = "Figure %s. Support for" % f["no"]
        for i, t in enumerate(text):
            if key in t:
                found["SIFIG_%s" % f["no"]] = i + 1
                break
    # A table's own title is word for word its contents line, so a search from page 1 lands on
    # the contents every time. Start past the first figure, which is always after the contents.
    start = min([v for k, v in found.items() if k.startswith("SIFIG_")] or [1]) - 1
    for t_ in getattr(SI, "SI_TABLES", []):
        key = "Table %s. %s" % (t_["no"], t_["title"][:40])
        for i in range(start, len(text)):
            if key in " ".join(text[i].split()):
                found["SITAB_%s" % t_["no"]] = i + 1
                break
    # The notes, the two headings and Table S2. Each anchor is the text as the document prints
    # it, and the search starts past the contents so a line in the list cannot match itself.
    flat = [" ".join(t.split()) for t in text]
    start = min([v for k, v in found.items() if k.startswith("SIFIG_")] or [1]) - 1
    anchors = {"SITAB_S2": "Table S2. Panel-level values",
               "SIEXT": "Analyses that support the main-text figures",
               "SINOTES": "Supplementary notes"}
    _keep = [e for e in SI.SI_NOTES if e[1] in ("keep", "rewrite") and len(e) > 3]
    for i, e in enumerate(_keep):
        anchors["SINOTE_N%d" % (i + 1)] = "N%d. %s" % (i + 1, e[2][:40])
    for bm, key in anchors.items():
        k = " ".join(key.split())
        for i in range(start, len(flat)):
            if k in flat[i]:
                found[bm] = i + 1
                break
    # References is the last heading in the file, so take the LAST match, never the first.
    for i in range(len(flat) - 1, start - 1, -1):
        if re.search(r"\bReferences\b", flat[i]):
            found["SIREFS"] = i + 1
            break
    d.close()
    return found, len(text)

if __name__ == "__main__":
    f, n = pages()
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(f, fh, indent=1, sort_keys=True)
    print("%d items located across %d pages" % (len(f), n))
    for k in sorted(f, key=lambda x: f[x]):
        print("   %-12s page %d" % (k, f[k]))
    print("wrote", os.path.abspath(OUT))
