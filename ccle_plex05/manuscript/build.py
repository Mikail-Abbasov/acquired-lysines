# -*- coding: utf-8 -*-
"""Assemble the single manuscript draft and the single Supporting Information draft.

Run from ccle_plex05/:   python3 manuscript/build.py

Outputs into ../results/:
  Manuscript_draft_<date>.docx
  Supporting_Information_draft_<date>.docx

House format: 11 pt Arial, single-spaced, no space before or after paragraphs, blank line between
paragraphs. Figures embedded at 6.5 in with the caption beneath. Figures carry no figure-level
title in the artwork.

References are cited in section text as [[key]] or [[key1,key2]] and numbered HERE, globally, in
order of first appearance — so a section can move without renumbering by hand.
"""
import os, re, sys, datetime, importlib.util
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
exec(open(os.path.expanduser("~/mk_section.py")).read())

import refs as REFMOD, sections as S, si as SI
spec = importlib.util.spec_from_file_location("drafted", os.path.join(HERE, "drafted.py"))
D = importlib.util.module_from_spec(spec); spec.loader.exec_module(D)

def _cap(c):
    """Caption blocks. Some are (text, is_title) tuples whose first entry repeats the figure
    title the build writes itself; some are plain strings with no title entry."""
    return [t for t, _ in c][1:] if c and isinstance(c[0], tuple) else list(c)

S.RESULTS[0]["body"], S.RESULTS[0]["cap"] = D.P_2_1, _cap(D.CAP_2_1)
S.RESULTS[1]["body"], S.RESULTS[1]["cap"] = D.P_2_2, _cap(D.CAP_2_2)
S.RESULTS[2]["body"], S.RESULTS[2]["cap"] = D.P_2_3, _cap(D.CAP_2_3)
S.RESULTS[3]["body"], S.RESULTS[3]["cap"] = D.P_2_4, _cap(D.CAP_2_4)
S.RESULTS[4]["body"], S.RESULTS[4]["cap"] = D.P_2_5, _cap(D.CAP_2_5)
S.RESULTS[5]["body"], S.RESULTS[5]["cap"] = D.P_2_6, _cap(D.CAP_2_6)
S.RESULTS[6]["body"], S.RESULTS[6]["cap"] = D.P_2_7, _cap(D.CAP_2_7)
S.RESULTS[7]["body"], S.RESULTS[7]["cap"] = D.P_2_8, _cap(D.CAP_2_8)

DATE = datetime.date.today().isoformat()
OUT = os.path.abspath(os.path.join(HERE, "..", "..", "results"))
FIGDIR = os.path.abspath(os.path.join(HERE, ".."))
PRINT = os.path.join(FIGDIR, "_build", "print")

# ---------------------------------------------------------------- reference numbering
ORDER, SEEN = [], set()
def _num(key):
    if key not in REFMOD.REFS:
        raise KeyError("unknown reference key: %s" % key)
    if key not in SEEN:
        SEEN.add(key); ORDER.append(key)
    return ORDER.index(key) + 1
def _fmt(nums):
    nums = sorted(nums); out, i = [], 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1: j += 1
        out.append("%d" % nums[i] if j == i else
                   ("%d,%d" % (nums[i], nums[j]) if j == i + 1 else "%d–%d" % (nums[i], nums[j])))
        i = j + 1
    # AMA style: superscript Arabic numerals, ranges collapsed. «» marks the run superscript.
    return "«" + ",".join(out) + "»"
def _hl_measured(items):
    """Highlights carry their own character count. Measure it here so the number can never
    drift from the text it describes; the bracket is not part of the bullet."""
    return ["%s  [%d]" % (t, len(t)) for t in items]

def cite(text):
    return re.sub(r"\[\[([a-z0-9,]+)\]\]",
                  lambda m: _fmt([_num(k) for k in m.group(1).split(",")]), text)

# ---------------------------------------------------------------- shrink figures for the draft
def print_copy(name):
    """A 2000 px-wide copy — about 300 dpi at 6.5 in — so the draft file stays manageable."""
    os.makedirs(PRINT, exist_ok=True)
    src, dst = os.path.join(FIGDIR, name), os.path.join(PRINT, name)
    if not os.path.exists(dst) or os.path.getmtime(dst) < os.path.getmtime(src):
        from PIL import Image
        im = Image.open(src); w, h = im.size
        im.resize((2000, int(2000 * h / w)), Image.LANCZOS).save(dst, optimize=True)
    return dst

def newdoc():
    doc = Document()
    for s in doc.sections:
        s.top_margin = s.bottom_margin = Inches(1); s.left_margin = s.right_margin = Inches(1)
    style(doc); return doc

def block(doc, paras, italic=False):
    for t in paras:
        para(doc, cite(t), italic=italic); blank(doc)

def todo(doc):
    para(doc, "[ NOT YET DRAFTED ]", bold=True); blank(doc)

# ================================================================ MANUSCRIPT
doc = newdoc()
para(doc, S.FRONT["title"], bold=True); blank(doc)
para(doc, S.FRONT["authors"]); blank(doc)
para(doc, S.FRONT["venue"], italic=True); blank(doc)
para(doc, "Draft assembled %s. %s" % (DATE, S.FRONT["title_note"]), italic=True); blank(doc)
para(doc, "Sections marked [ NOT YET DRAFTED ] carry their specification from the manuscript outline. "
          "Drafted sections carry prose, figure, caption and citations; the record of what was opened "
          "and closed while drafting them is in the section log, not here.", italic=True); blank(doc)

para(doc, "Highlights", bold=True); blank(doc)
para(doc, S.HIGHLIGHTS["spec"], italic=True); blank(doc)
for t in _hl_measured(S.HIGHLIGHTS["body"]):
    para(doc, "\u2022  " + t); blank(doc)

para(doc, "Graphical abstract", bold=True); blank(doc)
todo(doc); para(doc, S.GRAPHICAL_ABSTRACT["spec"], italic=True); blank(doc)

para(doc, "Abstract", bold=True); blank(doc)
block(doc, D.P_ABSTRACT)

para(doc, "Keywords", bold=True); blank(doc)
para(doc, "; ".join(S.KEYWORDS["body"])); blank(doc)
para(doc, S.KEYWORDS["note"], italic=True); blank(doc)

para(doc, "1.  Introduction", bold=True); blank(doc)
block(doc, D.P_INTRO)
para(doc, "2.  Results", bold=True); blank(doc)
for sec in S.RESULTS:
    para(doc, "%s  %s" % (sec["n"], sec["heading"]), bold=True); blank(doc)
    if sec["status"] == "todo":
        todo(doc); block(doc, sec["body"], italic=True)
    else:
        block(doc, sec["body"])
    figure(doc, print_copy(sec["fig"]), 6.5); blank(doc)
    para(doc, "Figure %d. %s." % (sec["figno"], sec["heading"]), bold=True); blank(doc)
    block(doc, sec["cap"], italic=(sec["status"] == "todo"))
para(doc, "3.  Discussion", bold=True); blank(doc)
block(doc, D.P_DISCUSSION)
para(doc, "4.  Conclusion", bold=True); blank(doc)
block(doc, D.P_CONCLUSION)
para(doc, "5.  Methods", bold=True); blank(doc)
block(doc, D.P_METHODS)
for name, st_, spec in S.DECLARATIONS:
    para(doc, name, bold=True); blank(doc)
    if st_ == "todo":
        todo(doc); para(doc, spec, italic=True); blank(doc)
    else:
        block(doc, S.DATA_AVAIL["body"])

para(doc, "References", bold=True); blank(doc)
para(doc, "Primary research papers only; no reviews. Every entry verified against PubMed or the "
          "article's own citation line. Numbered in order of first appearance.", italic=True); blank(doc)
for i, key in enumerate(ORDER, 1):
    para(doc, "%d.  %s" % (i, REFMOD.REFS[key]["text"])); blank(doc)
m_out = os.path.join(OUT, "Manuscript_draft_%s.docx" % DATE)
doc.save(m_out)
MAIN_ORDER = list(ORDER)
# The manuscript and the Supporting Information carry SEPARATE bibliographies, each numbered
# from 1. Resetting the counter here is what makes the SI number independently; a source cited
# in both documents therefore carries a different number in each, which is intended.
ORDER.clear(); SEEN.clear()

# ================================================================ SUPPORTING INFORMATION
sdoc = newdoc()
para(sdoc, "Supporting Information", bold=True); blank(sdoc)
para(sdoc, S.FRONT["title"], bold=True); blank(sdoc)
para(sdoc, S.FRONT["authors"]); blank(sdoc)
para(sdoc, "Draft assembled %s." % DATE, italic=True); blank(sdoc)
para(sdoc, "Supplementary figures are numbered sequentially in order of first citation and each "
           "names its parent figure. Figure 2 has two supplements, which is why one supplement per "
           "main figure was not used. Supplementary notes carry N-numbers and never S-numbers.",
     italic=True); blank(sdoc)
para(sdoc, "Contents", bold=True); blank(sdoc)
for f in SI.SI_FIGS:
    para(sdoc, "Figure %s — support for %s: %s" % (f["no"], f["parent"], f["title"]))
for t_ in getattr(SI, "SI_TABLES", []):
    para(sdoc, "Table %s — %s" % (t_["no"], t_["title"]))
blank(sdoc)
para(sdoc, "Supplementary notes follow the figures. Ten of the original seventeen topics survive "
           "the audit against the finalised figures; the rest are carried by a panel or rest on a "
           "withdrawn argument, and are listed with that status rather than silently dropped.")
blank(sdoc)
for f in SI.SI_FIGS:
    figure(sdoc, print_copy(f["fig"]), 6.5); blank(sdoc)
    para(sdoc, "Figure %s. Support for %s: %s" % (f["no"], f["parent"], f["title"]), bold=True); blank(sdoc)
    block(sdoc, f["cap"])
for t_ in getattr(SI, "SI_TABLES", []):
    para(sdoc, "Table %s. %s" % (t_["no"], t_["title"]), bold=True); blank(sdoc)
    table(sdoc, [t_["head"]] + [tuple(cite(c) for c in r) for r in t_["rows"]])
    blank(sdoc)
    block(sdoc, t_["foot"])

para(sdoc, "Supplementary notes", bold=True); blank(sdoc)
para(sdoc, SI.SI_NOTES_HEADER, italic=True); blank(sdoc)
kept    = [e[0] for e in SI.SI_NOTES if e[1] in ("keep", "rewrite")]
pending = [e[0] for e in SI.SI_NOTES if e[1] in ("keep", "rewrite") and len(e) < 4]
para(sdoc, "%d of the %d notes to be submitted %s written out below. Still to write: %s. The five "
           "marked cut and the two withdrawn need no text; they are listed so the decision is on the "
           "record and the numbering does not silently shift."
           % (len(kept) - len(pending), len(kept),
              "is" if len(kept) - len(pending) == 1 else "are", ", ".join(pending)),
     italic=True); blank(sdoc)
LABEL = {"keep": "", "cut": "  [CUT — carried by a panel]", "withdraw": "  [WITHDRAWN]",
         "rewrite": "  [REWRITTEN]"}
for e in SI.SI_NOTES:
    n, st_, t = e[0], e[1], e[2]
    para(sdoc, "%s.%s  %s%s" % (n, LABEL[st_], t, "" if len(e) > 3 else "   [ NOT YET WRITTEN ]"),
         bold=len(e) > 3); blank(sdoc)
    if len(e) > 3:
        block(sdoc, e[3])
para(sdoc, "Notes to appear in the submitted Supporting Information: %s (%d of 17). Renumber "
           "sequentially at submission." % (", ".join(kept), len(kept)), italic=True); blank(sdoc)
if ORDER:
    para(sdoc, "References", bold=True); blank(sdoc)
    para(sdoc, "The Supporting Information is numbered independently of the manuscript.", italic=True)
    blank(sdoc)
    for i, key in enumerate(ORDER, 1):
        para(sdoc, "%d. %s" % (i, REFMOD.REFS[key]["text"])); blank(sdoc)
s_out = os.path.join(OUT, "Supporting_Information_draft_%s.docx" % DATE)
sdoc.save(s_out)

print("manuscript : %s (%.1f MB)" % (m_out, os.path.getsize(m_out) / 1e6))
print("SI         : %s (%.1f MB)" % (s_out, os.path.getsize(s_out) / 1e6))
print("refs (ms)  : %d, in order: %s" % (len(MAIN_ORDER), ", ".join(MAIN_ORDER)))
print("refs (SI)  : %d, in order: %s" % (len(ORDER), ", ".join(ORDER)))
