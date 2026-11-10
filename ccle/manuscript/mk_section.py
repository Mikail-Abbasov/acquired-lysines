# -*- coding: utf-8 -*-
"""Build a manuscript section .docx: 11 pt Arial, single-spaced, no paragraph spacing,
with the section's figure embedded and captioned."""
import sys, re
from docx import Document
from docx.shared import Pt, Inches, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn

FONT, SIZE = "Arial", Pt(11)

# REAL SUPERSCRIPTS AND SUBSCRIPTS, not Unicode look-alikes. A glyph like U+00B2 carries its
# own design, so it does not match the surrounding Arial and a journal's typesetter cannot
# restyle it. Every one of them is written as the ordinary character in a run with Word's
# superscript or subscript property set, which is the button in Word's ribbon.
_SUP = {"\u2070": "0", "\u00b9": "1", "\u00b2": "2", "\u00b3": "3", "\u2074": "4",
        "\u2075": "5", "\u2076": "6", "\u2077": "7", "\u2078": "8", "\u2079": "9",
        "\u207a": "+", "\u207b": "-", "\u207c": "=", "\u207d": "(", "\u207e": ")",
        "\u207f": "n"}
_SUB = {"\u2080": "0", "\u2081": "1", "\u2082": "2", "\u2083": "3", "\u2084": "4",
        "\u2085": "5", "\u2086": "6", "\u2087": "7", "\u2088": "8", "\u2089": "9",
        "\u208a": "+", "\u208b": "-", "\u208c": "=", "\u208d": "(", "\u208e": ")"}
_JOIN = ",. "          # held inside a run only when another glyph of the same kind follows


def _script_runs(text):
    """(text, 'sup' | 'sub' | None) segments, with the glyphs mapped to plain characters.

    A separator is kept inside the segment only when the next non-separator character is
    another glyph of the same kind, so an affiliation marker reads as one superscript
    "1,2,3,4" rather than four superscripts with body-size commas between them.
    """
    out, i = [], 0
    while i < len(text):
        tbl = _SUP if text[i] in _SUP else (_SUB if text[i] in _SUB else None)
        if tbl is None:
            j = i
            while j < len(text) and text[j] not in _SUP and text[j] not in _SUB:
                j += 1
            out.append((text[i:j], None)); i = j
            continue
        kind = "sup" if tbl is _SUP else "sub"
        buf, j = [], i
        while j < len(text):
            if text[j] in tbl:
                buf.append(tbl[text[j]]); j += 1
            elif text[j] in _JOIN:
                k = j
                while k < len(text) and text[k] in _JOIN:
                    k += 1
                if k < len(text) and text[k] in tbl:
                    buf.append(text[j:k]); j = k
                else:
                    break
            else:
                break
        out.append(("".join(buf), kind)); i = j
    return [seg for seg in out if seg[0]]

def _rfonts(run):
    rPr = run._element.get_or_add_rPr(); rf = rPr.find(qn('w:rFonts'))
    if rf is None:
        rf = rPr.makeelement(qn('w:rFonts'), {}); rPr.insert(0, rf)
    for a in ('w:ascii', 'w:hAnsi', 'w:cs', 'w:eastAsia'):
        rf.set(qn(a), FONT)

def tab_stop(doc, inches=0.25):
    """Set the document's default tab interval.

    Word's default is 0.5 in, so a tab at the head of a paragraph indents half an inch. The
    value lives in settings.xml as w:defaultTabStop, in twips, 1440 to the inch.
    """
    st = doc.settings.element
    el = st.find(qn("w:defaultTabStop"))
    if el is None:
        el = st.makeelement(qn("w:defaultTabStop"), {}); st.append(el)
    el.set(qn("w:val"), str(int(round(inches * 1440))))


def style(doc):
    for name in ("Normal",):
        st = doc.styles[name]
        st.font.name = FONT; st.font.size = SIZE
        rPr = st.element.get_or_add_rPr(); rf = rPr.find(qn('w:rFonts'))
        if rf is None:
            rf = rPr.makeelement(qn('w:rFonts'), {}); rPr.insert(0, rf)
        for a in ('w:ascii', 'w:hAnsi', 'w:cs', 'w:eastAsia'):
            rf.set(qn(a), FONT)
        pf = st.paragraph_format
        pf.space_before = Pt(0); pf.space_after = Pt(0); pf.line_spacing = 1.0
        # JUSTIFIED, set on the style so anything that does not go through para() inherits it.
        pf.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY

def para(doc, text="", bold=False, italic=False, align=None, size=None):
    """size overrides the 11 pt body, which is what a figure legend uses at 8 pt."""
    p = doc.add_paragraph()
    pf = p.paragraph_format
    pf.space_before = Pt(0); pf.space_after = Pt(0); pf.line_spacing = 1.0
    p.alignment = align if align is not None else WD_ALIGN_PARAGRAPH.JUSTIFY
    # inline **bold**, *italic*, \u00absuperscript\u00bb
    for chunk, b, i, sup in _tokens(text):
        for seg, kind in _script_runs(chunk):
            r = p.add_run(seg)
            r.bold = bold or b; r.italic = italic or i
            if sup or kind == "sup":
                r.font.superscript = True
            elif kind == "sub":
                r.font.subscript = True
            r.font.name = FONT; r.font.size = size or SIZE; _rfonts(r)
    return p

def _tokens(t):
    """(text, bold, italic, superscript). \u00ab1,2\u00bb marks a superscript citation numeral."""
    out, i = [], 0
    for m in re.finditer(r"\*\*(.+?)\*\*|\*(.+?)\*|\u00ab(.+?)\u00bb", t, re.S):
        if m.start() > i: out.append((t[i:m.start()], False, False, False))
        if m.group(1) is not None: out.append((m.group(1), True, False, False))
        elif m.group(2) is not None: out.append((m.group(2), False, True, False))
        else: out.append((m.group(3), False, False, True))
        i = m.end()
    if i < len(t): out.append((t[i:], False, False, False))
    return out or [(t, False, False, False)]

def blank(doc): para(doc, "")

def table(doc, rows, header=True):
    t = doc.add_table(rows=0, cols=len(rows[0])); t.style = "Table Grid"
    for ri, row in enumerate(rows):
        cells = t.add_row().cells
        for ci, val in enumerate(row):
            cell = cells[ci]; cell.text = ""
            p = cell.paragraphs[0]
            pf = p.paragraph_format
            pf.space_before = Pt(0); pf.space_after = Pt(0); pf.line_spacing = 1.0
            for chunk, b, i, sup in _tokens(str(val)):
                for seg, kind in _script_runs(chunk):
                    r = p.add_run(seg)
                    r.bold = b or (header and ri == 0); r.italic = i
                    if sup or kind == "sup":
                        r.font.superscript = True
                    elif kind == "sub":
                        r.font.subscript = True
                    r.font.name = FONT; r.font.size = SIZE; _rfonts(r)
    return t

def figure(doc, path, width_in=6.5, align=None):
    """LEFT by default, not centred.

    Figures are placed at a width proportional to their own finalized width, so a narrow
    figure stays narrow. Centring it would split the leftover into two thin margins; left
    aligned, the space it does not need falls on the right in one piece, which is where the
    composite puts it too.
    """
    p = doc.add_paragraph()
    p.alignment = align if align is not None else WD_ALIGN_PARAGRAPH.LEFT
    pf = p.paragraph_format; pf.space_before = Pt(0); pf.space_after = Pt(0); pf.line_spacing = 1.0
    p.add_run().add_picture(path, width=Inches(width_in))
