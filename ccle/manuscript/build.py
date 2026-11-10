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
order of first appearance, so a section can move without renumbering by hand.
"""
import os, re, sys, datetime, importlib.util
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
# mk_section.py lived only in the sandbox home until 3 November 2026, outside the project
# and outside version control, so a reclaimed sandbox would have taken the build with it.
# It now sits beside this file and ships with the deposit. Home is kept as a fallback.
_mk = [p for p in (os.path.join(HERE, "mk_section.py"),
                   os.path.expanduser("~/mk_section.py")) if os.path.exists(p)]
if not _mk:
    raise SystemExit("mk_section.py not found beside build.py or in the home directory")
exec(open(_mk[0], encoding="utf-8").read())

import refs as REFMOD, sections as S, si as SI, spine as SP, panels as PN, dedup as DD
spec = importlib.util.spec_from_file_location("drafted", os.path.join(HERE, "drafted.py"))
D = importlib.util.module_from_spec(spec); spec.loader.exec_module(D)

def _cap(c):
    """Caption blocks. Some are (text, is_title) tuples whose first entry repeats the figure
    title the build writes itself; some are plain strings with no title entry."""
    return [t for t, _ in c][1:] if c and isinstance(c[0], tuple) else list(c)

# Each section's text, with the repeated-text pass applied (dedup.py) unless spine.py already
# carries a condensed version of it, which supersedes both.
def _cap1(x):
    return x[0] if isinstance(x, tuple) else x

def _only_cut(x, figno, legend):
    t = x[0] if isinstance(x, tuple) else x
    ls = PN.entry_letters(t) if legend else PN.body_letters(t, figno)
    return bool(ls) and all(l in PN.CUT.get(figno, []) for l in ls)

for _i, _sec in enumerate(S.RESULTS, 1):
    _f = _sec["figno"]
    _ob = DD.apply("P_2_%d" % _i, getattr(D, "P_2_%d" % _i))
    _oc = _cap(DD.apply("CAP_2_%d" % _i, getattr(D, "CAP_2_%d" % _i)))
    # A spine entry may cover only the panels that stay. The cut panels' original text is
    # appended so Extended Results carries it; it leaves the main text either way, so it costs
    # nothing against the word limit and nothing is lost from the paper.
    if _sec["n"] in SP.SPINE:
        _sec["body"] = list(SP.SPINE[_sec["n"]]) + [x for x in _ob if _only_cut(x, _f, False)]
        _sec["extended"] = list(SP.EXTENDED.get(_sec["n"], []))
    else:
        _sec["body"] = _ob
    if _f in SP.LEGEND:
        # A restored panel has no condensed entry and is no longer "only cut", so without this
        # it would reach the figure with no legend at all. Take the condensed entries, add the
        # ORIGINAL entry for any kept panel they do not cover, add the cut ones for Extended
        # Results, and order the result by panel letter with the panel-less entries held at the
        # ends where the scope note and the closing note belong.
        _have = {l for e in SP.LEGEND[_f] for l in PN.entry_letters(_cap1(e))}
        _restored = [x for x in _oc
                     if PN.entry_letters(_cap1(x))
                     and not set(PN.entry_letters(_cap1(x))) & set(PN.CUT.get(_f, []))
                     and not set(PN.entry_letters(_cap1(x))) & _have]
        _all = list(SP.LEGEND[_f]) + _restored + [x for x in _oc if _only_cut(x, _f, True)]
        def _key(x):
            ls = PN.entry_letters(_cap1(x))
            return (1, ls[0]) if ls else (0, "")
        _lead = [x for x in _all if not PN.entry_letters(_cap1(x))]
        _pan = sorted([x for x in _all if PN.entry_letters(_cap1(x))], key=_key)
        _sec["cap"] = (_lead[:1] + _pan + _lead[1:]) if _lead else _pan
    else:
        _sec["cap"] = _oc

# ---------------------------------------------------------------- seven display items
# Cell Chem Biol: "no more than 7 main figures and/or tables". The paper has 8, so Figures 5 and
# 6 merge. They are the same claim stated twice and the 148-word Summary already states them as
# one clause. Nothing is demoted to the supplement.
def _merge_items(a, b):
    sh_a, pan_a = PN.by_panel(a["body"], a["figno"], False)
    sh_b, pan_b = PN.by_panel(b["body"], b["figno"], False)
    cap_sh_a, cap_a = PN.by_panel(a["cap"], a["figno"], True)
    cap_sh_b, cap_b = PN.by_panel(b["cap"], b["figno"], True)
    body, cap, seen = [], [], set()
    def _t(x): return x[0] if isinstance(x, tuple) else x
    body += [_t(x) for x in sh_b[:1]]     # the opening question, from the second section
    for new, src, old in PN.MERGED_ORDER:
        src_pan = pan_a if src == a["figno"] else pan_b
        for x in src_pan.get(old, []):
            if id(x) not in seen:
                seen.add(id(x))
                # Citations stay in ORIGINAL numbering here. The one global fix() pass maps them
                # through PANEL_MAP, and rewriting them here as well remapped every merged-section
                # citation twice: 6a became 5a and then 5c.
                body.append(_t(x))
    body += [_t(x) for x in sh_a[-1:]]    # one closing bound, from the first section
    cap += [_t(x) for x in cap_sh_a[:1]]
    for new, src, old in PN.MERGED_ORDER:
        src_cap = cap_a if src == a["figno"] else cap_b
        for x in src_cap.get(old, []):
            if id(x) not in seen:
                seen.add(id(x))
                cap.append(PN.reletter_entry(_t(x), {old: new}))
    cap += cap_sh_a[-1:]
    # Everything not placed on the merged figure has to be carried out explicitly, or it is lost.
    kept = {(sf, old) for _, sf, old in PN.MERGED_ORDER}
    mv_b, mv_c = [], []
    for item, pans, caps in ((a, pan_a, cap_a), (b, pan_b, cap_b)):
        for l, xs in pans.items():
            if (item["figno"], l) not in kept:
                for x in xs:
                    if id(x) not in seen:
                        seen.add(id(x)); mv_b.append(x)
        for l, xs in caps.items():
            if (item["figno"], l) not in kept:
                for x in xs:
                    if id(x) not in seen:
                        seen.add(id(x)); mv_c.append(x)
    for x in list(sh_b[1:]) + list(sh_a[:-1]):
        if id(x) not in seen:
            seen.add(id(x)); mv_b.append(x)
    for x in list(cap_sh_b) + list(cap_sh_a[1:-1]):
        if id(x) not in seen:
            seen.add(id(x)); mv_c.append(x)
    n_in = len(a["body"]) + len(b["body"]) + len(a["cap"]) + len(b["cap"])
    n_out = len(body) + len(cap) + len(mv_b) + len(mv_c) - len(a.get("extended", [])) \
            - len(b.get("extended", []))
    assert n_in == n_out, "merge lost %d of %d items" % (n_in - n_out, n_in)
    return dict(n=a["n"], figno=a["figno"], fig=a["fig"], status=a["status"],
                heading=PN.MERGED_HEADING, body=body, cap=cap, merged=True,
                moved_body=mv_b + a.get("extended", []) + b.get("extended", []),
                moved_cap=mv_c)

MAIN = []
for _sec in S.RESULTS:
    if _sec["figno"] == PN.MERGE_FROM:
        continue
    if _sec["figno"] == PN.MERGE_INTO:
        MAIN.append(_merge_items(_sec, [r for r in S.RESULTS if r["figno"] == PN.MERGE_FROM][0]))
    else:
        MAIN.append(dict(_sec, merged=False))
# renumber: old 7 becomes 6, old 8 becomes 7
for _k, _sec in enumerate(MAIN, 1):
    _sec["src_figno"] = _sec["figno"]
    _sec["figno"] = _k
    _sec["n"] = "2.%d" % _k

DATE = datetime.date.today().isoformat()
# --preprint strips the front matter that exists for the journal and for us, not for a
# reader: the venue line (which records that the paper was retargeted off a special issue),
# the drafting note about square brackets, the Highlights block with its submission-page
# instructions and per-highlight character counts, and the graphical-abstract spec. A
# preprint is permanent, so none of that should be posted.
# The default build is the SUBMISSION document. --draft adds the working apparatus.
#
# This was the other way round until 3 November, and the naming hid what it meant: the flag
# was called --preprint, so the unflagged build looked like "the journal version" when it was
# the reading copy. Its first paragraph read "Retargeted 9 October 2026 from the Bioorganic &
# Medicinal Chemistry Tetrahedron Prize special issue honoring Benjamin F. Cravatt; that
# commitment still needs withdrawing." Thirteen paragraphs of that kind, notes to ourselves
# about titles to avoid, the section log, what the submission page says, and two reading-copy
# markers, all sat in the file named Manuscript_main_text.docx. None of it belongs in anything
# that leaves this folder, for either venue, so the clean document is now what you get by
# default and the annotated one has to be asked for by name.
DRAFT = "--draft" in sys.argv

# Word is authoritative for submission; PDFs are made once at finalization. results/ keeps
# the dated build history, CCB-submission/ holds the current file under a stable name, so
# nothing downstream has to guess which date is current.
SUBMIT = os.path.abspath(os.path.join(HERE, "..", "..", "CCB-submission"))


def _submit(path, stable):
    """Copy a built document into CCB-submission/ under a name that does not move.

    A draft run copies nothing: the submission folder holds only what can be uploaded.
    """
    if DRAFT or not os.path.isdir(SUBMIT):
        return
    import shutil
    shutil.copy2(path, os.path.join(SUBMIT, stable))
    print("  submission  %s" % stable)


AUTHOR = "Mikail E. Abbasov"


def stamp(path, title=None):
    """Put him in the file properties of a figure, as the documents already carry him.

    A figure left as matplotlib wrote it has an empty Author and names the library in
    Creator. Producer is left alone: it records which library wrote the bytes, which is a
    technical field and not a claim of authorship.
    """
    try:
        if path.lower().endswith(".pdf"):
            import pymupdf
            d = pymupdf.open(path)
            md = d.metadata or {}
            md.update(author=AUTHOR, title=title or os.path.basename(path),
                      subject=S.FRONT["title"], creator=AUTHOR,
                      keywords=", ".join(S.KEYWORDS["body"]))
            d.set_metadata(md); d.saveIncr(); d.close()
        elif path.lower().endswith(".png"):
            from PIL import Image
            from PIL.PngImagePlugin import PngInfo
            im = Image.open(path)
            dpi = im.info.get("dpi", (300, 300))
            meta = PngInfo()
            meta.add_text("Author", AUTHOR)
            meta.add_text("Title", title or os.path.basename(path))
            meta.add_text("Description", S.FRONT["title"])
            meta.add_text("Copyright", AUTHOR)
            im.save(path, dpi=dpi, pnginfo=meta)
        return True
    except Exception as e:
        print("  stamp failed on %s: %s" % (os.path.basename(path), e))
        return False


def _submit_figures():
    """Mirror the rendered figures into the submission folder as well.

    Cell Press takes figures as separate files, and anyone opening CCB-submission/ should not
    have to know where the renders live. PDFs only: that is the vector form a journal wants,
    and it keeps the whole set under 2 MB, where the PNGs would add 19 MB of review copies
    that nothing downstream needs.
    """
    if DRAFT or not os.path.isdir(SUBMIT):
        return
    import shutil
    dst = os.path.join(SUBMIT, "figures")
    os.makedirs(dst, exist_ok=True)
    n = 0
    # MAIN FIGURES GO UNDER THE PUBLISHED NUMBER, from the vector master each script writes.
    # This folder carried eight files called fig1 to fig8, which is the source numbering and
    # not what the paper cites. The map is FIG_MAP, built from panels.py, and the merged
    # Figure 5 has its own script because it is not a successor to either source.
    _vec = {v: "fig%d.pdf" % k for k, v in sorted(FIG_MAP.items())}
    _vec[[v for k, v in FIG_MAP.items() if k == PN.MERGE_INTO][0]] = "fig5pub.pdf"
    for _pub, _src in sorted(_vec.items()):
        _p = os.path.join(FIGDIR, _src)
        if not os.path.exists(_p):
            print("  submission  MISSING vector master %s for Figure %d" % (_src, _pub))
            continue
        shutil.copy2(_p, os.path.join(dst, "Figure%d.pdf" % _pub))
        n += 1
    # The supplementary set keeps its own S numbering, which the merge does not touch. It is
    # taken from the LIVE renders in FIGDIR, not from figures/si: that mirror was last written
    # on 4 November and every one of its ten files was older than the render beside it, so the
    # submission folder was shipping supplementary artwork a week behind the scripts.
    for f in sorted(os.listdir(FIGDIR)):
        if re.fullmatch(r"figS\d+\.pdf", f):
            shutil.copy2(os.path.join(FIGDIR, f), os.path.join(dst, f))
            n += 1
    # EVERY FILE IN THE PACKAGE CARRIES HIM, not just the documents.
    # The Title property is the figure's own title, not its filename.
    _titles = {}
    for _s in MAIN:
        _titles["Figure%d" % _s["figno"]] = us("Figure %d. %s." % (_s["figno"],
                                                                  _s["heading"].rstrip(".")))
    for _f in SI.SI_FIGS:
        _titles["fig%s" % _f["no"]] = fix("Figure %s. %s." % (_f["no"], _f["title"].rstrip(".")))
    _st = 0
    for f in sorted(os.listdir(dst)):
        if f.lower().endswith((".pdf", ".png")):
            _st += stamp(os.path.join(dst, f), _titles.get(os.path.splitext(f)[0]))
    _fin = os.path.join(SUBMIT, "figures-final")
    if os.path.isdir(_fin):
        for f in sorted(os.listdir(_fin)):
            if f.lower().endswith(".png"):
                _st += stamp(os.path.join(_fin, f),
                             title=_titles.get(re.sub(r"\.png$", "", f)))
    print("  submission  figures/ (%d PDFs, main figures under the published number), "
          "%d files stamped with the author" % (n, _st))
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
    """Highlights carry their own character count against the journal's 85.

    DRAFT only. Measuring it here is what keeps the number from drifting from the text it
    describes, but the bracket is ours and a submitted bullet ends at its last word.
    """
    return ["%s  [%d]" % (t, len(t)) for t in items] if DRAFT else list(items)

def cite(text):
    text = us(text)
    text = re.sub(r"\[\[([a-z0-9,]+)\]\]",
                  lambda m: _fmt([_num(k) for k in m.group(1).split(",")]), text)
    # A superscript numeral sits against the word it cites. The source writes " [[key]]" with
    # a space, which set the numeral a word-space away from its own word.
    return re.sub(r"[ \t]+(?=\u00ab)", "", text)

# ---------------------------------------------------------------- shrink figures for the draft
def print_copy(name):
    """A 2000 px-wide copy, about 300 dpi at 6.5 in, so the draft file stays manageable."""
    os.makedirs(PRINT, exist_ok=True)
    src, dst = os.path.join(FIGDIR, name), os.path.join(PRINT, name)
    if not os.path.exists(dst) or os.path.getmtime(dst) < os.path.getmtime(src):
        from PIL import Image
        im = Image.open(src); w, h = im.size
        im.resize((2000, int(2000 * h / w)), Image.LANCZOS).save(dst, optimize=True)
    return dst

def print_copy_src(src, name):
    """print_copy for a file that does not live in FIGDIR, cached under `name`.

    The finalized sheets live in CCB-submission/figures-final and are named by the PUBLISHED
    number, while FIGDIR's renders are named by the source number. Caching under the
    published name is what keeps the two from overwriting each other in _build/print.
    """
    os.makedirs(PRINT, exist_ok=True)
    dst = os.path.join(PRINT, name)
    if not os.path.exists(dst) or os.path.getmtime(dst) < os.path.getmtime(src):
        from PIL import Image
        im = Image.open(src); w, h = im.size
        im.resize((2000, int(2000 * h / w)), Image.LANCZOS).save(dst, optimize=True)
    return dst


CAPSIZE = Pt(8)                 # figure legends, against the 11 pt body


def legend_entry(e):
    """One legend entry, in the published form.

    The legends had grown three kinds of emphasis: a bold panel letter with its whole opening
    clause bold behind it, a bold warning sentence in the middle of an entry, and a closing
    caveat in italic. Three figures had none of it and two had all three. A published legend
    carries one convention, so every mark is stripped and the panel letter alone is set bold.
    """
    t = re.sub(r"\*\*(.+?)\*\*", r"\1", e, flags=re.S)
    t = re.sub(r"(?<!\*)\*(?!\*)(.+?)(?<!\*)\*(?!\*)", r"\1", t, flags=re.S)
    return re.sub(r"^(\([A-Za-z](?:\s*[,\u2013\u2014-]\s*[A-Za-z])*\))", r"**\1**", t.strip())


def native_in(path):
    """A finalized sheet's own width in inches. Every one is written at 300 dpi."""
    from PIL import Image
    return Image.open(path).size[0] / 300.0


def placer(paths, text_width=6.5):
    """The composite's rule, applied to a document.

    The composite brings every sheet to the size it will be FINALIZED at and then reduces
    them all by the SAME factor, so what the eye sees is the true relative width and the true
    relative type size. Placing each figure at the text width instead scales a narrow figure
    harder and prints its type larger, which is the one thing worth watching. Here the widest
    figure in the set fills the text width and every other one is placed at its own share of
    it, so a narrow sheet leaves its unused space on the right.
    """
    wide = max(native_in(p) for p in paths)
    k = text_width / wide
    return (lambda p: native_in(p) * k), k, wide


def newdoc():
    doc = Document()
    # House rule: the file properties carry him, not the generator's default.
    cp = doc.core_properties
    cp.author = cp.last_modified_by = "Mikail E. Abbasov"
    cp.title = S.FRONT["title"]
    for s in doc.sections:
        s.top_margin = s.bottom_margin = Inches(1); s.left_margin = s.right_margin = Inches(1)
    style(doc); return doc

AFFSIZE = Pt(8)                 # the affiliation block, against the 11 pt body


def affiliations(doc):
    """One department per line, at 8 pt, as a journal sets them.

    The four were one 11 pt paragraph running together, so the reader had to find the
    superscript markers inside a block of prose to tell the departments apart.
    """
    for a in re.split(r"(?=[\u00b9\u00b2\u00b3\u2074-\u2079])", S.FRONT["affiliations"]):
        if a.strip():
            para(doc, a.strip(), italic=True, size=AFFSIZE)
    blank(doc)


def block(doc, paras, italic=False):
    """A run of paragraphs with NO blank line between them.

    A blank line between every paragraph is how a working draft is set, not a manuscript.
    Blank lines are kept where they mark a boundary: between sections, and around a figure
    and its legend. Inside a section the paragraphs run on, and one blank closes the block.
    """
    for t in paras:
        # A TAB OPENS EACH ONE. With the blank lines gone, the indent is what marks where one
        # paragraph ends and the next begins.
        para(doc, "\t" + cite(t), italic=italic)
    if paras:
        blank(doc)

def todo(doc):
    para(doc, "[ NOT YET DRAFTED ]", bold=True); blank(doc)


# ---------------------------------------------------------------- US spelling
# Cell Press: "Spelling and use of periods and commas in numbers should conform to US usage."
# Applied at emission and counted, so the number is reportable rather than silent.
US = [("labelling","labeling"),("Labelling","Labeling"),("labelled","labeled"),
      ("Labelled","Labeled"),("tumours","tumors"),("tumour","tumor"),("Tumour","Tumor"),
      ("behaviours","behaviors"),("behaviour","behavior"),("modelled","modeled"),
      ("modelling","modeling"),("normalised","normalized"),("normalise","normalize"),
      ("ionised","ionized"),("oxidised","oxidized"),("localised","localized"),
      ("recognised","recognized"),("standardised","standardized"),("centred","centered"),
      ("artefact","artifact"),("fibre","fiber"),("defence","defense"),
      ("synthesised","synthesized"),
      # Longest form FIRST: plain str.replace turned "catalogued" into "catalogd" three
      # times in the manuscript before this was ordered. Capitalised "Catalogue" is left
      # alone on purpose, since it is part of COSMIC's own name.
      ("catalogued","cataloged"),("catalogue","catalog"),("neighbour","neighbor"),
      ("realised","realized"),("generalised","generalized"),("centre","center"),
      ("oesophagus","esophagus"),("analysed","analyzed"),("summarised","summarized")]
US_HITS = {}
def us(text):
    for a, b in US:
        if a in text:
            US_HITS[a] = US_HITS.get(a, 0) + text.count(a)
            text = text.replace(a, b)
    return text


# ---------------------------------------------------------------- legend arithmetic
# Cell Chemical Biology requires a legend to say what is shown, name the test, give the unit and
# the replicate structure, and define every measure and error bar. It does not require the legend
# to restate every value the panel plots, and with 40 panels that restatement is what puts the
# paper over the limit. Pure enumerations of values move to a supplemental table keyed by figure
# and panel; anything that defines how a quantity was built, or that bounds a claim, stays.
ARITH_ROWS = []
# Audit trail for the splitter: every entry it touched, what it kept, what it moved, and
# whether the cut fell inside a semicolon-delimited list. Written to panel_data so the
# relocation is reviewable instead of invisible.
ARITH_AUDIT = []
_KEEP = ("interval", "percentile", "bootstrap", "resampl", "shuffl", "permutation", "standardi",
         "defined", "definition", "replicat", "n =", "error bar", "scale bar", "see also",
         "mean divided", "unit is", "denominator", "not interchangeable", "test", "by construction",
         "lower bound", "upper bound", "rather than", "cannot", "is not", "so the", "because",
         # A qualifier is a bound, not a value. Table S2 states in print that nothing bounding a
         # claim was moved, and these are the forms that promise turned out not to cover: the
         # single exception, the structures that are not crystallographic, the one cell line that
         # carries the mutation, and which probe labelled what.
         "exception", "the only", "only by", "wild type", "cryo-em", "rrid", "accession")

# A plotted value is a bare number. An accession is not: 2SHP, 4DGP, 8GUA and CVCL_1344 each
# used to count as a value, which is how Figure 7d's whole PDB list, the RRID and the cryo-EM
# correction came to be classified as arithmetic and moved out of the legend. Require the digit
# run to stand free of letters and underscores, and strip the 10 of a x10^-n exponent, which is
# one value with its mantissa and not a second one.
_NUM = re.compile(r"(?<![A-Za-z0-9_])[-+\u2212]?\d[\d,.]*(?![A-Za-z0-9_])")
_EXP = re.compile(r"10(?=\u207b)")

def _values(s):
    return _NUM.findall(_EXP.sub("", s))

def split_arithmetic(entry, figno):
    """-> (legend text kept, [sentences moved]). The header sentence always stays.

    Splits on sentence boundaries only. A semicolon-delimited list is one unit and moves whole
    or stays whole: splitting on ';' cut twelve legends mid-list, leaving Figure 1c naming six
    allele-frequency bands and printing four, and Figure 6b saying two features differ and
    naming one."""
    parts = re.split(r"(?<=\.)\s+(?=[A-Z(])", entry)
    if not parts:
        return entry, []
    keep, moved = [parts[0]], []
    for s in parts[1:]:
        if len(_values(s)) >= 3 and not any(k in s.lower() for k in _KEEP):
            moved.append(s.strip())
        else:
            keep.append(s)
    if not moved:
        return entry, []
    ls = PN.entry_letters(entry)
    ARITH_ROWS.append(("Figure %d" % figno, ", ".join(ls) if ls else "general",
                       " ".join(moved)))
    _kept = " ".join(keep).rstrip() + " Values in Table S2."
    _idx = [i for i, s_ in enumerate(parts) if s_.strip() in moved]
    ARITH_AUDIT.append({"figno": figno,
                        "panels": ls,
                        "mid_list": any(parts[i].rstrip().endswith(";") for i in _idx)
                                 or any(parts[i - 1].rstrip().endswith(";") for i in _idx if i > 0),
                        "entry": entry, "kept": _kept, "moved": moved})
    return _kept, moved

# ---------------------------------------------------------------- the global figure/panel map
# The 5+6 merge renumbers old 7 to 6 and old 8 to 7, and panel cuts reletter what stays. Every
# citation anywhere in the paper, the SI included, has to follow. One map, one pass: successive
# replacement would cascade (c to b, then b to a).
def _letters_of(sec, legend):
    out = []
    for x in (sec["cap"] if legend else sec["body"]):
        t = x[0] if isinstance(x, tuple) else x
        out += PN.entry_letters(t) if legend else PN.body_letters(t, sec["src_figno"])
    return out

PANEL_MAP, FIG_MAP, CUT_OUT = {}, {}, set()
for sec in MAIN:
    srcs = [sec["src_figno"]] + ([PN.MERGE_FROM] if sec.get("merged") else [])
    FIG_MAP.update({sf: sec["figno"] for sf in srcs})
    if sec.get("merged"):
        for new, sf, old in PN.MERGED_ORDER:
            PANEL_MAP[(sf, old)] = (sec["figno"], new)
        for sf in srcs:
            for p in PN.CUT.get(sf, []):
                CUT_OUT.add((sf, p))
    else:
        sf = sec["src_figno"]
        allp = sorted({l for lg in (False, True) for l in _letters_of(sec, lg)})
        kept = [l for l in allp if l not in PN.CUT.get(sf, [])]
        for k, old in enumerate(kept):
            PANEL_MAP[(sf, old)] = (sec["figno"], chr(ord("a") + k))
        for p in PN.CUT.get(sf, []):
            CUT_OUT.add((sf, p))

# THE SUPPLEMENTARY NOTES ARE RENUMBERED, because the cut and withdrawn ones are not in the
# submitted document and a document that skips N6 and N9 reads as one with pieces missing.
# Three sentences cite a note by number and have to follow the renumber.
_SI_KEEP = [e for e in SI.SI_NOTES if e[1] in ("keep", "rewrite") and len(e) > 3]
NOTE_MAP = {e[0]: "N%d" % (i + 1) for i, e in enumerate(_SI_KEEP)}


def _note_ref(m):
    return m.group(1) + NOTE_MAP.get(m.group(2), m.group(2))


def _ref(m):
    """Rewrite one 'Figure N' or 'Figure Nx' citation.

    CASE. The paper cites panels as "Figure 8A" and this pattern read [a-z] only, so the
    match failed outright on the letter and the citation was left untouched: five "Figure 8"
    citations survived into a paper with seven figures, and every uppercase citation to
    Figures 5 through 8 kept its pre-merge number. Letters are looked up lowercase and
    written back in the case they were found in.
    """
    n, p = int(m.group(1)), m.group(2)
    if p is None:
        return "%s %d" % (m.group(0).split()[0], FIG_MAP.get(n, n))
    _lp = p.lower()
    if (n, _lp) in PANEL_MAP:
        nf, np_ = PANEL_MAP[(n, _lp)]
        return "%s %d%s" % (m.group(0).split()[0], nf,
                            np_.upper() if p.isupper() else np_)
    if (n, _lp) in CUT_OUT:
        return "Extended Results"
    # Neither mapped nor cut. Returning the match unchanged kept a stale figure number
    # alive through a renumber, which is how "Figure 8d" outlived Figure 8. The letter is
    # unknown, the number is not, so correct what is known.
    return "%s %d%s" % (m.group(0).split()[0], FIG_MAP.get(n, n), p)

def _range(m):
    a, b = int(m.group(2)), int(m.group(4))
    return "%s %d%s%d" % (m.group(1), FIG_MAP.get(a, a), m.group(3), FIG_MAP.get(b, b))

def fix(text):
    """Apply the renumber, the reletter and the US spelling to one string."""
    text = re.sub(r"(Figures?)\s+(\d+)(\s*(?:to|and|[\u2013\u2014-])\s*)(\d+)", _range, text)
    text = re.sub(r"Figures?\s+(\d+)([A-Za-z])?(?![a-zA-Z0-9])", _ref, text)
    text = re.sub(r"([Nn]otes?\s+(?:and\s+)?)(N\d+)\b", _note_ref, text)
    return us(text)

def fix_entry(x, rmap_legend=False):
    t = x[0] if isinstance(x, tuple) else x
    if rmap_legend:
        t = PN.reletter_entry(t, {old: new for (f, old), (nf, new) in PANEL_MAP.items()
                                  if f in _SRC_NOW})
    return fix(t)

def dropped(x, sec, legend=False):
    """True when every panel this item is attributable to has been cut."""
    t = x[0] if isinstance(x, tuple) else x
    srcs = [sec["src_figno"]] + ([PN.MERGE_FROM] if sec.get("merged") else [])
    ls = []
    for sf in srcs:
        got = PN.entry_letters(t) if legend else PN.body_letters(t, sf)
        ls += [(sf, l) for l in got]
    return bool(ls) and all(x_ in CUT_OUT for x_ in ls)

# ================================================================ MANUSCRIPT
_FINAL_DIR = os.path.join(SUBMIT, "figures-final")
_FINALS = [os.path.join(_FINAL_DIR, f) for f in sorted(os.listdir(_FINAL_DIR))
           if re.fullmatch(r"Figure\d+\.png", f)] if os.path.isdir(_FINAL_DIR) else []
if _FINALS:
    PLACE, _K_MAIN, _WIDE_MAIN = placer(_FINALS)
else:
    PLACE, _K_MAIN, _WIDE_MAIN = (lambda _p: 6.5), 1.0, 6.5

doc = newdoc()
para(doc, S.FRONT["title"], bold=True); blank(doc)
para(doc, S.FRONT["authors"]); blank(doc)
affiliations(doc)
# NO RUNNING TITLE. Cell Chemical Biology's components list does not ask for one, and a
# line of ours in a submitted manuscript is a line an editor has to decide what to do with.
# The string is kept in sections.py for the cover letter, where a short title is useful.
if DRAFT:
    para(doc, us(S.FRONT["venue"]), italic=True); blank(doc)
    para(doc, "Draft assembled %s. %s" % (DATE, S.FRONT["title_note"]), italic=True); blank(doc)
    # DRAFT only. The submitted file carries no build date and no note about its own title.
    para(doc, "Every section carries prose, figure, caption and citations; the record of what was "
              "opened and closed while drafting them is in the section log, not here. Square brackets "
              "mark the few items still to be supplied.", italic=True); blank(doc)

    para(doc, "Highlights", bold=True); blank(doc)
    para(doc, S.HIGHLIGHTS["spec"], italic=True); blank(doc)
    for t in _hl_measured(S.HIGHLIGHTS["body"]):
        para(doc, "\u2022  " + t); blank(doc)

para(doc, "Graphical abstract", bold=True); blank(doc)
if S.GRAPHICAL_ABSTRACT.get("status") == "drafted" and S.GRAPHICAL_ABSTRACT.get("fig"):
    figure(doc, print_copy(S.GRAPHICAL_ABSTRACT["fig"]),
           PLACE(os.path.join(FIGDIR, S.GRAPHICAL_ABSTRACT["fig"]))); blank(doc)
else:
    todo(doc)
if DRAFT:
    para(doc, S.GRAPHICAL_ABSTRACT["spec"], italic=True); blank(doc)

# Cell Press calls this the Summary, not the Abstract, and caps it at 150 words.
para(doc, "Summary", bold=True); blank(doc)
block(doc, D.P_ABSTRACT)

para(doc, "Keywords", bold=True); blank(doc)
# Cell Press: "up to 10 keywords with your paper, separated by commas".
para(doc, us(", ".join(S.KEYWORDS["body"]))); blank(doc)
if DRAFT:
    para(doc, S.KEYWORDS["note"], italic=True); blank(doc)

para(doc, "1.  Introduction", bold=True); blank(doc)
block(doc, [fix(p) for p in (getattr(SP, "INTRO", None) or D.P_INTRO)])
para(doc, "2.  Results", bold=True); blank(doc)
# The real failure mode: a KEPT panel with no spine text aimed at it. That happens when a
# condensed entry is labelled with the wrong letter, and the build then drops it silently, which
# has occurred three times while condensing. Cut-panel text in a spine entry is fine: the
# additive rule routes it to Extended Results.
SPINE_WARN = []
for _sec in S.RESULTS:
    if _sec["n"] not in SP.SPINE and _sec["figno"] not in SP.LEGEND:
        continue
    _f = _sec["figno"]
    _orig = _cap(DD.apply("CAP_2_%d" % (S.RESULTS.index(_sec) + 1),
                          getattr(D, "CAP_2_%d" % (S.RESULTS.index(_sec) + 1))))
    _all = sorted({l for e in _orig for l in PN.entry_letters(e)})
    _kept = [l for l in _all if l not in PN.CUT.get(_f, [])]
    _have_leg = {l for e in _sec["cap"] for l in PN.entry_letters(_cap1(e))}
    _have_bod = {l for t in _sec["body"] for l in PN.body_letters(_cap1(t), _f)}
    for _l in _kept:
        if _l not in _have_leg:
            SPINE_WARN.append("Figure %d panel %s is kept but has no legend entry" % (_f, _l))
        if _l not in _have_bod:
            SPINE_WARN.append("Figure %d panel %s is kept but no body text cites it"
                              % (_f, _l))

MOVED, EMIT = {}, {}
for sec in MAIN:
    _SRC_NOW = set([sec["src_figno"]] + ([PN.MERGE_FROM] if sec.get("merged") else []))
    if sec.get("merged"):
        # The merge already partitioned kept from moved, and relettered what it kept. Running
        # dropped() over it again would misjudge: the new letter f collides with the old 5f and
        # 6f that were cut, so the panel would be thrown away twice over.
        mb, ml = sec["body"], sec["cap"]
        movb, movl = sec["moved_body"], sec["moved_cap"]
    else:
        mb  = [p for p in sec["body"] if not dropped(p, sec)]
        ml  = [e for e in sec["cap"]  if not dropped(e, sec, legend=True)]
        movb = [p for p in sec["body"] if dropped(p, sec)] + sec.get("extended", [])
        movl = [e for e in sec["cap"]  if dropped(e, sec, legend=True)]
    rl = {old: new for (f, old), (nf, new) in PANEL_MAP.items() if f in _SRC_NOW}
    eb = [fix(p[0] if isinstance(p, tuple) else p) for p in mb]
    el = [fix(e[0] if isinstance(e, tuple) else e) if sec.get("merged")
          else fix(PN.reletter_entry(e[0] if isinstance(e, tuple) else e, rl)) for e in ml]
    el = [split_arithmetic(e, sec["figno"])[0] for e in el]
    MOVED[sec["figno"]] = (sec, movb, movl); EMIT[sec["figno"]] = (eb, el)
    para(doc, us("%s  %s" % (sec["n"], sec["heading"])), bold=True); blank(doc)
    block(doc, eb)
    # A composite is needed only where the published figure is not one rendered figure: where
    # panels were cut, or where two sources merged. That is Figures 1 and 5. For the other five
    # the current render IS the published figure, and preferring the composite "when it exists"
    # stood a 20 October crop in front of it, which is how the restyle of 3 November reached
    # every render and none of the main-text figures in the document.
    _srcs_fig = [sec["src_figno"]] + ([PN.MERGE_FROM] if sec.get("merged") else [])
    _needs_composite = bool(sec.get("merged")) or any(PN.CUT.get(_s) for _s in _srcs_fig)
    _rp = os.path.join(OUT, "_repanel", "fig%d_repanelled.png" % sec["figno"])
    # THE FINALIZED SHEET WINS, and it is named by the PUBLISHED number. Until 9 November
    # this line reached for FIGDIR's source-numbered render, so the document's Figure 6 was
    # source fig7 by luck of the mapping and its Figure 5 was a repanelled reading copy of a
    # merge that had never been rendered. figures-final now holds the seven published sheets
    # and this is where they enter the document.
    _fin = os.path.join(SUBMIT, "figures-final", "Figure%d.png" % sec["figno"])
    if os.path.exists(_fin):
        _art = print_copy_src(_fin, "Figure%d.png" % sec["figno"])
    else:
        _art = _rp if (_needs_composite and os.path.exists(_rp)) else print_copy(sec["fig"])
    figure(doc, _art, PLACE(_fin if os.path.exists(_fin) else _art)); blank(doc)
    # DRAFT MARKER, removed at submission, and now printed only for the two figures that are
    # genuinely composites. The panels are the current renderings but still carry the original
    # layout's proportions, while the legend below is the current specification.
    # Printing the map makes the mismatch impossible to read past, and doubles as the brief for
    # rebuilding the figure.
    _srcs_now = [sec["src_figno"]] + ([PN.MERGE_FROM] if sec.get("merged") else [])
    _map = sorted(((new_, "%d%s" % (sf, old_)) for (sf, old_), (nf, new_) in PANEL_MAP.items()
                   if sf in _srcs_now), key=lambda x: x[0])
    _gone = ", ".join("%d%s" % (sf, p) for sf in _srcs_now for p in PN.CUT.get(sf, []))
    if DRAFT and _art == _rp:
        para(doc, "[READING COPY, not final artwork. The panels above are the current "
                  "renderings, cropped and laid out in this figure's specified order with the "
                  "letters corrected: %s. Moved to Extended Results: %s. Panel proportions are "
                  "inherited from the original layout, and anything drawn on the figure rather than "
                  "on a panel is not inside a panel crop.]"
             % ("; ".join("%s from %s" % (n, o) for n, o in _map), _gone or "none"), italic=True)
    # ONE CONTINUOUS PARAGRAPH. The title was a paragraph of its own and every entry another,
    # with a blank line between each, which is not how a journal sets a legend. Bold run-in
    # title, then the entries running on, at 8 pt.
    para(doc, "**%s** %s" % (us("Figure %d. %s." % (sec["figno"], sec["heading"].rstrip("."))),
                             " ".join(legend_entry(cite(e)) for e in el)), size=CAPSIZE)
    blank(doc)
para(doc, "3.  Discussion", bold=True); blank(doc)
block(doc, [fix(p) for p in (getattr(SP, "DISCUSSION", None) or D.P_DISCUSSION)])
if not getattr(SP, "DISCUSSION", None):
    para(doc, "4.  Conclusion", bold=True); blank(doc)
    block(doc, [fix(p) for p in D.P_CONCLUSION])
para(doc, "Resource availability", bold=True); blank(doc)
block(doc, [fix(p) for p in S.RESOURCE_AVAIL["body"]])
para(doc, "5.  Methods", bold=True); blank(doc)
block(doc, [fix(p) for p in list(D.P_METHODS) + [DD.METHODS_ADDITION]
            + ([SP.METHODS_RELOCATED] if hasattr(SP, "METHODS_RELOCATED") else [])])
for name, st_, spec in S.DECLARATIONS:
    para(doc, name, bold=True); blank(doc)
    if st_ == "todo":
        todo(doc); para(doc, spec, italic=True); blank(doc)
    elif spec is None:                      # Data and code availability carries its own body
        block(doc, S.DATA_AVAIL["body"])
    else:
        block(doc, spec)

para(doc, "References", bold=True); blank(doc)
# The verification note is a standing rule of this lab, not a statement for a reader.
for i, key in enumerate(ORDER, 1):
    para(doc, "%d.  %s" % (i, REFMOD.REFS[key]["text"])); blank(doc)
m_out = os.path.join(OUT, ("Reading_copy_%s.docx" if DRAFT
                           else "Manuscript_%s.docx") % DATE)
doc.save(m_out)
# One clean main text serves both Cell Chemical Biology and bioRxiv: nothing in the apparatus
# was venue-specific, so there is no second version to keep in step.
_submit(m_out, "Manuscript_main_text.docx")
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
affiliations(sdoc)
if DRAFT:
    para(sdoc, "Draft assembled %s." % DATE, italic=True); blank(sdoc)
para(sdoc, "Supplementary figures are numbered sequentially in order of first citation and each "
           "names its parent figure. Figure 2 has two supplements, which is why one supplement per "
           "main figure was not used. Supplementary notes carry N-numbers and never S-numbers.",
     italic=True); blank(sdoc)
para(sdoc, "Contents", bold=True); blank(sdoc)
for f in SI.SI_FIGS:
    para(sdoc, fix("Figure %s, support for %s: %s" % (f["no"], f["parent"], f["title"])))
for t_ in getattr(SI, "SI_TABLES", []):
    para(sdoc, "Table %s. %s" % (t_["no"], t_["title"]))
blank(sdoc)
_WORD = {9:"Nine",10:"Ten",11:"Eleven",12:"Twelve",13:"Thirteen",16:"sixteen",17:"seventeen",18:"eighteen",19:"nineteen"}
_surv = [e for e in SI.SI_NOTES if e[1] not in ("cut", "withdraw")]
if DRAFT:
    para(sdoc, "Supplementary notes follow the figures. %s of the original %s topics survive "
               "the audit against the finalized figures; the rest are carried by a panel or rest "
               "on a withdrawn argument, and are listed with that status rather than silently "
               "dropped."
               % (_WORD.get(len(_surv), len(_surv)), _WORD.get(len(SI.SI_NOTES), len(SI.SI_NOTES))))
    blank(sdoc)
_SI_SRC = [os.path.join(FIGDIR, f["fig"]) for f in SI.SI_FIGS
           if os.path.exists(os.path.join(FIGDIR, f["fig"]))]
SPLACE = placer(_SI_SRC)[0] if _SI_SRC else (lambda _p: 6.5)
for f in SI.SI_FIGS:
    figure(sdoc, print_copy(f["fig"]), SPLACE(os.path.join(FIGDIR, f["fig"]))); blank(sdoc)
    para(sdoc, "**%s** %s" % (fix("Figure %s. Support for %s: %s." %
                                  (f["no"], f["parent"], f["title"].rstrip("."))),
                              " ".join(legend_entry(cite(fix(_c))) for _c in f["cap"])),
         size=CAPSIZE)
    blank(sdoc)
for t_ in getattr(SI, "SI_TABLES", []):
    para(sdoc, fix("Table %s. %s" % (t_["no"], t_["title"])), bold=True); blank(sdoc)
    table(sdoc, [t_["head"]] + [tuple(cite(c) for c in r) for r in t_["rows"]])
    blank(sdoc)
    block(sdoc, [fix(_b) for _b in t_["foot"]])

if ARITH_ROWS:
    para(sdoc, "Table S2. Panel-level values for the main figures.", bold=True); blank(sdoc)
    para(sdoc, "Every value plotted in a main figure panel, moved here so that each main legend "
               "states what is shown, the test, the unit, the replicate structure and the "
               "definition of every measure and error bar without also restating the panel. "
               "Nothing that bounds a claim was moved.", italic=True); blank(sdoc)
    table(sdoc, [("Figure", "Panel", "Values")] + [tuple(cite(c) for c in r) for r in ARITH_ROWS])
    blank(sdoc)

para(sdoc, "Extended Results", bold=True); blank(sdoc)
para(sdoc, "Analyses that support the main-text figures and are reported here in full. Panel "
           "letters are those of the source analysis; the main-text figures carry their own "
           "sequential lettering.", italic=True)
blank(sdoc)
for _f in sorted(MOVED):
    _sec, _mb, _ml = MOVED[_f]
    if not _mb and not _ml:
        continue
    _srcs = [_sec["src_figno"]] + ([PN.MERGE_FROM] if _sec.get("merged") else [])
    _ps = ", ".join("%d%s" % (sf, p) for sf in _srcs for p in PN.CUT.get(sf, []))
    # _ps is empty for a figure that lost no panels: these sections carry text moved out,
    # not panels, and the sentence used to read "panels  of the earlier draft".
    _head = ("Supporting main-text Figure %d, panels %s. %s." % (_f, _ps, _sec["heading"])
             if _ps else "Supporting main-text Figure %d. %s." % (_f, _sec["heading"]))
    para(sdoc, us(_head), bold=True); blank(sdoc)
    # fix(), not us(). The moved paragraphs were getting the US spelling pass and nothing
    # else, so a figure or note citation inside one kept its pre-merge number: "(note N18)"
    # survived into a document whose notes stop at N11. fix() is us() plus both renumbers.
    block(sdoc, [fix(p[0] if isinstance(p, tuple) else p) for p in _mb])
    if _ml:
        para(sdoc, "Legend, moved panels.", bold=True); blank(sdoc)
        block(sdoc, [fix(e[0] if isinstance(e, tuple) else e) for e in _ml], italic=True)

if getattr(SP, "EXTENDED_BACK", None):
    para(sdoc, "From the Discussion: the prospective screen this paper does not run.",
         bold=True); blank(sdoc)
    block(sdoc, [fix(p) for p in SP.EXTENDED_BACK])

para(sdoc, "Supplementary notes", bold=True); blank(sdoc)
# ONLY THE NOTES THAT ARE SUBMITTED, renumbered sequentially, with no status label and no
# account of which were cut. The audit that decided them is a working record and belongs in
# the deposit, not in the document. NOTE_MAP, built above, keeps the three citations to a
# note pointing at the right one after the renumber.
if DRAFT:
    para(sdoc, SI.SI_NOTES_HEADER, italic=True); blank(sdoc)
for e in _SI_KEEP:
    para(sdoc, fix("%s.  %s" % (NOTE_MAP[e[0]], e[2])), bold=True); blank(sdoc)
    block(sdoc, [fix(_b) for _b in e[3]])
if ORDER:
    para(sdoc, "References", bold=True); blank(sdoc)
    para(sdoc, "The Supporting Information is numbered independently of the manuscript.", italic=True)
    blank(sdoc)
    for i, key in enumerate(ORDER, 1):
        para(sdoc, "%d. %s" % (i, REFMOD.REFS[key]["text"])); blank(sdoc)
s_out = os.path.join(OUT, "Supporting_Information_draft_%s.docx" % DATE)
sdoc.save(s_out)
_submit(s_out, "Supporting_Information.docx")

# ================================================================ THE SEPARATE UPLOAD
# Cell Chemical Biology requires Highlights, the eTOC blurb and the Significance paragraph in one
# Word document separate from the main text. It is uploaded under the Editorial Manager option
# that includes "highlights and eTOC blurb", and none of it counts against the 7,000.
xdoc = newdoc()
para(xdoc, us(S.FRONT["title"]), bold=True); blank(xdoc)
# The three "spec" strings quote the journal's own instructions back at it, and the header
# describes the upload. Both are notes to us. DRAFT keeps them; the submitted file does not.
if DRAFT:
    para(xdoc, "Highlights, eTOC blurb and Significance. Uploaded as one document, separate "
               "from the main text.", italic=True); blank(xdoc)
para(xdoc, "Highlights", bold=True); blank(xdoc)
if DRAFT:
    para(xdoc, S.HIGHLIGHTS["spec"], italic=True); blank(xdoc)
for _h in _hl_measured(S.HIGHLIGHTS["body"]):
    para(xdoc, "\u2022  " + _h); blank(xdoc)
para(xdoc, "eTOC blurb", bold=True); blank(xdoc)
if DRAFT:
    para(xdoc, S.ETOC["spec"], italic=True); blank(xdoc)
para(xdoc, us(S.ETOC["body"])); blank(xdoc)
para(xdoc, "Significance", bold=True); blank(xdoc)
if DRAFT:
    para(xdoc, S.SIGNIFICANCE["spec"], italic=True); blank(xdoc)
para(xdoc, fix(S.SIGNIFICANCE["body"])); blank(xdoc)
x_out = os.path.join(OUT, "Highlights_eTOC_Significance_%s.docx" % DATE)
xdoc.save(x_out)
_submit(x_out, "Highlights_eTOC_Significance.docx")
_submit_figures()
print()
print("separate upload : %s" % x_out)
print("  highlights %d (limit 4), longest %d characters (limit 85)"
      % (len(S.HIGHLIGHTS["body"]), max(len(h) for h in S.HIGHLIGHTS["body"])))
print("  eTOC %d words (limit 50), Significance %d words (limit 300), neither counted"
      % (len(S.ETOC["body"].split()), len(S.SIGNIFICANCE["body"].split())))

# ---------------------------------------------------------------- word gate
# Verified from the Cell Chem Biol article-types page, 19 October 2026: a Resource has "a word
# count under 7,000 (including main figure legends but excluding STAR Methods text, supplemental
# item legends, and references section), and no more than 7 main figures and/or tables".
def _w(x):
    if isinstance(x, str): return len(re.sub(r"\[\[[a-z0-9,]+\]\]", "", x).split())
    if isinstance(x, (list, tuple)):
        return sum(_w(i[0] if isinstance(i, tuple) else i) for i in x)
    return 0
_INTRO_E = getattr(SP, "INTRO", None) or D.P_INTRO
_DISC_E  = getattr(SP, "DISCUSSION", None) or D.P_DISCUSSION
rows = [("Summary (limit 150)", _w(D.P_ABSTRACT)), ("1 Introduction", _w(_INTRO_E))]
for sec in MAIN:
    eb, el = EMIT[sec["figno"]]
    rows.append(("%s body%s" % (sec["n"], "" if sec["n"] in SP.SPINE else "   [FULL]"), _w(eb)))
    rows.append(("    Figure %d legend" % sec["figno"], _w(el)))
rows += [("3 Discussion", _w(_DISC_E)),
         ("4 Conclusion", 0 if getattr(SP, "DISCUSSION", None) else _w(D.P_CONCLUSION)),
         ("Resource availability", _w(S.RESOURCE_AVAIL["body"]))]
total = sum(n for _, n in rows)
print()
print("---- Cell Chem Biol countable budget ----")
for label, n in rows: print("  %-30s %6d" % (label, n))
print("  %-30s %6d" % ("TOTAL", total))
print("  %-30s %6d   %s" % ("limit", 7000,
      "PASS, %d spare" % (7000 - total) if total <= 7000 else "FAIL, %d over" % (total - 7000)))
print("  main display items %d of 7 allowed" % len(MAIN))
if SPINE_WARN:
    print()
    print("  !! a kept panel has no condensed text, so it will NOT be emitted:")
    for _w_ in SPINE_WARN:
        print("     " + _w_)
print("  excluded: Methods %d, references %d, supplemental"
      % (_w(list(D.P_METHODS) + [DD.METHODS_ADDITION]
            + ([SP.METHODS_RELOCATED] if hasattr(SP, "METHODS_RELOCATED") else [])),
         len(MAIN_ORDER)))
print()
print("  renumbering: %s" % ", ".join("old %d -> %d" % (k, v) for k, v in sorted(FIG_MAP.items())))
# The ledgers are keyed by SOURCE figure number, so panel_data/ holds fig8_head.json while the
# paper has 7 figures. Nothing said so, which is the same old-against-new trap that left nine
# "Figure 8" citations in the SI. Write the map out beside the ledgers so the deposit explains
# itself to anyone who opens it.
import json as _json
_led = os.path.join(OUT, "panel_data")
if os.path.isdir(_led):
    with open(os.path.join(_led, "figure_numbering.json"), "w", encoding="utf-8") as _fh:
        _json.dump({"note": "Ledger files and panel CSVs are named by SOURCE figure number, "
                            "which is the pre-merge numbering. The published paper uses the "
                            "values on the right of source_to_published. Panels listed under cut_to_extended "
                            "are not in the published figure; their text is in Extended Results "
                            "in the Supporting Information.",
                    "source_to_published": {str(k): v for k, v in sorted(FIG_MAP.items())},
                    "published_figure_count": len(sorted(set(FIG_MAP.values()))),
                    # Panel-level map, so a checker can test a cross-reference drawn INTO the
                    # artwork. The documents are renumbered by fix(); the figures are not, and
                    # figS6 drew "Figure 5e" and figS10 "Figure 8d" against a paper whose
                    # published names for those panels are 5f and 7d.
                    "panel_map": {"%d%s" % k: "%d%s" % v for k, v in sorted(PANEL_MAP.items())},
                    "cut_to_extended": {str(k): v for k, v in sorted(PN.CUT.items()) if v}},
                  _fh, indent=1)
    print("  wrote panel_data/figure_numbering.json")
    with open(os.path.join(_led, "arith_audit.json"), "w", encoding="utf-8") as _fh:
        _json.dump(ARITH_AUDIT, _fh, indent=1, ensure_ascii=False)
    print("  wrote panel_data/arith_audit.json (%d entries, %d cut inside a list)"
          % (len(ARITH_AUDIT), sum(1 for a in ARITH_AUDIT if a["mid_list"])))
print("  panels: %d kept, %d moved to Extended Results" % (len(PANEL_MAP), len(CUT_OUT)))

_BRIT = ("artefact", "tumour", "labelled", "behaviour", "catalogue", "ionised", "modelled",
         "neighbour", "normalised", "standardised", "centre", "realised", "oesophagus",
         "generalised", "analysed", "localised", "defence", "fibre", "recognised")
_MANGLED = ("catalogd", "modeld", "labeld", "normalizd", "standardizd", "realizd")
for _label, _path in (("manuscript", m_out), ("SI", s_out)):
    import docx as _dx
    _t = "\n".join(p.text for p in _dx.Document(_path).paragraphs)
    _b = sorted({w for w in _BRIT if re.search(r"\b%s" % w, _t)})
    _m = sorted({w for w in _MANGLED if w.strip() in _t})
    if _b or _m:
        print("  !! %s: British forms left %s; mangled by the spelling pass %s"
              % (_label, _b or "none", _m or "none"))
    else:
        print("  %s spelling: clean" % _label)
if US_HITS:
    print("  US spelling fixes: %s" % ", ".join("%s->%d" % kv for kv in sorted(US_HITS.items())))
print()
print("manuscript : %s (%.1f MB)" % (m_out, os.path.getsize(m_out) / 1e6))
print("SI         : %s (%.1f MB)" % (s_out, os.path.getsize(s_out) / 1e6))
print("refs (ms)  : %d, in order: %s" % (len(MAIN_ORDER), ", ".join(MAIN_ORDER)))
print("refs (SI)  : %d, in order: %s" % (len(ORDER), ", ".join(ORDER)))
