# -*- coding: utf-8 -*-
"""Move designated panels out of the main text into Extended Results, and reletter what stays.

Text only. The artwork is NOT touched: the figure PNGs still show every original panel, and
Mikail rebuilds them. Until he does, the main-text legends and the artwork disagree, which is
why the build prints that warning on every run.

Rules, chosen so no prose is ever silently mangled:
  * A legend entry belongs to the panels its leading "(a)" or "(b, c)" names. An entry covering
    several panels moves only if ALL of them are cut; otherwise it stays.
  * A body paragraph moves only if EVERY panel it cites is cut. A paragraph citing a mix of kept
    and cut panels stays in the main text, and its citation of the cut panel is redirected to
    Extended Results. Those paragraphs are listed by the build so they can be rewritten when the
    section is condensed, rather than being cut by a rule that cannot read them.
  * Retained panels reletter to a sequential run from a. Body citations are rewritten in ONE pass
    from a map, never by successive replacement, which would cascade c->b->a.
  * Moved panels keep their ORIGINAL letters in Extended Results, so each one still names the
    artwork panel it came from.
"""
import re

# Selected by role, never by cost: the expensive panels are expensive because they carry the
# argument. Every panel holding a spine claim, or a control that exists nowhere else, stays.
# Figure 1d keeps the driver-hotspot odds ratio the Summary quotes; Figure 4c is the gated half
# of the paired display and so is Figure 4's result; both were briefly cut and put back.
# Second pass, 19 Oct: 1f (the lineage control), 2f (the downstream-peptide route) and 3d
# (spectral quality, which carried no body text at all) move to Extended Results. Each is an
# audit or an extension rather than a spine claim, and none of the bounds travels with them.
# 20 Oct: the panel cuts were made to reach 7,000 words. Once the decision was taken to
# submit over the limit, they were buying something no longer needed, and they left the
# figures thin. Fourteen panels are restored, bringing every figure to six except old 7
# and old 8, which have five between them to give. Only these seven stay in Extended
# Results: 1e, whose signature enrichment the paper itself calls non-independent, and the
# six that would push the merged Figure 5 past the density it already has.
# 5 Nov: 1e is restored at the PI's direction. The figure returns to seven panels, keeps
# its source letters a to g, and is no longer a composite, so main-text Figure 1 is the
# render from fig1_v3.py rather than a repanelled reading copy. Six panels stay out, all
# of them from the merged Figure 5.
CUT = {1: [], 2: [], 3: [], 4: [],
       5: ["b", "d", "f"], 6: ["c", "e", "f"], 7: [], 8: []}

# ---------------------------------------------------------------- the 5+6 merge
# Cell Chem Biol allows no more than 7 main figures and/or tables, and the paper has 8. Rather
# than demote a result to the supplement, Figures 5 and 6 merge into one display item. They are
# the same claim stated twice: "Ligandability is reproducible, and poorly predicted by the site
# features tested" against "Engaged lysines share site properties across compounds, and no tested
# descriptor predicts which". The 148-word Summary already collapses them into a single clause.
MERGE_INTO, MERGE_FROM = 5, 6

# The merged figure's panels, in argument order: what is measured, it reproduces, structure fails
# to predict, no useful operating point, compounds fail to predict, breadth.
MERGED_ORDER = [("a", 6, "a"),   # what is measured, and what is done with it
                ("b", 6, "b"),   # 13 of 14 site features point the same way in every compound
                ("c", 5, "a"),   # the measured model against the reproducibility ceiling
                ("d", 5, "c"),   # no operating point at which the model is useful
                ("e", 6, "d"),   # no compound property predicts which lysines a compound prefers
                ("f", 5, "e")]   # breadth across chemotypes

MERGED_HEADING = ("Ligandability reproduces at the site level and resists prediction from site "
                  "structure or compound properties")

# After the merge the display items are 1,2,3,4,5(merged),6(was 7),7(was 8).
RENUMBER = {7: 6, 8: 7}

def new_figno(old):
    return RENUMBER.get(old, old)

def merged_map():
    """(source figure, source panel) -> new panel letter on the merged figure."""
    return {(f, p): new for new, f, p in MERGED_ORDER}

# CASE. Every legend in this paper is written "(A) The design.", and every one of these
# three parsers matched [a-z] only, so all three were blind to all of them. For a figure that
# neither moves nor loses a panel that is invisible. For the 5+6 merge it is not: by_panel
# found no letters, every entry was classed as shared, and published Figure 5 came out with
# a preamble and no panel entries at all. Letters are read case-insensitively and returned
# lowercase, which is the key PANEL_MAP and CUT use.
def entry_letters(entry):
    """Panel letters a legend entry covers, read from its leading parenthesis."""
    head = entry[:90]
    got = []
    for a, b in re.findall(r"\(([A-Za-z])(?:\s*[–—-]\s*([A-Za-z]))?\)", head):
        a = a.lower()
        got.append(a)
        if b:
            got += [chr(c) for c in range(ord(a) + 1, ord(b.lower()) + 1)]
    for a, b in re.findall(r"\(([A-Za-z]),\s*([A-Za-z])\)", head):
        got += [a.lower(), b.lower()]
    return sorted(set(got))

def body_letters(para, figno):
    return sorted(set(l.lower() for l in
                      re.findall(r"Figure\s+%d([A-Za-z])" % figno, para)))

def split_legend(cap, figno):
    """-> (main_entries, moved_entries) in original order."""
    cut = set(CUT.get(figno, []))
    main, moved = [], []
    for e in cap:
        ls = entry_letters(e)
        if ls and set(ls) <= cut:
            moved.append((ls, e))
        else:
            main.append((ls, e))
    return main, moved

def split_body(paras, figno):
    """-> (main_paras, moved_paras, mixed) where mixed lists (index, kept, cut) for the
    paragraphs that cite both and so cannot be moved by rule."""
    cut = set(CUT.get(figno, []))
    main, moved, mixed = [], [], []
    for i, p in enumerate(paras):
        ls = body_letters(p, figno)
        if ls and set(ls) <= cut:
            moved.append(p)
        else:
            if ls and set(ls) & cut:
                mixed.append((i, [l for l in ls if l not in cut], [l for l in ls if l in cut]))
            main.append(p)
    return main, moved, mixed

def reletter_map(figno, all_letters):
    """Retained original letter -> new sequential letter."""
    cut = set(CUT.get(figno, []))
    kept = [l for l in sorted(set(all_letters)) if l not in cut]
    return {old: chr(ord("a") + i) for i, old in enumerate(kept)}

def apply_reletter(text, figno, rmap):
    """One pass. Kept panels take their new letter; cut panels point at Extended Results."""
    cut = set(CUT.get(figno, []))
    def rep(m):
        old = m.group(1)
        if old in rmap:
            return "Figure %d%s" % (figno, rmap[old])
        if old in cut:
            return "Extended Results"
        return m.group(0)
    return re.sub(r"Figure\s+%d([a-z])" % figno, rep, text)

def reletter_entry(entry, rmap):
    """Rewrite a legend entry's own leading letters to the new run.

    Case-insensitive on the way in and case-preserving on the way out, because the legends
    are written "(A)" and rmap is keyed on lowercase letters.
    """
    def rep(m):
        inner = m.group(1)
        parts = re.split(r"(\s*[,–—-]\s*)", inner)
        out = []
        for seg in parts:
            s = seg.strip()
            if s.lower() in rmap:
                new = rmap[s.lower()]
                out.append(new.upper() if s.isupper() else new)
            else:
                out.append(seg)
        return "(" + "".join(out) + ")"
    return re.sub(r"\(([A-Za-z](?:\s*[,–—-]\s*[A-Za-z])*)\)", rep, entry, count=1)


def by_panel(items, figno, is_legend):
    """Split a body or legend into (shared, {letter: [entries]}) in original order.
    Used to build the merged display item panel by panel."""
    shared, per = [], {}
    for x in items:
        t = x[0] if isinstance(x, tuple) else x
        ls = entry_letters(t) if is_legend else body_letters(t, figno)
        if not ls:
            shared.append(x)
        else:
            per.setdefault(tuple(ls), []).append(x)
    flat = {}
    for ls, xs in per.items():
        for l in ls:
            flat.setdefault(l, []).extend(xs)
    return shared, flat
