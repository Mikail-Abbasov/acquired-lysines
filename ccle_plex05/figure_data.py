#!/usr/bin/env python3
"""
figure_data.py -- panel-level data deposit and a machine-readable ledger of headline numbers.

WHY THIS MODULE EXISTS
The Methods state that "every number in this paper is printed by the script that produces the
corresponding panel rather than transcribed from the artwork". That is true in substance -- the
figure scripts do print their numbers -- but as free text, so nothing can check it. Three of the
eighteen figure scripts kept a HEAD dict; the other fifteen did not. A claim that cannot be checked
is a claim a reviewer has to take on trust, and this paper does not ask for trust anywhere else.

Separately, no figure script wrote the data behind its own panels, so there was nothing to deposit
alongside the code: a reader could have the script and the artwork but not the table in between.

Both are the same fix. Every figure script now calls:

    fd.panel("fig2", "b", frame)     one tidy CSV per panel, written beside the figure
    fd.head("fig2", HEAD)            the panel's headline numbers as JSON

which gives the repository its deposited-data tier and gives verify.py something to diff against
the manuscript. The CSV is the data a reader would need to redraw the panel by hand; the JSON is
the set of values the prose quotes.

CONVENTIONS
  output          results/panel_data/<stem>_<panel>.csv and <stem>_head.json
  one row per     whatever the panel plots -- a bar, a point, a bin. Long/tidy, not wide.
  provenance      every CSV carries a header comment naming the script and the panel
  no rounding     values are written as computed; the figure rounds for display, the deposit
                  does not, so the two can differ and the CSV remains the source
"""
import os, json, csv, datetime

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.abspath(os.path.join(HERE, "..", "results", "panel_data"))


def _ensure():
    os.makedirs(OUT, exist_ok=True)


def panel(stem, letter, frame, note=""):
    """Write one panel's data. `frame` is a DataFrame or a list of dicts."""
    _ensure()
    path = os.path.join(OUT, "%s_%s.csv" % (stem, letter))
    try:
        import pandas as pd
        if not hasattr(frame, "to_csv"):
            frame = pd.DataFrame(list(frame))
        with open(path, "w", newline="") as fh:
            fh.write("# %s panel %s -- %s.py%s\n" % (stem, letter, stem,
                     ("  |  " + note) if note else ""))
            frame.to_csv(fh, index=False)
    except ImportError:
        rows = list(frame)
        if not rows:
            return path
        with open(path, "w", newline="") as fh:
            fh.write("# %s panel %s -- %s.py%s\n" % (stem, letter, stem,
                     ("  |  " + note) if note else ""))
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)
    print("    deposited %s (%d rows)" % (os.path.basename(path),
          len(frame) if hasattr(frame, "__len__") else -1))
    return path


def head(stem, HEAD):
    """Write the headline numbers this script computed, as JSON, for verify.py to diff."""
    _ensure()
    path = os.path.join(OUT, "%s_head.json" % stem)
    payload = {"script": stem + ".py",
               "written": datetime.date.today().isoformat(),
               "numbers": {k: HEAD[k] for k in sorted(HEAD)}}
    json.dump(payload, open(path, "w"), indent=1, default=str)
    print("    ledger    %s (%d values)" % (os.path.basename(path), len(HEAD)))
    return path
