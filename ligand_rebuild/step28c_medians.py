#!/usr/bin/env python3
"""
step28c_medians.py — per-population medians and IQRs for Figure 8 panel b, under the CB
anchoring used in the main analysis (arm A). Small output so the figure script does not have
to re-read 17 MB of union features.

INPUT   union_distribution.csv, measured_cb.csv, features_expanded.csv
OUTPUT  step28c_medians.csv
"""
import numpy as np, pandas as pd

FEATS = ["crowding_8A", "atoms_10A", "hydrophobic_8A", "aromatic_8A", "polar_8A",
         "anionic_8A", "cationic_8A", "net_charge_8A", "n_lys_8A"]
ek = set(map(tuple, pd.read_csv("step28_engaged_keys.csv").values))
feat = pd.read_csv("features_expanded.csv", low_memory=False).drop_duplicates(["uniprot", "position"])
key_eng = {(u, int(p)) for u, p in zip(feat.uniprot, feat.position) if (u, int(p)) in ek}

u = pd.read_csv("union_distribution.csv")
m = pd.read_csv("measured_cb.csv")
m["engaged"] = [(x, int(p)) in key_eng for x, p in zip(m.uniprot, m.position)]

rows = []
for pop, df, pre in (("acquired union", u, ""),
                     ("engaged", m[m.engaged], "cb_"),
                     ("labelled, not engaged", m[~m.engaged], "cb_")):
    for f in FEATS:
        v = pd.to_numeric(df[pre + f], errors="coerce").dropna()
        rows.append(dict(population=pop, feature=f, n=len(v), median=v.median(),
                         q25=v.quantile(.25), q75=v.quantile(.75), mean=v.mean(), sd=v.std()))
d = pd.DataFrame(rows)
d.to_csv("step28c_medians.csv", index=False)
print(d.pivot(index="feature", columns="population", values="median").to_string())
print("\nwrote step28c_medians.csv")
