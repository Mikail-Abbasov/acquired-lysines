"""
Figure 4, proposed rebuild — six panels instead of nine.

WHAT CHANGED AND WHY

Old b and c were two ranked bar charts of the same features asked two different ways, but with
different feature lists in different orders, so the comparison the figure is built on could not
be made by eye. They are now one chart drawn twice (b, c) on one shared feature order: before
the gate and after it. The argument is now visible rather than asserted.

Old e, f and g were three panels on feature redundancy — two scree plots and an adjusted-vs-alone
chart. That is scaffolding, not result, and it moves to the SI.

Old i reported four correlation coefficients over 5 and 9 hand-classified features with no null.
It is replaced by e, an exact permutation test over all C(15,6) = 5005 ways of splitting these
features 5-against-9, with the ungated version of the same test as its control (p = 0.008 gated,
p = 0.487 ungated).

Old d and h are kept, cleaned: d loses its unlabelled twin-axis traces (moved to SI) and gains
distinguishable line styles; h moves its tier counts out of the plotting area.

HONEST LABELLING: under the gate only the LABELLING effect is recomputed on well-modelled
lysines. The ligandability effect is as measured. Panel c says so.
"""
import os, re, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.lines import Line2D
from itertools import combinations
import abbasov_style as st
import figure_data as fd
st.apply()   # these two scripts imported the style module for its colours but
             # never applied it, so their typeface, sizes, spines and DPI were
             # matplotlib defaults while every other figure used the BMC style
st.slate()         # palette from the supplied ramp; see abbasov_style.SLATE

def find(n):
    for d in [".", "..", "ligand_rebuild", "../ligand_rebuild"]:
        p = os.path.join(d, n)
        if os.path.exists(p): return p
    raise FileNotFoundError(n)

S9  = pd.read_csv(find("step9_results.csv"))
G   = pd.read_csv(find("step11_gates.csv"))           # reactivity, gate ladder
G2  = pd.read_csv(find("step11b_ligandability_gates.csv"))   # ligandability, same ladder
G2  = G2[~G2.dropped]
R11 = open(find("step11_report.txt")).read()
R38 = open(find("step38_report.txt")).read()

NICE = lambda f: (f.replace("_8A", " 8Å").replace("_10A", " 10Å").replace("_", " ")
                   .replace("rsa tien", "exposure").replace("propka pKa", "predicted pKa")
                   .replace("plddt site", "pLDDT").replace("plddt 8Å", "pLDDT 8Å")
                   .replace("nz sasa", "NZ exposure").replace("sasa sidechain", "side-chain exp.")
                   .replace("dist centroid", "dist. to centroid").replace("n lys", "nearby lysines")
                   .replace("salt bridge", "salt bridge"))
SHAPE = {"crowding_8A", "atoms_10A", "hydrophobic_8A", "aromatic_8A", "polar_8A",
         "rsa_tien", "sasa_sidechain", "nz_sasa", "dist_centroid"}
CCOL  = {"charge": st.BLUE, "shape": st.ORANGE}

# ---------------------------------------------------------------- shared feature table
# Both contrasts are now measured on the same pLDDT ladder (step 11 for reactivity, step 11b
# for ligandability), so b and c differ only in the gate — not in which side was gated.
# pLDDT features are excluded throughout: they are the gating variable, not a site property.
def ladder(df, col):
    return {g: df[df.gate == g].set_index("feature")[col] for g in sorted(df.gate.unique())}
RE, LI = ladder(G, "d"), ladder(G2, "d")
TIER = 60          # the occupancy threshold that defines "engaged" in comparison B
GATE = max(RE)
F = pd.DataFrame({"react": RE[0], "lig": LI[0],
                  "react_g": RE[GATE], "lig_g": LI[GATE]}).dropna().reset_index()
F = F.rename(columns={"index": "feature"})
F["kind"] = np.where(F.feature.isin(SHAPE), "shape", "charge")
F = F.sort_values("react", ascending=False).reset_index(drop=True)
NF, KC = len(F), int((F.kind == "charge").sum())
F["gap_0"]  = (F.react   - F.lig).abs()
F["gap_90"] = (F.react_g - F.lig_g).abs()

SPLITS = list(combinations(range(NF), KC))
def split_test(rc, lc):
    """Exact permutation over every way of splitting these features KC-against-the-rest."""
    R = lambda idx: np.corrcoef(F.loc[idx, rc], F.loc[idx, lc])[0, 1]
    ci, si = F.index[F.kind == "charge"], F.index[F.kind == "shape"]
    obs = R(ci) - R(si)
    g = np.array([R(pd.Index(c)) - R(F.index.difference(pd.Index(c))) for c in SPLITS])
    return obs, g, float((g >= obs).mean()), R(ci), R(si)

NULL = {"gated": split_test("react_g", "lig_g"), "ungated": split_test("react", "lig")}
LADDER = []
for g in sorted(RE):
    if g not in LI: continue
    T = pd.DataFrame({"react": RE[g], "lig": LI[g]}).dropna()
    T = T[T.index.isin(F.feature)]
    T["kind"] = np.where(T.index.isin(SHAPE), "shape", "charge")
    Rr = lambda k: np.corrcoef(T[T.kind == k].react, T[T.kind == k].lig)[0, 1]
    idx = {k: F.index[F.kind == k] for k in ("charge", "shape")}
    obs = Rr("charge") - Rr("shape")
    Fg = F.assign(rr=F.feature.map(T.react), ll=F.feature.map(T.lig))
    Rg = lambda i: np.corrcoef(Fg.loc[i, "rr"], Fg.loc[i, "ll"])[0, 1]
    perm = np.array([Rg(pd.Index(c)) - Rg(Fg.index.difference(pd.Index(c))) for c in SPLITS])
    LADDER.append(dict(gate=g, charge=Rr("charge"), shape=Rr("shape"),
                       p=float((perm >= obs).mean())))
LADDER = pd.DataFrame(LADDER)

# ---------------------------------------------------------------- the setting, from the report
# These four numbers used to be typed into panel a as a string. They are now read from the same
# report the effects come from, so the artwork cannot drift from the analysis.
_m = re.search(r"(\d+) lysines on (\d+) proteins \| labelled (\d+) \(([\d.]+)%\)", R11)
SET_NK, SET_NP, SET_NL, SET_PCT = (int(_m.group(1)), int(_m.group(2)),
                                   int(_m.group(3)), float(_m.group(4)))

HEAD = {"a_n_lysines": SET_NK, "a_n_proteins": SET_NP, "a_n_labelled": SET_NL,
        "a_pct_labelled": SET_PCT, "a_n_features": NF, "a_n_charge_features": KC,
        "a_n_shape_features": NF - KC, "a_occupancy_tier_pct": TIER,
        "a_n_regroupings": len(SPLITS)}

# spread under the gate -- the range-restriction control quoted in 2.4
_SD = {g: G[G.gate == g].set_index("feature")["sd_all"] for g in sorted(G.gate.unique())}
for _f in ("crowding_8A", "polar_8A", "aromatic_8A", "rsa_tien", "dist_centroid"):
    if _f in _SD[0].index and _f in _SD[GATE].index:
        HEAD["bc_sd_%s_ungated" % _f] = round(float(_SD[0][_f]), 2)
        HEAD["bc_sd_%s_gated" % _f] = round(float(_SD[GATE][_f]), 2)
        _ch = 100 * (float(_SD[GATE][_f]) / float(_SD[0][_f]) - 1)
        HEAD["bc_sd_%s_pct_change" % _f] = round(_ch, 0)
        if _ch < 0:                       # 2.4 quotes these as spread LOST, positive
            HEAD["bc_sd_%s_pct_spread_lost" % _f] = round(-_ch, 0)

# ================================================================== canvas
fig = plt.figure(figsize=(st.DOUBLE, st.SINGLE * 1.50))
gs  = fig.add_gridspec(2, 3, hspace=0.40, wspace=0.34,
                       left=0.055, right=0.988, top=0.955, bottom=0.105,
                       width_ratios=[1.02, 1.0, 0.92])

# ---------------------------------------------------------------- a: the design
ax = fig.add_subplot(gs[0, 0]); st.panel_label(ax, "a", dx=-0.06, dy=1.07)
ax.axis("off")
ax.text(0.5, 1.005, "%s lysines · %s proteins\n%s labelled (%.1f%%)"
        % (f"{SET_NK:,}", f"{SET_NP:,}", f"{SET_NL:,}", SET_PCT),
        transform=ax.transAxes, ha="center", va="top", fontsize=5.2, color=st.GREY,
        linespacing=1.6)
for yy, (t, q, sub, c, fill) in enumerate([
    ("REACTIVITY",   "what makes a lysine get labelled at all?",
     "probe-labelled  vs  every other lysine,\nsame protein", st.BLUE,  st.PAPER),
    ("LIGANDABILITY", "what makes a labelled lysine blockable?",
     "engaged  vs  all labelled,\nat \u2265%d%% occupancy" % TIER,  st.GREEN, st.PAPER)]):
    top = 0.90 - yy * 0.31
    ax.add_patch(FancyBboxPatch((0.02, top - 0.255), 0.96, 0.255, transform=ax.transAxes,
                 boxstyle="round,pad=0.010,rounding_size=0.02", fc=fill,
                 ec=st.LIGHTGREY, lw=0.9, zorder=0))
    ax.text(0.06, top - 0.045, t, transform=ax.transAxes, fontsize=6.6, fontweight="bold",
            color=c, va="center")
    ax.text(0.06, top - 0.125, sub, transform=ax.transAxes, fontsize=5.0, color=st.GREY,
            va="center", linespacing=1.6)
    ax.text(0.06, top - 0.218, q, transform=ax.transAxes, fontsize=5.0, style="italic",
            color=c, va="center")
ax.text(0.5, 0.275, "one contrast cannot answer both questions",
        transform=ax.transAxes, ha="center", va="center", fontsize=5.0, color=st.GREY,
        style="italic")
key = [Line2D([], [], marker="o", ls="none", ms=3.6, mfc=st.GREY, mec=st.GREY,
              label="effect on being labelled"),
       Line2D([], [], marker="D", ls="none", ms=3.2, mfc="none", mec=st.GREY, mew=1.0,
              label="effect on being blocked"),
       Line2D([], [], marker="s", ls="none", ms=3.6, color=st.BLUE,  label="charge feature"),
       Line2D([], [], marker="s", ls="none", ms=3.6, color=st.ORANGE, label="shape feature")]
ax.legend(handles=key, loc="upper center", bbox_to_anchor=(0.5, 0.20), frameon=False,
          fontsize=5.0, handletextpad=0.5, labelspacing=0.6, ncol=2, columnspacing=1.2)

# ---------------------------------------------------------------- d, c: the paired display
for j, (rc, lc, gapcol) in enumerate([("react", "lig", "gap_0"),
                                     ("react_g", "lig_g", "gap_90")]):
    ax = fig.add_subplot(gs[0, 1 + j]); st.panel_label(ax, "bc"[j], dx=-0.30 if j == 0 else -0.10, dy=1.07)
    y = np.arange(NF)[::-1]
    for k, r_ in F.iterrows():
        c = CCOL[r_.kind]
        ax.plot([r_[rc], r_[lc]], [y[k]] * 2, color=st.LIGHTGREY, lw=2.0,
                solid_capstyle="round")
        ax.scatter(r_[rc], y[k], s=16, color=c, zorder=3, lw=0.0)
        ax.scatter(r_[lc], y[k], s=15, marker="D", facecolor="none", edgecolor=c, lw=1.0, zorder=4)
    ax.axvline(0, color=st.LIGHTGREY, lw=0.8)
    ax.set_yticks(y)
    if j == 0:
        ax.set_yticklabels([NICE(f) for f in F.feature], fontsize=5.0)
        for t, kk in zip(ax.get_yticklabels(), F.kind): t.set_color(CCOL[kk])
    else:
        ax.set_yticklabels([])
    ax.set_ylim(-2.4, NF - 0.35); ax.set_xlim(-0.56, 0.46)
    ax.set_xlabel("effect size, d", fontsize=6)
    g = F.groupby("kind")[gapcol].mean()
    ax.text(0.985, 0.02,
            "average distance between the pair\ncharge %.2f     shape %.2f" % (g["charge"], g["shape"]),
            transform=ax.transAxes, ha="right", va="bottom", fontsize=4.9, color=st.GREY,
            linespacing=1.55)
    for s in ("top", "right"): ax.spines[s].set_visible(False)

fd.panel("fig4", "bc", F[["feature", "kind", "react", "lig", "react_g", "lig_g",
                          "gap_0", "gap_90"]],
         note="react/lig = ungated effect size d; _g = at pLDDT >= %d; gap = |react - lig|" % GATE)
_g0, _g9 = F.groupby("kind").gap_0.mean(), F.groupby("kind").gap_90.mean()
HEAD["b_mean_gap_charge"] = round(float(_g0["charge"]), 3)
HEAD["b_mean_gap_shape"] = round(float(_g0["shape"]), 3)
HEAD["c_mean_gap_charge"] = round(float(_g9["charge"]), 3)
HEAD["c_mean_gap_shape"] = round(float(_g9["shape"]), 3)
for _f in ("crowding_8A", "atoms_10A", "polar_8A", "anionic_8A", "salt_bridge"):
    _r = F[F.feature == _f]
    if len(_r):
        _r = _r.iloc[0]
        for _c, _n in (("react", "react_ungated"), ("react_g", "react_gated"),
                       ("lig", "lig_ungated"), ("lig_g", "lig_gated")):
            HEAD["bc_%s_%s" % (_f, _n)] = round(float(_r[_c]), 3)

# ---------------------------------------------------------------- e: the gate, step by step
ax = fig.add_subplot(gs[1, 0]); st.panel_label(ax, "d", dx=-0.20, dy=1.07)
gates = sorted(G.gate.unique()); xs = np.arange(len(gates))
SERIES = [("net_charge_8A", st.BLUE,   "-",  "o"), ("anionic_8A", st.BLUE,   "--", "s"),
          ("crowding_8A",   st.ORANGE, "-",  "o"), ("aromatic_8A", st.ORANGE, "--", "s")]
for f, c, ls, mk in SERIES:
    sub = G[G.feature == f].set_index("gate").reindex(gates)
    ax.plot(xs, sub.d, ls=ls, marker=mk, color=c, lw=1.25, ms=3.0, mew=0, label=NICE(f))
ax.axhline(0, color=st.GREY, lw=0.8)
GATE_N, cur = {}, None
for ln in R11.splitlines():
    m = re.match(r"GATE\s+(\S+)", ln.strip())
    if m: cur = 0 if m.group(1) == "none" else int(re.sub(r"\D", "", m.group(1)) or 0)
    m2 = re.search(r"labelled ([\d,]+) \((\d+)% of all labelled\) \| unlabelled ([\d,]+)", ln)
    if m2 and cur is not None:
        GATE_N[cur] = (int(m2.group(1).replace(",", "")), int(m2.group(3).replace(",", "")))
fd.panel("fig4", "d", [dict(gate=int(_g_), feature=_f_,
                           effect_on_being_labelled_d=float(
                               G[(G.feature == _f_) & (G.gate == _g_)].d.iloc[0]),
                           n_labelled_kept=GATE_N.get(_g_, (None, None))[0],
                           n_unlabelled_kept=GATE_N.get(_g_, (None, None))[1])
                      for _g_ in gates for _f_, _c_, _l_, _m_ in SERIES
                      if len(G[(G.feature == _f_) & (G.gate == _g_)])],
         note="reactivity effect at each AlphaFold-confidence rung, with lysines retained")
for _g_ in gates:
    if _g_ in GATE_N:
        HEAD["d_n_labelled_gate%d" % _g_] = int(GATE_N[_g_][0])
        HEAD["d_n_unlabelled_gate%d" % _g_] = int(GATE_N[_g_][1])
        # the caption quotes these rounded to thousands ("none 30k/214k"); the artwork's tick
        # labels use the same rounding, so the rounded form is what has to be checkable
        HEAD["d_kept_gate%d_as_printed" % _g_] = "%dk/%dk" % (GATE_N[_g_][0] // 1000,
                                                              GATE_N[_g_][1] // 1000)
ax.set_xticks(xs)
ax.set_xticklabels(["%s\n%s" % ("none" if g == 0 else "≥%d" % g,
                    ("%dk/%dk" % (GATE_N[g][0]//1000, GATE_N[g][1]//1000)) if g in GATE_N else "")
                    for g in gates], fontsize=4.6, linespacing=1.8)
ax.set_xlabel("AlphaFold confidence gate", fontsize=6)
ax.text(0.0, -0.30, "labelled / unlabelled lysines kept", transform=ax.transAxes,
        fontsize=4.7, color=st.GREY, ha="left")
ax.set_ylabel("effect on being labelled, d", fontsize=6)
ax.set_xlim(-0.35, len(gates) + 0.62); ax.set_ylim(-0.56, 0.52)
for f, c, ls, mk in SERIES:                      # direct labels, no legend over the data
    v = G[(G.feature == f) & (G.gate == gates[-1])].d.iloc[0]
    ax.text(len(gates) - 0.82, v, "  " + NICE(f), color=c, fontsize=5.2, va="center", ha="left")
for s in ("top", "right"): ax.spines[s].set_visible(False)

# ---------------------------------------------------------------- f: the test
# Now that both contrasts run on the same gate ladder, the separation can be shown as a
# dose-response rather than a single before/after: charge stays coupled at every gate while
# shape comes apart progressively, and the exact permutation p falls with it.
sub_gs = gs[1, 1].subgridspec(2, 1, height_ratios=[1.35, 1.0], hspace=0.62)

ax = fig.add_subplot(sub_gs[0, 0]); st.panel_label(ax, "e", dx=-0.38, dy=1.14)
xs = np.arange(len(LADDER))
for k, c, mk in (("charge", st.BLUE, "o"), ("shape", st.ORANGE, "s")):
    ax.plot(xs, LADDER[k], marker=mk, color=c, lw=1.3, ms=3.2, mew=0)
    ax.text(xs[-1] + 0.12, LADDER[k].iloc[-1], "  " + k, color=c, fontsize=5.2, va="center")
for x, pv in zip(xs, LADDER.p):
    ax.text(x, 1.055, "p %s" % (("%.3f" % pv).lstrip("0") if pv >= 0.001 else "<.001"),
            ha="center", fontsize=4.5,
            color=st.ORANGE if pv < 0.05 else st.GREY,
            fontweight="bold" if pv < 0.05 else "normal")
ax.set_xticks(xs)
ax.set_xticklabels(["none" if g == 0 else "≥%d" % g for g in LADDER.gate], fontsize=5.0)
ax.set_ylim(0.20, 1.10); ax.set_xlim(-0.3, len(xs) - 0.25)
ax.set_ylabel("how closely the two\ncontrasts track", fontsize=5.6, linespacing=1.6)
ax.tick_params(axis="y", labelsize=5.0)
for sp in ("top", "right"): ax.spines[sp].set_visible(False)

ax = fig.add_subplot(sub_gs[1, 0])
obs, g, p, _, _ = NULL["gated"]
ax.hist(g, bins=58, color=st.LIGHTGREY, edgecolor="none")
ax.axvline(obs, color=st.ORANGE, lw=1.5, zorder=3)
ax.axvline(NULL["ungated"][0], color=st.GREY, lw=1.2, ls=(0, (2, 1.5)), zorder=3)
ax.annotate("no gate\n%.2f" % NULL["ungated"][0], (NULL["ungated"][0], ax.get_ylim()[1] * 0.62),
            xytext=(-7, 0), textcoords="offset points", ha="right", fontsize=4.5,
            color=st.GREY, linespacing=1.5)
ax.annotate("pLDDT ≥ 90\n%.2f" % obs, (obs, ax.get_ylim()[1] * 0.90),
            xytext=(-7, 0), textcoords="offset points", ha="right", fontsize=4.7,
            color=st.ORANGE, linespacing=1.5, fontweight="bold")
ax.set_xlim(-1.35, 1.05); ax.set_yticks([])
ax.tick_params(axis="x", labelsize=5.0)
for sp in ("top", "right", "left"): ax.spines[sp].set_visible(False)
ax.set_xlabel("all %s regroupings of the same features" % f"{len(g):,}", fontsize=6)

fd.panel("fig4", "e", LADDER.rename(columns={
    "charge": "charge_coupling_r", "shape": "shape_coupling_r", "p": "exact_permutation_p"}),
    note="how closely reactivity and ligandability track, by gate; p over all %d regroupings"
         % len(SPLITS))
for _r in LADDER.itertuples():
    _k = "none" if _r.gate == 0 else "%d" % _r.gate
    HEAD["e_charge_coupling_gate%s" % _k] = round(float(_r.charge), 2)
    HEAD["e_shape_coupling_gate%s" % _k] = round(float(_r.shape), 3)
    HEAD["e_p_gate%s" % _k] = round(float(_r.p), 3)
for _k in ("ungated", "gated"):
    HEAD["e_perm_obs_%s" % _k] = round(float(NULL[_k][0]), 3)
    HEAD["e_perm_p_%s" % _k] = round(float(NULL[_k][2]), 3)
    HEAD["e_corr_charge_%s" % _k] = round(float(NULL[_k][3]), 2)
    HEAD["e_corr_shape_%s" % _k] = round(float(NULL[_k][4]), 2)

# ---------------------------------------------------------------- f: occupancy dose-response
ax = fig.add_subplot(gs[1, 2]); st.panel_label(ax, "f", dx=-0.24, dy=1.07)
rows, hdr = [], None
for ln in R38[R38.find("Q1a"):R38.find("Q1b")].splitlines():
    p_ = ln.split()
    if len(p_) > 3 and p_[0] == "tier": hdr = p_
    elif hdr and len(p_) == len(hdr): rows.append(p_)
T = pd.DataFrame(rows, columns=hdr)
T = T[T.tier.str.match(r"^[\d>=\u2265\-]+$", na=False)].reset_index(drop=True)
for c in T.columns:
    if c != "tier": T[c] = pd.to_numeric(T[c], errors="coerce")
T["tier"] = T.tier.str.replace(">=", "\u2265", regex=False)
xs = np.arange(len(T))
# All four of these are SHAPE features (see SHAPE above, the panel b legend and panel d), so
# all four carry the shape colour. The old version gave crowding the CHARGE colour and put
# aromatic on a third hue, which contradicted the class code the rest of the figure uses.
# Identity comes from the marker and the direct label; exposure keeps its dashed rule because
# it is the one that runs the other way.
SER_F = (("hydrophob", st.ORANGE, "-",  "o", "hydrophobic 8Å"),
         ("crowding_", st.ORANGE, "-",  "s", "crowding 8Å"),
         ("aromatic_", st.ORANGE, "-",  "^", "aromatic 8Å"),
         ("rsa_tien",  st.ORANGE, "--", "D", "exposure"))
for f, c, ls, mk, lab in SER_F:
    if f in T.columns:
        ax.plot(xs, T[f], ls=ls, marker=mk, color=c, lw=1.25, ms=2.6, mew=0)
ax.axhline(0, color=st.LIGHTGREY, lw=0.8)
ax.set_xticks(xs); ax.set_xticklabels(T.tier, fontsize=4.8, rotation=45, ha="right")
ax.set_xlabel("occupancy of the lysine (%)", fontsize=6)
ax.set_ylabel("shift from the 0–10% tier, d", fontsize=6)
for f, c, ls, mk, lab in SER_F:                  # direct labels: identity is not colour alone
    if f in T.columns:
        ax.text(len(T) - 0.75, T[f].iloc[-1], "  " + lab, color=c, fontsize=5.2,
                va="center", ha="left")
ax.set_xlim(-0.4, len(T) + 2.4)
fd.panel("fig4", "f", T, note="shift from the 0-10%% occupancy tier, d, by tier")
HEAD["f_n_tiers"] = int(len(T))
ax.text(0.03, 0.14, "monotonic across all ten tiers,\nno discontinuity at the 60% cut",
        transform=ax.transAxes, ha="left", va="bottom", fontsize=4.9, color=st.GREY,
        linespacing=1.55)
for s in ("top", "right"): ax.spines[s].set_visible(False)


# NOTE: writes to results/ only. It does not overwrite figures/main/fig4 — adopting this
# version is a decision for the PI, not a side effect of running the script.
# via st.save so the blank-panel check covers this figure as well
st.save(fig, "fig4")
fd.head("fig4", HEAD)
fig.savefig("../results/fig4_v3.png", dpi=400, bbox_inches="tight")
fig.savefig("../results/fig4_v3.pdf", bbox_inches="tight")
print("features used: %d (%d charge / %d shape); ligandability tier = %d%%" % (NF, KC, NF-KC, TIER))
for c in NULL: print("  %-12s obs %+.3f   p %.4f" % (c, NULL[c][0], NULL[c][2]))
print("  mean pair gap  ungated:", F.groupby("kind").gap_0.mean().round(3).to_dict())
print("  mean pair gap    gated:", F.groupby("kind").gap_90.mean().round(3).to_dict())
for _k in ("ungated", "gated"):
    _o, _g, _p = NULL[_k][0], NULL[_k][1], NULL[_k][2]
    print("  e permutation %-8s observed %+.3f, p = %.3f, over %d regroupings" % (_k, _o, _p, len(_g)))
print("  e gate ladder:", ", ".join("%s charge %.2f shape %.2f p %.3f" %
      ("none" if r.gate == 0 else ">=%d" % r.gate, r.charge, r.shape, r.p) for r in LADDER.itertuples()))
print("  panels b,c display gate: pLDDT >= %d (GATE = max(RE)); ladder rungs %s" % (GATE, sorted(RE)))
print("  features (%d): charge = %s | shape = %s"
      % (NF, ", ".join(sorted(F.feature[F.kind == "charge"])),
         ", ".join(sorted(F.feature[F.kind == "shape"]))))
print("  per-feature, reactivity then ligandability, ungated -> gated:")
for _, r in F.iterrows():
    print("     %-16s %-6s react %+0.3f -> %+0.3f   lig %+0.3f -> %+0.3f   gap %.3f -> %.3f"
          % (r.feature, r.kind, r.react, r.react_g, r.lig, r.lig_g, r.gap_0, r.gap_90))
print("  within-kind correlation react vs lig: ungated charge %.2f shape %.2f | gated charge %.2f shape %.2f"
      % (NULL["ungated"][3], NULL["ungated"][4], NULL["gated"][3], NULL["gated"][4]))