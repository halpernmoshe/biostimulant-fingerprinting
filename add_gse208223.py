"""
Score GSE208223 (Ti nanoparticles, Arabidopsis shoot) on all 3 axes,
add to all_biostimulants_scores.csv, and regenerate combined figures.
"""
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

ROOT     = Path("C:/Users/moshe/Dropbox/ISF 2025")
SIG_DIR  = ROOT / "RNA_seq_axes" / "signatures"
CAL_FILE = ROOT / "RNA_seq_axes" / "calibration.csv"
SUPPL    = ROOT / "RNA_seq_data" / "suppl_files" / "GSE208223"
FIG_DIR  = ROOT / "state_space_figures"
FIG_DIR.mkdir(exist_ok=True)

# ── Load axis signatures ──────────────────────────────────────────────────────
sig_p   = pd.read_csv(SIG_DIR / "P_axis_logFC.csv",   index_col=0).iloc[:, 0]
sig_n   = pd.read_csv(SIG_DIR / "N_axis_logFC.csv",   index_col=0).iloc[:, 0]
sig_aba = pd.read_csv(SIG_DIR / "ABA_axis_logFC.csv", index_col=0).iloc[:, 0]
sigs    = {"P_axis": sig_p, "N_axis": sig_n, "ABA_axis": sig_aba}

# ── Load calibration ─────────────────────────────────────────────────────────
cal = pd.read_csv(CAL_FILE, index_col=0)
print("Calibration:\n", cal.to_string())

def to_pct(raw_score, axis):
    if axis not in cal.index:
        return raw_score
    c = cal.loc[axis, "ctrl_score"]
    t = cal.loc[axis, "treat_score"]
    rng = t - c
    if abs(rng) < 1e-9:
        return np.nan
    return 100.0 * (raw_score - c) / rng

# ── Score GSE208223 shoot 10-day samples ─────────────────────────────────────
print("\nScoring GSE208223 shoot samples...")
df_s = pd.read_csv(SUPPL / "GSE208223_Shoot_norm.txt.gz", sep="\t", index_col=0)
df_s.index = df_s.index.str.replace(r"\.\d+$", "", regex=True)
df_s.index.name = "GeneID"

# Sample layout:
# E1/E2/E3 = Control shoot 10d, F1/F2/F3 = Ti-treated shoot 10d
groups_208 = {
    "TiO2_control": (["E1", "E2", "E3"], "Control (TiO2)"),
    "TiO2_treated": (["F1", "F2", "F3"], "Ti nanoparticles"),
}

rows_208 = []
for grp_key, (cols, cat) in groups_208.items():
    for ax_name, sig in sigs.items():
        common = df_s.index.intersection(sig.index)
        score_per_sample = df_s.loc[common, cols].multiply(sig.loc[common], axis=0).sum(axis=0)
        mean_score = score_per_sample.mean()
        pct_score  = to_pct(mean_score, ax_name)
        rows_208.append({
            "gse":      "GSE208223",
            "category": cat,
            "group":    grp_key,
            ax_name:    pct_score,
        })
        print(f"  {grp_key} / {ax_name}: raw={mean_score:.1f}, pct={pct_score:.1f}%")

# Aggregate into one row per (gse, category, group) with three axis columns
df208 = pd.DataFrame(rows_208)
pivot_208 = (df208
             .groupby(["gse", "category", "group"])[["P_axis", "N_axis", "ABA_axis"]]
             .first()
             .reset_index())
print("\nGSE208223 pivot:\n", pivot_208.to_string(index=False))

# ── Merge with existing combined table ────────────────────────────────────────
existing = pd.read_csv(ROOT / "all_biostimulants_scores.csv")
print("\nExisting table:\n", existing.to_string(index=False))

# Remove old GSE208223 rows if any
existing = existing[existing["gse"] != "GSE208223"]
combined = pd.concat([existing, pivot_208], ignore_index=True)
combined.to_csv(ROOT / "all_biostimulants_scores.csv", index=False)
print("\nUpdated combined table saved:")
print(combined.to_string(index=False))

# ── Regenerate figures ────────────────────────────────────────────────────────
PALETTE = {
    "amino_acid":        ("#E69F00", "o", "Amino acids (AA)"),
    "humic_subst":       ("#56B4E9", "o", "Humic substances"),
    "GMV_treated":       ("#009E73", "s", "PGPR volatile (GMV)"),
    "diacetyl_treated":  ("#CC79A7", "s", "PGPR diacetyl"),
    "TiO2_treated":      ("#D55E00", "^", "Ti nanoparticles"),
    "control":           ("#999999", "x", "Control"),
    "GMV_control":       ("#bbbbbb", "x", "Control (PGPR GMV)"),
    "diacetyl_control":  ("#cccccc", "x", "Control (PGPR diacetyl)"),
    "TiO2_control":      ("#aaaaaa", "x", "Control (TiO2)"),
}
CTRL_GROUPS = {"control", "GMV_control", "diacetyl_control", "TiO2_control"}

# Pair-specific control groups:
# GSE138478 has GMV_control (for GMV_treated) and diacetyl_control (for diacetyl_treated)
# We match by group prefix: GMV_* and diacetyl_* within the same gse
CTRL_MATCH = {
    "GMV_treated":       "GMV_control",
    "diacetyl_treated":  "diacetyl_control",
    "TiO2_treated":      "TiO2_control",
    "amino_acid":        "control",
    "humic_subst":       "control",
}

axis_pairs = [("P_axis", "N_axis"), ("P_axis", "ABA_axis"), ("N_axis", "ABA_axis")]

for ax_x, ax_y in axis_pairs:
    if ax_x not in combined.columns or ax_y not in combined.columns:
        continue

    fig, ax = plt.subplots(figsize=(11, 8))

    # Reference landmarks (calibration endpoints)
    ax.scatter([0], [0],   s=150, color="#f0fff0", edgecolors="#88aa88",
               marker="s", zorder=1, label="_nolegend_")
    ax.scatter([100], [100], s=150, color="#fff0f0", edgecolors="#aa8888",
               marker="s", zorder=1, label="_nolegend_")

    # Build lookup: (gse, group) -> (x, y)
    pos = {}
    for _, row in combined.iterrows():
        if pd.notna(row.get(ax_x)) and pd.notna(row.get(ax_y)):
            pos[(row["gse"], row["group"])] = (row[ax_x], row[ax_y])

    # Draw arrows from matched control to each treatment
    for _, row in combined.iterrows():
        grp = row["group"]
        if grp in CTRL_GROUPS:
            continue
        if pd.isna(row.get(ax_x)) or pd.isna(row.get(ax_y)):
            continue
        ctrl_grp = CTRL_MATCH.get(grp)
        if ctrl_grp is None:
            continue
        ctrl_key = (row["gse"], ctrl_grp)
        if ctrl_key not in pos:
            continue

        cx, cy = pos[ctrl_key]
        tx, ty = row[ax_x], row[ax_y]
        color, marker, label = PALETTE.get(grp, ("#aaaaaa", "o", grp))

        ax.annotate("", xy=(tx, ty), xytext=(cx, cy),
                    arrowprops=dict(arrowstyle="-|>", color=color,
                                   lw=2.5, mutation_scale=18),
                    zorder=4)
        off_x = 6 if tx >= cx else -6
        off_y = 6 if ty >= cy else -6
        ax.annotate(label, (tx, ty),
                    textcoords="offset points", xytext=(off_x, off_y),
                    fontsize=8, color=color, fontweight="bold", zorder=5)

    # Plot all points
    for _, row in combined.iterrows():
        grp = row["group"]
        if pd.isna(row.get(ax_x)) or pd.isna(row.get(ax_y)):
            continue
        color, marker, label = PALETTE.get(grp, ("#aaaaaa", "o", grp))
        ax.scatter([row[ax_x]], [row[ax_y]], color=color, marker=marker,
                   s=80, alpha=0.85, edgecolors="black", linewidths=0.5, zorder=3)

    # Mark control centroids
    seen_ctrl = set()
    for _, row in combined.iterrows():
        grp = row["group"]
        if grp not in CTRL_GROUPS:
            continue
        if pd.isna(row.get(ax_x)) or pd.isna(row.get(ax_y)):
            continue
        key = (row["gse"], grp)
        if key in seen_ctrl:
            continue
        seen_ctrl.add(key)
        ax.scatter([row[ax_x]], [row[ax_y]], s=130, color="#444444",
                   edgecolors="black", marker="o", zorder=4)
        ax.annotate(f"{row['gse']}\n{grp}", (row[ax_x], row[ax_y]),
                    textcoords="offset points", xytext=(4, 6),
                    fontsize=6.5, color="#333333")

    # Reference lines
    ax.axhline(0,   color="#dddddd", lw=0.8, ls="--")
    ax.axvline(0,   color="#dddddd", lw=0.8, ls="--")
    ax.axhline(100, color="#ffcccc", lw=0.8, ls="--")
    ax.axvline(100, color="#ffcccc", lw=0.8, ls="--")

    # Legend (treatments only)
    legend_handles = []
    for grp_key, (c, m, lbl) in PALETTE.items():
        if grp_key in CTRL_GROUPS:
            continue
        if any(combined["group"] == grp_key):
            legend_handles.append(mpatches.Patch(color=c, label=lbl))
    legend_handles.append(mpatches.Patch(color="#444444", label="Control centroids"))
    ax.legend(handles=legend_handles, fontsize=8.5, loc="best",
              framealpha=0.9, edgecolor="#cccccc")

    xlabel = ax_x.replace("_", " ") + " (% of biological reference range)"
    ylabel = ax_y.replace("_", " ") + " (% of biological reference range)"
    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(
        f"Biostimulant physiological state-space: {ax_x.replace('_', ' ')} vs {ax_y.replace('_', ' ')}\n"
        f"GSE297649 (AA, HS) + GSE138478 (PGPR) + GSE208223 (Ti nanoparticles)",
        fontsize=11, fontweight="bold"
    )

    fname = FIG_DIR / f"combined_state_space_{ax_x}_vs_{ax_y}.png"
    fig.tight_layout()
    fig.savefig(fname, dpi=150)
    plt.close(fig)
    print(f"  Saved: {fname.name}")

print("\nDone.")
