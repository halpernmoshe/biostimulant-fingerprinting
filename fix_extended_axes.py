"""
Fix effect_vectors_all_axes.csv:
- GMV label was already fixed (Step 1 of extend_axis_scores.py)
- Replace CPM-scored diacetyl and Pmegaterium rows with raw count-scored versions
- Regenerate heatmap
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import pandas as pd
import numpy as np
import xlrd
import gzip
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path("C:/Users/moshe/Dropbox/ISF 2025")
SIG_DIR = ROOT / "RNA_seq_axes" / "signatures"
SCORE_DIR = ROOT / "RNA_seq_axes" / "scores"
FIG_DIR = ROOT / "state_space_figures"
EV_FILE = SCORE_DIR / "effect_vectors_all_axes.csv"

# Load calibration
cal = pd.read_csv(ROOT / "RNA_seq_axes/calibration.csv", index_col=0)

# Load all signatures
all_sigs = {}
for fname in SIG_DIR.glob("*_axis_logFC.csv"):
    name = fname.stem.replace("_axis_logFC", "")
    sig = pd.read_csv(fname, index_col=0).iloc[:, 0]
    all_sigs[name] = sig
GENESET_AXES = {n for n, s in all_sigs.items() if set(s.unique()).issubset({-1.0, 0.0, 1.0})}
print(f"Loaded {len(all_sigs)} signatures; {len(GENESET_AXES)} gene-set axes")


def score_all_axes_raw(count_df, sig_dict, geneset_axes):
    """Score using log2(count+1) — same scale as calibration."""
    log2e = np.log2(count_df.astype(float) + 1)
    results = {}
    for ax_name, sig in sig_dict.items():
        common = log2e.index.intersection(sig.index)
        if ax_name in geneset_axes:
            up_genes = sig[sig > 0].index.intersection(log2e.index)
            dn_genes = sig[sig < 0].index.intersection(log2e.index)
            s = pd.Series(0.0, index=log2e.columns)
            if len(up_genes):
                s += log2e.loc[up_genes].mean(axis=0)
            if len(dn_genes):
                s -= log2e.loc[dn_genes].mean(axis=0)
        else:
            s = log2e.loc[common].multiply(sig.loc[common], axis=0).sum(axis=0)
        results[ax_name] = s
    return results


def make_delta_rows(ds, ctrl_name, treat_name, ctrl_cols, treat_cols, scores):
    rows = []
    for ax_name, ax_scores in scores.items():
        ctrl_mean  = ax_scores[ctrl_cols].mean()
        treat_mean = ax_scores[treat_cols].mean()
        rows.append({
            "dataset": ds, "control": ctrl_name, "treatment": treat_name,
            "axis": ax_name, "delta": float(treat_mean - ctrl_mean),
            "ctrl_mean": float(ctrl_mean), "treat_mean": float(treat_mean),
        })
    return rows


# ── Score diacetyl (GSE138478 HMZ014_Dione) with log2(count+1) ───────────────
print("\nScoring diacetyl (HMZ014_Dione) with log2(count+1)...")
DIONE = ROOT / "RNA_seq_data/suppl_files/GSE138478/GSE138478_HMZ014_Dione_featureCounts_output.txt.gz"
dione_df = pd.read_csv(DIONE, sep="\t", comment="#", index_col=0)
anno = ["Chr", "Start", "End", "Strand", "Length"]
dione_df = dione_df[[c for c in dione_df.columns if c not in anno]]
dione_df.columns = [Path(c).stem for c in dione_df.columns]
dione_df = dione_df.astype(int)
ctrl_d  = [c for c in dione_df.columns if "_CK_" in c]
treat_d = [c for c in dione_df.columns if "_Dione_" in c]
print(f"  ctrl: {ctrl_d}")
print(f"  treat: {treat_d}")

scores_d = score_all_axes_raw(dione_df, all_sigs, GENESET_AXES)
dione_rows = make_delta_rows(
    "GSE138478", "diacetyl_control", "diacetyl_treated",
    ctrl_d, treat_d, scores_d
)
p_d = next(r for r in dione_rows if r["axis"] == "P")
lo, hi = cal.loc["P_axis", "ctrl_score"], cal.loc["P_axis", "treat_score"]
print(f"  P-axis: ctrl={p_d['ctrl_mean']:.0f} ({(p_d['ctrl_mean']-lo)/(hi-lo)*100:.1f}%)"
      f", treat={p_d['treat_mean']:.0f} ({(p_d['treat_mean']-lo)/(hi-lo)*100:.1f}%)"
      f", delta={p_d['delta']/(hi-lo)*100:.1f}%")
n_d = next(r for r in dione_rows if r["axis"] == "N")
lon, hin = cal.loc["N_axis", "ctrl_score"], cal.loc["N_axis", "treat_score"]
print(f"  N-axis delta: {n_d['delta']/(hin-lon)*100:.1f}%")


# ── Score P. megaterium (GSE199501) with log2(raw_count+1) ───────────────────
print("\nScoring P. megaterium (GSE199501) with log2(raw_count+1)...")
wb = xlrd.open_workbook("/tmp/gse199501_tmp.xls")
sh = wb.sheet_by_index(0)
headers = sh.row_values(0)
raw_cols = ["C-1953", "C-1954", "C-1955", "T-1956", "T-1957", "T-1958"]
gene_col_idx = headers.index("Gene ID ")
raw_col_idxs = {c: headers.index(c) for c in raw_cols}

genes, data = [], {c: [] for c in raw_cols}
for r in range(1, sh.nrows):
    row = sh.row_values(r)
    gid = str(row[gene_col_idx]).strip()
    if not gid.startswith("AT"):
        continue
    genes.append(gid)
    for c in raw_cols:
        val = row[raw_col_idxs[c]]
        data[c].append(int(val) if val != "" else 0)

pmeg_df = pd.DataFrame(data, index=genes)
pmeg_df.index.name = "GeneID"
ctrl_p  = [c for c in raw_cols if c.startswith("C-")]
treat_p = [c for c in raw_cols if c.startswith("T-")]
print(f"  ctrl: {ctrl_p}, treat: {treat_p}")
print(f"  Shape: {pmeg_df.shape}")

scores_p = score_all_axes_raw(pmeg_df, all_sigs, GENESET_AXES)
pmeg_rows = make_delta_rows(
    "GSE199501", "Pmegaterium_control", "Pmegaterium_treated",
    ctrl_p, treat_p, scores_p
)
p_p = next(r for r in pmeg_rows if r["axis"] == "P")
print(f"  P-axis: ctrl={p_p['ctrl_mean']:.0f} ({(p_p['ctrl_mean']-lo)/(hi-lo)*100:.1f}%)"
      f", treat={p_p['treat_mean']:.0f} ({(p_p['treat_mean']-lo)/(hi-lo)*100:.1f}%)"
      f", delta={p_p['delta']/(hi-lo)*100:.1f}%")
n_p = next(r for r in pmeg_rows if r["axis"] == "N")
print(f"  N-axis delta: {n_p['delta']/(hin-lon)*100:.1f}%")

# ── Rebuild effect_vectors_all_axes.csv ───────────────────────────────────────
print("\nRebuilding effect_vectors_all_axes.csv...")
ev = pd.read_csv(EV_FILE)
print(f"  Existing treatments: {ev['treatment'].unique().tolist()}")

# Drop CPM-scored versions of diacetyl and Pmegaterium (added by extend_axis_scores.py)
# Keep only: amino_acid, humic_subst, GMV_treated, TiO2_treated
keep_treatments = ["amino_acid", "humic_subst", "GMV_treated", "TiO2_treated"]
ev_clean = ev[ev["treatment"].isin(keep_treatments)].copy()
print(f"  After removing CPM-scored rows: {len(ev_clean)} rows")

# Append correctly scored diacetyl and Pmegaterium
new_rows = pd.DataFrame(dione_rows + pmeg_rows)
ev_final = pd.concat([ev_clean, new_rows], ignore_index=True)
ev_final.to_csv(EV_FILE, index=False)
print(f"  Saved: {EV_FILE.name} ({len(ev_final)} rows)")
print(f"  Treatments: {ev_final['treatment'].unique().tolist()}")

# ── Regenerate heatmap ────────────────────────────────────────────────────────
print("\nRegenerating heatmap...")
CLASS_MAP = {
    "amino_acid": "AA_biostimulant",
    "humic_subst": "HA_biostimulant",
    "GMV_treated": "PGPR",
    "diacetyl_treated": "PGPR",
    "Pmegaterium_treated": "PGPR",
    "TiO2_treated": "Nanoparticle",
}
ROW_LABELS = {
    "amino_acid": "Amino acids (GSE297649)",
    "humic_subst": "Humic substances (GSE297649)",
    "GMV_treated": "PGPR volatile GMV (GSE138478)",
    "diacetyl_treated": "PGPR diacetyl (GSE138478)",
    "TiO2_treated": "Ti nanoparticles (GSE208223)",
    "Pmegaterium_treated": "P. megaterium (GSE199501)",
}
ROW_COLORS = {
    "amino_acid": "#E69F00",
    "humic_subst": "#56B4E9",
    "GMV_treated": "#009E73",
    "diacetyl_treated": "#CC79A7",
    "TiO2_treated": "#D55E00",
    "Pmegaterium_treated": "#8B5CF6",
}

pivot = ev_final[ev_final["treatment"].isin(CLASS_MAP)].pivot_table(
    index="treatment", columns="axis", values="delta"
)
pivot_z = (pivot - pivot.mean()) / (pivot.std() + 1e-9)
print(f"  Pivot shape: {pivot_z.shape} ({pivot_z.shape[0]} treatments x {pivot_z.shape[1]} axes)")

fig, ax_h = plt.subplots(figsize=(16, 6))
im = ax_h.imshow(pivot_z.values, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
ax_h.set_xticks(range(pivot_z.shape[1]))
ax_h.set_xticklabels(pivot_z.columns, rotation=45, ha="right", fontsize=8)
ax_h.set_yticks(range(pivot_z.shape[0]))
rlabels = [ROW_LABELS.get(r, r) for r in pivot_z.index]
ax_h.set_yticklabels(rlabels, fontsize=9)
for i, t in enumerate(pivot_z.index):
    ax_h.get_yticklabels()[i].set_color(ROW_COLORS.get(t, "#333333"))
    ax_h.get_yticklabels()[i].set_fontweight("bold")
plt.colorbar(im, ax=ax_h, label="Effect (z-scored across axes)")
ax_h.set_title(
    "Biostimulant effect vectors across 15 axes — all 6 treatments\n"
    "(red=activation, blue=suppression; z-scored per axis; log2(count+1) scoring)",
    fontsize=11, fontweight="bold",
)
fig.tight_layout()
fname = FIG_DIR / "biostimulant_heatmap_all_axes.png"
fig.savefig(fname, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"  Saved: {fname.name}")

print("\n=== DONE ===")
