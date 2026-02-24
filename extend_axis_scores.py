"""
Extend effect_vectors_all_axes.csv with missing treatments:
1. Fix mislabeled GMV (was called 'diacetyl' in build_geneset_axes.py — column prefix bug)
2. Add true diacetyl (GSE138478 HMZ014_Dione file, CK vs Dione)
3. Add P. megaterium PGPR (GSE199501 CPM data)
Then regenerate heatmap and axis_selection.md.
"""
import gzip, sys
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT     = Path("C:/Users/moshe/Dropbox/ISF 2025")
SIG_DIR  = ROOT / "RNA_seq_axes" / "signatures"
SCORE_DIR = ROOT / "RNA_seq_axes" / "scores"
FIG_DIR  = ROOT / "state_space_figures"
EV_FILE  = SCORE_DIR / "effect_vectors_all_axes.csv"

# ── Load all 15 axis signatures ────────────────────────────────────────────────
print("Loading axis signatures...")
all_sigs = {}
for fname in SIG_DIR.glob("*_axis_logFC.csv"):
    name = fname.stem.replace("_axis_logFC", "")
    sig = pd.read_csv(fname, index_col=0).iloc[:, 0]
    sig.name = name
    all_sigs[name] = sig
print(f"  Loaded {len(all_sigs)} signatures: {sorted(all_sigs.keys())}")

# Gene-set axes (from GO annotations) — identified by presence of signatures with logFC = +/-1
GENESET_AXES = {name for name, sig in all_sigs.items() if set(sig.unique()).issubset({-1.0, 0.0, 1.0})}
print(f"  Gene-set axes ({len(GENESET_AXES)}): {sorted(GENESET_AXES)}")
DESEQ2_AXES = set(all_sigs.keys()) - GENESET_AXES
print(f"  DESeq2 axes ({len(DESEQ2_AXES)}): {sorted(DESEQ2_AXES)}")

# ── Scoring functions ──────────────────────────────────────────────────────────
def log2p1(df):
    return np.log2(df.astype(float) + 1)

def to_log2_expr(df, is_cpm=False):
    """Convert to log2 expression. If is_cpm, data is already CPM."""
    if is_cpm:
        return log2p1(df)
    # Raw counts: CPM normalize then log2(+1)
    lib = df.sum(axis=0)
    cpm = df.divide(lib, axis=1) * 1_000_000
    return log2p1(cpm)

def score_deseq2(log2expr, sig):
    """Dot-product: sum(log2expr * logFC) for common genes."""
    common = log2expr.index.intersection(sig.index)
    return log2expr.loc[common].multiply(sig.loc[common], axis=0).sum(axis=0)

def score_geneset(log2expr, sig):
    """Gene-set: mean(expr for +1 genes) - mean(expr for -1 genes)."""
    up_genes   = sig[sig > 0].index.intersection(log2expr.index)
    down_genes = sig[sig < 0].index.intersection(log2expr.index)
    score = pd.Series(0.0, index=log2expr.columns)
    if len(up_genes) > 0:
        score += log2expr.loc[up_genes].mean(axis=0)
    if len(down_genes) > 0:
        score -= log2expr.loc[down_genes].mean(axis=0)
    return score

def score_all_axes(log2expr, sig_dict, geneset_axes):
    """Score a log2 expression matrix on all axes. Returns dict of {axis: Series}."""
    results = {}
    for ax_name, sig in sig_dict.items():
        if ax_name in geneset_axes:
            results[ax_name] = score_geneset(log2expr, sig)
        else:
            results[ax_name] = score_deseq2(log2expr, sig)
    return results

def compute_delta_rows(ds_name, ctrl_name, treat_name, ctrl_cols, treat_cols, log2expr, sig_dict, geneset_axes):
    """Compute per-axis deltas for one treatment vs control."""
    scores = score_all_axes(log2expr, sig_dict, geneset_axes)
    rows = []
    for ax_name, ax_scores in scores.items():
        ctrl_mean  = ax_scores[ctrl_cols].mean()
        treat_mean = ax_scores[treat_cols].mean()
        delta = treat_mean - ctrl_mean
        rows.append({
            'dataset': ds_name,
            'control': ctrl_name,
            'treatment': treat_name,
            'axis': ax_name,
            'delta': float(delta),
            'ctrl_mean': float(ctrl_mean),
            'treat_mean': float(treat_mean),
        })
    return rows

# ── Step 1: Fix mislabeled GMV in existing effect_vectors_all_axes.csv ─────────
print("\nStep 1: Fixing GMV mislabeling...")
ev = pd.read_csv(EV_FILE)
print(f"  Before: treatments = {ev['treatment'].unique().tolist()}")

# In build_geneset_axes.py, HMZ008 columns were misclassified as 'diacetyl_*'
# because the code checked column prefix for 'GMV' but columns start with 'HMZ008'.
# The 'diacetyl_treated' rows in effect_vectors_all_axes.csv are actually GMV volatile data.
ev['control']   = ev['control'].replace('diacetyl_control', 'GMV_control')
ev['treatment'] = ev['treatment'].replace('diacetyl_treated', 'GMV_treated')
print(f"  After:  treatments = {ev['treatment'].unique().tolist()}")

# ── Step 2: Score true diacetyl (GSE138478 HMZ014_Dione) ─────────────────────
print("\nStep 2: Scoring GSE138478 diacetyl (HMZ014_Dione)...")
DIONE_FILE = ROOT / "RNA_seq_data/suppl_files/GSE138478/GSE138478_HMZ014_Dione_featureCounts_output.txt.gz"

dione_df = pd.read_csv(DIONE_FILE, sep='\t', comment='#', index_col=0)
anno_cols = ['Chr', 'Start', 'End', 'Strand', 'Length']
dione_data_cols = [c for c in dione_df.columns if c not in anno_cols]
dione_df = dione_df[dione_data_cols]
# Clean column names (strip path and .bam)
dione_df.columns = [Path(c).stem for c in dione_df.columns]
print(f"  Columns: {dione_df.columns.tolist()}")
print(f"  Shape: {dione_df.shape}")

# Group assignment from column names:
# HMZ014_CK_2DAT_rep1_th2 → control; HMZ014_Dione_2DAT_rep1_th2 → treated
dione_ctrl  = [c for c in dione_df.columns if '_CK_' in c]
dione_treat = [c for c in dione_df.columns if '_Dione_' in c]
print(f"  Control: {dione_ctrl}")
print(f"  Treated: {dione_treat}")

# CPM normalize (raw counts)
dione_log2 = to_log2_expr(dione_df, is_cpm=False)

dione_rows = compute_delta_rows(
    ds_name='GSE138478', ctrl_name='diacetyl_control', treat_name='diacetyl_treated',
    ctrl_cols=dione_ctrl, treat_cols=dione_treat,
    log2expr=dione_log2, sig_dict=all_sigs, geneset_axes=GENESET_AXES
)
print(f"  Computed {len(dione_rows)} axis scores")
# Spot-check P-axis
p_row = next(r for r in dione_rows if r['axis'] == 'P')
print(f"  P-axis delta (diacetyl): ctrl={p_row['ctrl_mean']:.0f}, treat={p_row['treat_mean']:.0f}, delta={p_row['delta']:.0f}")

# ── Step 3: Score GSE199501 P. megaterium (CPM data) ─────────────────────────
print("\nStep 3: Scoring GSE199501 P. megaterium...")
try:
    import xlrd
    XLS_GZ = ROOT / "RNA_seq_data/suppl_files/GSE199501/GSE199501_Processed_data_DEG_total_CPM.xls.gz"
    with gzip.open(XLS_GZ, 'rb') as f:
        xls_bytes = f.read()
    with open('/tmp/gse199501_tmp.xls', 'wb') as f:
        f.write(xls_bytes)

    wb = xlrd.open_workbook('/tmp/gse199501_tmp.xls')
    sh = wb.sheet_by_index(0)
    headers = sh.row_values(0)
    print(f"  Headers: {headers}")

    # CPM columns are 'C-1953 CPM', 'C-1954 CPM', 'C-1955 CPM', 'T-1956 CPM', ...
    cpm_cols = [h for h in headers if 'CPM' in str(h)]
    gene_col = headers.index('Gene ID ')
    col_indices = {h: headers.index(h) for h in cpm_cols + ['Gene ID ']}

    genes, data = [], {c: [] for c in cpm_cols}
    for row_idx in range(1, sh.nrows):
        row = sh.row_values(row_idx)
        gene_id = str(row[col_indices['Gene ID ']]).strip()
        if not gene_id.startswith('AT'):
            continue
        genes.append(gene_id)
        for c in cpm_cols:
            data[c].append(float(row[col_indices[c]]) if row[col_indices[c]] != '' else 0.0)

    gse199_cpm = pd.DataFrame(data, index=genes)
    gse199_cpm.index.name = 'GeneID'
    print(f"  Shape: {gse199_cpm.shape}")
    print(f"  Columns: {gse199_cpm.columns.tolist()}")

    # Convert CPM → log2(CPM+1)
    gse199_log2 = to_log2_expr(gse199_cpm, is_cpm=True)

    ctrl_cpm_cols  = [c for c in cpm_cols if c.startswith('C-')]
    treat_cpm_cols = [c for c in cpm_cols if c.startswith('T-')]
    print(f"  Control CPM cols: {ctrl_cpm_cols}")
    print(f"  Treated CPM cols: {treat_cpm_cols}")

    pmeg_rows = compute_delta_rows(
        ds_name='GSE199501', ctrl_name='Pmegaterium_control', treat_name='Pmegaterium_treated',
        ctrl_cols=ctrl_cpm_cols, treat_cols=treat_cpm_cols,
        log2expr=gse199_log2, sig_dict=all_sigs, geneset_axes=GENESET_AXES
    )
    print(f"  Computed {len(pmeg_rows)} axis scores")
    # Spot-check P-axis
    p_row = next(r for r in pmeg_rows if r['axis'] == 'P')
    print(f"  P-axis delta (P. megaterium): ctrl={p_row['ctrl_mean']:.0f}, treat={p_row['treat_mean']:.0f}, delta={p_row['delta']:.0f}")

except ImportError:
    print("  ERROR: xlrd not installed. Cannot read XLS file.")
    pmeg_rows = []

# ── Step 4: Combine and save updated effect_vectors_all_axes.csv ──────────────
print("\nStep 4: Combining and saving...")
new_df = pd.DataFrame(dione_rows + pmeg_rows)
ev_updated = pd.concat([ev, new_df], ignore_index=True)
ev_updated.to_csv(EV_FILE, index=False)
print(f"  Saved: {EV_FILE.name} ({len(ev_updated)} rows)")
print(f"  Treatments: {ev_updated['treatment'].unique().tolist()}")

# ── Step 5: Update axis_selection.md with CLASS_MAP including all treatments ───
print("\nStep 5: Recomputing axis F-statistics...")
from scipy import stats

# Load per-sample scores (for F-statistics we need per-sample not just deltas)
# We'll use the deltas as proxy treatment means (limited but consistent with before)
CLASS_MAP = {
    'amino_acid': 'AA_biostimulant',
    'humic_subst': 'HA_biostimulant',
    'GMV_treated': 'PGPR',
    'diacetyl_treated': 'PGPR',
    'Pmegaterium_treated': 'PGPR',
    'TiO2_treated': 'Nanoparticle',
}

# Use treat_mean as proxy for within-class values (n=1 per treatment at this level)
# F-stat across treatment groups using delta values
treat_data = ev_updated[ev_updated['treatment'].isin(CLASS_MAP.keys())].copy()
treat_data['class'] = treat_data['treatment'].map(CLASS_MAP)

f_stats = {}
for ax_name in all_sigs.keys():
    ax_data = treat_data[treat_data['axis'] == ax_name]
    if ax_data.empty:
        continue
    classes = ax_data['class'].unique()
    if len(classes) < 2:
        continue
    groups_list = [ax_data[ax_data['class'] == c]['delta'].values for c in classes]
    groups_list = [g for g in groups_list if len(g) > 0]
    if len(groups_list) < 2:
        continue
    try:
        f, p = stats.f_oneway(*groups_list)
        f_stats[ax_name] = {'F': float(f), 'p': float(p), 'n_classes': len(classes)}
    except Exception:
        pass

f_df = pd.DataFrame(f_stats).T.sort_values('F', ascending=False)
print("\n  Axis F-statistics (updated):")
print(f_df.round(4).to_string())

axis_md_lines = [
    "# Axis Selection Analysis\n\n",
    "## F-statistic Ranking (updated with all 6 biostimulant treatments)\n\n",
    "Classes: AA biostimulant, HA biostimulant, PGPR (GMV+diacetyl+P.megaterium), Nanoparticle\n\n",
    f_df.round(4).to_markdown(), "\n\n",
    "## Interpretation\n",
    "High F = axis separates biostimulant classes well.\n\n",
    "### Top discriminative axes:\n",
]
for ax in f_df.index[:5]:
    row = f_df.loc[ax]
    axis_md_lines.append(f"- **{ax}**: F={row['F']:.2f}, p={row['p']:.4f}\n")

axis_md_lines.append("\n### Low-priority axes:\n")
for ax in f_df.index[-5:]:
    row = f_df.loc[ax]
    axis_md_lines.append(f"- **{ax}**: F={row['F']:.2f} (low discriminative power)\n")

with open(ROOT / "agent_coordination" / "axis_selection.md", 'w') as f:
    f.writelines(axis_md_lines)
print("  Saved: agent_coordination/axis_selection.md")

# ── Step 6: Regenerate heatmap ────────────────────────────────────────────────
print("\nStep 6: Regenerating heatmap...")
try:
    pivot = ev_updated[ev_updated['treatment'].isin(CLASS_MAP.keys())].pivot_table(
        index='treatment', columns='axis', values='delta'
    )
    # Normalize columns to z-scores for visualization
    pivot_z = (pivot - pivot.mean()) / (pivot.std() + 1e-9)

    # Nice row labels
    ROW_LABELS = {
        'amino_acid': 'Amino acids\n(GSE297649)',
        'humic_subst': 'Humic substances\n(GSE297649)',
        'GMV_treated': 'PGPR volatile GMV\n(GSE138478)',
        'diacetyl_treated': 'PGPR diacetyl\n(GSE138478)',
        'TiO2_treated': 'Ti nanoparticles\n(GSE208223)',
        'Pmegaterium_treated': 'P. megaterium\n(GSE199501)',
    }
    row_labels = [ROW_LABELS.get(r, r) for r in pivot_z.index]

    # Row color bands for biostimulant class
    ROW_COLORS = {
        'amino_acid': '#E69F00',
        'humic_subst': '#56B4E9',
        'GMV_treated': '#009E73',
        'diacetyl_treated': '#CC79A7',
        'TiO2_treated': '#D55E00',
        'Pmegaterium_treated': '#8B5CF6',
    }

    n_rows = pivot_z.shape[0]
    n_cols = pivot_z.shape[1]
    fig_w = max(14, n_cols * 0.85)
    fig_h = max(6, n_rows * 0.9)

    fig, (ax_heat, ax_cbar_ax) = plt.subplots(
        1, 2, figsize=(fig_w, fig_h),
        gridspec_kw={'width_ratios': [20, 1]}
    )

    im = ax_heat.imshow(pivot_z.values, aspect='auto', cmap='RdBu_r', vmin=-2.5, vmax=2.5)
    ax_heat.set_xticks(range(n_cols))
    ax_heat.set_xticklabels(pivot_z.columns, rotation=45, ha='right', fontsize=8)
    ax_heat.set_yticks(range(n_rows))
    ax_heat.set_yticklabels(row_labels, fontsize=9)

    # Colored row labels
    for i, treatment in enumerate(pivot_z.index):
        color = ROW_COLORS.get(treatment, '#333333')
        ax_heat.get_yticklabels()[i].set_color(color)
        ax_heat.get_yticklabels()[i].set_fontweight('bold')

    plt.colorbar(im, cax=ax_cbar_ax, label='Effect size (z-scored across axes)')
    ax_heat.set_title(
        "Biostimulant effect vectors across 15 axes\n"
        "(red = positive activation, blue = suppression; z-scored per axis)\n"
        "All 6 treatments: CPM-normalized counts",
        fontsize=11, fontweight='bold'
    )

    fig.tight_layout()
    fname = FIG_DIR / "biostimulant_heatmap_all_axes.png"
    fig.savefig(fname, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {fname.name}")

except Exception as e:
    print(f"  Heatmap error: {e}")
    import traceback
    traceback.print_exc()

print("\n=== ALL DONE ===")
print(f"Updated treatments: {ev_updated['treatment'].unique().tolist()}")
