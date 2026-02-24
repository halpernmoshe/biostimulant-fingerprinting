"""
Build JA-axis and SA-axis from GSE90077 (Arabidopsis MeJA / BTH treatment).
JA: MeJA 8h vs Mock 8h (peak JA response timepoint)
SA: BTH 8h vs Mock 8h (BTH = benzothiadiazole, functional SA analog)
Score biostimulant datasets on both axes.
"""
import sys
sys.stdout.reconfigure(encoding='utf-8')
import gzip
import pandas as pd
import numpy as np
from pathlib import Path

ROOT     = Path("C:/Users/moshe/Dropbox/ISF 2025")
SIG_DIR  = ROOT / "RNA_seq_axes" / "signatures"
SCORE_DIR = ROOT / "RNA_seq_axes" / "scores"
CAL_FILE = ROOT / "RNA_seq_axes" / "calibration.csv"
GSE90077 = ROOT / "RNA_seq_data/suppl_files/GSE90077/GSE90077_JBM_counts.txt.gz"

# ── Load count matrix ─────────────────────────────────────────────────────────
print("Loading GSE90077...")
with gzip.open(GSE90077, 'rt') as fh:
    labels = fh.readline().strip().split('\t')
    df = pd.read_csv(fh, sep='\t', index_col=0)

# Map sample IDs to treatment labels
label_map = dict(zip(df.columns, labels))
print(f"  {df.shape[0]} genes x {df.shape[1]} samples")
print(f"  Treatments: {sorted(set(l.rsplit('_', 1)[0] for l in labels))}")

# Check library sizes
lib = df.sum(axis=0)
print(f"  Library size range: {lib.min():,} - {lib.max():,}")

# ── Run DESeq2 via pydeseq2 ───────────────────────────────────────────────────
try:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    print("\nRunning DESeq2 via pydeseq2...")

    def run_deseq2(ctrl_cols, treat_cols, contrast_name):
        """Run DESeq2 on a subset of samples."""
        all_cols = ctrl_cols + treat_cols
        counts_sub = df[all_cols].T  # samples x genes
        # Remove genes with all zeros
        counts_sub = counts_sub.loc[:, counts_sub.sum() > 0]

        # Build metadata
        meta = pd.DataFrame({
            'sample': all_cols,
            'condition': ['ctrl'] * len(ctrl_cols) + ['treat'] * len(treat_cols)
        }).set_index('sample')

        print(f"\n  {contrast_name}: ctrl={ctrl_cols}, treat={treat_cols}")
        print(f"  Non-zero genes: {counts_sub.shape[1]}")

        dds = DeseqDataSet(
            counts=counts_sub,
            metadata=meta,
            design_factors="condition",
            ref_level=["condition", "ctrl"],
            quiet=True,
        )
        dds.deseq2()
        stat_res = DeseqStats(dds, contrast=["condition", "treat", "ctrl"], quiet=True)
        stat_res.summary()
        res = stat_res.results_df
        res = res.dropna(subset=['log2FoldChange'])
        print(f"  DESeq2 results: {len(res)} genes with logFC")
        print(f"  Sig (padj<0.05): {(res.padj < 0.05).sum()}")
        return res

    # Identify column IDs for each condition
    def get_cols(pattern):
        return [col for col, lab in label_map.items() if lab.startswith(pattern)]

    # JA-axis: MeJA_8h vs Mock_8h
    ja_res = run_deseq2(
        ctrl_cols=get_cols('Mock_8h'),
        treat_cols=get_cols('MeJA_8h'),
        contrast_name='JA-axis (MeJA_8h vs Mock_8h)'
    )
    ja_sig = ja_res['log2FoldChange'].rename('logFC')
    ja_sig.index.name = 'GeneID'
    ja_sig.to_csv(SIG_DIR / "JA_axis_logFC.csv", header=True)
    print(f"  Saved: JA_axis_logFC.csv ({len(ja_sig)} genes)")

    # SA-axis: BTH_8h vs Mock_8h
    sa_res = run_deseq2(
        ctrl_cols=get_cols('Mock_8h'),
        treat_cols=get_cols('BTH_8h'),
        contrast_name='SA-axis (BTH_8h vs Mock_8h)'
    )
    sa_sig = sa_res['log2FoldChange'].rename('logFC')
    sa_sig.index.name = 'GeneID'
    sa_sig.to_csv(SIG_DIR / "SA_axis_logFC.csv", header=True)
    print(f"  Saved: SA_axis_logFC.csv ({len(sa_sig)} genes)")

except ImportError:
    print("  pydeseq2 not available — using manual logFC (log2(CPM_treat+1) - log2(CPM_ctrl+1))")

    def manual_logfc(ctrl_cols, treat_cols):
        """Simple log2(mean CPM treat / mean CPM ctrl) as logFC proxy."""
        lib_loc = df[ctrl_cols + treat_cols].sum(axis=0)
        cpm = df[ctrl_cols + treat_cols].divide(lib_loc, axis=1) * 1e6
        ctrl_mean = np.log2(cpm[ctrl_cols].mean(axis=1) + 1)
        treat_mean = np.log2(cpm[treat_cols].mean(axis=1) + 1)
        return treat_mean - ctrl_mean

    def get_cols(pattern):
        return [col for col, lab in label_map.items() if lab.startswith(pattern)]

    ja_logfc = manual_logfc(get_cols('Mock_8h'), get_cols('MeJA_8h'))
    ja_logfc.name = 'logFC'; ja_logfc.index.name = 'GeneID'
    ja_logfc.to_csv(SIG_DIR / "JA_axis_logFC.csv", header=True)
    print(f"  Saved: JA_axis_logFC.csv (manual logFC, {len(ja_logfc)} genes)")

    sa_logfc = manual_logfc(get_cols('Mock_8h'), get_cols('BTH_8h'))
    sa_logfc.name = 'logFC'; sa_logfc.index.name = 'GeneID'
    sa_logfc.to_csv(SIG_DIR / "SA_axis_logFC.csv", header=True)
    print(f"  Saved: SA_axis_logFC.csv (manual logFC, {len(sa_logfc)} genes)")

    ja_sig = ja_logfc
    sa_sig = sa_logfc

# ── Score biostimulant datasets on JA and SA axes ────────────────────────────
print("\nScoring biostimulant datasets on JA and SA axes...")

def score_on_sig(log2expr, sig):
    """Dot-product score: sum(log2(count+1) * logFC) for common genes."""
    common = log2expr.index.intersection(sig.index)
    return log2expr.loc[common].multiply(sig.loc[common], axis=0).sum(axis=0)

def cpm_log2(count_df):
    lib = count_df.sum(axis=0)
    return np.log2(count_df.divide(lib, axis=1) * 1e6 + 1)

# Calibration: use Mock_8h vs MeJA_8h scores ON the GSE90077 data itself
# (same approach as P-axis: ctrl=ctrl condition score, treat=treat condition score)
mock8h_cols = [col for col, lab in label_map.items() if lab.startswith('Mock_8h')]
meja8h_cols = [col for col, lab in label_map.items() if lab.startswith('MeJA_8h')]
bth8h_cols  = [col for col, lab in label_map.items() if lab.startswith('BTH_8h')]

df_log2 = np.log2(df.astype(float) + 1)  # all GSE90077 samples, raw counts

ja_scores_cal = score_on_sig(df_log2, ja_sig)
sa_scores_cal = score_on_sig(df_log2, sa_sig)

ja_ctrl_score  = ja_scores_cal[mock8h_cols].mean()
ja_treat_score = ja_scores_cal[meja8h_cols].mean()
sa_ctrl_score  = sa_scores_cal[mock8h_cols].mean()
sa_treat_score = sa_scores_cal[bth8h_cols].mean()

print(f"\nCalibration scores:")
print(f"  JA-axis: Mock_8h={ja_ctrl_score:.1f}  MeJA_8h={ja_treat_score:.1f}  range={ja_treat_score-ja_ctrl_score:.1f}")
print(f"  SA-axis: Mock_8h={sa_ctrl_score:.1f}  BTH_8h={sa_treat_score:.1f}   range={sa_treat_score-sa_ctrl_score:.1f}")

# Update calibration.csv
cal = pd.read_csv(CAL_FILE, index_col=0)
for ax_name, ctrl_s, treat_s, ctrl_lab, treat_lab in [
    ('JA_axis', ja_ctrl_score, ja_treat_score, 'Mock_8h', 'MeJA_8h'),
    ('SA_axis', sa_ctrl_score, sa_treat_score, 'Mock_8h', 'BTH_8h'),
]:
    if ax_name not in cal.index:
        cal.loc[ax_name] = {
            'ctrl_score': float(ctrl_s), 'treat_score': float(treat_s),
            'ctrl_label': ctrl_lab, 'treat_label': treat_lab, 'source': 'GSE90077'
        }
    else:
        cal.loc[ax_name, ['ctrl_score','treat_score']] = [float(ctrl_s), float(treat_s)]
cal.to_csv(CAL_FILE)
print(f"\nUpdated calibration.csv with JA_axis and SA_axis")

def calibrate(score, ax_name):
    lo = cal.loc[ax_name, 'ctrl_score']
    hi = cal.loc[ax_name, 'treat_score']
    return (float(score) - lo) / (hi - lo) * 100

# ── Score each biostimulant dataset ──────────────────────────────────────────
import gzip as gz

def load_htseq_dir(gse_dir, pattern="*.txt.gz"):
    from pathlib import Path as P
    files = list(P(gse_dir).glob(pattern))
    dfs = []
    for f in files:
        try:
            with gz.open(f, 'rt') as fh:
                lines = [l.strip().split('\t') for l in fh if not l.startswith('__')]
            genes  = [l[0] for l in lines if len(l) == 2]
            counts = [int(l[1]) for l in lines if len(l) == 2]
            if genes:
                dfs.append(pd.Series(counts, index=genes, name=f.stem.split('.')[0]))
        except Exception:
            continue
    if not dfs:
        return None
    return pd.concat(dfs, axis=1).fillna(0).astype(int)

def load_featurecounts(fpath):
    df_fc = pd.read_csv(fpath, sep='\t', comment='#', index_col=0)
    anno = ['Chr','Start','End','Strand','Length']
    dc = [c for c in df_fc.columns if c not in anno]
    df_fc = df_fc[dc]
    df_fc.columns = [Path(c).stem for c in df_fc.columns]
    return df_fc.fillna(0).astype(int)

datasets = {
    'GSE297649': {
        'load': lambda: load_htseq_dir(ROOT / "RNA_seq_data/suppl_files/GSE297649"),
        'groups': lambda df: {c: ('amino_acid' if ('wtaa' in c.lower()) else
                                  ('humic_subst' if ('wtsh' in c.lower() or 'wths' in c.lower()) else 'control'))
                              for c in df.columns}
    },
    'GSE138478_GMV': {
        'load': lambda: load_featurecounts(ROOT / "RNA_seq_data/suppl_files/GSE138478/GSE138478_HMZ008_featureCounts_output.txt.gz"),
        'groups': lambda df: {c: ('GMV_control' if '0.05C' in c or '0.5C' in c else 'GMV_treated')
                              for c in df.columns}
    },
    'GSE138478_diacetyl': {
        'load': lambda: load_featurecounts(ROOT / "RNA_seq_data/suppl_files/GSE138478/GSE138478_HMZ014_Dione_featureCounts_output.txt.gz"),
        'groups': lambda df: {c: ('diacetyl_control' if '_CK_' in c else 'diacetyl_treated')
                              for c in df.columns}
    },
}

print("\n=== BIOSTIMULANT SCORES ON JA AND SA AXES ===")
results = []
for ds_name, cfg in datasets.items():
    df_ds = cfg['load']()
    if df_ds is None:
        print(f"  {ds_name}: NOT LOADED")
        continue
    groups = cfg['groups'](df_ds)
    # Use log2(count+1) for consistency with calibration
    df_log2_ds = np.log2(df_ds.astype(float) + 1)
    for sig, ax_name, cal_name in [(ja_sig, 'JA_axis', 'JA_axis'), (sa_sig, 'SA_axis', 'SA_axis')]:
        scores = score_on_sig(df_log2_ds, sig)
        # Group means
        group_means = {}
        for col, score in scores.items():
            grp = groups.get(col, 'unknown')
            group_means.setdefault(grp, []).append(score)
        for grp, vals in group_means.items():
            mean_score = np.mean(vals)
            pct = calibrate(mean_score, cal_name)
            results.append({'dataset': ds_name, 'group': grp, 'axis': ax_name,
                            'score': mean_score, 'pct': pct, 'n': len(vals)})

res_df = pd.DataFrame(results)

print(f"\n{'Dataset':<22} {'Group':<22} {'JA (%)':>8} {'SA (%)':>8}")
print("-" * 65)
datasets_shown = set()
for _, r in res_df.iterrows():
    key = (r.dataset, r.group)
    ja_row = res_df[(res_df.dataset==r.dataset) & (res_df.group==r.group) & (res_df.axis=='JA_axis')]
    sa_row = res_df[(res_df.dataset==r.dataset) & (res_df.group==r.group) & (res_df.axis=='SA_axis')]
    if key not in datasets_shown and not ja_row.empty and not sa_row.empty:
        datasets_shown.add(key)
        ja_pct = ja_row.iloc[0]['pct']
        sa_pct = sa_row.iloc[0]['pct']
        print(f"  {r.dataset:<20} {r.group:<22} {ja_pct:>+8.1f}% {sa_pct:>+8.1f}%")

# Compute effect vectors (delta = treated - control)
print("\n=== EFFECT VECTORS (treated - control) ===")
CTRL_MAP = {
    'GSE297649': 'control', 'GSE138478_GMV': 'GMV_control', 'GSE138478_diacetyl': 'diacetyl_control'
}
TREAT_MAP = {
    'GSE297649_amino_acid': ('GSE297649', 'amino_acid'),
    'GSE297649_humic_subst': ('GSE297649', 'humic_subst'),
    'GSE138478_GMV_treated': ('GSE138478_GMV', 'GMV_treated'),
    'GSE138478_diacetyl_treated': ('GSE138478_diacetyl', 'diacetyl_treated'),
}

print(f"\n{'Treatment':<30} {'JA delta (%)':>14} {'SA delta (%)':>14}")
print("-" * 65)
for treat_label, (ds_name, treat_grp) in TREAT_MAP.items():
    ctrl_grp = CTRL_MAP.get(ds_name)
    for ax_name in ['JA_axis', 'SA_axis']:
        ctrl_row  = res_df[(res_df.dataset==ds_name) & (res_df.group==ctrl_grp) & (res_df.axis==ax_name)]
        treat_row = res_df[(res_df.dataset==ds_name) & (res_df.group==treat_grp) & (res_df.axis==ax_name)]
        if ctrl_row.empty or treat_row.empty:
            continue

    ja_ctrl  = res_df[(res_df.dataset==ds_name) & (res_df.group==ctrl_grp) & (res_df.axis=='JA_axis')].iloc[0]['pct']
    ja_treat = res_df[(res_df.dataset==ds_name) & (res_df.group==treat_grp) & (res_df.axis=='JA_axis')].iloc[0]['pct']
    sa_ctrl  = res_df[(res_df.dataset==ds_name) & (res_df.group==ctrl_grp) & (res_df.axis=='SA_axis')].iloc[0]['pct']
    sa_treat = res_df[(res_df.dataset==ds_name) & (res_df.group==treat_grp) & (res_df.axis=='SA_axis')].iloc[0]['pct']

    ja_delta = ja_treat - ja_ctrl
    sa_delta = sa_treat - sa_ctrl
    print(f"  {treat_label:<30} {ja_delta:>+14.1f}% {sa_delta:>+14.1f}%")

res_df.to_csv(SCORE_DIR / "GSE90077_JA_SA_axis_scores.csv", index=False)
print(f"\nSaved: GSE90077_JA_SA_axis_scores.csv")
print("\n=== DONE ===")
