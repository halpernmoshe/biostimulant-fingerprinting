"""
Build ET-axis (ethylene) from GSE220957.
WT_ET_rep1/2/3 vs WT_mock_rep1/2/3 (ACC treatment, 1h).
Also build a validation SA-axis (1h) and ABA-axis (1h) from same dataset.
Score all biostimulant datasets on ET-axis.
"""
import sys, gzip
sys.stdout.reconfigure(encoding='utf-8')
import pandas as pd
import numpy as np
from pathlib import Path

ROOT     = Path("C:/Users/moshe/Dropbox/ISF 2025")
SIG_DIR  = ROOT / "RNA_seq_axes/signatures"
SCORE_DIR = ROOT / "RNA_seq_axes/scores"
CAL_FILE  = ROOT / "RNA_seq_axes/calibration.csv"
DATA_FILE = ROOT / "RNA_seq_data/suppl_files/GSE220957/GSE220957_WT_mpk6_count_matrix_from_salmon.csv.gz"

# ── Load count matrix ─────────────────────────────────────────────────────────
print("Loading GSE220957 count matrix...")
df = pd.read_csv(DATA_FILE, index_col=0, compression='gzip')
# Round fractional Salmon counts to integers for DESeq2
df = df.round().astype(int)
print(f"  Shape: {df.shape}")
print(f"  Columns: {list(df.columns)}")

# Strip gene ID artifacts (AtRTD2 IDs should be standard TAIR AT*G*)
df.index = df.index.str.strip()
# Keep only standard Arabidopsis gene IDs
df = df[df.index.str.match(r'^AT[1-5MC]G\d+$')]
print(f"  After filtering to AT*G* genes: {df.shape[0]} genes")

# Sample groups
CTRL   = ['WT_mock_rep1', 'WT_mock_rep2', 'WT_mock_rep3']
ET_TR  = ['WT_ET_rep1', 'WT_ET_rep2', 'WT_ET_rep3']
SA_TR  = ['WT_SA_rep1', 'WT_SA_rep2', 'WT_SA_rep3']
ABA_TR = ['WT_ABA_rep1', 'WT_ABA_rep2', 'WT_ABA_rep3']
JA_TR  = ['WT_JA_rep1', 'WT_JA_rep2', 'WT_JA_rep3']

# ── DESeq2 helper ─────────────────────────────────────────────────────────────
try:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    def run_deseq2(ctrl_cols, treat_cols, name):
        counts_sub = df[ctrl_cols + treat_cols].T
        counts_sub = counts_sub.loc[:, counts_sub.sum() > 0]
        meta = pd.DataFrame({
            'sample': ctrl_cols + treat_cols,
            'condition': ['ctrl'] * len(ctrl_cols) + ['treat'] * len(treat_cols)
        }).set_index('sample')
        print(f"\n  {name}: {len(ctrl_cols)} ctrl + {len(treat_cols)} treat")
        print(f"  Non-zero genes: {counts_sub.shape[1]}")
        dds = DeseqDataSet(counts=counts_sub, metadata=meta,
                           design_factors="condition",
                           ref_level=["condition", "ctrl"], quiet=True)
        dds.deseq2()
        stat_res = DeseqStats(dds, contrast=["condition", "treat", "ctrl"], quiet=True)
        stat_res.summary()
        res = stat_res.results_df.dropna(subset=['log2FoldChange'])
        print(f"  Genes with logFC: {len(res)}, sig (padj<0.05): {(res.padj < 0.05).sum()}")
        return res

    print("\nRunning DESeq2 for ET-axis...")
    et_res = run_deseq2(CTRL, ET_TR, "ET-axis (ACC 1h)")

    # Validate with known ET markers
    et_markers = {'ERF1': 'AT3G23240', 'EIN3': 'AT3G20770', 'HLS1': 'AT4G37580', 'ACS5': 'AT5G65800'}
    print("\nET marker genes:")
    for name, gid in et_markers.items():
        if gid in et_res.index:
            r = et_res.loc[gid]
            print(f"  {name} ({gid}): logFC={r['log2FoldChange']:+.3f}, padj={r.get('padj', 'NA')}")

    et_sig = et_res['log2FoldChange'].rename('logFC')
    et_sig.index.name = 'GeneID'
    et_sig.to_csv(SIG_DIR / "ET_axis_logFC.csv", header=True)
    print(f"\nSaved ET_axis_logFC.csv ({len(et_sig)} genes)")

    # Optionally build SA-axis v2 for cross-validation
    print("\nRunning DESeq2 for SA-axis v2 (1h)...")
    sa2_res = run_deseq2(CTRL, SA_TR, "SA-axis v2 (1h)")
    sa2_markers = {'PR1': 'AT2G14610', 'WRKY70': 'AT3G56400', 'NPR1': 'AT1G64280'}
    print("\nSA marker genes (v2):")
    for name, gid in sa2_markers.items():
        if gid in sa2_res.index:
            r = sa2_res.loc[gid]
            print(f"  {name} ({gid}): logFC={r['log2FoldChange']:+.3f}, padj={r.get('padj', 'NA')}")
    sa2_sig = sa2_res['log2FoldChange'].rename('logFC')
    sa2_sig.index.name = 'GeneID'
    sa2_sig.to_csv(SIG_DIR / "SA_axis_v2_logFC.csv", header=True)
    print(f"Saved SA_axis_v2_logFC.csv ({len(sa2_sig)} genes)")

    # Also check JA at 1h (compare with 8h axis we built)
    print("\nRunning DESeq2 for JA-axis v2 (1h, cross-validation)...")
    ja2_res = run_deseq2(CTRL, JA_TR, "JA v2 (1h)")
    print(f"  JA 1h: {(ja2_res.padj < 0.05).sum()} sig genes")

except ImportError:
    print("pydeseq2 not available — using manual logFC")
    lib = df.sum(axis=0)
    cpm = df.divide(lib, axis=1) * 1e6

    def manual_logfc(ctrl_cols, treat_cols):
        ctrl_mean = np.log2(cpm[ctrl_cols].mean(axis=1) + 1)
        treat_mean = np.log2(cpm[treat_cols].mean(axis=1) + 1)
        return (treat_mean - ctrl_mean).rename('logFC')

    et_sig = manual_logfc(CTRL, ET_TR); et_sig.index.name = 'GeneID'
    et_sig.to_csv(SIG_DIR / "ET_axis_logFC.csv", header=True)
    print(f"Saved ET_axis_logFC.csv (manual, {len(et_sig)} genes)")
    sa2_sig = manual_logfc(CTRL, SA_TR); sa2_sig.index.name = 'GeneID'
    sa2_sig.to_csv(SIG_DIR / "SA_axis_v2_logFC.csv", header=True)

# ── Calibration ───────────────────────────────────────────────────────────────
print("\nComputing calibration scores...")

def score_on_sig(log2expr, sig):
    common = log2expr.index.intersection(sig.index)
    return log2expr.loc[common].multiply(sig.loc[common], axis=0).sum(axis=0)

df_log2 = np.log2(df.astype(float) + 1)

et_scores = score_on_sig(df_log2, et_sig)
et_ctrl  = et_scores[CTRL].mean()
et_treat = et_scores[ET_TR].mean()
print(f"ET-axis: ctrl={et_ctrl:.2f}, treat={et_treat:.2f}, range={et_treat-et_ctrl:.2f}")

sa2_scores = score_on_sig(df_log2, sa2_sig)
sa2_ctrl  = sa2_scores[CTRL].mean()
sa2_treat = sa2_scores[SA_TR].mean()
print(f"SA-axis-v2: ctrl={sa2_ctrl:.2f}, treat={sa2_treat:.2f}, range={sa2_treat-sa2_ctrl:.2f}")

# Update calibration.csv
cal = pd.read_csv(CAL_FILE, index_col=0)
cal.loc['ET_axis'] = {
    'ctrl_score': float(et_ctrl), 'treat_score': float(et_treat),
    'ctrl_label': 'WT_mock_1h', 'treat_label': 'WT_ACC_1h', 'source': 'GSE220957'
}
cal.loc['SA_axis_v2'] = {
    'ctrl_score': float(sa2_ctrl), 'treat_score': float(sa2_treat),
    'ctrl_label': 'WT_mock_1h', 'treat_label': 'WT_SA_1h', 'source': 'GSE220957'
}
cal.to_csv(CAL_FILE)
print("\nUpdated calibration.csv with ET_axis and SA_axis_v2")

def calibrate_et(score): return (float(score) - et_ctrl) / (et_treat - et_ctrl) * 100
def calibrate_sa2(score): return (float(score) - sa2_ctrl) / (sa2_treat - sa2_ctrl) * 100

# ── Score biostimulant datasets ───────────────────────────────────────────────
print("\n=== BIOSTIMULANT ET-AXIS SCORES ===")
import gzip as gz

def load_htseq_dir(gse_dir):
    files = list(Path(gse_dir).glob("*.txt.gz"))
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
    return pd.concat(dfs, axis=1).fillna(0).astype(int) if dfs else None

def load_featurecounts(fpath):
    df_fc = pd.read_csv(fpath, sep='\t', comment='#', index_col=0)
    anno = ['Chr','Start','End','Strand','Length']
    dc = [c for c in df_fc.columns if c not in anno]
    df_fc = df_fc[dc]
    df_fc.columns = [Path(c).stem for c in df_fc.columns]
    return df_fc.fillna(0).astype(int)

results_et = {}
results_sa2 = {}

def score_ds(name, count_df, treat_cols_ds, ctrl_cols_ds, results_dict, calibrate_fn, ax_label):
    log2e = np.log2(count_df.astype(float) + 1)
    scores = score_on_sig(log2e, et_sig if 'ET' in ax_label else sa2_sig)
    t = scores[treat_cols_ds].mean()
    c = scores[ctrl_cols_ds].mean()
    dt = calibrate_fn(t) - calibrate_fn(c)
    results_dict[name] = dt
    print(f"  {name}: ctrl={calibrate_fn(c):+.1f}% treat={calibrate_fn(t):+.1f}% delta={dt:+.2f}%")

# GSE297649
ds = load_htseq_dir(ROOT / "RNA_seq_data/suppl_files/GSE297649")
if ds is not None:
    aa_treat = [c for c in ds.columns if 'wtaa' in c.lower()]
    hs_treat = [c for c in ds.columns if 'wtsh' in c.lower() or 'wths' in c.lower()]
    ctrl = [c for c in ds.columns if c not in aa_treat + hs_treat]
    score_ds('amino_acid', ds, aa_treat, ctrl, results_et, calibrate_et, 'ET')
    score_ds('humic_subst', ds, hs_treat, ctrl, results_et, calibrate_et, 'ET')
    score_ds('amino_acid', ds, aa_treat, ctrl, results_sa2, calibrate_sa2, 'SA2')
    score_ds('humic_subst', ds, hs_treat, ctrl, results_sa2, calibrate_sa2, 'SA2')

# GMV volatile
gmv_df = load_featurecounts(ROOT / "RNA_seq_data/suppl_files/GSE138478/GSE138478_HMZ008_featureCounts_output.txt.gz")
if gmv_df is not None:
    gmv_ctrl  = [c for c in gmv_df.columns if '0.05C' in c or '0.5C' in c]
    gmv_treat = [c for c in gmv_df.columns if '0.05T' in c]
    if not gmv_treat:
        gmv_treat = [c for c in gmv_df.columns if c not in gmv_ctrl]
    score_ds('GMV_treated', gmv_df, gmv_treat, gmv_ctrl, results_et, calibrate_et, 'ET')
    score_ds('GMV_treated', gmv_df, gmv_treat, gmv_ctrl, results_sa2, calibrate_sa2, 'SA2')

# Diacetyl
dia_df = load_featurecounts(ROOT / "RNA_seq_data/suppl_files/GSE138478/GSE138478_HMZ014_Dione_featureCounts_output.txt.gz")
if dia_df is not None:
    dia_ctrl  = [c for c in dia_df.columns if '_CK_' in c]
    dia_treat = [c for c in dia_df.columns if '_Dione_' in c]
    score_ds('diacetyl_treated', dia_df, dia_treat, dia_ctrl, results_et, calibrate_et, 'ET')
    score_ds('diacetyl_treated', dia_df, dia_treat, dia_ctrl, results_sa2, calibrate_sa2, 'SA2')

# TiO2 (rLog)
tio2_df = pd.read_csv(ROOT / "RNA_seq_data/suppl_files/GSE208223/GSE208223_Shoot_norm.txt.gz",
                       sep='\t', index_col=0, compression='gzip')
tio2_df.index = tio2_df.index.str.replace(r'\.\d+$', '', regex=True).str.strip('"')

for ax_label, sig_use, calibrate_fn, results_dict in [
    ('ET', et_sig, calibrate_et, results_et),
    ('SA2', sa2_sig, calibrate_sa2, results_sa2)
]:
    common = tio2_df.index.intersection(sig_use.index)
    tio2_s = tio2_df.loc[common].multiply(sig_use.loc[common], axis=0).sum(axis=0)
    t2c = tio2_s[['E1','E2','E3']].mean()
    t2t = tio2_s[['F1','F2','F3']].mean()
    results_dict['TiO2_treated'] = calibrate_fn(t2t) - calibrate_fn(t2c)
    print(f"  TiO2 {ax_label}: delta={results_dict['TiO2_treated']:+.2f}%")

print("\n=== SUMMARY ===")
print(f"\n{'Treatment':<25} {'ET delta':>10} {'SA-v2 delta':>12}")
print("-" * 50)
for t in results_et:
    et_d  = results_et.get(t, np.nan)
    sa2_d = results_sa2.get(t, np.nan)
    print(f"  {t:<23} {et_d:>+10.1f}% {sa2_d:>+12.1f}%")

# Update biostimulant_6axis_summary.csv
summary_file = SCORE_DIR / "biostimulant_6axis_summary.csv"
summary = pd.read_csv(summary_file)
summary = summary[~summary['axis'].isin(['ET_axis', 'SA_axis_v2'])]

new_rows = []
for t, d in results_et.items():
    new_rows.append({'treatment': t, 'axis': 'ET_axis', 'delta_pct': d})
for t, d in results_sa2.items():
    new_rows.append({'treatment': t, 'axis': 'SA_axis_v2', 'delta_pct': d})
summary = pd.concat([summary, pd.DataFrame(new_rows)], ignore_index=True)
summary.to_csv(summary_file, index=False)
print(f"\nSaved updated summary with ET and SA-v2 axes: {summary_file}")
print("Done.")
