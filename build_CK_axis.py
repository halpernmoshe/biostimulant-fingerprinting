"""
Build CK-axis (cytokinin) from GSE149641.
Col-0 BA (benzyl adenine, 0.1 µM, 1h) vs Col-0 NaOH control.
Salmon .sf.txt.gz files with NumReads at transcript level.
Aggregate to gene level (strip .1/.2 isoform suffix), then run DESeq2.
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
DATA_DIR  = ROOT / "RNA_seq_data/suppl_files/GSE149641"

# ── Sample manifest (Col-0 only) ──────────────────────────────────────────────
SAMPLES = {
    'GSM4506435_COL.NAOH.1.sf.txt.gz': 'ctrl',
    'GSM4506436_COL.NAOH.2.sf.txt.gz': 'ctrl',
    'GSM4506437_COL.BA.1.sf.txt.gz':   'treat',
    'GSM4506438_COL.BA.2.sf.txt.gz':   'treat',
}

# ── Load and aggregate to gene level ─────────────────────────────────────────
print("Loading Salmon NumReads files...")
sample_data = {}
for fname, condition in SAMPLES.items():
    fpath = DATA_DIR / fname
    with gzip.open(fpath, 'rt') as fh:
        df_tx = pd.read_csv(fh, sep='\t', index_col=0)
    # Aggregate transcript → gene (strip isoform: AT1G01010.1 → AT1G01010)
    df_tx.index = df_tx.index.str.replace(r'\.\d+$', '', regex=True)
    gene_counts = df_tx['NumReads'].groupby(df_tx.index).sum()
    gene_counts = gene_counts.round().astype(int)
    sample_name = fname.replace('.sf.txt.gz', '').replace('GSM4506435_', '').replace('GSM4506436_', '').replace('GSM4506437_', '').replace('GSM4506438_', '')
    sample_data[fname.split('_', 1)[1].replace('.sf.txt.gz', '')] = gene_counts
    print(f"  {fname}: {len(gene_counts)} genes, total reads={gene_counts.sum():,}")

counts_df = pd.DataFrame(sample_data)
print(f"\nCount matrix: {counts_df.shape}")
print(f"Sample names: {list(counts_df.columns)}")
print(f"Library sizes: {counts_df.sum().to_dict()}")

# Remove all-zero genes
counts_df = counts_df.loc[counts_df.sum(axis=1) > 0]
print(f"After removing zero genes: {counts_df.shape}")

# ── Run DESeq2 via pydeseq2 ───────────────────────────────────────────────────
ctrl_cols  = [c for c in counts_df.columns if 'NAOH' in c]
treat_cols = [c for c in counts_df.columns if 'BA' in c]
print(f"\nControl samples: {ctrl_cols}")
print(f"Treatment samples: {treat_cols}")

all_cols = ctrl_cols + treat_cols
counts_sub = counts_df[all_cols].T  # samples x genes

meta = pd.DataFrame({
    'sample': all_cols,
    'condition': ['ctrl'] * len(ctrl_cols) + ['treat'] * len(treat_cols)
}).set_index('sample')

print(f"\nRunning DESeq2 (n={len(ctrl_cols)} ctrl, {len(treat_cols)} treat)...")

try:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

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
    print(f"DESeq2 results: {len(res)} genes with logFC")
    print(f"Significant (padj<0.05): {(res.padj < 0.05).sum()}")

    # Check key CK marker genes
    ck_markers = {
        'ARR5':  'AT3G48100',
        'ARR7':  'AT1G19050',
        'ARR15': 'AT1G74890',
        'LOG4':  'AT2G33380',
    }
    print("\nKey CK marker genes:")
    for name, gene_id in ck_markers.items():
        if gene_id in res.index:
            row = res.loc[gene_id]
            print(f"  {name} ({gene_id}): logFC={row['log2FoldChange']:+.3f}, padj={row.get('padj', 'NA')}")
        else:
            print(f"  {name} ({gene_id}): NOT IN RESULTS")

    ck_sig = res['log2FoldChange'].rename('logFC')
    ck_sig.index.name = 'GeneID'

except ImportError:
    print("pydeseq2 not available — using manual logFC")
    lib = counts_df.sum(axis=0)
    cpm = counts_df.divide(lib, axis=1) * 1e6
    ctrl_mean = np.log2(cpm[ctrl_cols].mean(axis=1) + 1)
    treat_mean = np.log2(cpm[treat_cols].mean(axis=1) + 1)
    ck_sig = (treat_mean - ctrl_mean).rename('logFC')
    ck_sig.index.name = 'GeneID'
    print(f"Manual logFC computed: {len(ck_sig)} genes")

# Save signature
sig_path = SIG_DIR / "CK_axis_logFC.csv"
ck_sig.to_csv(sig_path, header=True)
print(f"\nSaved: {sig_path} ({len(ck_sig)} genes)")

# ── Calibration ───────────────────────────────────────────────────────────────
print("\nComputing calibration scores...")

def score_on_sig(log2expr, sig):
    common = log2expr.index.intersection(sig.index)
    return log2expr.loc[common].multiply(sig.loc[common], axis=0).sum(axis=0)

df_log2 = np.log2(counts_df.astype(float) + 1)

scores = score_on_sig(df_log2, ck_sig)
ctrl_score  = scores[ctrl_cols].mean()
treat_score = scores[treat_cols].mean()
cal_range   = treat_score - ctrl_score

print(f"  Control (NaOH): {ctrl_score:.2f}")
print(f"  Treatment (BA):  {treat_score:.2f}")
print(f"  Range: {cal_range:.2f}")

# Update calibration.csv
cal = pd.read_csv(CAL_FILE, index_col=0)
cal.loc['CK_axis'] = {
    'ctrl_score': float(ctrl_score),
    'treat_score': float(treat_score),
    'ctrl_label': 'Col0_NaOH_1h',
    'treat_label': 'Col0_BA_1h',
    'source': 'GSE149641'
}
cal.to_csv(CAL_FILE)
print(f"\nUpdated calibration.csv with CK_axis")
print(cal.to_string())

# ── Score biostimulant datasets on CK-axis ────────────────────────────────────
print("\n=== BIOSTIMULANT CK-AXIS SCORES ===")

import gzip as gz

def load_htseq_dir(gse_dir, pattern="*.txt.gz"):
    files = list(Path(gse_dir).glob(pattern))
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

def calibrate_ck(score):
    return (float(score) - ctrl_score) / cal_range * 100

def score_dataset_ck(ds_name, count_df, treat_cols_ds, ctrl_cols_ds):
    log2e = np.log2(count_df.astype(float) + 1)
    scores = score_on_sig(log2e, ck_sig)
    t = scores[treat_cols_ds].mean()
    c = scores[ctrl_cols_ds].mean()
    dt = calibrate_ck(t) - calibrate_ck(c)
    print(f"  {ds_name}: ctrl={calibrate_ck(c):+.1f}% treat={calibrate_ck(t):+.1f}% delta={dt:+.2f}%")
    return dt

results_ck = {}

# GSE297649 (amino acid + humic)
ds = load_htseq_dir(ROOT / "RNA_seq_data/suppl_files/GSE297649")
if ds is not None:
    aa_treat = [c for c in ds.columns if 'wtaa' in c.lower()]
    hs_treat = [c for c in ds.columns if 'wtsh' in c.lower() or 'wths' in c.lower()]
    ctrl = [c for c in ds.columns if c not in aa_treat + hs_treat]
    results_ck['amino_acid']  = score_dataset_ck('amino_acid (GSE297649)', ds, aa_treat, ctrl)
    results_ck['humic_subst'] = score_dataset_ck('humic_subst (GSE297649)', ds, hs_treat, ctrl)

# GSE138478 GMV volatile (low P condition)
gmv_df = load_featurecounts(ROOT / "RNA_seq_data/suppl_files/GSE138478/GSE138478_HMZ008_featureCounts_output.txt.gz")
if gmv_df is not None:
    gmv_ctrl  = [c for c in gmv_df.columns if '0.05C' in c or '0.5C' in c]
    gmv_treat = [c for c in gmv_df.columns if '0.05T' in c]  # low P only
    if not gmv_treat:
        gmv_treat = [c for c in gmv_df.columns if c not in gmv_ctrl]
    results_ck['GMV_treated'] = score_dataset_ck('GMV_volatile (GSE138478)', gmv_df, gmv_treat, gmv_ctrl)

# GSE138478 diacetyl
dia_df = load_featurecounts(ROOT / "RNA_seq_data/suppl_files/GSE138478/GSE138478_HMZ014_Dione_featureCounts_output.txt.gz")
if dia_df is not None:
    dia_ctrl  = [c for c in dia_df.columns if '_CK_' in c]
    dia_treat = [c for c in dia_df.columns if '_Dione_' in c]
    results_ck['diacetyl_treated'] = score_dataset_ck('diacetyl (GSE138478)', dia_df, dia_treat, dia_ctrl)

# GSE208223 TiO2 (rLog normalized — use directly)
tio2_df = pd.read_csv(ROOT / "RNA_seq_data/suppl_files/GSE208223/GSE208223_Shoot_norm.txt.gz",
                       sep='\t', index_col=0, compression='gzip')
tio2_df.index = tio2_df.index.str.replace(r'\.\d+$', '', regex=True).str.strip('"')
# For TiO2, use rLog directly (already log2-scale)
common_tio2 = tio2_df.index.intersection(ck_sig.index)
tio2_scores_raw = tio2_df.loc[common_tio2].multiply(ck_sig.loc[common_tio2], axis=0).sum(axis=0)
t2_ctrl  = tio2_scores_raw[['E1','E2','E3']].mean()
t2_treat = tio2_scores_raw[['F1','F2','F3']].mean()
results_ck['TiO2_treated'] = calibrate_ck(t2_treat) - calibrate_ck(t2_ctrl)
print(f"  TiO2 (GSE208223): ctrl={calibrate_ck(t2_ctrl):+.1f}% treat={calibrate_ck(t2_treat):+.1f}% delta={results_ck['TiO2_treated']:+.2f}%")

print("\n=== SUMMARY: CK-axis deltas ===")
for t, d in results_ck.items():
    print(f"  {t:<25}: {d:+.2f}%")

# Update biostimulant_6axis_summary.csv
summary_file = SCORE_DIR / "biostimulant_6axis_summary.csv"
summary = pd.read_csv(summary_file)
summary = summary[summary['axis'] != 'CK_axis']  # remove old CK rows if any

new_rows = [{'treatment': t, 'axis': 'CK_axis', 'delta_pct': d} for t, d in results_ck.items()]
summary = pd.concat([summary, pd.DataFrame(new_rows)], ignore_index=True)
summary.to_csv(summary_file, index=False)
print(f"\nSaved updated 6-axis summary (now with CK): {summary_file}")
print("Done.")
