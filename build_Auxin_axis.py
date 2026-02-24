"""
Build Auxin-axis from GSE104385.
Col-0 IAA 10µM 2h vs Col-0 EtOH mock (n=4 each).
Files: 9-column format with Count column (col index 6, header 'Count').
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
DATA_DIR  = ROOT / "RNA_seq_data/suppl_files/GSE104385"

# ── Sample manifest (Col-0 only) ──────────────────────────────────────────────
CTRL_FILES  = [
    "GSM2796785_col_1.ACTGTCG.gene_counts.txt.gz",
    "GSM2796786_col_2.GTTCTCA.gene_counts.txt.gz",
    "GSM2796787_col_3.TTGCCCC.gene_counts.txt.gz",
    "GSM2796788_col_4.TAGGTAA.gene_counts.txt.gz",
]
TREAT_FILES = [
    "GSM2796789_col_pos_1.TGCTGGG.gene_counts.txt.gz",
    "GSM2796790_col_pos_2.TCAAACA.gene_counts.txt.gz",
    "GSM2796791_col_pos_3.AGCGCCG.gene_counts.txt.gz",
    "GSM2796792_col_pos_4.ATCGAGC.gene_counts.txt.gz",
]

def load_file(fpath):
    """Load 9-col file: Feature Chr Start End Strand Length Count CPM RPKM"""
    with gzip.open(fpath, 'rt') as fh:
        df = pd.read_csv(fh, sep='\t', index_col=0)
    # Column 'Count' (7th, 0-indexed 6th after index)
    if 'Count' in df.columns:
        return df['Count'].astype(int)
    # Fallback: last numeric column before CPM
    print(f"  WARNING: 'Count' not found in {fpath.name}, columns: {list(df.columns)}")
    return df.iloc[:, 5].astype(int)

# ── Load all samples ──────────────────────────────────────────────────────────
print("Loading GSE104385 HTSeq count files...")
sample_data = {}
for flist, label in [(CTRL_FILES, 'ctrl'), (TREAT_FILES, 'treat')]:
    for fname in flist:
        fpath = DATA_DIR / fname
        sname = fname.split('_', 2)[2].split('.')[0]  # e.g. "col_1"
        s = load_file(fpath)
        sample_data[sname] = s
        print(f"  {sname} ({label}): {len(s)} genes, total={s.sum():,}")

counts_df = pd.DataFrame(sample_data)
counts_df = counts_df.loc[counts_df.sum(axis=1) > 0]
# Keep only standard Arabidopsis gene IDs
counts_df = counts_df[counts_df.index.str.match(r'^AT[1-5MC]G\d+')]
print(f"\nCount matrix: {counts_df.shape}")

ctrl_cols  = [f.split('_', 2)[2].split('.')[0] for f in CTRL_FILES]
treat_cols = [f.split('_', 2)[2].split('.')[0] for f in TREAT_FILES]
print(f"Control: {ctrl_cols}")
print(f"Treated: {treat_cols}")

# Library size check
lib = counts_df.sum(axis=0)
print(f"\nLibrary sizes:\n{lib.to_string()}")
ctrl_lib  = lib[ctrl_cols].mean()
treat_lib = lib[treat_cols].mean()
print(f"Mean ctrl: {ctrl_lib/1e6:.2f}M, Mean treat: {treat_lib/1e6:.2f}M, ratio: {treat_lib/ctrl_lib:.3f}x")

# ── DESeq2 ────────────────────────────────────────────────────────────────────
counts_sub = counts_df[ctrl_cols + treat_cols].T
meta = pd.DataFrame({
    'sample': ctrl_cols + treat_cols,
    'condition': ['ctrl'] * len(ctrl_cols) + ['treat'] * len(treat_cols)
}).set_index('sample')

print(f"\nRunning DESeq2 (n={len(ctrl_cols)} ctrl, {len(treat_cols)} treat)...")
try:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    dds = DeseqDataSet(counts=counts_sub, metadata=meta,
                       design_factors="condition",
                       ref_level=["condition", "ctrl"], quiet=True)
    dds.deseq2()
    stat_res = DeseqStats(dds, contrast=["condition", "treat", "ctrl"], quiet=True)
    stat_res.summary()
    res = stat_res.results_df.dropna(subset=['log2FoldChange'])
    print(f"DESeq2: {len(res)} genes with logFC, {(res.padj < 0.05).sum()} sig (padj<0.05)")

    # Validate auxin markers
    auxin_markers = {
        'IAA1':  'AT4G14560', 'IAA3':  'AT1G04240',
        'SAUR6': 'AT2G45210', 'GH3.3': 'AT2G23170',
        'ARF7':  'AT5G20730',
    }
    print("\nAuxin marker genes:")
    for name, gid in auxin_markers.items():
        if gid in res.index:
            r = res.loc[gid]
            print(f"  {name} ({gid}): logFC={r['log2FoldChange']:+.3f}, padj={r.get('padj', 'NA'):.4g}")
        else:
            print(f"  {name} ({gid}): NOT IN RESULTS")

    aux_sig = res['log2FoldChange'].rename('logFC')
    aux_sig.index.name = 'GeneID'

except ImportError:
    print("pydeseq2 not available — manual logFC")
    lib = counts_df.sum(axis=0)
    cpm = counts_df.divide(lib, axis=1) * 1e6
    aux_sig = (np.log2(cpm[treat_cols].mean(axis=1) + 1) - np.log2(cpm[ctrl_cols].mean(axis=1) + 1)).rename('logFC')
    aux_sig.index.name = 'GeneID'

# Save signature
sig_path = SIG_DIR / "Auxin_axis_logFC.csv"
aux_sig.to_csv(sig_path, header=True)
print(f"\nSaved: {sig_path} ({len(aux_sig)} genes)")

# ── Calibration ───────────────────────────────────────────────────────────────
def score_on_sig(log2expr, sig):
    common = log2expr.index.intersection(sig.index)
    return log2expr.loc[common].multiply(sig.loc[common], axis=0).sum(axis=0)

df_log2 = np.log2(counts_df.astype(float) + 1)
scores = score_on_sig(df_log2, aux_sig)
ctrl_score  = scores[ctrl_cols].mean()
treat_score = scores[treat_cols].mean()
cal_range   = treat_score - ctrl_score
print(f"\nCalibration: ctrl={ctrl_score:.2f}, treat={treat_score:.2f}, range={cal_range:.2f}")

cal = pd.read_csv(CAL_FILE, index_col=0)
cal.loc['Auxin_axis'] = {
    'ctrl_score': float(ctrl_score), 'treat_score': float(treat_score),
    'ctrl_label': 'Col0_EtOH_2h', 'treat_label': 'Col0_IAA_2h', 'source': 'GSE104385'
}
cal.to_csv(CAL_FILE)
print("Updated calibration.csv with Auxin_axis")

def calibrate_aux(score): return (float(score) - ctrl_score) / cal_range * 100

# ── Score biostimulants ───────────────────────────────────────────────────────
print("\n=== BIOSTIMULANT AUXIN-AXIS SCORES ===")
import gzip as gz

def load_htseq_dir(gse_dir):
    files = list(Path(gse_dir).glob("*.txt.gz"))
    dfs = []
    for f in files:
        try:
            with gz.open(f, 'rt') as fh:
                lines = [l.strip().split('\t') for l in fh if not l.startswith('__')]
            genes  = [l[0] for l in lines if len(l) == 2]
            counts_ = [int(l[1]) for l in lines if len(l) == 2]
            if genes:
                dfs.append(pd.Series(counts_, index=genes, name=f.stem.split('.')[0]))
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

# Also try loading the 9-col format from this same dataset dir
def load_9col_dir(gse_dir, ctrl_pattern, treat_pattern):
    """Load 9-col HTSeq files (Feature Chr Start End Strand Length Count CPM RPKM)."""
    ctrl_dfs, treat_dfs = [], []
    for f in Path(gse_dir).glob("*.txt.gz"):
        try:
            with gz.open(f, 'rt') as fh:
                df_f = pd.read_csv(fh, sep='\t', index_col=0)
            if 'Count' in df_f.columns:
                series = df_f['Count'].astype(int)
            else:
                continue
            sname = f.stem.split('.')[0]
            if ctrl_pattern in f.name:
                ctrl_dfs.append(series.rename(sname))
            elif treat_pattern in f.name:
                treat_dfs.append(series.rename(sname))
        except Exception:
            continue
    if not ctrl_dfs or not treat_dfs:
        return None, [], []
    df_all = pd.concat(ctrl_dfs + treat_dfs, axis=1).fillna(0).astype(int)
    ctrl_names = [d.name for d in ctrl_dfs]
    treat_names = [d.name for d in treat_dfs]
    return df_all, ctrl_names, treat_names

results_aux = {}

def score_dataset_aux(name, count_df, treat_cols_ds, ctrl_cols_ds):
    log2e = np.log2(count_df.astype(float) + 1)
    s = score_on_sig(log2e, aux_sig)
    t = s[treat_cols_ds].mean()
    c = s[ctrl_cols_ds].mean()
    dt = calibrate_aux(t) - calibrate_aux(c)
    results_aux[name] = dt
    print(f"  {name}: ctrl={calibrate_aux(c):+.1f}% treat={calibrate_aux(t):+.1f}% delta={dt:+.2f}%")

# GSE297649 (simple 2-col HTSeq files)
ds = load_htseq_dir(ROOT / "RNA_seq_data/suppl_files/GSE297649")
if ds is not None:
    aa_t  = [c for c in ds.columns if 'wtaa' in c.lower()]
    hs_t  = [c for c in ds.columns if 'wtsh' in c.lower() or 'wths' in c.lower()]
    ctrl = [c for c in ds.columns if c not in aa_t + hs_t]
    score_dataset_aux('amino_acid',  ds, aa_t,  ctrl)
    score_dataset_aux('humic_subst', ds, hs_t,  ctrl)

# GMV volatile
gmv_df = load_featurecounts(ROOT / "RNA_seq_data/suppl_files/GSE138478/GSE138478_HMZ008_featureCounts_output.txt.gz")
if gmv_df is not None:
    gmv_ctrl  = [c for c in gmv_df.columns if '0.05C' in c or '0.5C' in c]
    gmv_treat = [c for c in gmv_df.columns if '0.05T' in c]
    if not gmv_treat: gmv_treat = [c for c in gmv_df.columns if c not in gmv_ctrl]
    score_dataset_aux('GMV_treated', gmv_df, gmv_treat, gmv_ctrl)

# Diacetyl
dia_df = load_featurecounts(ROOT / "RNA_seq_data/suppl_files/GSE138478/GSE138478_HMZ014_Dione_featureCounts_output.txt.gz")
if dia_df is not None:
    dia_ctrl  = [c for c in dia_df.columns if '_CK_' in c]
    dia_treat = [c for c in dia_df.columns if '_Dione_' in c]
    score_dataset_aux('diacetyl_treated', dia_df, dia_treat, dia_ctrl)

# TiO2 (rLog)
tio2_df = pd.read_csv(ROOT / "RNA_seq_data/suppl_files/GSE208223/GSE208223_Shoot_norm.txt.gz",
                       sep='\t', index_col=0, compression='gzip')
tio2_df.index = tio2_df.index.str.replace(r'\.\d+$', '', regex=True).str.strip('"')
common = tio2_df.index.intersection(aux_sig.index)
tio2_s = tio2_df.loc[common].multiply(aux_sig.loc[common], axis=0).sum(axis=0)
t2c = tio2_s[['E1','E2','E3']].mean()
t2t = tio2_s[['F1','F2','F3']].mean()
results_aux['TiO2_treated'] = calibrate_aux(t2t) - calibrate_aux(t2c)
print(f"  TiO2: delta={results_aux['TiO2_treated']:+.2f}%")

print("\n=== SUMMARY: Auxin-axis deltas ===")
for t, d in results_aux.items():
    print(f"  {t:<25}: {d:+.2f}%")

# Update summary CSV
summary_file = SCORE_DIR / "biostimulant_6axis_summary.csv"
summary = pd.read_csv(summary_file)
summary = summary[summary['axis'] != 'Auxin_axis']
new_rows = [{'treatment': t, 'axis': 'Auxin_axis', 'delta_pct': d} for t, d in results_aux.items()]
summary = pd.concat([summary, pd.DataFrame(new_rows)], ignore_index=True)
summary.to_csv(summary_file, index=False)
print(f"\nSaved updated summary with Auxin axis: {summary_file}")
print("Done.")
