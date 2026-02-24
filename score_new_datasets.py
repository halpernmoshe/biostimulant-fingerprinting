#!/usr/bin/env python3
import sys, warnings, gzip, os
warnings.filterwarnings("ignore")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import numpy as np
import pandas as pd
from pathlib import Path

ROOT = Path("C:/Users/moshe/Dropbox/ISF 2025")
SIG_DIR = ROOT / "RNA_seq_axes" / "signatures"

CALIB = {
    'P_axis':    (16652.763941, -872.395312),
    'N_axis':    (18594.784332, 923.997646),
    'ABA_axis':  (9504.783016, -5668.164126),
    'JA_axis':   (6463.471145, -2892.189941),
    'SA_axis':   (11155.108413, -5008.787578),
    'Fe_axis':   (2751.369024, -1508.685827),
    'CK_axis':   (2790.343455, -3025.991931),
    'ET_axis':   (2605.979079, -1853.526591),
    'Auxin_axis':(13106.379233, -4090.714534),
}

def calib_range(axis_name):
    t, c = CALIB[axis_name]
    return t - c

def read_htseq_counts(directory):
    files = sorted(Path(directory).glob("*.gz"))
    dfs = []
    for f in files:
        with gzip.open(f, "rt") as fh:
            lines = [l.strip().split("\t") for l in fh]
        genes = [l[0] for l in lines if len(l) == 2 and not l[0].startswith("__")]
        counts_vals = [int(l[1]) for l in lines if len(l) == 2 and not l[0].startswith("__")]
        if genes:
            dfs.append(pd.Series(counts_vals, index=genes, name=f.name))
    if not dfs:
        return pd.DataFrame()
    return pd.concat(dfs, axis=1).fillna(0).astype(int)

def read_featurecounts_aggr(directory):
    files = sorted(Path(directory).glob("*.gz"))
    gene_counts = {}
    for f in files:
        sample_name = f.name
        with gzip.open(f, "rt") as fh:
            lines = fh.readlines()
        data_lines = [l.strip() for l in lines[2:] if l.strip()]
        gene_to_count = {}
        for line in data_lines:
            parts = line.split("\t")
            if len(parts) >= 8:
                gene_id = parts[6]
                count = int(parts[7])
                gene_to_count[gene_id] = gene_to_count.get(gene_id, 0) + count
        gene_counts[sample_name] = gene_to_count
    df = pd.DataFrame(gene_counts).fillna(0).astype(int)
    return df

def score_contrast(counts, treat_cols, ctrl_cols, axis_name):
    sig_path = SIG_DIR / (axis_name + "_logFC.csv")
    if not sig_path.exists():
        return None, None, None
    sig = pd.read_csv(sig_path, index_col=0)
    sig_col = 'log2FoldChange' if 'log2FoldChange' in sig.columns else sig.columns[0]
    sig_v = sig[sig_col]
    all_cols = treat_cols + ctrl_cols
    common = counts.index.intersection(sig_v.index)
    sig_v = sig_v.loc[common]
    lib_sizes = counts[all_cols].sum(axis=0)
    cpm = counts.loc[common, all_cols].div(lib_sizes) * 1e6
    expr = np.log2(cpm + 1)
    scores = expr.T.dot(sig_v)
    treat_score = scores[treat_cols].mean()
    ctrl_score = scores[ctrl_cols].mean()
    delta = treat_score - ctrl_score
    cr = calib_range(axis_name)
    pct = delta / cr * 100
    return pct, treat_score, ctrl_score

axes_to_score = ['P_axis', 'N_axis', 'ABA_axis', 'JA_axis', 'SA_axis',
                 'Fe_axis', 'CK_axis', 'ET_axis', 'Auxin_axis']

print("="*60)
print("GSE222356: Serendipita indica on Arabidopsis roots")
print("="*60)

d222 = ROOT / "RNA_seq_data" / "suppl_files" / "GSE222356"
c222 = read_htseq_counts(d222)
print(f"Shape: {c222.shape}")
print(f"Sample reads:\n{c222.sum().to_string()}")

file_meta = {
    'GSM6921484_WTCHG_123535_02_htseq_filter.txt': ('mock', '3dai'),
    'GSM6921485_WTCHG_123535_04_htseq_filter.txt': ('mock', '10dai'),
    'GSM6921486_WTCHG_123535_05_htseq_filter.txt': ('treated', '3dai'),
    'GSM6921487_WTCHG_123535_06_htseq_filter.txt': ('treated', '10dai'),
    'GSM6921488_WTCHG_124311_07_htseq_filter.txt': ('mock', '3dai'),
    'GSM6921489_WTCHG_124311_12_htseq_filter.txt': ('mock', '10dai'),
    'GSM6921490_WTCHG_124311_16_htseq_filter.txt': ('treated', '3dai'),
    'GSM6921491_WTCHG_124311_19_htseq_filter.txt': ('treated', '10dai'),
}

results_222 = {}
for timepoint in ['3dai', '10dai']:
    treat_gsm = [k for k,(c,t) in file_meta.items() if c=='treated' and t==timepoint]
    ctrl_gsm = [k for k,(c,t) in file_meta.items() if c=='mock' and t==timepoint]
    treat_cols = [c for c in c222.columns if any(g.replace('.txt','') in c for g in treat_gsm)]
    ctrl_cols = [c for c in c222.columns if any(g.replace('.txt','') in c for g in ctrl_gsm)]
    print(f"\n{timepoint}: treat={treat_cols}, ctrl={ctrl_cols}")
    row = {}
    for axis in axes_to_score:
        pct, ts, cs = score_contrast(c222, treat_cols, ctrl_cols, axis)
        if pct is not None:
            row[axis] = pct
            print(f"  {axis}: {pct:+.1f}%")
    results_222[timepoint] = row

print("\n" + "="*60)
print("GSE225213: Rhizophagus AMF on Arabidopsis Col-0 roots")
print("="*60)

d225 = ROOT / "RNA_seq_data" / "suppl_files" / "GSE225213"
c225 = read_featurecounts_aggr(d225)
print(f"Shape: {c225.shape}")
print(f"Sample reads:\n{c225.sum().to_string()}")

treat_am = [c for c in c225.columns if any(x in c for x in ['GSM7041574','GSM7041575','GSM7041576','GSM7041577'])]
ctrl_am = [c for c in c225.columns if any(x in c for x in ['GSM7041585','GSM7041586','GSM7041587','GSM7041588'])]
print(f"\n+AM={treat_am}")
print(f"-AM={ctrl_am}")

results_225 = {}
for axis in axes_to_score:
    pct, ts, cs = score_contrast(c225, treat_am, ctrl_am, axis)
    if pct is not None:
        results_225[axis] = pct
        print(f"  {axis}: {pct:+.1f}%")

print("\n" + "="*60)
print("SUMMARY")
print("="*60)
header = f"{'Axis':12s} | {'Serendipita 3dai':16s} | {'Serendipita 10dai':16s} | {'AMF Rhizoph':12s}"
print(header)
print("-"*65)
for axis in axes_to_score:
    v1 = results_222.get('3dai',{}).get(axis, float('nan'))
    v2 = results_222.get('10dai',{}).get(axis, float('nan'))
    v3 = results_225.get(axis, float('nan'))
    print(f"  {axis:10s} | {v1:+12.1f}%   | {v2:+12.1f}%   | {v3:+12.1f}%")
