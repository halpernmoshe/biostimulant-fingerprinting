import pandas as pd
import numpy as np
import gzip
import os
import json

PYTHON_ROOT = "C:/Users/moshe/Dropbox/ISF 2025"
SIG_DIR = f"{PYTHON_ROOT}/RNA_seq_axes/signatures"
SCORES_DIR = f"{PYTHON_ROOT}/RNA_seq_axes/scores"
DATA_DIR = f"{PYTHON_ROOT}/RNA_seq_data/suppl_files"

AXES = {
    'P_axis':   ('P_axis_logFC.csv',   26286),
    'N_axis':   ('N_axis_logFC.csv',   22617),
    'ABA_axis': ('ABA_axis_logFC.csv', 19453),
    'JA_axis':  ('JA_axis_logFC.csv',  12847),
    'SA_axis':  ('SA_axis_logFC.csv',  20189),
    'Fe_axis':  ('Fe_axis_logFC.csv',   4260),
}

def load_signature(axis_name):
    fname, calib = AXES[axis_name]
    sig = pd.read_csv(f"{SIG_DIR}/{fname}", index_col=0)
    sig.index = sig.index.str.split('.').str[0]
    return sig['logFC'], calib

def compute_scores(counts_df, sig_v):
    common = counts_df.index.intersection(sig_v.index)
    counts_sub = counts_df.loc[common]
    sig_sub = sig_v.loc[common]
    lib_sizes = counts_sub.sum(axis=0)
    cpm = counts_sub.div(lib_sizes) * 1e6
    expr = np.log2(cpm + 1)
    scores = expr.T.dot(sig_sub)
    return scores, len(common)

def score_groups(counts_df, treat_cols, ctrl_cols, sig_v, calib, axis_name, label):
    scores, n_genes = compute_scores(counts_df, sig_v)
    treat_found = [c for c in treat_cols if c in scores.index]
    ctrl_found  = [c for c in ctrl_cols  if c in scores.index]
    if not treat_found:
        print(f"  WARNING [{axis_name}] no treat cols for {label}")
        print(f"  Available: {scores.index.tolist()[:10]}")
        return None, None, None
    if not ctrl_found:
        print(f"  WARNING [{axis_name}] no ctrl cols for {label}")
        return None, None, None
    treat_mean = scores[treat_found].mean()
    ctrl_mean  = scores[ctrl_found].mean()
    delta = treat_mean - ctrl_mean
    pct   = delta / calib * 100
    print(f"  {label} | {axis_name}: n_genes={n_genes}, delta={delta:.2f}, {pct:.2f}%")
    return delta, pct, n_genes


# === DATASET 1: GSE87338 ===
print("=" * 70)
print("DATASET 1: GSE87338 (Castrillo 2017)")
print("=" * 70)

gz_path1 = f"{DATA_DIR}/GSE87338/GSE87338_PBI_counts.txt.gz"

with gzip.open(gz_path1, 'rt') as f:
    sample_names_line = f.readline().strip().split('\t')
    pbi_header_line   = f.readline().strip().split('\t')
    data_lines = []
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) > 1:
            data_lines.append(parts)

gene_ids = [row[0] for row in data_lines]
count_matrix = []
for row in data_lines:
    try:
        count_matrix.append([int(x) for x in row[1:]])
    except ValueError:
        count_matrix.append([0] * len(sample_names_line))

counts_gse87338 = pd.DataFrame(count_matrix, index=gene_ids, columns=sample_names_line)
counts_gse87338.index = counts_gse87338.index.str.split('.').str[0]

print(f"GSE87338 loaded: {counts_gse87338.shape[0]} genes x {counts_gse87338.shape[1]} samples")

all_cols = counts_gse87338.columns.tolist()
col0_highP_noBact = [c for c in all_cols if 'Col-0_HighP_noBact' in c]
col0_lowP_noBact  = [c for c in all_cols if 'Col-0_LowP_noBact'  in c]
col0_highP_SynCom = [c for c in all_cols if 'Col-0_HighP_SynCom' in c]
col0_lowP_SynCom  = [c for c in all_cols if 'Col-0_LowP_SynCom'  in c]

print(f"Col-0 groups: HighP_noBact={len(col0_highP_noBact)}, LowP_noBact={len(col0_lowP_noBact)}, HighP_SynCom={len(col0_highP_SynCom)}, LowP_SynCom={len(col0_lowP_SynCom)}")

comparisons_87338 = [
    ("SynCom_LowP_vs_noBact_LowP",   col0_lowP_SynCom,  col0_lowP_noBact,  "SynCom effect at Low Pi"),
    ("SynCom_HighP_vs_noBact_HighP",  col0_highP_SynCom, col0_highP_noBact, "SynCom effect at High Pi"),
    ("LowP_noBact_vs_HighP_noBact",   col0_lowP_noBact,  col0_highP_noBact, "Pi starvation (no bacteria)"),
    ("LowP_SynCom_vs_HighP_SynCom",   col0_lowP_SynCom,  col0_highP_SynCom, "Pi starvation (with SynCom)"),
]

results_gse87338 = {}
print("\nScoring GSE87338:")
for label, treat_cols, ctrl_cols, description in comparisons_87338:
    print(f"\n  [{label}] {description}")
    row = {'description': description, 'n_treat': len(treat_cols), 'n_ctrl': len(ctrl_cols)}
    for axis_name in AXES:
        sig_v, calib = load_signature(axis_name)
        delta, pct, n_genes = score_groups(counts_gse87338, treat_cols, ctrl_cols, sig_v, calib, axis_name, label)
        row[axis_name + '_delta'] = round(float(delta), 3) if delta is not None else None
        row[axis_name + '_pct']   = round(float(pct), 3)   if pct   is not None else None
        row['n_genes_' + axis_name] = int(n_genes) if n_genes is not None else None
    results_gse87338[label] = row


# === DATASET 2: GSE157128 ===
print("\n" + "=" * 70)
print("DATASET 2: GSE157128 (Harbort 2020)")
print("=" * 70)

gz_path2 = f"{DATA_DIR}/GSE157128/GSE157128_Col-0_RNA-Seq_raw_counts.txt.gz"

with gzip.open(gz_path2, 'rt') as f:
    counts_raw = pd.read_csv(f, sep='\t', index_col=0, comment='#')

annot_cols = ['Gene', 'Chromosome', 'Start', 'End', 'Strand', 'Length', 'BioType', 'Description']
count_cols = [c for c in counts_raw.columns if c not in annot_cols]
counts_gse157128 = counts_raw[count_cols].copy()
counts_gse157128.index = counts_gse157128.index.str.split('.').str[0]
counts_gse157128 = counts_gse157128.apply(pd.to_numeric, errors="coerce").fillna(0)
counts_gse157128 = counts_gse157128.groupby(counts_gse157128.index).sum()

print(f"GSE157128 loaded: {counts_gse157128.shape[0]} genes x {counts_gse157128.shape[1]} samples")
print("Columns:", counts_gse157128.columns.tolist())

syncom33   = ['SynCom33_1', 'SynCom33_2', 'SynCom33_3', 'SynCom33_4']
syncom35   = ['SynCom35_1', 'SynCom35_2', 'SynCom35_3', 'SynCom35_4']
mock       = ['Mock_1', 'Mock_2', 'Mock_3', 'Mock_4']
hk33       = ['HKSynCom33_1', 'HKSynCom33_2', 'HKSynCom33_3', 'HKSynCom33_4']
hk35       = ['HKSynCom35_1', 'HKSynCom35_2', 'HKSynCom35_3', 'HKSynCom35_4']
mock_flg22 = ['Mock_Flg22_1', 'Mock_Flg22_2', 'Mock_Flg22_3', 'Mock_Flg22_4']
sc33_flg22 = ['SynCom33_Flg22_1', 'SynCom33_Flg22_2', 'SynCom33_Flg22_3', 'SynCom33_Flg22_4']
sc35_flg22 = ['SynCom35_Flg22_1', 'SynCom35_Flg22_2', 'SynCom35_Flg22_3', 'SynCom35_Flg22_4']

comparisons_157128 = [
    ("SynCom33_vs_Mock",             syncom33,    mock,       "Living non-suppressive SynCom33 vs axenic"),
    ("SynCom35_vs_Mock",             syncom35,    mock,       "Living suppressive SynCom35 vs axenic"),
    ("HKSynCom33_vs_Mock",           hk33,        mock,       "Heat-killed SynCom33 vs axenic (MAMP-only)"),
    ("HKSynCom35_vs_Mock",           hk35,        mock,       "Heat-killed SynCom35 vs axenic (MAMP-only)"),
    ("flg22_vs_Mock",                mock_flg22,  mock,       "flg22 MAMP alone vs axenic"),
    ("SynCom33_Flg22_vs_Mock_Flg22", sc33_flg22,  mock_flg22, "SynCom33 effect in flg22 context"),
    ("SynCom35_Flg22_vs_Mock_Flg22", sc35_flg22,  mock_flg22, "SynCom35 effect in flg22 context"),
]

results_gse157128 = {}
print("\nScoring GSE157128:")
for label, treat_cols, ctrl_cols, description in comparisons_157128:
    print(f"\n  [{label}] {description}")
    row = {'description': description, 'n_treat': len(treat_cols), 'n_ctrl': len(ctrl_cols)}
    for axis_name in AXES:
        sig_v, calib = load_signature(axis_name)
        delta, pct, n_genes = score_groups(counts_gse157128, treat_cols, ctrl_cols, sig_v, calib, axis_name, label)
        row[axis_name + '_delta'] = round(float(delta), 3) if delta is not None else None
        row[axis_name + '_pct']   = round(float(pct), 3)   if pct   is not None else None
        row['n_genes_' + axis_name] = int(n_genes) if n_genes is not None else None
    results_gse157128[label] = row


# === SAVE RESULTS ===
print("\n" + "=" * 70)
print("SAVING RESULTS")
print("=" * 70)

rows = []
for label, row in results_gse87338.items():
    r = {'dataset': 'GSE87338', 'comparison': label}
    r.update(row)
    rows.append(r)
for label, row in results_gse157128.items():
    r = {'dataset': 'GSE157128', 'comparison': label}
    r.update(row)
    rows.append(r)

df_out = pd.DataFrame(rows)
out_csv = f"{SCORES_DIR}/PGPR_validation_scores.csv"
df_out.to_csv(out_csv, index=False)
print(f"Saved CSV: {out_csv}")

pct_axes = ['P_axis', 'N_axis', 'ABA_axis', 'JA_axis', 'SA_axis', 'Fe_axis']
pct_cols = [ax + '_pct' for ax in pct_axes]

print("\nSUMMARY TABLE (% calibration range):")
print(f"{'Dataset':<12} {'Comparison':<44} {'P':>7} {'N':>7} {'ABA':>7} {'JA':>7} {'SA':>7} {'Fe':>7}")
print("-" * 105)
for _, r in df_out.iterrows():
    vals = [r.get(c, None) for c in pct_cols]
    vals_str = [f"{v:7.1f}" if v is not None else "    N/A" for v in vals]
    print(f"{r['dataset']:<12} {r['comparison']:<44} {''.join(vals_str)}")

json_out = f"{SCORES_DIR}/PGPR_validation_scores.json"
all_results = {'GSE87338': results_gse87338, 'GSE157128': results_gse157128}

def make_serial(obj):
    if isinstance(obj, dict):
        return {k: make_serial(v) for k, v in obj.items()}
    elif isinstance(obj, float) and np.isnan(obj):
        return None
    elif isinstance(obj, (np.float64, np.float32)):
        return float(obj)
    elif isinstance(obj, (np.int64, np.int32)):
        return int(obj)
    return obj

with open(json_out, 'w') as f:
    json.dump(make_serial(all_results), f, indent=2)
print(f"Saved JSON: {json_out}")
print("\nDone.")

