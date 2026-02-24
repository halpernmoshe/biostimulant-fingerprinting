#!/usr/bin/env python3
"""
Cross-species projection: PGPR and humic substance biostimulant effects
in Arabidopsis scored on a tomato phosphate starvation response (PSR) axis.

Approach: Build tomato-P-axis-in-TAIR-space by mapping Solyc IDs -> TAIR IDs
via the ortholog table, then score Arabidopsis count data on that signature.
Compare delta scores on conserved tomato gene set vs full Arabidopsis P-axis.
"""

import gzip
import re
import sys
import os
import pandas as pd
import numpy as np
import xlrd

ROOT = 'C:/Users/moshe/Dropbox/ISF 2025'

# ============================================================================
# 1. Load signatures
# ============================================================================
print('Loading signatures...')
P_ath = pd.read_csv(ROOT + '/RNA_seq_axes/signatures/P_axis_logFC.csv', index_col=0)['logFC']
P_tom_raw = pd.read_csv(ROOT + '/RNA_seq_axes/signatures/P_axis_tomato_logFC.csv', index_col=0)['logFC']
P_tom_raw.index = [re.sub(r'\.\d+$', '', idx) for idx in P_tom_raw.index]

print(f'  Arabidopsis P-axis: {len(P_ath)} genes')
print(f'  Tomato P-axis:      {len(P_tom_raw)} genes (Solyc IDs, version stripped)')

# ============================================================================
# 2. Build ortholog map (Solyc -> TAIR) and tomato-axis-in-TAIR-space
# ============================================================================
print('\nBuilding Solyc->TAIR ortholog map...')
orth = pd.read_csv(ROOT + '/RNA_seq_axes/tair_to_solyc_orthologs.tsv', sep='\t')
orth.columns = ['tair_id', 'solyc_id', 'solyc_name', 'orthology_type']
orth['solyc_id_base'] = orth['solyc_id'].str.replace(r'\.\d+$', '', regex=True)
orth_valid = orth.dropna(subset=['solyc_id_base']).copy()

# Priority: one2one > one2many > many2many (best orthology type gets TAIR assignment)
orth_valid['priority'] = orth_valid['orthology_type'].map(
    {'ortholog_one2one': 0, 'ortholog_one2many': 1, 'ortholog_many2many': 2}).fillna(3)
solyc_to_tair = (orth_valid.sort_values('priority')
                  .groupby('solyc_id_base')['tair_id'].first()
                  .to_dict())

# For each tomato P-axis gene, map to TAIR (first/best hit; skip duplicates)
tomato_axis_tair = {}
n_mapped = 0
for solyc_id, logfc in P_tom_raw.items():
    if solyc_id in solyc_to_tair:
        tair_id = solyc_to_tair[solyc_id]
        if tair_id not in tomato_axis_tair:
            tomato_axis_tair[tair_id] = logfc
            n_mapped += 1

T_sig = pd.Series(tomato_axis_tair)
pct_tomato_mapped = n_mapped / len(P_tom_raw) * 100
n_tom_sig = len(T_sig)
print(f'  Tomato genes with TAIR ortholog: {n_mapped} / {len(P_tom_raw)} ({pct_tomato_mapped:.1f}%)')
print(f'  Unique TAIR genes in tomato-axis-in-TAIR-space: {n_tom_sig}')

# ============================================================================
# 3. Scoring function
# ============================================================================
def score_sample(counts_s, sig):
    """Dot product score: sum(log2(count+1) * logFC) over shared genes."""
    common = sig.index.intersection(counts_s.index)
    if len(common) < 50:
        return np.nan, len(common)
    expr = np.log2(counts_s[common].astype(float) + 1)
    return (sig[common] * expr).sum(), len(common)


# ============================================================================
# 4. GSE297649: Humic substances + Amino acids (Arabidopsis)
# ============================================================================
print('\n=== GSE297649: Humic substances + Amino acids (Arabidopsis) ===')
suppl297 = ROOT + '/RNA_seq_data/suppl_files/GSE297649'
htseq_map = {
    'WTaa_1': 'GSM8996165_htseqcount_WTaa_1.txt.gz',
    'WTaa_2': 'GSM8996166_htseqcount_WTaa_2.txt.gz',
    'WTaa_3': 'GSM8996167_htseqcount_WTaa_3.txt.gz',
    'WTcontrol_1': 'GSM8996168_htseqcount_WTcontrol_1.txt.gz',
    'WTcontrol_2': 'GSM8996169_htseqcount_WTcontrol_2.txt.gz',
    'WTcontrol_3': 'GSM8996170_htseqcount_WTcontrol_3.txt.gz',
    'WTsh_1': 'GSM8996171_htseqcount_WTsh_1.txt.gz',
    'WTsh_2': 'GSM8996172_htseqcount_WTsh_2.txt.gz',
    'WTsh_3': 'GSM8996173_htseqcount_WTsh_3.txt.gz',
}

def get_grp_297649(sname):
    if 'WTaa' in sname: return 'amino_acid'
    if 'WTsh' in sname: return 'humic_subst'
    return 'control'

scores_297649 = {}
for sname, fname in htseq_map.items():
    fpath = suppl297 + '/' + fname
    with gzip.open(fpath, 'rt') as fh:
        lines = [l.strip().split('\t') for l in fh if not l.startswith('__')]
    genes = [l[0] for l in lines if len(l) == 2 and l[0].startswith('AT')]
    vals  = [int(l[1]) for l in lines if len(l) == 2 and l[0].startswith('AT')]
    series = pd.Series(vals, index=genes)
    s_ath, n_ath = score_sample(series, P_ath)
    s_tom, n_tom = score_sample(series, T_sig)
    scores_297649[sname] = {
        'group': get_grp_297649(sname),
        's_ath': s_ath, 's_tom': s_tom,
        'n_ath': n_ath, 'n_tom': n_tom
    }

df_297649 = pd.DataFrame(scores_297649).T
print('\nPer-sample scores (GSE297649):')
print(df_297649[['group', 's_ath', 's_tom', 'n_ath', 'n_tom']].to_string())

ctrl_ath  = df_297649.loc[df_297649.group == 'control', 's_ath'].mean()
ctrl_tom  = df_297649.loc[df_297649.group == 'control', 's_tom'].mean()
n_ath_gse297 = int(df_297649.loc[df_297649.group == 'control', 'n_ath'].mean())
n_tom_gse297 = int(df_297649.loc[df_297649.group == 'control', 'n_tom'].mean())

results_297649 = []
for group, biostim in [('amino_acid', 'Amino acid biostimulant'),
                        ('humic_subst', 'Humic substance biostimulant')]:
    grp_ath = df_297649.loc[df_297649.group == group, 's_ath'].mean()
    grp_tom = df_297649.loc[df_297649.group == group, 's_tom'].mean()
    delta_ath = grp_ath - ctrl_ath
    delta_tom = grp_tom - ctrl_tom
    ratio = delta_tom / delta_ath if delta_ath != 0 else float('nan')
    pg_ath = delta_ath / n_ath_gse297
    pg_tom = delta_tom / n_tom_gse297
    print(f'\n  {group}:')
    print(f'    delta_ath = {delta_ath:+.2f}  (n_genes = {n_ath_gse297})')
    print(f'    delta_tom = {delta_tom:+.2f}  (n_genes = {n_tom_gse297})')
    print(f'    ratio (tom/ath) = {ratio:.3f}')
    print(f'    per-gene ath = {pg_ath:.6f},  per-gene tom = {pg_tom:.6f}')
    results_297649.append({
        'dataset': 'GSE297649', 'treatment': group, 'biostim_class': biostim,
        'delta_ath': delta_ath, 'n_ath': n_ath_gse297,
        'delta_tom': delta_tom, 'n_tom': n_tom_gse297,
        'ratio_tom_ath': ratio,
        'per_gene_ath': pg_ath, 'per_gene_tom': pg_tom
    })


# ============================================================================
# 5. GSE138478: Diacetyl PGPR volatiles (Arabidopsis)
# ============================================================================
print('\n=== GSE138478: Diacetyl PGPR volatiles (Arabidopsis) ===')
suppl138 = ROOT + '/RNA_seq_data/suppl_files/GSE138478'

def load_featurecounts(gz_path):
    with gzip.open(gz_path, 'rt') as fh:
        lines = list(fh)
    header = lines[1].strip().split('\t')
    sample_cols = header[6:]
    data = {}
    for line in lines[2:]:
        parts = line.strip().split('\t')
        if not parts[0].startswith('AT'):
            continue
        gene = parts[0]
        for i, sc in enumerate(sample_cols):
            data.setdefault(sc, {})[gene] = int(parts[6 + i])
    return {k: pd.Series(v) for k, v in data.items()}

counts_hmz008 = load_featurecounts(suppl138 + '/GSE138478_HMZ008_featureCounts_output.txt.gz')
counts_hmz014 = load_featurecounts(suppl138 + '/GSE138478_HMZ014_Dione_featureCounts_output.txt.gz')

grp_hmz008 = {
    'HMZ008_0.05C_5DAT_rep1.bam': 'GMV_control',
    'HMZ008_0.05C_5DAT_rep2.bam': 'GMV_control',
    'HMZ008_0.5C_5DAT_rep1.bam':  'GMV_control',
    'HMZ008_0.5C_5DAT_rep2.bam':  'GMV_control',
    'HMZ008_0.05T_5DAT_rep1.bam': 'GMV_treated',
    'HMZ008_0.05T_5DAT_rep2.bam': 'GMV_treated',
    'HMZ008_0.5T_5DAT_rep1.bam':  'GMV_treated',
    'HMZ008_0.5T_5DAT_rep2.bam':  'GMV_treated',
}
grp_hmz014 = {
    'HMZ014_CK_2DAT_rep1_th2.bam':    'diacetyl_control',
    'HMZ014_CK_2DAT_rep2_th2.bam':    'diacetyl_control',
    'HMZ014_Dione_2DAT_rep1_th2.bam': 'diacetyl_treated',
    'HMZ014_Dione_2DAT_rep2_th2.bam': 'diacetyl_treated',
}

def score_fc_groups(counts_dict, group_map, ctrl_grp, trt_grp, label, biostim, dataset):
    ctrl_ath, ctrl_tom, trt_ath, trt_tom = [], [], [], []
    n_ath_v, n_tom_v = [], []
    for sname, series in counts_dict.items():
        grp = group_map.get(sname)
        s_ath, n_ath = score_sample(series, P_ath)
        s_tom, n_tom = score_sample(series, T_sig)
        if grp == ctrl_grp:
            ctrl_ath.append(s_ath)
            ctrl_tom.append(s_tom)
            n_ath_v.append(n_ath)
            n_tom_v.append(n_tom)
        elif grp == trt_grp:
            trt_ath.append(s_ath)
            trt_tom.append(s_tom)
    delta_ath = float(np.mean(trt_ath)) - float(np.mean(ctrl_ath))
    delta_tom = float(np.mean(trt_tom)) - float(np.mean(ctrl_tom))
    n_ath = int(np.mean(n_ath_v))
    n_tom = int(np.mean(n_tom_v))
    ratio = delta_tom / delta_ath if delta_ath != 0 else float('nan')
    pg_ath = delta_ath / n_ath
    pg_tom = delta_tom / n_tom
    print(f'\n  {label}:')
    print(f'    delta_ath = {delta_ath:+.2f}  (n_genes = {n_ath})')
    print(f'    delta_tom = {delta_tom:+.2f}  (n_genes = {n_tom})')
    print(f'    ratio (tom/ath) = {ratio:.3f}')
    print(f'    per-gene ath = {pg_ath:.6f},  per-gene tom = {pg_tom:.6f}')
    return {
        'dataset': dataset, 'treatment': label, 'biostim_class': biostim,
        'delta_ath': delta_ath, 'n_ath': n_ath,
        'delta_tom': delta_tom, 'n_tom': n_tom,
        'ratio_tom_ath': ratio,
        'per_gene_ath': pg_ath, 'per_gene_tom': pg_tom
    }

r1_138 = score_fc_groups(counts_hmz008, grp_hmz008, 'GMV_control', 'GMV_treated',
                          'diacetyl_GMV', 'PGPR volatile (GMV diacetyl)', 'GSE138478')
r2_138 = score_fc_groups(counts_hmz014, grp_hmz014, 'diacetyl_control', 'diacetyl_treated',
                          'pure_diacetyl', 'PGPR volatile (pure diacetyl)', 'GSE138478')


# ============================================================================
# 6. GSE139154: Bacillus GB03 volatiles (Arabidopsis)
# ============================================================================
print('\n=== GSE139154: Bacillus GB03 volatiles (Arabidopsis) ===')
suppl139 = ROOT + '/RNA_seq_data/suppl_files/GSE139154'
with gzip.open(suppl139 + '/GSE139154_all_sample_cnt.tsv.gz', 'rt') as fh:
    df139 = pd.read_csv(fh, sep='\t', index_col=0)
ctrl139 = [c for c in df139.columns if 'Col-0_C' in c]
trt139  = [c for c in df139.columns if 'Col-0_T' in c]
print(f'  Ctrl: {ctrl139}, Trt: {trt139}')

ctrl_ath, ctrl_tom, trt_ath, trt_tom = [], [], [], []
n_ath_v, n_tom_v = [], []
for c in ctrl139:
    s, n = score_sample(df139[c], P_ath)
    ctrl_ath.append(s); n_ath_v.append(n)
    s2, n2 = score_sample(df139[c], T_sig)
    ctrl_tom.append(s2); n_tom_v.append(n2)
for c in trt139:
    s, _ = score_sample(df139[c], P_ath); trt_ath.append(s)
    s2, _ = score_sample(df139[c], T_sig); trt_tom.append(s2)

delta_ath_139 = np.mean(trt_ath) - np.mean(ctrl_ath)
delta_tom_139 = np.mean(trt_tom) - np.mean(ctrl_tom)
n_ath_139 = int(np.mean(n_ath_v))
n_tom_139 = int(np.mean(n_tom_v))
ratio_139 = delta_tom_139 / delta_ath_139 if delta_ath_139 != 0 else float('nan')
pg_ath_139 = delta_ath_139 / n_ath_139
pg_tom_139 = delta_tom_139 / n_tom_139
print(f'\n  Bacillus GB03:')
print(f'    delta_ath = {delta_ath_139:+.2f}  (n_genes = {n_ath_139})')
print(f'    delta_tom = {delta_tom_139:+.2f}  (n_genes = {n_tom_139})')
print(f'    ratio (tom/ath) = {ratio_139:.3f}')
print(f'    per-gene ath = {pg_ath_139:.6f},  per-gene tom = {pg_tom_139:.6f}')

r_139 = {
    'dataset': 'GSE139154', 'treatment': 'Bacillus_GB03',
    'biostim_class': 'PGPR volatile (Bacillus GB03)',
    'delta_ath': delta_ath_139, 'n_ath': n_ath_139,
    'delta_tom': delta_tom_139, 'n_tom': n_tom_139,
    'ratio_tom_ath': ratio_139,
    'per_gene_ath': pg_ath_139, 'per_gene_tom': pg_tom_139
}


# ============================================================================
# 7. GSE87337: SynCom bacterial community (Arabidopsis)
# ============================================================================
print('\n=== GSE87337: SynCom bacterial community (Arabidopsis) ===')
suppl873 = ROOT + '/RNA_seq_data/suppl_files/GSE87337'
with gzip.open(suppl873 + '/GSE87337_PBI_preTreat_counts.txt.gz', 'rt') as fh:
    lines873 = list(fh)

exp_header = lines873[0].strip().split('\t')
pbi_ids = lines873[1].strip().split('\t')[1:]  # row 1: Gene, PBI115, PBI116, ...

syncom_idx = [i for i, c in enumerate(exp_header)
              if 'LowP_SynCom' in c and 'Col-0' in c and 'preTreat(-P)' in c]
nobact_idx  = [i for i, c in enumerate(exp_header)
               if 'LowP_noBact' in c and 'Col-0' in c and 'preTreat(-P)' in c]
syncom_pbi  = [pbi_ids[i] for i in syncom_idx if i < len(pbi_ids)]
nobact_pbi  = [pbi_ids[i] for i in nobact_idx if i < len(pbi_ids)]
print(f'  SynCom PBI IDs ({len(syncom_pbi)}): {syncom_pbi}')
print(f'  NoBact PBI IDs ({len(nobact_pbi)}): {nobact_pbi}')

target_pbi = set(syncom_pbi + nobact_pbi)
pbi_to_col = {pbi: i for i, pbi in enumerate(pbi_ids)}
data873 = {pbi: {} for pbi in target_pbi}

for line in lines873[2:]:
    parts = line.strip().split('\t')
    if not parts[0].startswith('AT'):
        continue
    gene = parts[0]
    for pbi in target_pbi:
        col_idx = pbi_to_col[pbi]
        data873[pbi][gene] = int(parts[col_idx + 1])  # +1 for gene column

series873 = {pbi: pd.Series(d) for pbi, d in data873.items()}

ctrl_ath, ctrl_tom, trt_ath, trt_tom = [], [], [], []
n_ath_v, n_tom_v = [], []
for pbi in nobact_pbi:
    s, n = score_sample(series873[pbi], P_ath)
    ctrl_ath.append(s); n_ath_v.append(n)
    s2, n2 = score_sample(series873[pbi], T_sig)
    ctrl_tom.append(s2); n_tom_v.append(n2)
for pbi in syncom_pbi:
    s, _ = score_sample(series873[pbi], P_ath); trt_ath.append(s)
    s2, _ = score_sample(series873[pbi], T_sig); trt_tom.append(s2)

delta_ath_873 = np.mean(trt_ath) - np.mean(ctrl_ath)
delta_tom_873 = np.mean(trt_tom) - np.mean(ctrl_tom)
n_ath_873 = int(np.mean(n_ath_v))
n_tom_873 = int(np.mean(n_tom_v))
ratio_873 = delta_tom_873 / delta_ath_873 if delta_ath_873 != 0 else float('nan')
pg_ath_873 = delta_ath_873 / n_ath_873
pg_tom_873 = delta_tom_873 / n_tom_873
print(f'\n  SynCom:')
print(f'    delta_ath = {delta_ath_873:+.2f}  (n_genes = {n_ath_873})')
print(f'    delta_tom = {delta_tom_873:+.2f}  (n_genes = {n_tom_873})')
print(f'    ratio (tom/ath) = {ratio_873:.3f}')
print(f'    per-gene ath = {pg_ath_873:.6f},  per-gene tom = {pg_tom_873:.6f}')

r_873 = {
    'dataset': 'GSE87337', 'treatment': 'SynCom',
    'biostim_class': 'SynCom bacterial community',
    'delta_ath': delta_ath_873, 'n_ath': n_ath_873,
    'delta_tom': delta_tom_873, 'n_tom': n_tom_873,
    'ratio_tom_ath': ratio_873,
    'per_gene_ath': pg_ath_873, 'per_gene_tom': pg_tom_873
}


# ============================================================================
# 8. GSE199501: Priestia megaterium PGPR (Arabidopsis)
# ============================================================================
print('\n=== GSE199501: Priestia megaterium PGPR (Arabidopsis) ===')
xls_path = ROOT + '/RNA_seq_data/suppl_files/GSE199501/GSE199501_Processed_data_DEG_total_CPM.xls_tmp.xls'
wb = xlrd.open_workbook(xls_path)
ws = wb.sheets()[0]
headers = ws.row_values(0)
print(f'  Headers: {headers[:10]}')

# Ctrl = C-1953, C-1954, C-1955 (col idx 1,2,3); Trt = T-1956, T-1957, T-1958 (col idx 4,5,6)
ctrl_idx = [1, 2, 3]
trt_idx  = [4, 5, 6]

ctrl_data = {headers[i]: {} for i in ctrl_idx}
trt_data  = {headers[i]: {} for i in trt_idx}

for row in range(1, ws.nrows):
    rv = ws.row_values(row)
    gene = str(rv[0]).strip()
    if not gene.startswith('AT'):
        continue
    for i in ctrl_idx:
        ctrl_data[headers[i]][gene] = int(rv[i]) if rv[i] != '' else 0
    for i in trt_idx:
        trt_data[headers[i]][gene] = int(rv[i]) if rv[i] != '' else 0

ctrl_s = {k: pd.Series(v) for k, v in ctrl_data.items()}
trt_s  = {k: pd.Series(v) for k, v in trt_data.items()}

ctrl_ath, ctrl_tom, trt_ath, trt_tom = [], [], [], []
n_ath_v, n_tom_v = [], []
for k, series in ctrl_s.items():
    s, n = score_sample(series, P_ath)
    ctrl_ath.append(s); n_ath_v.append(n)
    s2, n2 = score_sample(series, T_sig)
    ctrl_tom.append(s2); n_tom_v.append(n2)
for k, series in trt_s.items():
    s, _ = score_sample(series, P_ath); trt_ath.append(s)
    s2, _ = score_sample(series, T_sig); trt_tom.append(s2)

delta_ath_199 = np.mean(trt_ath) - np.mean(ctrl_ath)
delta_tom_199 = np.mean(trt_tom) - np.mean(ctrl_tom)
n_ath_199 = int(np.mean(n_ath_v))
n_tom_199 = int(np.mean(n_tom_v))
ratio_199 = delta_tom_199 / delta_ath_199 if delta_ath_199 != 0 else float('nan')
pg_ath_199 = delta_ath_199 / n_ath_199
pg_tom_199 = delta_tom_199 / n_tom_199
print(f'\n  P. megaterium:')
print(f'    delta_ath = {delta_ath_199:+.2f}  (n_genes = {n_ath_199})')
print(f'    delta_tom = {delta_tom_199:+.2f}  (n_genes = {n_tom_199})')
print(f'    ratio (tom/ath) = {ratio_199:.3f}')
print(f'    per-gene ath = {pg_ath_199:.6f},  per-gene tom = {pg_tom_199:.6f}')

r_199 = {
    'dataset': 'GSE199501', 'treatment': 'P_megaterium',
    'biostim_class': 'PGPR (Priestia megaterium)',
    'delta_ath': delta_ath_199, 'n_ath': n_ath_199,
    'delta_tom': delta_tom_199, 'n_tom': n_tom_199,
    'ratio_tom_ath': ratio_199,
    'per_gene_ath': pg_ath_199, 'per_gene_tom': pg_tom_199
}


# ============================================================================
# 9. Compile all results
# ============================================================================
all_results = results_297649 + [r1_138, r2_138, r_139, r_873, r_199]
df_all = pd.DataFrame(all_results)

print('\n' + '=' * 80)
print('FULL RESULTS TABLE')
print('=' * 80)
print(df_all.to_string(index=False))

# Save
out_path = ROOT + '/RNA_seq_axes/scores/cross_species_tomato_paxis_scores.csv'
df_all.to_csv(out_path, index=False)
print(f'\nSaved: {out_path}')

print('\nDone.')
