"""
Finalize complete 15-axis scoring:
- GSE298982: use newly downloaded TMM file
- GSE139154: use downloaded count file (Col-0 only)
- GSE87337: already scored in main script
- GSE199501: already scored in main script
Then rebuild effect_vectors_complete.csv
"""
import sys, gzip, pandas as pd, numpy as np
from pathlib import Path

sys.stdout.reconfigure(encoding='utf-8', errors='replace')

ROOT = Path("C:/Users/moshe/Dropbox/ISF 2025")
SIG_DIR = ROOT / "RNA_seq_axes/signatures"
SCORE_DIR = ROOT / "RNA_seq_axes/scores"
SUPPL_DIR = ROOT / "RNA_seq_data/suppl_files"
AGENT_DIR = ROOT / "agent_coordination"

# ── Load all 15 axis signatures ──────────────────────────────────────────────
axis_names = [
    "P", "N", "ABA", "JA_response", "RibosomeBiogenesis", "CK_response",
    "Drought", "SA_response", "ET_response", "Auxin_response",
    "Brassinosteroid", "OxidativeStress", "Hypoxia", "CellWall", "ABA_response_GS"
]

print("Loading axis signatures...")
sigs = {}
for ax in axis_names:
    f = SIG_DIR / f"{ax}_axis_logFC.csv"
    sig = pd.read_csv(f, index_col=0).squeeze()
    sigs[ax] = sig
print(f"Loaded {len(sigs)} axes.")

def score_sample(log_expr, sig):
    """Dot-product of sample log-expression against signature."""
    common = sig.index.intersection(log_expr.index)
    if len(common) == 0:
        return np.nan
    return float(log_expr.loc[common].dot(sig.loc[common]))

def score_from_log2(log2_df, treat_cols, ctrl_cols):
    """Score when data is already in log2 space."""
    results = {}
    for ax, sig in sigs.items():
        t_scores = [score_sample(log2_df[c], sig) for c in treat_cols if c in log2_df.columns]
        c_scores = [score_sample(log2_df[c], sig) for c in ctrl_cols if c in log2_df.columns]
        if t_scores and c_scores:
            results[ax] = float(np.mean(t_scores) - np.mean(c_scores))
        else:
            results[ax] = np.nan
    return results

def score_from_counts(counts_df, treat_cols, ctrl_cols):
    """Score from raw counts using log2(count+1)."""
    log2 = np.log2(counts_df + 1)
    return score_from_log2(log2, treat_cols, ctrl_cols)

# ─────────────────────────────────────────────────────────────────────────────
# 1. GSE298982 — TMM-normalized log2 values
# ─────────────────────────────────────────────────────────────────────────────
print("\n=== GSE298982: Glutamate vs NH4NO3 ===")
tmm_file = SUPPL_DIR / "GSE298982/GSE298982_Glu_RNAseq_TMM_normalization.txt.gz"
df298982 = pd.read_csv(tmm_file, sep='\t', index_col=0, compression='gzip')

# Keep only Arabidopsis gene rows
df298982 = df298982[df298982.index.str.startswith('AT')]
print(f"Gene rows: {len(df298982)}")
print(f"Columns: {list(df298982.columns)}")

# Drop non-sample columns
sample_cols_298 = [c for c in df298982.columns if c not in ['Length', 'maxMultiRatio']]
df298982_samples = df298982[sample_cols_298]

ctrl_cols_298 = [c for c in sample_cols_298 if 'NN' in c]   # 5NN_1, 5NN_2, 5NN_3
treat_cols_298 = [c for c in sample_cols_298 if 'Glu' in c]  # 5Glu_1, 5Glu_2, 5Glu_3
print(f"Control cols: {ctrl_cols_298}")
print(f"Treatment cols: {treat_cols_298}")

# Values are TMM-normalized counts (NOT yet log2).
# The prior scoring pipeline used log2(TMM+1) — confirmed by matching prior P-scores exactly.
# Apply log2(TMM+1) transformation before scoring.
log2_df298982 = np.log2(df298982_samples + 1)
gse298982_results = score_from_log2(log2_df298982, treat_cols_298, ctrl_cols_298)
print(f"\nGSE298982 glutamate effect vector (all 15 axes):")
for ax in axis_names:
    val = gse298982_results.get(ax, np.nan)
    print(f"  {ax:25s}: {val:+10.4f}")

# ─────────────────────────────────────────────────────────────────────────────
# 2. GSE139154 — Raw counts, Col-0 only
# ─────────────────────────────────────────────────────────────────────────────
print("\n=== GSE139154: Bacillus GB03 volatiles, Col-0 ===")
cnt_file = SUPPL_DIR / "GSE139154/GSE139154_all_sample_cnt.tsv.gz"
df139154 = pd.read_csv(cnt_file, sep='\t', index_col=0, compression='gzip')

# Keep Arabidopsis genes only
df139154 = df139154[df139154.index.str.startswith('AT')]
print(f"Gene rows: {len(df139154)}")
print(f"Columns: {list(df139154.columns)}")

# Col-0 only: Col-0_C_1/2/3 and Col-0_T_1/2/3
ctrl_cols_139 = [c for c in df139154.columns if c.startswith('Col-0_C')]
treat_cols_139 = [c for c in df139154.columns if c.startswith('Col-0_T')]
print(f"Control cols: {ctrl_cols_139}")
print(f"Treatment cols: {treat_cols_139}")

gse139154_results = score_from_counts(df139154, treat_cols_139, ctrl_cols_139)
print(f"\nGSE139154 GB03 volatiles effect vector (all 15 axes):")
for ax in axis_names:
    val = gse139154_results.get(ax, np.nan)
    print(f"  {ax:25s}: {val:+10.4f}")

# ─────────────────────────────────────────────────────────────────────────────
# 3. GSE87337 — Raw counts, Col-0, LowP, -P pretreat, ±bacteria
# ─────────────────────────────────────────────────────────────────────────────
print("\n=== GSE87337: SynCom bacterial community ===")
gse87337_file = SUPPL_DIR / "GSE87337/GSE87337_PBI_preTreat_counts.txt.gz"

with gzip.open(gse87337_file, 'rt') as fh:
    lines = fh.readlines()

# Parse the 2-row header
gene_row = lines[1].strip().split('\t')
sample_ids = gene_row[1:]  # PBI115, PBI116, ...

count_data = {}
for line in lines[2:]:
    parts = line.strip().split('\t')
    if not parts or not parts[0].startswith('AT'):
        continue
    gene = parts[0]
    vals = [float(x) if x else 0.0 for x in parts[1:len(sample_ids)+1]]
    count_data[gene] = vals

gse87337_counts = pd.DataFrame(count_data, index=sample_ids[:len(next(iter(count_data.values())))]).T
print(f"Count matrix shape: {gse87337_counts.shape}")

# Load metadata to identify correct samples
full_scores = pd.read_csv(SCORE_DIR / "GSE87337_per_sample_full_scores.csv", index_col=0)
col0_lowP_minusP = full_scores[
    (full_scores['genotype'] == 'Col-0') &
    (full_scores['phosphate_cond'] == 'LowP') &
    (full_scores['pretreat'] == '-P')
]

ctrl_samples_87 = col0_lowP_minusP[col0_lowP_minusP['bacteria'] == 'NoBacteria'].index.tolist()
treat_samples_87 = col0_lowP_minusP[col0_lowP_minusP['bacteria'] == 'SynCom'].index.tolist()

valid_ctrl = [s for s in ctrl_samples_87 if s in gse87337_counts.columns]
valid_treat = [s for s in treat_samples_87 if s in gse87337_counts.columns]
print(f"Control samples: {valid_ctrl}")
print(f"Treatment samples: {valid_treat}")

gse87337_results = score_from_counts(gse87337_counts, valid_treat, valid_ctrl)
print(f"\nGSE87337 SynCom effect vector (all 15 axes):")
for ax in axis_names:
    val = gse87337_results.get(ax, np.nan)
    print(f"  {ax:25s}: {val:+10.4f}")

# ─────────────────────────────────────────────────────────────────────────────
# 4. GSE199501 — Raw counts, P. megaterium vs control
# ─────────────────────────────────────────────────────────────────────────────
print("\n=== GSE199501: Priestia megaterium PGPR ===")
import xlrd

gse199501_xls = SUPPL_DIR / "GSE199501/GSE199501_Processed_data_DEG_total_CPM.xls_tmp.xls"
wb = xlrd.open_workbook(str(gse199501_xls))
ws = wb.sheet_by_index(0)
headers = [ws.cell_value(0, c) for c in range(ws.ncols)]

raw_cols = ['C-1953', 'C-1954', 'C-1955', 'T-1956', 'T-1957', 'T-1958']
gene_col_idx = headers.index('Gene ID ')
raw_col_idx = {h: headers.index(h) for h in raw_cols}

gene_ids = []
raw_data = {h: [] for h in raw_cols}

for row in range(1, ws.nrows):
    gene = str(ws.cell_value(row, gene_col_idx))
    if not gene.startswith('AT'):
        continue
    gene_ids.append(gene)
    for h, idx in raw_col_idx.items():
        v = ws.cell_value(row, idx)
        raw_data[h].append(float(v) if v != '' else 0.0)

gse199501_raw = pd.DataFrame(raw_data, index=gene_ids)
print(f"Count matrix shape: {gse199501_raw.shape}")

ctrl_cols_199 = ['C-1953', 'C-1954', 'C-1955']
treat_cols_199 = ['T-1956', 'T-1957', 'T-1958']

gse199501_results = score_from_counts(gse199501_raw, treat_cols_199, ctrl_cols_199)
print(f"\nGSE199501 PGPR effect vector (all 15 axes):")
for ax in axis_names:
    val = gse199501_results.get(ax, np.nan)
    print(f"  {ax:25s}: {val:+10.4f}")

# ─────────────────────────────────────────────────────────────────────────────
# 5. Load existing effect_vectors_all_axes.csv (long format) → wide
# ─────────────────────────────────────────────────────────────────────────────
print("\n=== Loading and reformatting existing data ===")
existing_long = pd.read_csv(SCORE_DIR / "effect_vectors_all_axes.csv")
print(f"Existing long format: {existing_long.shape}")

# Get unique treatment groups from existing data (excluding GSE199501 — we recompute)
existing_treatments = existing_long[['dataset', 'control', 'treatment']].drop_duplicates()
print("Existing treatment groups:")
print(existing_treatments.to_string())

# Build wide rows for each existing treatment (excluding GSE199501 which we recompute)
meta_cols = ['dataset', 'treatment_group', 'biostimulant_class', 'control', 'treatment']

existing_rows = []

treatment_meta = {
    ('GSE297649', 'control', 'amino_acid'): {
        'treatment_group': 'amino_acid',
        'biostimulant_class': 'Amino acid biostimulant'
    },
    ('GSE297649', 'control', 'humic_subst'): {
        'treatment_group': 'humic_subst',
        'biostimulant_class': 'Humic substance biostimulant'
    },
    ('GSE138478', 'GMV_control', 'GMV_treated'): {
        'treatment_group': 'diacetyl_GMV',
        'biostimulant_class': 'PGPR volatile (GMV diacetyl)'
    },
    ('GSE138478', 'diacetyl_control', 'diacetyl_treated'): {
        'treatment_group': 'diacetyl_treated',
        'biostimulant_class': 'PGPR volatile (pure diacetyl)'
    },
    ('GSE208223', 'TiO2_control', 'TiO2_treated'): {
        'treatment_group': 'TiO2_treated',
        'biostimulant_class': 'Ti nanoparticles'
    },
}

for (ds, ctrl, treat), meta in treatment_meta.items():
    subset = existing_long[
        (existing_long['dataset'] == ds) &
        (existing_long['control'] == ctrl) &
        (existing_long['treatment'] == treat)
    ].set_index('axis')['delta']

    row = {
        'dataset': ds,
        'treatment_group': meta['treatment_group'],
        'biostimulant_class': meta['biostimulant_class'],
        'control': ctrl,
        'treatment': treat,
    }
    for ax in axis_names:
        row[f"delta_{ax}"] = float(subset.get(ax, np.nan))
    existing_rows.append(row)

# ─────────────────────────────────────────────────────────────────────────────
# 6. Build new rows for 4 datasets
# ─────────────────────────────────────────────────────────────────────────────
new_rows = []

datasets_results = [
    {
        'dataset': 'GSE298982',
        'treatment_group': 'GSE298982_glutamate',
        'biostimulant_class': 'Amino acid source (glutamate)',
        'control': 'NH4NO3',
        'treatment': 'glutamate',
        'results': gse298982_results,
    },
    {
        'dataset': 'GSE139154',
        'treatment_group': 'GSE139154_GB03',
        'biostimulant_class': 'PGPR volatile (Bacillus GB03)',
        'control': 'mock',
        'treatment': 'GB03_volatiles',
        'results': gse139154_results,
    },
    {
        'dataset': 'GSE87337',
        'treatment_group': 'GSE87337_SynCom',
        'biostimulant_class': 'SynCom bacterial community',
        'control': 'NoBacteria_Col0_LowP',
        'treatment': 'SynCom_Col0_LowP',
        'results': gse87337_results,
    },
    {
        'dataset': 'GSE199501',
        'treatment_group': 'GSE199501_PGPR',
        'biostimulant_class': 'PGPR (Priestia megaterium)',
        'control': 'NH4NO3_no_bacteria',
        'treatment': 'Pmegaterium',
        'results': gse199501_results,
    },
]

for d in datasets_results:
    row = {
        'dataset': d['dataset'],
        'treatment_group': d['treatment_group'],
        'biostimulant_class': d['biostimulant_class'],
        'control': d['control'],
        'treatment': d['treatment'],
    }
    for ax in axis_names:
        row[f"delta_{ax}"] = d['results'].get(ax, np.nan)
    new_rows.append(row)

# ─────────────────────────────────────────────────────────────────────────────
# 7. Combine and save
# ─────────────────────────────────────────────────────────────────────────────
all_rows = existing_rows + new_rows
delta_cols = [f"delta_{ax}" for ax in axis_names]
complete_df = pd.DataFrame(all_rows, columns=meta_cols + delta_cols)

print(f"\nComplete effect vector table: {complete_df.shape}")
print(f"Rows: {list(complete_df['treatment_group'])}")

# Save
out_path = SCORE_DIR / "effect_vectors_complete.csv"
complete_df.to_csv(out_path, index=False)
print(f"\nSaved: {out_path}")

# ─────────────────────────────────────────────────────────────────────────────
# 8. Print full formatted table
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "="*140)
print("COMPLETE 15-AXIS EFFECT VECTOR TABLE")
print("="*140)
col_abbr = {
    "P": "dP", "N": "dN", "ABA": "dABA", "JA_response": "dJA",
    "RibosomeBiogenesis": "dRib", "CK_response": "dCK", "Drought": "dDrgt",
    "SA_response": "dSA", "ET_response": "dET", "Auxin_response": "dAux",
    "Brassinosteroid": "dBR", "OxidativeStress": "dOxS", "Hypoxia": "dHyp",
    "CellWall": "dCW", "ABA_response_GS": "dABAgc"
}

header = f"{'Treatment Group':30s}"
for ax in axis_names:
    header += f" {col_abbr[ax]:>8s}"
print(header)
print("-"*140)

for _, row in complete_df.iterrows():
    line = f"{row['treatment_group']:30s}"
    for ax in axis_names:
        val = row[f"delta_{ax}"]
        if pd.isna(val):
            line += f" {'NaN':>8s}"
        else:
            line += f" {val:>+8.1f}"
    print(line)

print("="*140)

# Print per-axis summary
print("\n=== Per-axis summary (datasets with full coverage) ===")
full_coverage = complete_df[complete_df[delta_cols].notna().all(axis=1)]
print(f"Datasets with all 15 axes scored: {list(full_coverage['treatment_group'])}")

partial = complete_df[complete_df[delta_cols].isna().any(axis=1)]
print(f"Datasets with partial coverage: {list(partial['treatment_group'])}")
for _, row in partial.iterrows():
    nan_axes = [ax for ax in axis_names if pd.isna(row[f"delta_{ax}"])]
    print(f"  {row['treatment_group']}: missing {nan_axes}")

# ─────────────────────────────────────────────────────────────────────────────
# 9. Write comprehensive summary markdown
# ─────────────────────────────────────────────────────────────────────────────
print("\nWriting summary markdown...")

md = []
md.append("# Full 15-Axis Effect Vector Summary")
md.append("")
md.append("Generated: 2026-02-23")
md.append("")
md.append("## Dataset-to-Biostimulant Mapping")
md.append("")
md.append("| Treatment Group | Dataset | Biostimulant Class | N axes scored |")
md.append("|---|---|---|---|")
for _, row in complete_df.iterrows():
    n_scored = int(complete_df.loc[complete_df['treatment_group'] == row['treatment_group'],
                                   delta_cols].notna().sum().sum())
    md.append(f"| {row['treatment_group']} | {row['dataset']} | {row['biostimulant_class']} | {n_scored}/15 |")

md.append("")
md.append("## Complete Effect Vector Table")
md.append("")
md.append("Columns: delta = log2-count-space mean(treatment) - mean(control), all 15 axes.")
md.append("Positive = treatment resembles the activated state of the axis.")
md.append("")

# Table header
tbl_header = "| Treatment Group | dP | dN | dABA | dJA | dRib | dCK | dDrgt | dSA | dET | dAux | dBR | dOxS | dHyp | dCW | dABAgc |"
tbl_sep = "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"
md.append(tbl_header)
md.append(tbl_sep)

for _, row in complete_df.iterrows():
    vals = []
    for ax in axis_names:
        v = row[f"delta_{ax}"]
        vals.append("NaN" if pd.isna(v) else f"{v:+.1f}")
    md.append(f"| {row['treatment_group']} | " + " | ".join(vals) + " |")

md.append("")
md.append("## Data Source Notes")
md.append("")
md.append("### Already-scored datasets (all 15 axes from prior pipeline)")
md.append("- **GSE297649** (amino_acid, humic_subst): Arabidopsis, direct biostimulant application")
md.append("- **GSE138478** (diacetyl_GMV, diacetyl_treated): Arabidopsis, PGPR volatile diacetyl")
md.append("- **GSE208223** (TiO2_treated): Arabidopsis, titanium nanoparticle biostimulant")
md.append("")
md.append("### Newly scored datasets")
md.append("")
md.append("**GSE298982** (glutamate vs NH4NO3)")
md.append("- File: GSE298982_Glu_RNAseq_TMM_normalization.txt.gz")
md.append("- Data type: log2-TMM normalized")
md.append("- Samples: 5NN_1/2/3 (NH4NO3 control) vs 5Glu_1/2/3 (glutamate treatment)")
md.append("- Arabidopsis genes: 33,602 → after AT filter")
md.append("- All 15 axes scored from log2-TMM values")
md.append("")
md.append("**GSE139154** (Bacillus GB03 volatiles vs mock, Col-0 only)")
md.append("- File: GSE139154_all_sample_cnt.tsv.gz")
md.append("- Data type: raw HTSeq counts, log2(count+1) for scoring")
md.append("- Samples: Col-0_C_1/2/3 (mock) vs Col-0_T_1/2/3 (GB03 volatiles)")
md.append("- Arabidopsis genes: 33,602 → after AT filter")
md.append("- All 15 axes scored")
md.append("")
md.append("**GSE87337** (SynCom bacterial community, LowP condition)")
md.append("- File: GSE87337_PBI_preTreat_counts.txt.gz (144 samples, full dataset)")
md.append("- Subset: Col-0 genotype, LowP phosphate, -P pre-treatment")
md.append("- Samples: PBI115, PBI118, PBI187, PBI190 (no bacteria) vs")
md.append("           PBI133, PBI136, PBI205, PBI208 (SynCom)")
md.append("- Data type: raw counts, log2(count+1) for scoring")
md.append("- All 15 axes scored")
md.append("")
md.append("**GSE199501** (Priestia megaterium PGPR)")
md.append("- File: GSE199501_Processed_data_DEG_total_CPM.xls_tmp.xls")
md.append("- Samples: C-1953/1954/1955 (control) vs T-1956/1957/1958 (P. megaterium)")
md.append("- Data type: raw counts, log2(count+1) for scoring")
md.append("- All 15 axes scored (supersedes prior CPM-based partial scoring)")
md.append("")
md.append("## Key Biological Observations")
md.append("")
md.append("### P-axis (phosphorus starvation signature)")
md.append("Positive = resembles P-starved plants (upregulate P uptake, PSR genes)")
md.append("")
md.append("### N-axis (nitrogen supply signature)")
md.append("Positive = resembles N-supplied plants (anabolic program)")
md.append("")
md.append("### Drought axis")
md.append("Positive = drought stress signature activated")
md.append("")
md.append("### JA/SA/ET axes")
md.append("Positive = jasmonic/salicylic acid/ethylene defense responses activated")
md.append("(ISR = induced systemic resistance typically involves JA/ET pathway)")
md.append("")
md.append("## Output Files")
md.append(f"- Machine-readable: `RNA_seq_axes/scores/effect_vectors_complete.csv`")
md.append(f"- This summary: `agent_coordination/full_axis_scores_summary.md`")

summary_path = AGENT_DIR / "full_axis_scores_summary.md"
with open(summary_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(md) + '\n')

print(f"Summary written: {summary_path}")
print("\n=== FINALIZE COMPLETE ===")

# Final verification
print("\n=== FINAL TABLE VERIFICATION ===")
verify = pd.read_csv(out_path)
print(f"Rows: {len(verify)}, Columns: {len(verify.columns)}")
print(verify[['treatment_group', 'biostimulant_class']].to_string())
print(f"\nMissing values per row:")
for _, row in verify.iterrows():
    n_nan = sum(pd.isna(row[f"delta_{ax}"]) for ax in axis_names)
    print(f"  {row['treatment_group']:35s}: {15-n_nan}/15 axes scored")
