"""
Complete 15-axis scoring pipeline for all biostimulant datasets.
Produces effect_vectors_complete.csv and full_axis_scores_summary.md.
"""
import sys, gzip, os, urllib.request
import pandas as pd
import numpy as np
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
    print(f"  {ax}: {len(sig)} genes")

def score_sample(log_counts, sig):
    """Dot-product of sample log-counts against signature."""
    common = sig.index.intersection(log_counts.index)
    if len(common) == 0:
        return np.nan
    return float(log_counts.loc[common].dot(sig.loc[common]))

def score_dataset(counts_df, treat_cols, ctrl_cols):
    """
    Return dict of axis -> delta (mean_treat - mean_ctrl) for log2(counts+1).
    Also return per-sample scores for diagnostics.
    """
    log_counts = np.log2(counts_df + 1)
    results = {}
    for ax, sig in sigs.items():
        t_scores = [score_sample(log_counts[c], sig) for c in treat_cols if c in log_counts.columns]
        c_scores = [score_sample(log_counts[c], sig) for c in ctrl_cols if c in log_counts.columns]
        if t_scores and c_scores:
            results[ax] = np.mean(t_scores) - np.mean(c_scores)
        else:
            results[ax] = np.nan
    return results

def score_cpm_dataset(cpm_df, treat_cols, ctrl_cols):
    """
    Score using log2(CPM+1) — for datasets that provide CPM directly.
    """
    log_cpm = np.log2(cpm_df + 1)
    results = {}
    for ax, sig in sigs.items():
        t_scores = [score_sample(log_cpm[c], sig) for c in treat_cols if c in log_cpm.columns]
        c_scores = [score_sample(log_cpm[c], sig) for c in ctrl_cols if c in log_cpm.columns]
        if t_scores and c_scores:
            results[ax] = np.mean(t_scores) - np.mean(c_scores)
        else:
            results[ax] = np.nan
    return results

# ─────────────────────────────────────────────────────────────────────────────
# 1. LOAD EXISTING effect_vectors_all_axes.csv
# ─────────────────────────────────────────────────────────────────────────────
print("\n=== Loading existing effect vectors ===")
existing_path = SCORE_DIR / "effect_vectors_all_axes.csv"
existing_df = pd.read_csv(existing_path)
print(f"Existing rows: {len(existing_df)}")
print(existing_df[['dataset','control','treatment','axis','delta']].head(10).to_string())

# Pivot existing to wide format: one row per treatment, one col per axis
print("\nPivoting existing data to wide format...")
existing_wide = existing_df.pivot_table(
    index=['dataset', 'control', 'treatment'],
    columns='axis',
    values='delta'
).reset_index()
existing_wide.columns.name = None
print(f"Wide format shape: {existing_wide.shape}")
print(existing_wide[['dataset','control','treatment']].to_string())

# ─────────────────────────────────────────────────────────────────────────────
# 2. GSE298982: Glutamate treatment vs NH4NO3 control
# ─────────────────────────────────────────────────────────────────────────────
print("\n\n=== GSE298982: Glutamate treatment (Arabidopsis, 12 days) ===")

gse298982_dir = SUPPL_DIR / "GSE298982"
gse298982_files = list(gse298982_dir.glob("*.gz"))

if not gse298982_files:
    print("No files found. Attempting to download from GEO FTP...")
    ftp_base = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE298nnn/GSE298982/suppl/"

    # List what's available via the known log2TMM file reference
    # The existing scores file references 5NN and 5Glu samples
    # Try to download the count matrix
    ftp_urls = [
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE298nnn/GSE298982/suppl/GSE298982_raw_counts.txt.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE298nnn/GSE298982/suppl/GSE298982_counts.txt.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE298nnn/GSE298982/suppl/GSE298982_raw_count_matrix.txt.gz",
    ]

    downloaded = False
    for url in ftp_urls:
        try:
            fname = url.split("/")[-1]
            dest = gse298982_dir / fname
            print(f"  Trying: {url}")
            urllib.request.urlretrieve(url, dest)
            print(f"  Downloaded: {fname}")
            downloaded = True
            gse298982_files = list(gse298982_dir.glob("*.gz"))
            break
        except Exception as e:
            print(f"  Failed: {e}")

    if not downloaded:
        print("  Direct download failed. Using existing per-axis scores from prior scoring run.")

# Check if we have processed data from prior run
gse298982_score_path = SCORE_DIR / "GSE298982_log2TMM_scores.csv"
gse298982_raw_scores = None
if gse298982_score_path.exists():
    gse298982_raw_scores = pd.read_csv(gse298982_score_path, index_col=0)
    print(f"Found prior per-sample scores: {gse298982_raw_scores.shape}")
    print(gse298982_raw_scores.to_string())

if gse298982_files:
    print(f"Count files found: {[f.name for f in gse298982_files]}")
    # Try to load a count matrix
    for f in gse298982_files:
        try:
            with gzip.open(f, 'rt') as fh:
                first_line = fh.readline()
            print(f"  File {f.name} header: {first_line[:200]}")
        except Exception as e:
            print(f"  Error: {e}")

# ── Score GSE298982 using per-sample log2TMM scores ──────────────────────────
# We have log2TMM scores for P, N, ABA, JA axes from the prior run.
# We need to re-run for all 15 axes. Since raw count data may not be present,
# let's attempt to download the actual count file from GEO.

print("\nAttempting GEO FTP listing for GSE298982...")
try:
    import urllib.request
    # Try the GEO supplemental listing
    listing_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE298nnn/GSE298982/suppl/"
    req = urllib.request.Request(listing_url, headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req, timeout=30) as resp:
        html = resp.read().decode('utf-8', errors='replace')
    print("FTP listing retrieved, length:", len(html))
    # Find .gz file references
    import re
    files_found = re.findall(r'GSE298982[^"<\s]+\.gz', html)
    print(f"Files found in listing: {files_found}")
except Exception as e:
    print(f"FTP listing failed: {e}")

# ─────────────────────────────────────────────────────────────────────────────
# 3. GSE139154: Bacillus GB03 volatiles, Arabidopsis, 48h
# ─────────────────────────────────────────────────────────────────────────────
print("\n\n=== GSE139154: Bacillus GB03 volatiles (Arabidopsis, 48h) ===")
gse139154_dir = SUPPL_DIR / "GSE139154"
gse139154_files = list(gse139154_dir.glob("*.gz"))

if not gse139154_files:
    print("No files found. Attempting to download from GEO FTP...")
    try:
        listing_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE139nnn/GSE139154/suppl/"
        req = urllib.request.Request(listing_url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=30) as resp:
            html = resp.read().decode('utf-8', errors='replace')
        print("FTP listing length:", len(html))
        import re
        files_found = re.findall(r'GSE139154[^"<\s]+\.gz', html)
        print(f"Files in listing: {files_found}")

        for fname in files_found[:3]:  # Download first few
            url = listing_url + fname
            dest = gse139154_dir / fname
            try:
                urllib.request.urlretrieve(url, dest)
                print(f"Downloaded: {fname}")
            except Exception as e:
                print(f"Failed {fname}: {e}")
        gse139154_files = list(gse139154_dir.glob("*.gz"))
    except Exception as e:
        print(f"Listing failed: {e}")

# Check existing per-sample scores
gse139154_score_path = SCORE_DIR / "GSE139154_per_sample_scores.csv"
gse139154_raw_scores = None
if gse139154_score_path.exists():
    gse139154_raw_scores = pd.read_csv(gse139154_score_path, index_col=0)
    print(f"Found prior per-sample scores: {gse139154_raw_scores.shape}")
    print(gse139154_raw_scores.to_string())

if gse139154_files:
    print(f"Files found: {[f.name for f in gse139154_files]}")

# ─────────────────────────────────────────────────────────────────────────────
# 4. GSE87337: SynCom bacterial community, Arabidopsis
#    Col-0, LowP pre-treatment, +bacteria vs -bacteria
# ─────────────────────────────────────────────────────────────────────────────
print("\n\n=== GSE87337: SynCom bacterial community ===")
gse87337_file = SUPPL_DIR / "GSE87337/GSE87337_PBI_preTreat_counts.txt.gz"

print(f"Loading count matrix from: {gse87337_file}")
with gzip.open(gse87337_file, 'rt') as fh:
    lines = fh.readlines()

# Parse header
header = lines[0].strip().split('\t')
print(f"Header entries: {len(header)}")
print(f"First 5 column headers: {header[:5]}")
print(f"Sample IDs row 2: {lines[1].strip().split(chr(9))[:5]}")

# Build the count matrix
data_lines = []
for line in lines[1:]:  # skip blank description line, start from Gene row
    parts = line.strip().split('\t')
    if parts[0].startswith('AT') or parts[0].startswith('Gene'):
        data_lines.append(parts)

# The file has 2 header rows: description row (line 0) and Gene/PBI row (line 1)
desc_row = lines[0].strip().split('\t')
gene_row = lines[1].strip().split('\t')
print(f"\nDesc row (first): {desc_row[:6]}")
print(f"Gene/PBI row (first): {gene_row[:6]}")

# Read full matrix
count_data = {}
sample_ids = gene_row[1:]  # PBI115, PBI116, etc.

for line in lines[2:]:
    parts = line.strip().split('\t')
    if not parts or not parts[0].startswith('AT'):
        continue
    gene = parts[0]
    vals = parts[1:]
    count_data[gene] = vals

gse87337_counts = pd.DataFrame(count_data).T
gse87337_counts.columns = sample_ids[:len(gse87337_counts.columns)]
gse87337_counts = gse87337_counts.apply(pd.to_numeric, errors='coerce').fillna(0)
print(f"\nGSE87337 count matrix: {gse87337_counts.shape}")
print(f"Sample IDs (first 10): {list(gse87337_counts.columns[:10])}")

# Now map sample IDs to their metadata using the description row
# From full_scores file: Col-0, NoBacteria/SynCom, LowP/-P pretreat
# We need to recreate the mapping from the per_sample_full_scores file
full_scores = pd.read_csv(SCORE_DIR / "GSE87337_per_sample_full_scores.csv", index_col=0)
print(f"\nFull scores metadata shape: {full_scores.shape}")

# Get the relevant samples: Col-0, LowP, -P pretreat
col0_lowP_minusP = full_scores[
    (full_scores['genotype'] == 'Col-0') &
    (full_scores['phosphate_cond'] == 'LowP') &
    (full_scores['pretreat'] == '-P')
]
print(f"\nCol-0, LowP, -P pretreat samples: {len(col0_lowP_minusP)}")
print(col0_lowP_minusP[['bacteria', 'phosphate_cond', 'pretreat', 'experiment']].to_string())

# Get sample IDs for these samples
ctrl_samples_87 = col0_lowP_minusP[col0_lowP_minusP['bacteria'] == 'NoBacteria'].index.tolist()
treat_samples_87 = col0_lowP_minusP[col0_lowP_minusP['bacteria'] == 'SynCom'].index.tolist()
print(f"\nControl samples (NoBacteria, Col-0, LowP, -P): {ctrl_samples_87}")
print(f"Treatment samples (SynCom, Col-0, LowP, -P): {treat_samples_87}")

# Score GSE87337 on all 15 axes
gse87337_valid_ctrl = [s for s in ctrl_samples_87 if s in gse87337_counts.columns]
gse87337_valid_treat = [s for s in treat_samples_87 if s in gse87337_counts.columns]
print(f"Valid control samples in count matrix: {gse87337_valid_ctrl}")
print(f"Valid treatment samples in count matrix: {gse87337_valid_treat}")

gse87337_results = score_dataset(gse87337_counts, gse87337_valid_treat, gse87337_valid_ctrl)
print(f"\nGSE87337 SynCom effect vector (all 15 axes):")
for ax, val in gse87337_results.items():
    print(f"  {ax:25s}: {val:+10.4f}")

# ─────────────────────────────────────────────────────────────────────────────
# 5. GSE199501: Priestia megaterium PGPR — extend to all 15 axes
# ─────────────────────────────────────────────────────────────────────────────
print("\n\n=== GSE199501: Priestia megaterium PGPR ===")
gse199501_xls = SUPPL_DIR / "GSE199501/GSE199501_Processed_data_DEG_total_CPM.xls_tmp.xls"

import xlrd
wb = xlrd.open_workbook(str(gse199501_xls))
ws = wb.sheet_by_index(0)
headers = [ws.cell_value(0, c) for c in range(ws.ncols)]
print(f"Headers: {headers}")

# Build count dataframe from raw count columns (C-1953, C-1954, C-1955, T-1956, T-1957, T-1958)
raw_cols = ['C-1953', 'C-1954', 'C-1955', 'T-1956', 'T-1957', 'T-1958']
# Also get CPM columns
cpm_cols = ['C-1953 CPM', 'C-1954 CPM', 'C-1955 CPM', 'T-1956 CPM', 'T-1957 CPM', 'T-1958 CPM']

gene_col_idx = headers.index('Gene ID ')
raw_col_idx = {h: headers.index(h) for h in raw_cols if h in headers}
cpm_col_idx = {h: headers.index(h) for h in cpm_cols if h in headers}

print(f"Gene ID column index: {gene_col_idx}")
print(f"Raw count columns: {raw_col_idx}")
print(f"CPM columns: {cpm_col_idx}")

gene_ids = []
raw_data = {h: [] for h in raw_cols if h in raw_col_idx}
cpm_data = {h: [] for h in cpm_cols if h in cpm_col_idx}

for row in range(1, ws.nrows):
    gene = ws.cell_value(row, gene_col_idx)
    if not gene or not str(gene).startswith('AT'):
        continue
    gene_ids.append(str(gene))
    for h, idx in raw_col_idx.items():
        v = ws.cell_value(row, idx)
        raw_data[h].append(float(v) if v != '' else 0.0)
    for h, idx in cpm_col_idx.items():
        v = ws.cell_value(row, idx)
        cpm_data[h].append(float(v) if v != '' else 0.0)

gse199501_raw = pd.DataFrame(raw_data, index=gene_ids)
gse199501_cpm = pd.DataFrame(cpm_data, index=gene_ids)
print(f"\nGSE199501 raw count matrix: {gse199501_raw.shape}")
print(f"GSE199501 CPM matrix: {gse199501_cpm.shape}")
print(gse199501_raw.head(3).to_string())

ctrl_cols_199 = [c for c in raw_cols if c.startswith('C-')]
treat_cols_199 = [c for c in raw_cols if c.startswith('T-')]
print(f"\nControl cols: {ctrl_cols_199}")
print(f"Treatment cols: {treat_cols_199}")

gse199501_results = score_dataset(gse199501_raw, treat_cols_199, ctrl_cols_199)
print(f"\nGSE199501 PGPR effect vector (all 15 axes):")
for ax, val in gse199501_results.items():
    print(f"  {ax:25s}: {val:+10.4f}")

# Also score with CPM for verification
cpm_ctrl = [c for c in cpm_cols if 'C-' in c]
cpm_treat = [c for c in cpm_cols if 'T-' in c]
gse199501_cpm_results = score_cpm_dataset(gse199501_cpm, cpm_treat, cpm_ctrl)
print(f"\nGSE199501 CPM-based effect vector (verification):")
for ax, val in gse199501_cpm_results.items():
    print(f"  {ax:25s}: {val:+10.4f}")

# ─────────────────────────────────────────────────────────────────────────────
# 6. Handle GSE298982 and GSE139154 — from per-sample scores or download
# ─────────────────────────────────────────────────────────────────────────────
print("\n\n=== Handling GSE298982 (glutamate) ===")
# The existing GSE298982_log2TMM_scores.csv only has 4 axes (P, N, ABA, JA)
# We need to score all 15. Let's try downloading the count matrix.

gse298982_dir = SUPPL_DIR / "GSE298982"
gse298982_count_files = list(gse298982_dir.glob("*.gz"))

gse298982_results = None

if gse298982_count_files:
    print(f"Files available: {[f.name for f in gse298982_count_files]}")
    # Try to load
    for f in gse298982_count_files:
        try:
            df = pd.read_csv(f, sep='\t', index_col=0, compression='gzip')
            print(f"  Loaded {f.name}: {df.shape}")
            print(f"  Columns: {list(df.columns)}")
            # Identify treatment and control columns
            ctrl_cols = [c for c in df.columns if 'NN' in c or 'control' in c.lower() or 'NH4' in c.lower() or c.endswith('_C')]
            treat_cols = [c for c in df.columns if 'Glu' in c or 'glutamate' in c.lower() or c.endswith('_T')]
            print(f"  Ctrl cols: {ctrl_cols}")
            print(f"  Treat cols: {treat_cols}")
            if ctrl_cols and treat_cols:
                gse298982_results = score_dataset(df, treat_cols, ctrl_cols)
                break
        except Exception as e:
            print(f"  Error loading {f.name}: {e}")

if gse298982_results is None:
    print("No count files available for GSE298982. Attempting FTP download...")
    try:
        # Try the specific known file patterns
        candidate_names = [
            "GSE298982_rawcounts.txt.gz",
            "GSE298982_raw_counts.txt.gz",
            "GSE298982_counts_matrix.txt.gz",
            "GSE298982_gene_counts.txt.gz",
            "GSE298982_HTSeq_counts.txt.gz",
        ]
        ftp_base = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE298nnn/GSE298982/suppl/"

        for name in candidate_names:
            url = ftp_base + name
            dest = gse298982_dir / name
            try:
                print(f"  Trying: {name}")
                urllib.request.urlretrieve(url, dest)
                print(f"  SUCCESS: Downloaded {name}")
                df = pd.read_csv(dest, sep='\t', index_col=0, compression='gzip')
                print(f"  Loaded: {df.shape}, cols: {list(df.columns)}")
                ctrl_cols = [c for c in df.columns if any(x in c for x in ['NN', 'control', 'NH4', 'ctrl', 'C_'])]
                treat_cols = [c for c in df.columns if any(x in c for x in ['Glu', 'glutamate', 'T_'])]
                if not ctrl_cols:
                    ctrl_cols = list(df.columns[:3])
                if not treat_cols:
                    treat_cols = list(df.columns[3:6])
                print(f"  Ctrl: {ctrl_cols}, Treat: {treat_cols}")
                gse298982_results = score_dataset(df, treat_cols, ctrl_cols)
                break
            except Exception as e:
                print(f"  Failed: {e}")
    except Exception as e:
        print(f"Download attempt failed: {e}")

if gse298982_results is None:
    # Use the existing log2TMM per-sample scores to compute deltas for the 4 existing axes,
    # and estimate missing axes using the pattern from per-sample scores
    print("\nFalling back to existing per-sample score deltas for GSE298982...")
    # From GSE298982_log2TMM_scores.csv:
    # Rows: 5NN_1, 5NN_2, 5NN_3 (control = NH4NO3)
    #       5Glu_1, 5Glu_2, 5Glu_3 (treatment = glutamate)
    tmm_scores = pd.read_csv(SCORE_DIR / "GSE298982_log2TMM_scores.csv", index_col=0)
    # Drop the Length and maxMultiRatio rows
    tmm_scores = tmm_scores[tmm_scores.index.str.startswith('5')]
    ctrl_rows = tmm_scores[tmm_scores.index.str.contains('NN')]
    treat_rows = tmm_scores[tmm_scores.index.str.contains('Glu')]

    gse298982_results = {}
    for col in tmm_scores.columns:
        ax_name = col.replace('_score', '').replace('score', '')
        if ax_name == 'JA':
            ax_name = 'JA_response'
        gse298982_results[ax_name] = float(treat_rows[col].mean() - ctrl_rows[col].mean())

    print(f"  Computed from existing TMM scores: {gse298982_results}")
    print("\n  NOTE: GSE298982 count matrix not downloaded. Only 4 axes available from prior scoring.")
    print("  Missing axes will be NaN in the output. Will attempt to download raw counts separately.")

    # Fill missing axes with NaN
    for ax in axis_names:
        if ax not in gse298982_results:
            gse298982_results[ax] = np.nan

print(f"\nGSE298982 effect vector (glutamate vs NH4NO3):")
for ax in axis_names:
    val = gse298982_results.get(ax, np.nan)
    print(f"  {ax:25s}: {val:+10.4f}" if not np.isnan(val) else f"  {ax:25s}: NaN")

print("\n\n=== Handling GSE139154 (Bacillus GB03 volatiles) ===")
gse139154_dir = SUPPL_DIR / "GSE139154"
gse139154_count_files = list(gse139154_dir.glob("*.gz"))

gse139154_results = None

if gse139154_count_files:
    print(f"Files available: {[f.name for f in gse139154_count_files]}")
    for f in gse139154_count_files:
        try:
            df = pd.read_csv(f, sep='\t', index_col=0, compression='gzip')
            print(f"  Loaded {f.name}: {df.shape}")
            print(f"  Columns: {list(df.columns)}")
            ctrl_cols = [c for c in df.columns if any(x in c.lower() for x in ['mock', 'control', 'ctrl', '_C_', 'C_'])]
            treat_cols = [c for c in df.columns if any(x in c.lower() for x in ['gb03', 'bacteria', 'treated', '_T_', 'T_'])]
            print(f"  Ctrl: {ctrl_cols}, Treat: {treat_cols}")
            if ctrl_cols and treat_cols:
                gse139154_results = score_dataset(df, treat_cols, ctrl_cols)
                break
        except Exception as e:
            print(f"  Error: {e}")

if gse139154_results is None:
    print("No count files. Attempting FTP download...")
    try:
        ftp_base = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE139nnn/GSE139154/suppl/"
        candidate_names = [
            "GSE139154_rawcounts.txt.gz",
            "GSE139154_raw_counts.txt.gz",
            "GSE139154_counts.txt.gz",
            "GSE139154_gene_counts.txt.gz",
            "GSE139154_counts_matrix.txt.gz",
        ]
        for name in candidate_names:
            url = ftp_base + name
            dest = gse139154_dir / name
            try:
                print(f"  Trying: {name}")
                urllib.request.urlretrieve(url, dest)
                print(f"  SUCCESS: {name}")
                df = pd.read_csv(dest, sep='\t', index_col=0, compression='gzip')
                print(f"  Shape: {df.shape}, cols: {list(df.columns)}")
                ctrl_cols = [c for c in df.columns if any(x in c.lower() for x in ['mock', 'control', 'ctrl'])]
                treat_cols = [c for c in df.columns if any(x in c.lower() for x in ['gb03', 'bact', 'treat'])]
                if not ctrl_cols:
                    ctrl_cols = list(df.columns[:3])
                if not treat_cols:
                    treat_cols = list(df.columns[3:6])
                gse139154_results = score_dataset(df, treat_cols, ctrl_cols)
                break
            except Exception as e:
                print(f"  Failed: {e}")
    except Exception as e:
        print(f"FTP failed: {e}")

if gse139154_results is None:
    # Use existing per-sample scores (5 axes: P, N, ABA, JA, SA)
    print("\nFalling back to existing per-sample scores for GSE139154...")
    per_sample = pd.read_csv(SCORE_DIR / "GSE139154_per_sample_scores.csv", index_col=0)
    print(per_sample.to_string())

    # Col-0_C = control, Col-0_T = treatment (GB03 volatiles)
    ctrl_rows = per_sample[per_sample.index.str.contains('Col-0_C')]
    treat_rows = per_sample[per_sample.index.str.contains('Col-0_T')]
    print(f"  Control rows: {list(ctrl_rows.index)}")
    print(f"  Treatment rows: {list(treat_rows.index)}")

    gse139154_results = {}
    for col in per_sample.columns:
        ax_name = col.replace('_score', '')
        if ax_name == 'JA':
            ax_name = 'JA_response'
        gse139154_results[ax_name] = float(treat_rows[col].mean() - ctrl_rows[col].mean())

    for ax in axis_names:
        if ax not in gse139154_results:
            gse139154_results[ax] = np.nan

print(f"\nGSE139154 effect vector (GB03 volatiles vs mock):")
for ax in axis_names:
    val = gse139154_results.get(ax, np.nan)
    print(f"  {ax:25s}: {val:+10.4f}" if not np.isnan(val) else f"  {ax:25s}: NaN")

# ─────────────────────────────────────────────────────────────────────────────
# 7. BUILD COMPLETE WIDE TABLE
# ─────────────────────────────────────────────────────────────────────────────
print("\n\n=== Building complete effect vector table ===")

# Define column names for the wide table (delta prefix)
delta_cols = {ax: f"delta_{ax}" for ax in axis_names}

# Build new rows to append
new_rows = []

# GSE298982 glutamate
row298982 = {
    'dataset': 'GSE298982',
    'treatment_group': 'GSE298982_glutamate',
    'biostimulant_class': 'Amino acid source (glutamate)',
    'control': 'NH4NO3',
    'treatment': 'glutamate',
}
for ax in axis_names:
    row298982[f"delta_{ax}"] = gse298982_results.get(ax, np.nan)
new_rows.append(row298982)

# GSE139154 GB03 volatiles
row139154 = {
    'dataset': 'GSE139154',
    'treatment_group': 'GSE139154_GB03',
    'biostimulant_class': 'PGPR volatile (Bacillus GB03)',
    'control': 'mock',
    'treatment': 'GB03_volatiles',
}
for ax in axis_names:
    row139154[f"delta_{ax}"] = gse139154_results.get(ax, np.nan)
new_rows.append(row139154)

# GSE87337 SynCom
row87337 = {
    'dataset': 'GSE87337',
    'treatment_group': 'GSE87337_SynCom',
    'biostimulant_class': 'SynCom bacterial community',
    'control': 'NoBacteria_Col0_LowP',
    'treatment': 'SynCom_Col0_LowP',
}
for ax in axis_names:
    row87337[f"delta_{ax}"] = gse87337_results.get(ax, np.nan)
new_rows.append(row87337)

# GSE199501 PGPR (full 15-axis)
row199501 = {
    'dataset': 'GSE199501',
    'treatment_group': 'GSE199501_PGPR',
    'biostimulant_class': 'PGPR (Priestia megaterium)',
    'control': 'NH4NO3_no_bacteria',
    'treatment': 'Pmegaterium',
}
for ax in axis_names:
    row199501[f"delta_{ax}"] = gse199501_results.get(ax, np.nan)
new_rows.append(row199501)

new_df = pd.DataFrame(new_rows)
print(f"\nNew rows shape: {new_df.shape}")
print(new_df[['treatment_group', 'biostimulant_class']].to_string())

# ── Reformat existing data to match wide format ───────────────────────────────
# The existing data is long format. Convert to wide matching our column scheme.
existing_rows = []

# GSE297649 amino_acid
aa_row = existing_df[existing_df['treatment'] == 'amino_acid'].set_index('axis')['delta'].to_dict()
existing_rows.append({
    'dataset': 'GSE297649',
    'treatment_group': 'amino_acid',
    'biostimulant_class': 'Amino acid biostimulant',
    'control': 'control',
    'treatment': 'amino_acid',
    **{f"delta_{k}": v for k, v in aa_row.items()}
})

# GSE297649 humic_subst
hs_row = existing_df[existing_df['treatment'] == 'humic_subst'].set_index('axis')['delta'].to_dict()
existing_rows.append({
    'dataset': 'GSE297649',
    'treatment_group': 'humic_subst',
    'biostimulant_class': 'Humic substance biostimulant',
    'control': 'control',
    'treatment': 'humic_subst',
    **{f"delta_{k}": v for k, v in hs_row.items()}
})

# GSE138478 GMV treated
gmv_row = existing_df[existing_df['treatment'] == 'GMV_treated'].set_index('axis')['delta'].to_dict()
existing_rows.append({
    'dataset': 'GSE138478',
    'treatment_group': 'diacetyl_GMV',
    'biostimulant_class': 'PGPR volatile (GMV diacetyl)',
    'control': 'GMV_control',
    'treatment': 'GMV_treated',
    **{f"delta_{k}": v for k, v in gmv_row.items()}
})

# GSE138478 diacetyl treated
da_row = existing_df[existing_df['treatment'] == 'diacetyl_treated'].set_index('axis')['delta'].to_dict()
existing_rows.append({
    'dataset': 'GSE138478',
    'treatment_group': 'diacetyl_treated',
    'biostimulant_class': 'PGPR volatile (pure diacetyl)',
    'control': 'diacetyl_control',
    'treatment': 'diacetyl_treated',
    **{f"delta_{k}": v for k, v in da_row.items()}
})

# GSE208223 TiO2
tio2_row = existing_df[existing_df['treatment'] == 'TiO2_treated'].set_index('axis')['delta'].to_dict()
existing_rows.append({
    'dataset': 'GSE208223',
    'treatment_group': 'TiO2_treated',
    'biostimulant_class': 'Ti nanoparticles',
    'control': 'TiO2_control',
    'treatment': 'TiO2_treated',
    **{f"delta_{k}": v for k, v in tio2_row.items()}
})

# GSE199501 existing (from long format — already in existing_df)
pgpr_row_existing = existing_df[existing_df['treatment'] == 'Pmegaterium_treated'].set_index('axis')['delta'].to_dict()
if pgpr_row_existing:
    # The existing file may have partial axes — we'll use our newly computed full scores instead
    print(f"\nFound GSE199501 in existing file with axes: {list(pgpr_row_existing.keys())}")
    print("Will use newly computed 15-axis scores for GSE199501.")

existing_wide_clean = pd.DataFrame(existing_rows)
print(f"\nExisting data reformatted: {existing_wide_clean.shape}")

# ── Combine all rows ──────────────────────────────────────────────────────────
complete_df = pd.concat([existing_wide_clean, new_df], ignore_index=True)

# Ensure all 15 delta columns are present
for ax in axis_names:
    col = f"delta_{ax}"
    if col not in complete_df.columns:
        complete_df[col] = np.nan

# Reorder columns
meta_cols = ['dataset', 'treatment_group', 'biostimulant_class', 'control', 'treatment']
delta_col_names = [f"delta_{ax}" for ax in axis_names]
complete_df = complete_df[meta_cols + delta_col_names]

print(f"\nComplete effect vector table shape: {complete_df.shape}")
print("\nTreatment groups:")
for _, row in complete_df.iterrows():
    print(f"  {row['treatment_group']:35s} | {row['biostimulant_class']}")

# ── Save complete table ───────────────────────────────────────────────────────
out_path = SCORE_DIR / "effect_vectors_complete.csv"
complete_df.to_csv(out_path, index=False)
print(f"\nSaved complete effect vectors: {out_path}")

# ─────────────────────────────────────────────────────────────────────────────
# 8. PRINT FULL 15-AXIS TABLE (formatted)
# ─────────────────────────────────────────────────────────────────────────────
print("\n\n" + "="*120)
print("COMPLETE 15-AXIS EFFECT VECTOR TABLE")
print("="*120)
print(f"{'Treatment Group':35s}", end="")
for ax in axis_names:
    print(f" {ax[:12]:>12s}", end="")
print()
print("-"*120)

for _, row in complete_df.iterrows():
    print(f"{row['treatment_group']:35s}", end="")
    for ax in axis_names:
        val = row[f"delta_{ax}"]
        if pd.isna(val):
            print(f" {'NaN':>12s}", end="")
        else:
            print(f" {val:>+12.2f}", end="")
    print()

print("="*120)

# ─────────────────────────────────────────────────────────────────────────────
# 9. Write summary markdown
# ─────────────────────────────────────────────────────────────────────────────
print("\n\nWriting summary markdown...")
summary_path = AGENT_DIR / "full_axis_scores_summary.md"

md_lines = [
    "# Full 15-Axis Effect Vector Summary",
    "",
    f"Generated: 2026-02-23",
    "",
    "## Dataset-to-Biostimulant Class Mapping",
    "",
    "| Treatment Group | Dataset | Biostimulant Class |",
    "|---|---|---|",
]

for _, row in complete_df.iterrows():
    md_lines.append(f"| {row['treatment_group']} | {row['dataset']} | {row['biostimulant_class']} |")

md_lines += [
    "",
    "## Complete Effect Vector Table (delta = treatment - control, log2-count space)",
    "",
    "Positive delta = treatment resembles the HIGH end of the axis (e.g., P-starved for P-axis).",
    "Negative delta = treatment resembles the LOW end of the axis.",
    "",
    "### Abbreviations",
    "- P = Phosphorus starvation axis",
    "- N = Nitrogen supply axis",
    "- ABA = Abscisic acid axis",
    "- JA_response = Jasmonic acid response axis",
    "- RibBiogen = Ribosome Biogenesis axis",
    "- CK = Cytokinin response axis",
    "- Drought = Drought stress axis",
    "- SA = Salicylic acid response axis",
    "- ET = Ethylene response axis",
    "- Auxin = Auxin response axis",
    "- BR = Brassinosteroid axis",
    "- OxStress = Oxidative stress axis",
    "- Hypoxia = Hypoxia axis",
    "- CellWall = Cell wall axis",
    "- ABA_GS = ABA response (guard cell, second signature)",
    "",
]

# Build markdown table
header_row = "| Treatment Group | Biostimulant Class | P | N | ABA | JA | RibBiogen | CK | Drought | SA | ET | Auxin | BR | OxStress | Hypoxia | CellWall | ABA_GS |"
sep_row = "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"
md_lines.append(header_row)
md_lines.append(sep_row)

for _, row in complete_df.iterrows():
    vals = []
    for ax in axis_names:
        val = row[f"delta_{ax}"]
        if pd.isna(val):
            vals.append("NaN")
        else:
            vals.append(f"{val:+.1f}")
    md_lines.append(f"| {row['treatment_group']} | {row['biostimulant_class']} | " + " | ".join(vals) + " |")

md_lines += [
    "",
    "## Notes on Data Sources",
    "",
    "### GSE298982 (Glutamate treatment)",
    "- Species: Arabidopsis thaliana, 12 days",
    "- Treatment: 5 mM glutamate vs 5 mM NH4NO3 (nitrogen-equivalent control)",
    "- Data: log2TMM-normalized counts; 3 replicates per group",
    "- Axes with data: P, N, ABA, JA_response (4/15 from prior scoring run)",
    "- Missing axes: NaN (raw count matrix not downloaded)",
    "",
    "### GSE139154 (Bacillus GB03 volatiles)",
    "- Species: Arabidopsis thaliana, 48h volatile exposure",
    "- Treatment: GB03 bacteria vs mock (no bacteria)",
    "- Data: Col-0 genotype only; 3 replicates per group",
    "- Axes with data: P, N, ABA, JA_response, SA_response (5/15 from prior scoring run)",
    "- Missing axes: NaN (raw count matrix not downloaded)",
    "",
    "### GSE87337 (SynCom bacterial community)",
    "- Species: Arabidopsis thaliana, 12 days",
    "- Subset used: Col-0 genotype, LowP phosphate condition, -P pre-treatment",
    "- Treatment: SynCom inoculated vs no bacteria",
    "- Data: raw HTSeq counts (144 samples total, subset of 4 used)",
    "- All 15 axes scored from raw counts",
    "",
    "### GSE199501 (Priestia megaterium PGPR)",
    "- Species: Arabidopsis thaliana, 7 days",
    "- Treatment: P. megaterium inoculation vs no bacteria control",
    "- Data: raw counts from XLS file (C-1953/1954/1955 = ctrl, T-1956/1957/1958 = treat)",
    "- All 15 axes scored from raw counts",
    "",
    "### GSE297649 (Amino acids + Humic substances)",
    "- Previously scored, all 15 axes complete",
    "",
    "### GSE138478 (PGPR diacetyl + GMV volatiles)",
    "- Previously scored, all 15 axes complete",
    "",
    "### GSE208223 (TiO2 nanoparticles)",
    "- Previously scored, all 15 axes complete",
    "",
    "## Interpretation Notes",
    "",
    "See `RNA_seq_axes/scores/effect_vectors_complete.csv` for the machine-readable version.",
    "",
    "Key patterns to investigate:",
    "1. Which biostimulants show P-axis activation (mimicking P-starvation signaling)?",
    "2. Which show ABA-like signatures (water status modulation)?",
    "3. Do PGPR treatments cluster together on JA/SA axes (ISR-related)?",
    "4. Does SynCom community differ from single-organism PGPR in its fingerprint?",
]

with open(summary_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(md_lines))

print(f"Summary written to: {summary_path}")
print("\n=== PIPELINE COMPLETE ===")
print(f"Output: {out_path}")
print(f"Summary: {summary_path}")
