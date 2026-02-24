"""
Re-score GSE138478 (PGPR) using log2(CPM+1) to check library size effects.
featureCounts format: Geneid, Chr, Start, End, Strand, Length, sample1, sample2, ...
"""
import gzip, io
import pandas as pd
import numpy as np
from pathlib import Path

ROOT = Path("C:/Users/moshe/Dropbox/ISF 2025")
SUPPL = ROOT / "RNA_seq_data/suppl_files/GSE138478"
SIG_DIR = ROOT / "RNA_seq_axes/signatures"
CAL_FILE = ROOT / "RNA_seq_axes/calibration.csv"

# ── Read featureCounts files ──────────────────────────────────────────────────
def read_feature_counts(gz_path):
    """Read featureCounts file. Skip comment lines, extract gene_id + counts."""
    with gzip.open(gz_path, "rt") as f:
        lines = f.readlines()
    # Skip lines starting with '#'
    data_lines = [l for l in lines if not l.startswith("#")]
    from io import StringIO
    df = pd.read_csv(StringIO("".join(data_lines)), sep="\t")
    print(f"  Columns: {df.columns.tolist()}")
    # First col = Geneid, cols 1-5 = Chr/Start/End/Strand/Length, rest = samples
    sample_cols = df.columns[6:].tolist()
    print(f"  Samples ({len(sample_cols)}): {sample_cols}")
    df = df.set_index("Geneid")[sample_cols]
    return df

print("=== Reading HMZ008 (PGPR volatile GMV) ===")
hmz008 = read_feature_counts(SUPPL / "GSE138478_HMZ008_featureCounts_output.txt.gz")
print(f"  Shape: {hmz008.shape}")

print("\n=== Reading HMZ014 (PGPR diacetyl) ===")
hmz014 = read_feature_counts(SUPPL / "GSE138478_HMZ014_Dione_featureCounts_output.txt.gz")
print(f"  Shape: {hmz014.shape}")

# ── Library sizes ─────────────────────────────────────────────────────────────
print("\n=== Library Sizes ===")
for name, df in [("HMZ008", hmz008), ("HMZ014_Dione", hmz014)]:
    lib = df.sum(axis=0)
    median = lib.median()
    print(f"\n{name}:")
    for col in df.columns:
        ratio = lib[col] / median
        flag = " *** DEPTH OUTLIER ***" if ratio > 1.5 or ratio < 0.67 else ""
        print(f"  {col}: {lib[col]:,}  (ratio={ratio:.2f}){flag}")

# ── Group assignments ─────────────────────────────────────────────────────────
# HMZ008: column names are like "GSE138478_HMZ008_featureCounts_output__HMZ008_0.05C_5DAT_rep1"
# We need to shorten these and assign groups

def assign_groups_hmz008(col):
    # Extract the part after the last __
    short = col.split("__")[-1] if "__" in col else col
    if "_C_" in short or "C_5DAT" in short:
        return "GMV_control", short
    else:
        return "GMV_treated", short

def assign_groups_hmz014(col):
    short = col.split("__")[-1] if "__" in col else col
    if "CK_" in short:
        return "diacetyl_control", short
    else:
        return "diacetyl_treated", short

groups_008 = {col: assign_groups_hmz008(col) for col in hmz008.columns}
groups_014 = {col: assign_groups_hmz014(col) for col in hmz014.columns}

print("\n=== Group assignments ===")
for col, (grp, short) in groups_008.items():
    print(f"  HMZ008: {short} -> {grp}")
for col, (grp, short) in groups_014.items():
    print(f"  HMZ014: {short} -> {grp}")

# ── Calibration ───────────────────────────────────────────────────────────────
cal = pd.read_csv(CAL_FILE, index_col=0)

def calibrate(score, axis):
    lo = cal.loc[axis, "ctrl_score"]
    hi = cal.loc[axis, "treat_score"]
    return (score - lo) / (hi - lo) * 100

# ── Score on each axis ────────────────────────────────────────────────────────
axes = ["P_axis", "N_axis", "ABA_axis"]

def score_df(count_df, groups_map, dataset_label):
    """Score a count matrix (raw and CPM) on all 3 axes."""
    lib = count_df.sum(axis=0)
    cpm = count_df.divide(lib, axis=1) * 1_000_000
    log_raw = np.log2(count_df + 1)
    log_cpm = np.log2(cpm + 1)

    results = []
    for axis in axes:
        sig = pd.read_csv(SIG_DIR / f"{axis}_logFC.csv", index_col=0)
        sig.columns = ["logFC"]
        common = sig.index.intersection(count_df.index)
        sig_c = sig.loc[common, "logFC"]
        print(f"\n  {axis}: {len(sig_c)} genes in common with {dataset_label}")

        for col in count_df.columns:
            grp, short = groups_map[col]
            raw_score = (log_raw.loc[common, col] * sig_c).sum()
            cpm_score = (log_cpm.loc[common, col] * sig_c).sum()
            results.append({
                "dataset": dataset_label,
                "sample": short,
                "group": grp,
                "axis": axis,
                "lib_size": lib[col],
                "raw_score": raw_score,
                "raw_pct": calibrate(raw_score, axis),
                "cpm_score": cpm_score,
                "cpm_pct": calibrate(cpm_score, axis),
            })
    return pd.DataFrame(results)

print("\n=== Scoring HMZ008 ===")
res008 = score_df(hmz008, groups_008, "HMZ008_GMV")
print("\n=== Scoring HMZ014 ===")
res014 = score_df(hmz014, groups_014, "HMZ014_diacetyl")

res_all = pd.concat([res008, res014], ignore_index=True)

# ── Compare raw vs CPM effect vectors ────────────────────────────────────────
print("\n\n=== EFFECT VECTORS: RAW vs CPM ===")
comparisons = [
    ("HMZ008_GMV",      "GMV_treated",      "GMV_control",      "PGPR GMV volatile"),
    ("HMZ014_diacetyl", "diacetyl_treated", "diacetyl_control", "PGPR diacetyl"),
]

for dataset, treat_grp, ctrl_grp, label in comparisons:
    d = res_all[res_all.dataset == dataset]
    print(f"\n--- {label} ---")
    print(f"  Library sizes (control): ", end="")
    ctrl_libs = d[d.group == ctrl_grp]["lib_size"].unique()
    for l in ctrl_libs:
        print(f"{l:,} ", end="")
    print(f"\n  Library sizes (treated): ", end="")
    treat_libs = d[d.group == treat_grp]["lib_size"].unique()
    for l in treat_libs:
        print(f"{l:,} ", end="")
    print()

    for axis in axes:
        raw_ctrl  = d[(d.axis==axis) & (d.group==ctrl_grp)]["raw_pct"].mean()
        raw_treat = d[(d.axis==axis) & (d.group==treat_grp)]["raw_pct"].mean()
        cpm_ctrl  = d[(d.axis==axis) & (d.group==ctrl_grp)]["cpm_pct"].mean()
        cpm_treat = d[(d.axis==axis) & (d.group==treat_grp)]["cpm_pct"].mean()
        raw_delta = raw_treat - raw_ctrl
        cpm_delta = cpm_treat - cpm_ctrl
        flag = ""
        if (raw_delta * cpm_delta < 0):
            flag = " *** DIRECTION FLIP! ***"
        elif abs(cpm_delta) < abs(raw_delta) * 0.3 and abs(raw_delta) > 1:
            flag = " ** >70% reduction after CPM **"
        print(f"  {axis:12}: raw_delta={raw_delta:+.2f}%  cpm_delta={cpm_delta:+.2f}%{flag}")

# ── Save results ──────────────────────────────────────────────────────────────
out = ROOT / "RNA_seq_axes/scores/GSE138478_per_sample_scores_CPM.csv"
res_all.to_csv(out, index=False)
print(f"\nSaved: {out}")
