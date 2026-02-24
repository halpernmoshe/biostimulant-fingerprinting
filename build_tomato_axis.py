#!/usr/bin/env python3
"""Build native tomato P-axis from GSE284240 and score GSE217631 on it."""
import sys, warnings, gzip
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from pathlib import Path
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

ROOT      = Path("C:/Users/moshe/Dropbox/ISF 2025")
SIG_DIR   = ROOT / "RNA_seq_axes" / "signatures"
SCORE_DIR = ROOT / "RNA_seq_axes" / "scores"
SIG_DIR.mkdir(exist_ok=True)
SCORE_DIR.mkdir(exist_ok=True)

print("=" * 60)
print("Building native tomato P-axis from GSE284240")
print("=" * 60)

# 1. Load count data
print("\n1. Loading GSE284240 count data...")
xl_path = ROOT / "RNA_seq_data/suppl_files/GSE284240/GSE284240_Raw_gene_count_annot.xlsx"
df_raw = pd.read_excel(xl_path, sheet_name=0, index_col=0)
if "Description" in df_raw.columns:
    df_raw = df_raw.drop(columns=["Description"])
print(f"   Shape: {df_raw.shape}")
print(f"   Columns: {list(df_raw.columns)}")
print(f"   First 3 gene IDs: {list(df_raw.index[:3])}")

df_raw.index = df_raw.index.str.replace(r"\.\d+$", "", regex=True)

wt_pf_cols = [c for c in df_raw.columns if c.startswith("WT_PF")]
wt_p_cols  = [c for c in df_raw.columns if c.startswith("WT_P") and "PF" not in c]
print(f"   WT_PF (LowPi) cols: {wt_pf_cols}")
print(f"   WT_P  (HighPi) cols: {wt_p_cols}")

all_cols = wt_p_cols + wt_pf_cols
counts = df_raw[all_cols].astype(int)

keep = counts.sum(axis=1) >= 10
counts = counts[keep]
print(f"   After low-count filter: {counts.shape[0]} genes")

# 2. DESeq2
print("\n2. Running DESeq2 (WT_PF vs WT_P)...")
mat = counts.T
metadata = pd.DataFrame({
    "condition": ["ctrl"] * len(wt_p_cols) + ["treat"] * len(wt_pf_cols)
}, index=all_cols)

keep2 = mat.std(axis=0) > 0
mat = mat.loc[:, keep2]
print(f"   After variance filter: {mat.shape[1]} genes")

dds = DeseqDataSet(counts=mat, metadata=metadata, design_factors="condition")
dds.deseq2()
stats_obj = DeseqStats(dds, contrast=["condition", "treat", "ctrl"])
stats_obj.summary()
res = stats_obj.results_df
lfc = res["log2FoldChange"].dropna()
padj = res["padj"].dropna()
print(f"   DESeq2 complete: {len(lfc)} genes with logFC")
print(f"   Significant (padj<0.05): {(padj < 0.05).sum()} genes")

marker_genes = {
    "Solyc11g069450": "SlPT4 (phosphate transporter, expect UP)",
    "Solyc10g078260": "SlSPX1 (expect UP)",
    "Solyc03g026040": "SlPHO1 (expect UP)",
    "Solyc09g010090": "SlPAP (expect UP)",
}
print("\n   Marker gene check:")
for gid, desc in marker_genes.items():
    if gid in lfc.index:
        print(f"   {gid}: logFC={lfc[gid]:.3f}  -- {desc}")
    else:
        print(f"   {gid}: NOT FOUND -- {desc}")

sig_out = SIG_DIR / "P_axis_tomato_logFC.csv"
lfc.to_frame("logFC").rename_axis("GeneID").to_csv(sig_out)
print(f"\n   Saved: {sig_out}")

# 3. Calibration
print("\n3. Computing tomato P-axis calibration...")

def score_sample(expr_log2, lfc_sig):
    """Dot product over shared genes."""
    common = lfc_sig.index.intersection(expr_log2.index)
    return (lfc_sig[common] * expr_log2[common]).sum()

log2_counts = np.log2(counts.astype(float) + 1)

ctrl_scores = [score_sample(log2_counts[c], lfc) for c in wt_p_cols]
trt_scores  = [score_sample(log2_counts[c], lfc) for c in wt_pf_cols]

ctrl_mean = np.mean(ctrl_scores)
trt_mean  = np.mean(trt_scores)
ctrl_list = ", ".join(f"{s:.1f}" for s in ctrl_scores)
trt_list  = ", ".join(f"{s:.1f}" for s in trt_scores)
print(f"   WT_P  (HighPi/ctrl) mean score: {ctrl_mean:.1f}  [{ctrl_list}]")
print(f"   WT_PF (LowPi/treat) mean score: {trt_mean:.1f}  [{trt_list}]")
print(f"   Range: {trt_mean - ctrl_mean:.1f} units")

# 4. Score GSE217631
print("\n4. Scoring GSE217631 on native tomato P-axis...")
fpkm_path = ROOT / "RNA_seq_data/suppl_files/GSE217631/GSE217631_fpkm_sample.xls.gz"
with gzip.open(fpkm_path, "rt") as fh:
    fpkm = pd.read_csv(fh, sep="\t", index_col=0)

fpkm.index = fpkm.index.str.replace(r"\.\d+$", "", regex=True)
print(f"   GSE217631 shape: {fpkm.shape}")

log2_fpkm = np.log2(fpkm.astype(float) + 1)

cols = list(fpkm.columns)
print(f"   Columns: {cols}")

leaf_ctrl = [c for c in cols if "LeafCtr" in c or ("leaf" in c.lower() and "ctr" in c.lower())]
leaf_trt  = [c for c in cols if "LeafTrt" in c or ("leaf" in c.lower() and "trt" in c.lower())]
root_ctrl = [c for c in cols if "RootCtr" in c or ("root" in c.lower() and "ctr" in c.lower())]
root_trt  = [c for c in cols if "RootTrt" in c or ("root" in c.lower() and "trt" in c.lower())]

if not leaf_ctrl:
    leaf_ctrl = cols[0:3]; leaf_trt = cols[3:6]
    root_ctrl = cols[6:9]; root_trt = cols[9:12]
    print("   Used positional column assignment")

groups = {"LeafCtrl": leaf_ctrl, "LeafTrt": leaf_trt,
          "RootCtrl": root_ctrl, "RootTrt":  root_trt}

common_genes = lfc.index.intersection(log2_fpkm.index)
print(f"   Genes in tomato P-axis: {len(lfc)}")
print(f"   Genes in GSE217631 FPKM: {len(log2_fpkm)}")
pct_overlap = len(common_genes) / len(lfc) * 100
print(f"   Overlap (no orthology needed!): {len(common_genes)} ({pct_overlap:.1f}% of axis)")

results = []
for grp, sample_cols in groups.items():
    avail = [c for c in sample_cols if c in log2_fpkm.columns]
    grp_scores = [score_sample(log2_fpkm[c], lfc) for c in avail]
    grp_mean = np.mean(grp_scores)
    grp_pct  = (grp_mean - ctrl_mean) / (trt_mean - ctrl_mean) * 100
    results.append({"group": grp, "raw_score": grp_mean, "pct_of_range": grp_pct})
    print(f"   {grp}: raw={grp_mean:.1f}  pct={grp_pct:.2f}%")

results_df = pd.DataFrame(results)

# 5. Within-experiment deltas
print("\n5. Treatment effects (PSI-475 vs Control):")
for tissue in ["Leaf", "Root"]:
    mask_c = results_df["group"] == (tissue + "Ctrl")
    mask_t = results_df["group"] == (tissue + "Trt")
    ctrl_pct = results_df.loc[mask_c, "pct_of_range"].values[0]
    trt_pct  = results_df.loc[mask_t, "pct_of_range"].values[0]
    delta    = trt_pct - ctrl_pct
    print(f"   {tissue}: Ctrl={ctrl_pct:.2f}%  Trt={trt_pct:.2f}%  Delta P-axis={delta:+.2f}%")

print("\n   Context:")
print("   Arabidopsis amino acids (GSE297649): Delta P-axis = +2.4% (Arabidopsis calibration)")
print("   Tomato PSI-475 leaf (native tomato calibration): see above")
print("   Tomato PSI-475 root (native tomato calibration): see above")

# 6. Save
out = SCORE_DIR / "GSE217631_native_tomato_paxis_scores.csv"
results_df.to_csv(out, index=False)
print(f"\n6. Saved: {out}")

calib_row = pd.DataFrame([{
    "axis": "P_axis_tomato",
    "ctrl_score": ctrl_mean, "treat_score": trt_mean,
    "ctrl_label": "HighPi_WT", "treat_label": "LowPi_WT",
    "source": "GSE284240"
}])
calib_path = ROOT / "RNA_seq_axes" / "calibration_tomato.csv"
calib_row.to_csv(calib_path, index=False)
print(f"   Saved tomato calibration: {calib_path}")

print("\n=== DONE ===")
