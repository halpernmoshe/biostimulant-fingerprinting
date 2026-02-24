#!/usr/bin/env python3
"""build_K_axis.py - Build K-axis from GSE180004 and 10-axis heatmap."""
import sys, warnings, gzip, io, tarfile, xlrd
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
from io import StringIO

ROOT = Path("C:/Users/moshe/Dropbox/ISF 2025")
SUPPL = ROOT / "RNA_seq_data" / "suppl_files"
AXES_DIR = ROOT / "RNA_seq_axes"
SIG_DIR = AXES_DIR / "signatures"
SCORE_DIR = AXES_DIR / "scores"
PLOT_DIR = ROOT / "state_space_figures"
for d in [SIG_DIR, SCORE_DIR, PLOT_DIR]:
    d.mkdir(parents=True, exist_ok=True)

LOG = []
def log(msg):
    ts = datetime.now().strftime("%H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line)
    LOG.append(line)

def score_on_axis(expr_df, lfc_sig):
    common = lfc_sig.index.intersection(expr_df.index)
    if len(common) < 100:
        log(f"  WARNING: only {len(common)} genes")
        return pd.Series(dtype=float)
    sig = lfc_sig[common]
    expr = expr_df.loc[common].astype(float)
    return expr.multiply(sig, axis=0).sum(axis=0)

def score_counts_on_axis(counts_df, lfc_sig):
    return score_on_axis(np.log2(counts_df.astype(float) + 1), lfc_sig)

def compute_delta_pct(treat_scores, ctrl_scores, calib_range):
    delta = float(np.mean(treat_scores)) - float(np.mean(ctrl_scores))
    pct = delta / calib_range * 100
    return delta, pct

# ============================================================
log("=" * 60)
log("STEP 1: Build K-axis from GSE180004 (4n vs 2n root)")
log("=" * 60)

K_SIG_FILE = SIG_DIR / "K_axis_logFC.csv"
EXTRACTED_DIR = SUPPL / "GSE180004_extracted"

n_extracted = len(list(EXTRACTED_DIR.glob("*.gz"))) if EXTRACTED_DIR.exists() else 0
if n_extracted < 12:
    EXTRACTED_DIR.mkdir(exist_ok=True)
    tar_path = SUPPL / "GSE180004" / "GSE180004_RAW.tar"
    log(f"  Extracting {tar_path.name} ...")
    with tarfile.open(tar_path) as tar:
        tar.extractall(EXTRACTED_DIR)
    log(f"  Extracted OK")
else:
    log(f"  Already extracted: {n_extracted} files")

def load_fpkm_counts(directory):
    files = sorted(Path(directory).glob("*.gene.FPKM.txt.gz"))
    dfs = []
    for f in files:
        with gzip.open(f, "rt") as fh:
            df_f = pd.read_csv(fh, sep="	")
        gene_col = df_f.columns[0]
        cnt = df_f["expected_count"].round().astype(int)
        parts = f.name.replace(".gene.FPKM.txt.gz", "").split("_", 1)
        short = parts[1] if len(parts) > 1 else parts[0]
        dfs.append(pd.Series(cnt.values, index=df_f[gene_col].values, name=short))
    df = pd.concat(dfs, axis=1).fillna(0).astype(int)
    log(f"  Loaded {df.shape[0]} genes x {df.shape[1]} samples")
    return df

all_counts = load_fpkm_counts(EXTRACTED_DIR)

if K_SIG_FILE.exists():
    K_sig = pd.read_csv(K_SIG_FILE, index_col=0)["logFC"]
    log(f"  Loaded existing K-axis: {len(K_sig)} genes")
else:
    log("  Computing K-axis with pydeseq2...")
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    ctrl_cols = [c for c in all_counts.columns if "2n_Control_root" in c]
    treat_cols = [c for c in all_counts.columns if "4n_Control_root" in c]
    log(f"  ctrl (2n): {ctrl_cols}")
    log(f"  treat (4n): {treat_cols}")
    all_cols = ctrl_cols + treat_cols
    mat = all_counts[all_cols].T
    metadata = pd.DataFrame(
        {"condition": ["ctrl"] * len(ctrl_cols) + ["treat"] * len(treat_cols)},
        index=all_cols)
    keep = (mat > 0).sum(axis=0) >= 3
    mat = mat.loc[:, keep]
    log(f"  After filtering: {mat.shape[1]} genes")
    dds = DeseqDataSet(counts=mat, metadata=metadata, design_factors="condition")
    dds.deseq2()
    stats = DeseqStats(dds, contrast=["condition", "treat", "ctrl"])
    stats.summary()
    res = stats.results_df
    K_sig = res["log2FoldChange"].dropna()
    K_sig.index.name = "gene_id"
    K_sig.to_frame(name="logFC").to_csv(K_SIG_FILE)
    log(f"  Saved K-axis: {len(K_sig)} genes")

hak5 = "AT4G13420"
if hak5 in K_sig.index:
    log(f"  HAK5 logFC = {K_sig[hak5]:.3f}")
else:
    log("  HAK5 NOT in K-axis")

# ============================================================
log("=" * 60)
log("STEP 2: Calibrate K-axis")
log("=" * 60)

calib_rows = []
for col in all_counts.columns:
    counts_col = all_counts[[col]]
    sc = score_counts_on_axis(counts_col, K_sig)
    if len(sc) == 0: continue
    score_val = float(sc.iloc[0])
    if "2n_Control_root" in col: group = "2n_ctrl"
    elif "4n_Control_root" in col: group = "4n_lowK"
    else: continue
    calib_rows.append({"sample": col, "group": group, "raw_score": score_val})

calib_df = pd.DataFrame(calib_rows)
ctrl_mean_K = calib_df[calib_df.group == "2n_ctrl"]["raw_score"].mean()
treat_mean_K = calib_df[calib_df.group == "4n_lowK"]["raw_score"].mean()
K_RANGE = treat_mean_K - ctrl_mean_K
log(f"  2n_ctrl mean: {ctrl_mean_K:.1f}")
log(f"  4n_lowK mean: {treat_mean_K:.1f}")
log(f"  K_RANGE: {K_RANGE:.1f}")
calib_df["pct_calibration"] = (calib_df["raw_score"] - ctrl_mean_K) / K_RANGE * 100
calib_df["axis"] = "K"
calib_df.to_csv(AXES_DIR / "calibration_K.csv", index=False)
log("  Saved calibration_K.csv")


# ============================================================
log("=" * 60)
log("STEP 3: Score 9 biostimulant treatments on K-axis")
log("=" * 60)

results = []

def append_result(dataset, tgroup, bclass, ctrl_sc, treat_sc):
    delta, pct = compute_delta_pct(treat_sc.values, ctrl_sc.values, K_RANGE)
    results.append({"dataset": dataset, "treatment_group": tgroup,
                    "biostimulant_class": bclass,
                    "ctrl_K_mean": float(ctrl_sc.mean()),
                    "treat_K_mean": float(treat_sc.mean()),
                    "delta_K": delta, "delta_K_pct": pct})
    log(f"  {tgroup}: delta_K={delta:.1f}, pct={pct:.1f}%")


# ============================================================
log("=" * 60)
log("STEP 3: Score 9 biostimulant treatments")
log("=" * 60)
results = []

def append_result(dataset, tgroup, bclass, ctrl_sc, treat_sc):
    delta, pct = compute_delta_pct(treat_sc.values, ctrl_sc.values, K_RANGE)
    results.append({"dataset": dataset, "treatment_group": tgroup,
                    "biostimulant_class": bclass,
                    "ctrl_K_mean": float(ctrl_sc.mean()),
                    "treat_K_mean": float(treat_sc.mean()),
                    "delta_K": delta, "delta_K_pct": pct})
    log(f"  {tgroup}: {pct:.1f}%")

# 1. GSE297649 (HTSeq counts)
log("  Loading GSE297649...")
gse_dir = SUPPL / "GSE297649"
counts_ht = {}
for f_297 in sorted(gse_dir.glob("*.txt.gz")):
    with gzip.open(f_297, "rt") as fh:
        lines_h = [l.strip().split("	") for l in fh if not l.startswith("__")]
    genes_h = [l[0] for l in lines_h if len(l) == 2]
    cnts_h = [int(l[1]) for l in lines_h if len(l) == 2]
    counts_ht[f_297.stem] = pd.Series(cnts_h, index=genes_h)
df_297649 = pd.DataFrame(counts_ht)
aa_c = [c for c in df_297649.columns if "WTaa" in c]
hs_c = [c for c in df_297649.columns if "WThumic" in c]
ct_c = [c for c in df_297649.columns if "WTctrl" in c]
log(f"  aa={len(aa_c)}, humic={len(hs_c)}, ctrl={len(ct_c)}")
if aa_c and ct_c:
    append_result("GSE297649", "amino_acid", "Amino acid biostimulant",
                  score_counts_on_axis(df_297649[ct_c], K_sig),
                  score_counts_on_axis(df_297649[aa_c], K_sig))
if hs_c and ct_c:
    append_result("GSE297649", "humic_subst", "Humic substance biostimulant",
                  score_counts_on_axis(df_297649[ct_c], K_sig),
                  score_counts_on_axis(df_297649[hs_c], K_sig))

# 2. GSE138478 (featureCounts)
log("  Loading GSE138478...")
gse_dir = SUPPL / "GSE138478"
def load_fc(fpath):
    with gzip.open(fpath, "rt") as fhx:
        lines_fc = fhx.readlines()
    data_fc = [l for l in lines_fc if not l.startswith("#")]
    df_fc = pd.read_csv(StringIO("".join(data_fc)), sep="	", index_col=0)
    drop_fc = ["Chr", "Start", "End", "Strand", "Length"]
    return df_fc.drop(columns=[c for c in drop_fc if c in df_fc.columns])
gmv_fc = load_fc(gse_dir / "GSE138478_HMZ008_featureCounts_output.txt.gz")
gmv_c = [c for c in gmv_fc.columns if "_C_" in c]
gmv_t = [c for c in gmv_fc.columns if "_T_" in c]
if gmv_c and gmv_t:
    append_result("GSE138478", "diacetyl_GMV", "PGPR volatile (GMV diacetyl)",
                  score_counts_on_axis(gmv_fc[gmv_c], K_sig),
                  score_counts_on_axis(gmv_fc[gmv_t], K_sig))
dia_fc = load_fc(gse_dir / "GSE138478_HMZ014_Dione_featureCounts_output.txt.gz")
dia_c = [c for c in dia_fc.columns if "CK" in c]
dia_t = [c for c in dia_fc.columns if "Dione" in c]
if dia_c and dia_t:
    append_result("GSE138478", "diacetyl_treated", "PGPR volatile (pure diacetyl)",
                  score_counts_on_axis(dia_fc[dia_c], K_sig),
                  score_counts_on_axis(dia_fc[dia_t], K_sig))

# 3. GSE208223 (log2 normalized, already log2 scale)
log("  Loading GSE208223...")
with gzip.open(SUPPL / "GSE208223" / "GSE208223_Root_norm.txt.gz", "rt") as fh:
    df_208 = pd.read_csv(fh, sep="	", index_col=0)
df_208.index = df_208.index.str.replace(r"\.\d+$", "", regex=True)
ctrl_sc_208 = score_on_axis(df_208[["A1","A2","A3"]], K_sig)
treat_sc_208 = score_on_axis(df_208[["B1","B2","B3"]], K_sig)
if len(ctrl_sc_208) > 0:
    append_result("GSE208223", "TiO2_treated", "Ti nanoparticles", ctrl_sc_208, treat_sc_208)

# 4. GSE298982 (log2 TMM normalized)
log("  Loading GSE298982...")
with gzip.open(SUPPL / "GSE298982" / "GSE298982_Glu_RNAseq_TMM_normalization.txt.gz", "rt") as fh:
    lines298 = fh.readlines()
header298 = lines298[0].lstrip("#").strip().split("	")
rows298 = [l.strip().split("	") for l in lines298[1:] if l.strip()]
df_298 = pd.DataFrame(rows298, columns=header298).set_index("geneID")
ctrl_298 = ["5NN_1", "5NN_2", "5NN_3"]
treat_298 = ["5Glu_1", "5Glu_2", "5Glu_3"]
df_298_vals = df_298[ctrl_298 + treat_298].astype(float)
ctrl_sc_298 = score_on_axis(df_298_vals[ctrl_298], K_sig)
treat_sc_298 = score_on_axis(df_298_vals[treat_298], K_sig)
if len(ctrl_sc_298) > 0:
    append_result("GSE298982", "GSE298982_glutamate", "Amino acid source (glutamate)",
                  ctrl_sc_298, treat_sc_298)

# 5. GSE139154 (raw counts)
log("  Loading GSE139154...")
with gzip.open(SUPPL / "GSE139154" / "GSE139154_all_sample_cnt.tsv.gz", "rt") as fh:
    header139 = fh.readline().strip().split("	")
    df_139 = pd.read_csv(fh, sep="	", index_col=0, header=None)
df_139.columns = header139
ctrl_139 = [c for c in header139 if "Col-0_C_" in c]
treat_139 = [c for c in header139 if "Col-0_T_" in c]
if ctrl_139 and treat_139:
    append_result("GSE139154", "GSE139154_GB03", "PGPR volatile (Bacillus GB03)",
                  score_counts_on_axis(df_139[ctrl_139].astype(int), K_sig),
                  score_counts_on_axis(df_139[treat_139].astype(int), K_sig))

# 6. GSE87337 (raw counts, NoBact vs SynCom, Col-0, LowP)
log("  Loading GSE87337...")
with gzip.open(SUPPL / "GSE87337" / "GSE87337_PBI_preTreat_counts.txt.gz", "rt") as fh:
    desc_names_87 = fh.readline().strip().split("	")
    desc_codes_87 = fh.readline().strip().split("	")[1:]
    df_87337 = pd.read_csv(fh, sep="	", index_col=0, header=None)
df_87337.columns = desc_codes_87
noBact_c = [code for code,name in zip(desc_codes_87,desc_names_87)
             if "noBact" in name and "Col-0" in name
             and "LowP" in name and "preTreat(-P)" in name]
synCom_t = [code for code,name in zip(desc_codes_87,desc_names_87)
             if "SynCom" in name and "Col-0" in name
             and "LowP" in name and "preTreat(-P)" in name]
log(f"  noBact={noBact_c}")
log(f"  synCom={synCom_t}")
if noBact_c and synCom_t:
    append_result("GSE87337", "GSE87337_SynCom", "SynCom bacterial community",
                  score_counts_on_axis(df_87337[noBact_c].astype(int), K_sig),
                  score_counts_on_axis(df_87337[synCom_t].astype(int), K_sig))

# 7. GSE199501 (CPM from old Excel xls.gz)
log("  Loading GSE199501...")
with gzip.open(SUPPL / "GSE199501" / "GSE199501_Processed_data_DEG_total_CPM.xls.gz", "rb") as fh:
    xls_data = fh.read()
wb_xls = xlrd.open_workbook(file_contents=xls_data)
sh_xls = wb_xls.sheets()[0]
hdr_xls = [sh_xls.cell_value(0, j) for j in range(sh_xls.ncols)]
rows_xls = [[sh_xls.cell_value(i, j) for j in range(sh_xls.ncols)] for i in range(1, sh_xls.nrows)]
df_199 = pd.DataFrame(rows_xls, columns=hdr_xls).set_index("Gene ID ")
ctrl_cpm = ["C-1953 CPM", "C-1954 CPM", "C-1955 CPM"]
treat_cpm = ["T-1956 CPM", "T-1957 CPM", "T-1958 CPM"]
df_199v = df_199[ctrl_cpm + treat_cpm].astype(float)
ctrl_sc_199 = score_on_axis(np.log2(df_199v[ctrl_cpm] + 1), K_sig)
treat_sc_199 = score_on_axis(np.log2(df_199v[treat_cpm] + 1), K_sig)
if len(ctrl_sc_199) > 0:
    append_result("GSE199501", "GSE199501_PGPR", "PGPR (Priestia megaterium)",
                  ctrl_sc_199, treat_sc_199)

# Save K-axis results
df_res = pd.DataFrame(results)
df_res.to_csv(SCORE_DIR / "K_axis_biostimulant_scores.csv", index=False)
log(f"  Saved K_axis_biostimulant_scores.csv ({len(df_res)} rows)")
for _, r in df_res.iterrows():
    log(f"  {r.treatment_group:35s}  {r.delta_K_pct:.1f}%")

# ============================================================
log("=" * 60)
log("STEP 4: Build comprehensive 10-axis heatmap")
log("=" * 60)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

# Load calibration ranges for all axes
calib_main = pd.read_csv(AXES_DIR / "calibration.csv")
calib_ranges = {row["axis"]: row["treat_score"] - row["ctrl_score"]
                for _, row in calib_main.iterrows()}
log(f"  Axes with calibration: {list(calib_ranges.keys())}")

# Load 6-axis summary (5 treatments x 10 axes, already in % of calib range)
df6 = pd.read_csv(SCORE_DIR / "biostimulant_6axis_summary.csv")
matrix_6 = df6.pivot(index="treatment", columns="axis", values="delta_pct")
log(f"  6-axis matrix: {matrix_6.shape}")

# Load 9-treatment effect vectors
ev = pd.read_csv(SCORE_DIR / "effect_vectors_complete.csv")

# Map from effect_vectors treatment names to 6-axis treatment names
map_to_6 = {
    "amino_acid": "amino_acid",
    "humic_subst": "humic_subst",
    "diacetyl_GMV": "GMV_treated",
    "diacetyl_treated": "diacetyl_treated",
    "TiO2_treated": "TiO2_treated",
    "GSE298982_glutamate": None,
    "GSE139154_GB03": None,
    "GSE87337_SynCom": None,
    "GSE199501_PGPR": None,
}

display_names = {
    "amino_acid": "Amino acid (GSE297649)",
    "humic_subst": "Humic substances (GSE297649)",
    "diacetyl_GMV": "PGPR mix VOC (GMV, GSE138478)",
    "diacetyl_treated": "Pure diacetyl (GSE138478)",
    "TiO2_treated": "TiO2 nanoparticles (GSE208223)",
    "GSE298982_glutamate": "Glutamate N-source (GSE298982)",
    "GSE139154_GB03": "Bacillus GB03 VOC (GSE139154)",
    "GSE87337_SynCom": "SynCom community (GSE87337)",
    "GSE199501_PGPR": "P. megaterium (GSE199501)",
}

axes_ordered = ["P_axis", "N_axis", "ABA_axis", "JA_axis", "SA_axis",
                "Fe_axis", "CK_axis", "ET_axis", "Auxin_axis", "K_axis"]
axis_display = {
    "P_axis": "P-starvation", "N_axis": "N-supply", "ABA_axis": "ABA",
    "JA_axis": "JA", "SA_axis": "SA (BTH)", "Fe_axis": "Fe-deficiency",
    "CK_axis": "Cytokinin", "ET_axis": "Ethylene",
    "Auxin_axis": "Auxin (IAA)", "K_axis": "K-demand (4n)"
}

all_treatments = list(ev["treatment_group"])
matrix_data = {}
for tname in all_treatments:
    rd = {}
    ev_row = ev[ev["treatment_group"] == tname].iloc[0]
    rd["P_axis"]   = ev_row["delta_P"] / calib_ranges["P_axis"] * 100
    rd["N_axis"]   = ev_row["delta_N"] / calib_ranges["N_axis"] * 100
    rd["ABA_axis"] = ev_row["delta_ABA"] / calib_ranges["ABA_axis"] * 100
    six_name = map_to_6.get(tname)
    if six_name and six_name in matrix_6.index:
        row6 = matrix_6.loc[six_name]
        for ax in ["JA_axis", "SA_axis", "Fe_axis", "CK_axis", "ET_axis", "Auxin_axis"]:
            rd[ax] = float(row6.get(ax, float("nan")))
    else:
        for ax in ["JA_axis", "SA_axis", "Fe_axis", "CK_axis", "ET_axis", "Auxin_axis"]:
            rd[ax] = float("nan")
    k_row_df = df_res[df_res["treatment_group"] == tname]
    rd["K_axis"] = float(k_row_df["delta_K_pct"].iloc[0]) if len(k_row_df) > 0 else float("nan")
    matrix_data[tname] = rd

full_matrix = pd.DataFrame(matrix_data).T[axes_ordered]
full_matrix.index = [display_names.get(t, t) for t in full_matrix.index]
full_matrix.columns = [axis_display.get(a, a) for a in full_matrix.columns]
log(f"  Full matrix: {full_matrix.shape}")
full_matrix.to_csv(SCORE_DIR / "comprehensive_10axis_matrix.csv")
log("  Saved comprehensive_10axis_matrix.csv")

# Build heatmap
plot_m = full_matrix.clip(-50, 80).fillna(0)
annot_m = full_matrix.round(0)

fig, ax = plt.subplots(figsize=(14, 8))
sns.heatmap(
    plot_m, center=0, cmap="RdBu_r",
    annot=annot_m, fmt=".0f",
    linewidths=0.5, linecolor="lightgray",
    ax=ax, vmin=-50, vmax=80,
    cbar_kws={"label": "% of calibration range", "shrink": 0.8}
)
ax.set_title(
    "Biostimulant transcriptomic fingerprints (10 physiological axes)" + chr(10) +
    "(% of calibration range; NaN=white/0; clamped [-50, +80])",
    fontsize=12, pad=12)
ax.set_xlabel("Physiological axis", fontsize=11)
ax.set_ylabel("Biostimulant treatment", fontsize=11)
plt.xticks(rotation=40, ha="right", fontsize=10)
plt.yticks(rotation=0, fontsize=9)
plt.tight_layout()
out_fig = PLOT_DIR / "comprehensive_heatmap_10axes.png"
plt.savefig(out_fig, dpi=150, bbox_inches="tight")
plt.close()
log(f"  Saved heatmap: {out_fig}")

# ============================================================
# ============================================================
# Report
# ============================================================
# ============================================================
# Write report
# ============================================================
hak5_val = float(K_sig.get("AT4G13420", float("nan")))
report_path = ROOT / "agent_coordination" / "K_axis_report.md"
ENDL = chr(10)
report_lines = [
    "# K-Axis Build Report" + ENDL + ENDL,
    f"Generated: {datetime.now().strftime(chr(37)+chr(89)+chr(45)+chr(109)+chr(45)+chr(100))}" + ENDL + ENDL,
    "## Experimental design" + ENDL + ENDL,
    "K-axis from GSE180004: 4n vs 2n Arabidopsis root, control (Fischer et al. 2022)" + ENDL + ENDL,
    f"2n_ctrl mean: {ctrl_mean_K:.1f}, 4n_lowK mean: {treat_mean_K:.1f}, K_RANGE: {K_RANGE:.1f}" + ENDL + ENDL,
    f"HAK5 (AT4G13420) logFC = {hak5_val:.3f}" + ENDL + ENDL,
    "## K-axis scores" + ENDL + ENDL,
    "| Treatment | delta_K | % calib |" + ENDL,
    "|-----------|---------|---------|" + ENDL,
]
for _, r2 in df_res.iterrows():
    report_lines.append(f"| {r2.treatment_group} | {r2.delta_K:.1f} | {r2.delta_K_pct:.1f}% |" + ENDL)
report_lines.append(ENDL + "## Log" + ENDL + ENDL)
for logline in LOG:
    report_lines.append(logline + ENDL)
with open(report_path, "w", encoding="utf-8") as rf:
    rf.writelines(report_lines)
log(f"  Report: {report_path}")

log("=" * 60)
log("ALL DONE")
log("=" * 60)
