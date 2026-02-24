#!/usr/bin/env python3
import sys, warnings, gzip, xlrd
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

def log(msg):
    print(f"[{datetime.now().strftime(chr(37)+chr(72)+chr(58)+chr(37)+chr(77)+chr(58)+chr(37)+chr(83))}] {msg}")

K_sig = pd.read_csv(SIG_DIR / "K_axis_logFC.csv", index_col=0)["logFC"]
log(f"Loaded K-axis: {len(K_sig)} genes")

calib_K = pd.read_csv(AXES_DIR / "calibration_K.csv")
ctrl_K = calib_K[calib_K.group=="2n_ctrl"]["raw_score"].mean()
treat_K = calib_K[calib_K.group=="4n_lowK"]["raw_score"].mean()
K_RANGE = treat_K - ctrl_K
log(f"K_RANGE = {K_RANGE:.1f}")

def score_log2counts(df, sig):
    common = sig.index.intersection(df.index)
    return np.log2(df.loc[common].astype(float)+1).multiply(sig[common], axis=0).sum(axis=0)

def score_log2expr(df, sig):
    common = sig.index.intersection(df.index)
    return df.loc[common].astype(float).multiply(sig[common], axis=0).sum(axis=0)

def delta_pct(treat_sc, ctrl_sc):
    d = float(treat_sc.mean()) - float(ctrl_sc.mean())
    return d, d / K_RANGE * 100

results = []


# 1. GSE297649: WTcontrol=ctrl, WTaa=amino acid, WTsh=humic substance
log("Loading GSE297649...")
gse_dir = SUPPL / "GSE297649"
cnts = {}
for f in sorted(gse_dir.glob("*.txt.gz")):
    with gzip.open(f, "rt") as fh:
        lns = [l.strip().split("	") for l in fh if not l.startswith("__")]
    g = [l[0] for l in lns if len(l)==2]
    c = [int(l[1]) for l in lns if len(l)==2]
    cnts[f.stem] = pd.Series(c, index=g)
df297 = pd.DataFrame(cnts)
aa_c = [c for c in df297.columns if "WTaa" in c]
sh_c = [c for c in df297.columns if "WTsh" in c]
ct_c = [c for c in df297.columns if "WTcontrol" in c]
log(f"  aa={len(aa_c)}, humic={len(sh_c)}, ctrl={len(ct_c)}")
if aa_c and ct_c:
    d,p = delta_pct(score_log2counts(df297[aa_c],K_sig),score_log2counts(df297[ct_c],K_sig))
    results.append({"dataset":"GSE297649","treatment_group":"amino_acid",
                    "biostimulant_class":"Amino acid biostimulant","delta_K":d,"delta_K_pct":p})
    log(f"  amino_acid: {p:.1f}%")
if sh_c and ct_c:
    d,p = delta_pct(score_log2counts(df297[sh_c],K_sig),score_log2counts(df297[ct_c],K_sig))
    results.append({"dataset":"GSE297649","treatment_group":"humic_subst",
                    "biostimulant_class":"Humic substance biostimulant","delta_K":d,"delta_K_pct":p})
    log(f"  humic_subst: {p:.1f}%")

# 2. GSE138478: featureCounts, GMV experiment (C_5DAT=ctrl, T_5DAT=treat)
log("Loading GSE138478...")
gse_dir = SUPPL / "GSE138478"
def load_fc(fpath):
    with gzip.open(fpath, "rt") as fh:
        lns = fh.readlines()
    data = [l for l in lns if not l.startswith("#")]
    df = pd.read_csv(StringIO("".join(data)), sep="	", index_col=0)
    drop = ["Chr","Start","End","Strand","Length"]
    return df.drop(columns=[c for c in drop if c in df.columns])
gmv = load_fc(gse_dir / "GSE138478_HMZ008_featureCounts_output.txt.gz")
gmv_c = [c for c in gmv.columns if "C_5DAT" in c]
gmv_t = [c for c in gmv.columns if "T_5DAT" in c]
log(f"  GMV ctrl={len(gmv_c)}, treat={len(gmv_t)}")
if gmv_c and gmv_t:
    d,p = delta_pct(score_log2counts(gmv[gmv_t],K_sig),score_log2counts(gmv[gmv_c],K_sig))
    results.append({"dataset":"GSE138478","treatment_group":"diacetyl_GMV",
                    "biostimulant_class":"PGPR volatile (GMV diacetyl)","delta_K":d,"delta_K_pct":p})
    log(f"  diacetyl_GMV: {p:.1f}%")
dia = load_fc(gse_dir / "GSE138478_HMZ014_Dione_featureCounts_output.txt.gz")
dia_c = [c for c in dia.columns if "CK" in c]
dia_t = [c for c in dia.columns if "Dione" in c]
if dia_c and dia_t:
    d,p = delta_pct(score_log2counts(dia[dia_t],K_sig),score_log2counts(dia[dia_c],K_sig))
    results.append({"dataset":"GSE138478","treatment_group":"diacetyl_treated",
                    "biostimulant_class":"PGPR volatile (pure diacetyl)","delta_K":d,"delta_K_pct":p})
    log(f"  diacetyl_treated: {p:.1f}%")

# 3. GSE208223: log2-normalized, A1-3=ctrl, B1-3=TiO2
log("Loading GSE208223...")
with gzip.open(SUPPL/"GSE208223"/"GSE208223_Root_norm.txt.gz","rt") as fh:
    df208 = pd.read_csv(fh, sep="	", index_col=0)
df208.index = df208.index.str.replace(r"\.\d+$","",regex=True)
d,p = delta_pct(score_log2expr(df208[["B1","B2","B3"]],K_sig),
                score_log2expr(df208[["A1","A2","A3"]],K_sig))
results.append({"dataset":"GSE208223","treatment_group":"TiO2_treated",
                "biostimulant_class":"Ti nanoparticles","delta_K":d,"delta_K_pct":p})
log(f"  TiO2_treated: {p:.1f}%")

# 4. GSE298982: raw TMM (linear), apply log2(TMM+1)
log("Loading GSE298982...")
with gzip.open(SUPPL/"GSE298982"/"GSE298982_Glu_RNAseq_TMM_normalization.txt.gz","rt") as fh:
    lns298 = fh.readlines()
hdr298 = lns298[0].lstrip("#").strip().split("	")
df298 = pd.DataFrame([l.strip().split("	") for l in lns298[1:] if l.strip()],
                     columns=hdr298).set_index("geneID")
ctrl298 = ["5NN_1","5NN_2","5NN_3"]
treat298 = ["5Glu_1","5Glu_2","5Glu_3"]
# Apply log2(raw_TMM+1) since values are linear TMM normalized counts
d,p = delta_pct(score_log2counts(df298[treat298].astype(float),K_sig),
                score_log2counts(df298[ctrl298].astype(float),K_sig))
results.append({"dataset":"GSE298982","treatment_group":"GSE298982_glutamate",
                "biostimulant_class":"Amino acid source (glutamate)","delta_K":d,"delta_K_pct":p})
log(f"  GSE298982_glutamate: {p:.1f}%")

# 5. GSE139154: raw counts
log("Loading GSE139154...")
with gzip.open(SUPPL/"GSE139154"/"GSE139154_all_sample_cnt.tsv.gz","rt") as fh:
    hdr139 = fh.readline().strip().split("	")
    df139 = pd.read_csv(fh, sep="	", index_col=0, header=None)
df139.columns = hdr139
ctrl139 = [c for c in hdr139 if "Col-0_C_" in c]
treat139 = [c for c in hdr139 if "Col-0_T_" in c]
d,p = delta_pct(score_log2counts(df139[treat139].astype(int),K_sig),
                score_log2counts(df139[ctrl139].astype(int),K_sig))
results.append({"dataset":"GSE139154","treatment_group":"GSE139154_GB03",
                "biostimulant_class":"PGPR volatile (Bacillus GB03)","delta_K":d,"delta_K_pct":p})
log(f"  GSE139154_GB03: {p:.1f}%")


# 6. GSE87337: NoBact vs SynCom, Col-0, LowP, preTreat with -P condition
log("Loading GSE87337...")
with gzip.open(SUPPL/"GSE87337"/"GSE87337_PBI_preTreat_counts.txt.gz","rt") as fh:
    desc87 = fh.readline().strip().split("	")
    codes87 = fh.readline().strip().split("	")[1:]
    df87 = pd.read_csv(fh, sep="	", index_col=0, header=None)
df87.columns = codes87
ptag = "preTreat(-P)"  # the pretreatment condition tag
noBact = [code for code,nm in zip(codes87,desc87)
          if "noBact" in nm and "Col-0" in nm and "LowP" in nm and ptag in nm]
synCom = [code for code,nm in zip(codes87,desc87)
          if "SynCom" in nm and "Col-0" in nm and "LowP" in nm and ptag in nm]
log(f"  noBact={noBact}")
log(f"  synCom={synCom}")
d,p = delta_pct(score_log2counts(df87[synCom].astype(int),K_sig),
                score_log2counts(df87[noBact].astype(int),K_sig))
results.append({"dataset":"GSE87337","treatment_group":"GSE87337_SynCom",
                "biostimulant_class":"SynCom bacterial community","delta_K":d,"delta_K_pct":p})
log(f"  GSE87337_SynCom: {p:.1f}%")

# 7. GSE199501: CPM (linear), apply log2(CPM+1)
log("Loading GSE199501...")
with gzip.open(SUPPL/"GSE199501"/"GSE199501_Processed_data_DEG_total_CPM.xls.gz","rb") as fh:
    xls_data = fh.read()
wb = xlrd.open_workbook(file_contents=xls_data)
sh = wb.sheets()[0]
hdr = [sh.cell_value(0,j) for j in range(sh.ncols)]
rows_xls = [[sh.cell_value(i,j) for j in range(sh.ncols)] for i in range(1,sh.nrows)]
df199 = pd.DataFrame(rows_xls, columns=hdr).set_index("Gene ID ")
ctrl_cpm = ["C-1953 CPM","C-1954 CPM","C-1955 CPM"]
treat_cpm = ["T-1956 CPM","T-1957 CPM","T-1958 CPM"]
df199v = df199[ctrl_cpm+treat_cpm].astype(float)
d,p = delta_pct(score_log2counts(df199v[treat_cpm],K_sig),
                score_log2counts(df199v[ctrl_cpm],K_sig))
results.append({"dataset":"GSE199501","treatment_group":"GSE199501_PGPR",
                "biostimulant_class":"PGPR (Priestia megaterium)","delta_K":d,"delta_K_pct":p})
log(f"  GSE199501_PGPR: {p:.1f}%")

# Save all results
df_res = pd.DataFrame(results)
df_res.to_csv(SCORE_DIR/"K_axis_biostimulant_scores.csv", index=False)
log(f"Saved K_axis_biostimulant_scores.csv ({len(df_res)} rows)")
log("")
log("Summary:")
for _, r in df_res.iterrows():
    log(f"  {r.treatment_group:35s}  {r.delta_K_pct:.1f}%")
