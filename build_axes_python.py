#!/usr/bin/env python3
"""
build_axes_python.py  —  RNA-seq axis building + biostimulant scoring
======================================================================
Axes:
  P-axis:   GSE90061  (Arabidopsis +/-Pi root tips, DESeq2, HTSeq counts)
  N-axis:   GSE97500  (Arabidopsis N-supply time-series, Varala 2018, DESeq2)
  ABA-axis: GSE127910 (Arabidopsis + ABA vs DMSO, TPM fold-change)

Biostimulants scored:
  GSE297649  (Arabidopsis + amino acids + humic substances)

Outputs:
  RNA_seq_axes/signatures/  — logFC CSVs per axis
  RNA_seq_axes/scores/      — per-sample and summary scores
  all_axes_scores_rnaseq.csv
  state_space_figures/      — PNG plots
"""

import sys, warnings, gzip, urllib.request
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime

ROOT      = Path(__file__).parent
SUPPL     = ROOT / "RNA_seq_data" / "suppl_files"
MATRIX    = ROOT / "RNA_seq_data" / "series_matrices"
AXES_DIR  = ROOT / "RNA_seq_axes"
SIG_DIR   = AXES_DIR / "signatures"
SCORE_DIR = AXES_DIR / "scores"
PLOT_DIR  = ROOT / "state_space_figures"
for d in [AXES_DIR, SIG_DIR, SCORE_DIR, PLOT_DIR]:
    d.mkdir(parents=True, exist_ok=True)

LOG = []
def log(msg):
    ts   = datetime.now().strftime("%H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line)
    LOG.append(line)


# ─────────────────────────────────────────────────────────────
# Utilities
# ─────────────────────────────────────────────────────────────

def read_htseq_counts(directory, pattern="*.gz"):
    """Read all HTSeq count .gz files in a directory into a DataFrame."""
    files = sorted(Path(directory).glob(pattern))
    dfs   = []
    for f in files:
        try:
            with gzip.open(f, "rt") as fh:
                lines = [l.strip().split("\t") for l in fh if not l.startswith("__")]
        except Exception:
            continue
        genes  = [l[0] for l in lines if len(l) == 2]
        counts = [int(l[1]) for l in lines if len(l) == 2]
        if genes:
            dfs.append(pd.Series(counts, index=genes, name=f.stem.split(".")[0]))
    if not dfs:
        return pd.DataFrame()
    df = pd.concat(dfs, axis=1).fillna(0).astype(int)
    log(f"  Loaded {df.shape[0]} genes x {df.shape[1]} samples from {Path(directory).name}")
    return df


def parse_series_matrix(path):
    """Parse a GEO series matrix .gz file. Returns dict of metadata lists."""
    meta = {}
    with gzip.open(path, "rb") as f:
        raw = f.read()
    text = raw.decode("latin-1")
    for line in text.split("\n"):
        if line.startswith("!Sample_"):
            parts = line.strip().split("\t")
            key   = parts[0].lstrip("!")
            vals  = [v.strip('"') for v in parts[1:]]
            if key not in meta:
                meta[key] = vals
            else:
                # duplicate key = repeat characteristics
                meta.setdefault(key + "_extra", []).append(vals)
    return meta


def run_deseq2(counts_df, treat_cols, ctrl_cols, name="contrast"):
    """DESeq2 via pydeseq2. Returns logFC Series (treat vs ctrl)."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds  import DeseqStats

    all_cols = ctrl_cols + treat_cols
    mat      = counts_df[all_cols].T          # samples x genes
    metadata = pd.DataFrame({
        "condition": ["ctrl"] * len(ctrl_cols) + ["treat"] * len(treat_cols)
    }, index=all_cols)

    keep = mat.std(axis=0) > 0
    mat  = mat.loc[:, keep]
    log(f"  DESeq2 ({name}): {mat.shape[1]} genes, treat={len(treat_cols)} ctrl={len(ctrl_cols)}")

    dds = DeseqDataSet(counts=mat, metadata=metadata, design_factors="condition")
    dds.deseq2()
    stats = DeseqStats(dds, contrast=["condition", "treat", "ctrl"])
    stats.summary()
    res = stats.results_df
    lfc = res["log2FoldChange"].dropna()
    log(f"  Done: {len(lfc)} genes with logFC")
    return lfc


def score_on_axis(counts_df, lfc_sig):
    """Dot product: for each sample, sum(log2(count+1) * logFC)."""
    common = lfc_sig.index.intersection(counts_df.index)
    if len(common) < 100:
        log(f"  WARNING: only {len(common)} genes in common — skipping")
        return pd.Series(dtype=float)
    sig    = lfc_sig[common]
    expr   = np.log2(counts_df.loc[common].astype(float) + 1)
    return expr.multiply(sig, axis=0).sum(axis=0)


def score_on_axis_tpm(tpm_df, lfc_sig):
    """Dot product using TPM (log2-transformed) instead of counts."""
    common = lfc_sig.index.intersection(tpm_df.index)
    if len(common) < 100:
        log(f"  WARNING: only {len(common)} genes in common for TPM scoring — skipping")
        return pd.Series(dtype=float)
    sig  = lfc_sig[common]
    expr = np.log2(tpm_df.loc[common].astype(float) + 1)
    return expr.multiply(sig, axis=0).sum(axis=0)


def download_file(url, dest):
    """Download a file from URL to dest path (skip if exists)."""
    dest = Path(dest)
    if dest.exists() and dest.stat().st_size > 1000:
        log(f"  Already present: {dest.name}")
        return True
    log(f"  Downloading: {url}")
    try:
        urllib.request.urlretrieve(url, dest)
        log(f"  OK: {dest.name} ({dest.stat().st_size // 1024} KB)")
        return True
    except Exception as e:
        log(f"  FAIL: {e}")
        return False


# ══════════════════════════════════════════════════════════════
# STEP 1a: P-axis — load from file or recompute
# ══════════════════════════════════════════════════════════════
log("=" * 55)
log("STEP 1a: P-axis signature (GSE90061)")
log("=" * 55)

p_sig_file = SIG_DIR / "P_axis_logFC.csv"

if p_sig_file.exists():
    P_sig = pd.read_csv(p_sig_file, index_col=0)["logFC"]
    log(f"  Loaded existing P-axis signature: {len(P_sig)} genes")
else:
    log("  Computing P-axis from GSE90061 count files...")
    counts_P = read_htseq_counts(SUPPL / "GSE90061")
    hp_cols  = [c for c in counts_P.columns if "hp" in c]
    lp_cols  = [c for c in counts_P.columns if "lp" in c
                and "lfe" not in c and "lpM" not in c]
    log(f"  High-Pi: {hp_cols}")
    log(f"  Low-Pi:  {lp_cols}")
    if hp_cols and lp_cols:
        P_sig = run_deseq2(counts_P, lp_cols, hp_cols, "LowPi_vs_HighPi")
        P_sig.to_frame("logFC").rename_axis("GeneID").to_csv(p_sig_file)
        log(f"  Saved: {p_sig_file.name}")
    else:
        P_sig = pd.Series(dtype=float)
        log("  ERROR: could not identify hp/lp columns")


# ══════════════════════════════════════════════════════════════
# STEP 1b: N-axis — GSE97500 (Varala 2018, shoot t=120 N vs KCl)
# ══════════════════════════════════════════════════════════════
log("\n" + "=" * 55)
log("STEP 1b: N-axis signature (GSE97500)")
log("=" * 55)

n_sig_file = SIG_DIR / "N_axis_logFC.csv"
N_sig = pd.Series(dtype=float)

if n_sig_file.exists():
    N_sig = pd.read_csv(n_sig_file, index_col=0)["logFC"]
    log(f"  Loaded existing N-axis signature: {len(N_sig)} genes")
else:
    sm_path = MATRIX / "GSE97500_series_matrix.txt.gz"
    if not sm_path.exists():
        # Try to download it
        url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE97nnn/GSE97500/matrix/GSE97500_series_matrix.txt.gz"
        download_file(url, sm_path)

    if sm_path.exists():
        log("  Parsing GSE97500 series matrix...")
        meta = parse_series_matrix(sm_path)

        gsm_ids = meta.get("Sample_geo_accession", [])
        titles  = meta.get("Sample_title", [])
        log(f"  Found {len(gsm_ids)} samples")
        if titles:
            log(f"  Sample titles (first 12): {titles[:12]}")

        # Title format: S0-1 = Shoot t=0 rep1
        #               S5C-1 = Shoot t=5 KCl-control rep1
        #               S5N-1 = Shoot t=5 Nitrogen rep1
        #               R120N-3 = Root t=120 Nitrogen rep3
        # Strategy: use Shoot, t=120min, N vs C (3 reps each)
        gsm_to_title = dict(zip(gsm_ids, titles))

        shoot_N_120 = []   # High-N at 120 min
        shoot_C_120 = []   # KCl control at 120 min

        for gsm, title in gsm_to_title.items():
            t = title.strip()
            if t.startswith("S") and "120" in t:
                if t.startswith("S120N"):
                    shoot_N_120.append(gsm)
                elif t.startswith("S120C"):
                    shoot_C_120.append(gsm)

        log(f"  Shoot t=120 N-treated GSMs: {shoot_N_120}")
        log(f"  Shoot t=120 KCl-ctrl GSMs:  {shoot_C_120}")

        if len(shoot_N_120) >= 2 and len(shoot_C_120) >= 2:
            # Load count data
            counts_N = read_htseq_counts(SUPPL / "GSE97500")

            # Map GSM to column name (column names start with GSM prefix)
            def find_col(gsm, columns):
                for c in columns:
                    if c.startswith(gsm):
                        return c
                return None

            n_treat_cols = [find_col(g, counts_N.columns) for g in shoot_N_120]
            n_ctrl_cols  = [find_col(g, counts_N.columns) for g in shoot_C_120]
            n_treat_cols = [c for c in n_treat_cols if c]
            n_ctrl_cols  = [c for c in n_ctrl_cols  if c]

            log(f"  Matched N-treat columns: {n_treat_cols}")
            log(f"  Matched KCl-ctrl columns: {n_ctrl_cols}")

            if len(n_treat_cols) >= 2 and len(n_ctrl_cols) >= 2:
                N_sig = run_deseq2(
                    counts_N, n_treat_cols, n_ctrl_cols,
                    "NitrogenSupply_vs_KCl_t120_shoot"
                )
                N_sig.to_frame("logFC").rename_axis("GeneID").to_csv(n_sig_file)
                log(f"  N-axis signature saved: {n_sig_file.name}")
            else:
                log("  ERROR: could not match GSMs to count columns")
        else:
            log("  Could not find shoot t=120 samples in series matrix")
            log("  Falling back: using t=60 shoot samples...")
            shoot_N_60 = [g for g, t in gsm_to_title.items() if t.startswith("S60N")]
            shoot_C_60 = [g for g, t in gsm_to_title.items() if t.startswith("S60C")]
            log(f"  t=60 N GSMs: {shoot_N_60}, C GSMs: {shoot_C_60}")
            if len(shoot_N_60) >= 2 and len(shoot_C_60) >= 2:
                counts_N = read_htseq_counts(SUPPL / "GSE97500")
                n_treat_cols = [c for c in counts_N.columns
                                if any(c.startswith(g) for g in shoot_N_60)]
                n_ctrl_cols  = [c for c in counts_N.columns
                                if any(c.startswith(g) for g in shoot_C_60)]
                if len(n_treat_cols) >= 2 and len(n_ctrl_cols) >= 2:
                    N_sig = run_deseq2(
                        counts_N, n_treat_cols, n_ctrl_cols,
                        "NitrogenSupply_vs_KCl_t60_shoot"
                    )
                    N_sig.to_frame("logFC").rename_axis("GeneID").to_csv(n_sig_file)
                    log(f"  N-axis (t=60) saved: {n_sig_file.name}")
    else:
        log("  GSE97500 series matrix not available — N-axis skipped")


# ══════════════════════════════════════════════════════════════
# STEP 1c: ABA-axis — GSE127910 (ABA vs DMSO, TPM fold-change)
# ══════════════════════════════════════════════════════════════
log("\n" + "=" * 55)
log("STEP 1c: ABA-axis signature (GSE127910)")
log("=" * 55)

aba_sig_file = SIG_DIR / "ABA_axis_logFC.csv"
ABA_sig = pd.Series(dtype=float)

GEO_FTP_ABA = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE127nnn/GSE127910/suppl"

if aba_sig_file.exists():
    ABA_sig = pd.read_csv(aba_sig_file, index_col=0)["logFC"]
    log(f"  Loaded existing ABA-axis signature: {len(ABA_sig)} genes")
else:
    aba_dir = SUPPL / "GSE127910"
    aba_dir.mkdir(exist_ok=True)

    ctrl_file = aba_dir / "GSE127910_Control.txt.gz"
    aba_file  = aba_dir / "GSE127910_Treatment_A.txt.gz"   # natural ABA

    ok_ctrl = download_file(f"{GEO_FTP_ABA}/GSE127910_Control.txt.gz", ctrl_file)
    ok_aba  = download_file(f"{GEO_FTP_ABA}/GSE127910_Treatment_A.txt.gz", aba_file)

    if ok_ctrl and ok_aba:
        def read_tpm_gz(path):
            with gzip.open(path, "rb") as f:
                raw = f.read()
            text = raw.decode("latin-1")
            rows = [l.split("\t") for l in text.strip().split("\n") if l]
            header = rows[0]
            data   = {r[0]: [float(x) for x in r[1:]] for r in rows[1:] if len(r) > 1}
            n_cols = len(rows[1]) - 1 if len(rows) > 1 else 1
            df = pd.DataFrame(data, index=range(n_cols)).T
            df.columns = header[1:] if len(header) > 1 else range(n_cols)
            return df

        try:
            ctrl_tpm = read_tpm_gz(ctrl_file)
            aba_tpm  = read_tpm_gz(aba_file)
            log(f"  Control TPM: {ctrl_tpm.shape}, ABA TPM: {aba_tpm.shape}")

            # Align on gene IDs
            common_genes = ctrl_tpm.index.intersection(aba_tpm.index)
            log(f"  Common genes: {len(common_genes)}")

            if len(common_genes) >= 1000:
                ctrl_mean = ctrl_tpm.loc[common_genes].astype(float).mean(axis=1)
                aba_mean  = aba_tpm.loc[common_genes].astype(float).mean(axis=1)
                # Avoid log(0)
                lfc = np.log2((aba_mean + 0.1) / (ctrl_mean + 0.1))
                ABA_sig = lfc.rename("logFC")
                ABA_sig.index.name = "GeneID"
                ABA_sig.to_frame().to_csv(aba_sig_file)
                log(f"  ABA-axis signature saved: {aba_sig_file.name} ({len(ABA_sig)} genes)")
            else:
                log("  Not enough common genes — check gene ID format")
        except Exception as e:
            log(f"  ABA TPM parsing error: {e}")
    else:
        log("  Download failed — ABA axis skipped")


# ══════════════════════════════════════════════════════════════
# STEP 1d: Calibration anchors
# Score reference samples on their own axes → absolute scale
#   P-axis  0% = HighPi samples   (GSE90061)
#   P-axis 100% = LowPi samples   (GSE90061)
#   N-axis  0% = KCl-ctrl t=120   (GSE97500)
#   N-axis 100% = N-supply t=120  (GSE97500)
#   ABA-axis 0% = DMSO control    (GSE127910)
#   ABA-axis 100% = ABA treatment (GSE127910)
# ══════════════════════════════════════════════════════════════
log("\n" + "=" * 55)
log("STEP 1d: Calibration anchors (reference experiments)")
log("=" * 55)

CAL_FILE = AXES_DIR / "calibration.csv"

def load_calibration():
    empty = pd.DataFrame(columns=["ctrl_score", "treat_score",
                                   "ctrl_label", "treat_label", "source"])
    empty.index.name = "axis"
    if not CAL_FILE.exists():
        return empty
    try:
        df = pd.read_csv(CAL_FILE)
        if "axis" in df.columns:
            df = df.set_index("axis")
        return df
    except Exception:
        return empty

def save_calibration(cal_df):
    cal_df.to_csv(CAL_FILE, index=True)

cal = load_calibration()

# --- P-axis calibration (GSE90061) ---
if "P_axis" not in cal.index and not P_sig.empty:
    try:
        log("  Computing P-axis calibration (GSE90061)...")
        counts_P_cal = read_htseq_counts(SUPPL / "GSE90061")
        hp_cal = [c for c in counts_P_cal.columns if "hp" in c]
        lp_cal = [c for c in counts_P_cal.columns if "lp" in c
                  and "lfe" not in c and "lpM" not in c]
        p_scores_hp = score_on_axis(counts_P_cal[hp_cal], P_sig)
        p_scores_lp = score_on_axis(counts_P_cal[lp_cal], P_sig)
        ctrl_p  = float(p_scores_hp.mean())
        treat_p = float(p_scores_lp.mean())
        cal.loc["P_axis"] = [ctrl_p, treat_p, "HighPi", "LowPi", "GSE90061"]
        save_calibration(cal)
        log(f"  P-axis: HighPi={ctrl_p:.0f}  LowPi={treat_p:.0f}  "
            f"range={treat_p - ctrl_p:.0f}")
    except Exception as e:
        log(f"  P calibration failed: {e}")
elif "P_axis" in cal.index:
    log(f"  P-axis calibration loaded: "
        f"ctrl={cal.loc['P_axis','ctrl_score']:.0f}  "
        f"treat={cal.loc['P_axis','treat_score']:.0f}")

# --- N-axis calibration (GSE97500, reuse already-parsed metadata) ---
if "N_axis" not in cal.index and not N_sig.empty:
    try:
        log("  Computing N-axis calibration (GSE97500)...")
        sm_path_n = MATRIX / "GSE97500_series_matrix.txt.gz"
        if sm_path_n.exists():
            meta_n = parse_series_matrix(sm_path_n)
            gsm_ids_n = meta_n.get("Sample_geo_accession", [])
            titles_n  = meta_n.get("Sample_title", [])
            gsm_to_t  = dict(zip(gsm_ids_n, titles_n))
            n120_N = [g for g, t in gsm_to_t.items() if t.startswith("S120N")]
            n120_C = [g for g, t in gsm_to_t.items() if t.startswith("S120C")]
            counts_N_cal = read_htseq_counts(SUPPL / "GSE97500")
            nc_cols = [c for c in counts_N_cal.columns
                       if any(c.startswith(g) for g in n120_C)]
            nt_cols = [c for c in counts_N_cal.columns
                       if any(c.startswith(g) for g in n120_N)]
            if nc_cols and nt_cols:
                n_scores_c = score_on_axis(counts_N_cal[nc_cols], N_sig)
                n_scores_t = score_on_axis(counts_N_cal[nt_cols], N_sig)
                ctrl_n  = float(n_scores_c.mean())
                treat_n = float(n_scores_t.mean())
                cal.loc["N_axis"] = [ctrl_n, treat_n,
                                     "KCl_t120", "N_supply_t120", "GSE97500"]
                save_calibration(cal)
                log(f"  N-axis: KCl={ctrl_n:.0f}  N-supply={treat_n:.0f}  "
                    f"range={treat_n - ctrl_n:.0f}")
    except Exception as e:
        log(f"  N calibration failed: {e}")
elif "N_axis" in cal.index:
    log(f"  N-axis calibration loaded: "
        f"ctrl={cal.loc['N_axis','ctrl_score']:.0f}  "
        f"treat={cal.loc['N_axis','treat_score']:.0f}")

# --- ABA-axis calibration (GSE127910, TPM) ---
if "ABA_axis" not in cal.index and not ABA_sig.empty:
    try:
        log("  Computing ABA-axis calibration (GSE127910)...")
        aba_dir_cal = SUPPL / "GSE127910"
        ctrl_f = aba_dir_cal / "GSE127910_Control.txt.gz"
        aba_f  = aba_dir_cal / "GSE127910_Treatment_A.txt.gz"
        if ctrl_f.exists() and aba_f.exists():
            def read_tpm_gz_cal(path):
                with gzip.open(path, "rb") as f:
                    raw = f.read()
                text = raw.decode("latin-1")
                rows = [l.split("\t") for l in text.strip().split("\n") if l]
                header = rows[0]
                data = {r[0]: [float(x) for x in r[1:]]
                        for r in rows[1:] if len(r) > 1}
                n_c = len(rows[1]) - 1 if len(rows) > 1 else 1
                df = pd.DataFrame(data, index=range(n_c)).T
                df.columns = header[1:] if len(header) > 1 else range(n_c)
                return df
            ct = read_tpm_gz_cal(ctrl_f)
            at = read_tpm_gz_cal(aba_f)
            aba_c_scores = score_on_axis_tpm(ct, ABA_sig)
            aba_t_scores = score_on_axis_tpm(at, ABA_sig)
            ctrl_a  = float(aba_c_scores.mean())
            treat_a = float(aba_t_scores.mean())
            cal.loc["ABA_axis"] = [ctrl_a, treat_a, "DMSO", "ABA", "GSE127910"]
            save_calibration(cal)
            log(f"  ABA-axis: DMSO={ctrl_a:.0f}  ABA={treat_a:.0f}  "
                f"range={treat_a - ctrl_a:.0f}")
    except Exception as e:
        log(f"  ABA calibration failed: {e}")
elif "ABA_axis" in cal.index:
    log(f"  ABA-axis calibration loaded: "
        f"ctrl={cal.loc['ABA_axis','ctrl_score']:.0f}  "
        f"treat={cal.loc['ABA_axis','treat_score']:.0f}")

log(f"\n  Calibration table:\n{cal.to_string()}")


# ══════════════════════════════════════════════════════════════
# STEP 2: Load biostimulant count data (GSE297649)
# ══════════════════════════════════════════════════════════════
log("\n" + "=" * 55)
log("STEP 2: Score biostimulant GSE297649")
log("  (Arabidopsis + amino acids + humic substances)")
log("=" * 55)

counts_B  = read_htseq_counts(SUPPL / "GSE297649")
aa_cols   = [c for c in counts_B.columns if "WTaa" in c]
sh_cols   = [c for c in counts_B.columns if "WTsh" in c]
ctrl_cols = [c for c in counts_B.columns if "WTcontrol" in c]
log(f"  Amino acid (WTaa): {aa_cols}")
log(f"  Humic subst (WTsh): {sh_cols}")
log(f"  Control: {ctrl_cols}")

# Axes available
axes = {
    "P_axis":   P_sig,
    "N_axis":   N_sig,
    "ABA_axis": ABA_sig,
}

all_rows = []
for group_label, group_cols, category in [
    ("amino_acid",   aa_cols,   "Amino acid biostimulant"),
    ("humic_subst",  sh_cols,   "Humic substance biostimulant"),
    ("control",      ctrl_cols, "Control"),
]:
    if not group_cols:
        continue
    for ax_name, sig in axes.items():
        if sig.empty:
            continue
        scores = score_on_axis(counts_B[group_cols], sig)
        for sample, score in scores.items():
            all_rows.append({
                "gse":      "GSE297649",
                "gsm":      sample,
                "group":    group_label,
                "category": category,
                "axis":     ax_name,
                "score":    float(score),
            })

per_sample_df = pd.DataFrame(all_rows)
if not per_sample_df.empty:
    per_sample_df.to_csv(SCORE_DIR / "GSE297649_per_sample_scores.csv", index=False)
    log(f"  Per-sample scores saved")

    summary = (per_sample_df
               .groupby(["gse", "category", "group", "axis"])["score"]
               .agg(["mean", "std", "count"])
               .reset_index())
    log("\n  Score summary:")
    log(summary.to_string(index=False))
    summary.to_csv(SCORE_DIR / "GSE297649_summary_scores.csv", index=False)


# ══════════════════════════════════════════════════════════════
# STEP 3: Unified output table (wide format, % of spread)
# ══════════════════════════════════════════════════════════════
log("\n" + "=" * 55)
log("STEP 3: Unified scores table")
log("=" * 55)

if not per_sample_df.empty:
    pivot = (per_sample_df
             .groupby(["gse", "category", "group", "axis"])["score"]
             .mean()
             .unstack("axis")
             .reset_index())

    # % calibrated against reference experiment anchors:
    #   0% = reference control state (HighPi / KCl / DMSO)
    #  100% = reference stressed/treated state (LowPi / N-supply / ABA)
    # Values outside 0-100 are valid — they mean "more extreme than the reference"
    for ax in [c for c in pivot.columns if c.endswith("_axis")]:
        if ax not in pivot.columns:
            continue
        if ax in cal.index:
            ctrl_anchor  = float(cal.loc[ax, "ctrl_score"])
            treat_anchor = float(cal.loc[ax, "treat_score"])
            rng = treat_anchor - ctrl_anchor
            if abs(rng) > 0:
                pivot[f"pct_{ax}"] = 100 * (pivot[ax] - ctrl_anchor) / rng
            else:
                pivot[f"pct_{ax}"] = np.nan
            log(f"  {ax}: anchors ctrl={ctrl_anchor:.0f}, treat={treat_anchor:.0f} "
                f"(from {cal.loc[ax,'source']})")
        else:
            # Fallback: within-experiment range if calibration missing
            log(f"  {ax}: no calibration — using within-experiment range (fallback)")
            ctrl_row  = pivot[pivot["group"] == "control"]
            ctrl_mean = float(ctrl_row[ax].values[0]) if not ctrl_row.empty else pivot[ax].min()
            rng = pivot[ax].max() - pivot[ax].min()
            pivot[f"pct_{ax}"] = 100 * (pivot[ax] - ctrl_mean) / rng if rng > 0 else 0.0

    out_file = ROOT / "all_axes_scores_rnaseq.csv"
    pivot.to_csv(out_file, index=False)
    log(f"  Saved: {out_file.name}")

    # Show pct columns
    pct_cols = ["category", "group"] + [c for c in pivot.columns if c.startswith("pct_")]
    log(f"\n{pivot[pct_cols].to_string(index=False)}")


# ══════════════════════════════════════════════════════════════
# STEP 4: Visualization
# ══════════════════════════════════════════════════════════════
log("\n" + "=" * 55)
log("STEP 4: Figures")
log("=" * 55)

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import matplotlib.colors as mcolors

    COLORS = {
        "amino_acid":  "#E69F00",
        "humic_subst": "#56B4E9",
        "control":     "#666666",
    }
    LABELS = {
        "amino_acid":  "Amino acids",
        "humic_subst": "Humic substances",
        "control":     "Control",
    }
    BIOSTIM_GROUPS = ["amino_acid", "humic_subst"]
    groups_order   = ["control", "amino_acid", "humic_subst"]

    # Build per-group pct lookup from pivot
    pct_lookup = {}   # group -> {axis -> pct_value}
    if 'pivot' in dir() and not pivot.empty:
        for _, row in pivot.iterrows():
            grp = row["group"]
            pct_lookup[grp] = {
                col.replace("pct_", ""): float(row[col])
                for col in pivot.columns if col.startswith("pct_")
            }

    pct_axes = sorted(pct_lookup.get("control", {}).keys())   # e.g. ["ABA_axis","N_axis","P_axis"]
    has_pct  = bool(pct_axes and "control" in pct_lookup)

    # ── Figure 1: biostimulant effect VECTORS (delta from control) ──────────────
    # Each biostimulant is a signed bar: how far it shifts each axis vs control
    if has_pct and len(pct_axes) >= 2:
        ctrl_vals = pct_lookup["control"]
        fig, ax = plt.subplots(figsize=(9, 5))
        x     = np.arange(len(pct_axes))
        width = 0.35
        ax_labels = {
            "P_axis":   "P-axis\n(← P-sufficient  P-starved →)",
            "N_axis":   "N-axis\n(← N-starved  N-supplied →)",
            "ABA_axis": "ABA-axis\n(← non-stressed  ABA-stressed →)",
        }
        for i, grp in enumerate(BIOSTIM_GROUPS):
            if grp not in pct_lookup:
                continue
            deltas = [pct_lookup[grp].get(ax_n, 0) - ctrl_vals.get(ax_n, 0)
                      for ax_n in pct_axes]
            # per-sample SD of delta (from per_sample_df)
            errs = []
            for ax_n in pct_axes:
                if ax_n in cal.index:
                    rng = float(cal.loc[ax_n, "treat_score"]) - float(cal.loc[ax_n, "ctrl_score"])
                    grp_scores = per_sample_df[
                        (per_sample_df["group"] == grp) &
                        (per_sample_df["axis"] == ax_n)]["score"]
                    ctrl_scores = per_sample_df[
                        (per_sample_df["group"] == "control") &
                        (per_sample_df["axis"] == ax_n)]["score"]
                    se = np.sqrt(grp_scores.std()**2 / len(grp_scores) +
                                 ctrl_scores.std()**2 / len(ctrl_scores)) / abs(rng) * 100
                    errs.append(se)
                else:
                    errs.append(0)
            bars = ax.bar(x + (i - 0.5) * width, deltas, width,
                          yerr=errs, capsize=4,
                          color=COLORS[grp], alpha=0.85,
                          label=LABELS[grp], zorder=3)
            # colour negative bars darker to show direction clearly
            for bar, d in zip(bars, deltas):
                if d < 0:
                    bar.set_alpha(0.6)

        ax.axhline(0, color="black", linewidth=1.0, zorder=4)
        ax.set_xticks(x)
        ax.set_xticklabels([ax_labels.get(a, a) for a in pct_axes], fontsize=9)
        ax.set_ylabel("Shift from control (% of reference range)", fontsize=10)
        ax.set_title("Biostimulant effect vectors — GSE297649, Arabidopsis\n"
                     "Direction shows which physiological state each treatment promotes",
                     fontsize=11, fontweight="bold")
        ax.legend(fontsize=10)
        ax.grid(axis="y", alpha=0.3, zorder=0)
        plt.tight_layout()
        out = PLOT_DIR / "GSE297649_effect_vectors.png"
        plt.savefig(out, dpi=200)
        plt.close()
        log(f"  Figure 1 (effect vectors) saved: {out.name}")

    # ── Figure 2: 2D state-space with arrows from control ──────────────────────
    # One panel per axis pair; control = start of arrow, treatment = tip
    axis_pairs = [("P_axis", "N_axis"), ("P_axis", "ABA_axis"), ("N_axis", "ABA_axis")]

    # Reference landmarks in calibrated % space (always at 0 or 100)
    landmarks = {
        "P_axis":   {"ctrl_label": "P-sufficient\n(HighPi)", "treat_label": "P-starved\n(LowPi)"},
        "N_axis":   {"ctrl_label": "N-limited\n(KCl ctrl)", "treat_label": "N-supplied\n(+KNO3)"},
        "ABA_axis": {"ctrl_label": "Unstressed\n(DMSO)",     "treat_label": "ABA-stressed\n(+ABA)"},
    }

    for ax_x, ax_y in axis_pairs:
        if not has_pct:
            continue
        if ax_x.replace("pct_","") not in pct_axes or ax_y.replace("pct_","") not in pct_axes:
            # strip pct_ prefix for lookup
            pass
        x_key = ax_x  # e.g. "P_axis"
        y_key = ax_y

        if x_key not in (pct_lookup.get("control", {})) or \
           y_key not in (pct_lookup.get("control", {})):
            continue

        ctrl_x = pct_lookup["control"][x_key]
        ctrl_y = pct_lookup["control"][y_key]

        fig, ax = plt.subplots(figsize=(8, 7))

        # Draw reference corners (0,0) and (100,100) as grey landmarks
        lm = landmarks.get(ax_x, {})
        lm2 = landmarks.get(ax_y, {})
        ax.scatter([0], [0], s=120, color="#dddddd", edgecolors="#aaaaaa",
                   zorder=2, marker="s")
        ax.annotate(f"{lm.get('ctrl_label','')} /\n{lm2.get('ctrl_label','')}",
                    (0, 0), textcoords="offset points", xytext=(-8, -28),
                    fontsize=7, color="#888888", ha="center")
        ax.scatter([100], [100], s=120, color="#ffdddd", edgecolors="#cc8888",
                   zorder=2, marker="s")
        ax.annotate(f"{lm.get('treat_label','')} /\n{lm2.get('treat_label','')}",
                    (100, 100), textcoords="offset points", xytext=(8, 8),
                    fontsize=7, color="#aa4444", ha="left")

        # Plot individual replicates as small dots
        for grp in groups_order:
            # Get per-sample scores in pct space
            x_samp = per_sample_df[
                (per_sample_df["group"] == grp) &
                (per_sample_df["axis"] == ax_x)
            ].set_index("gsm")["score"]
            y_samp = per_sample_df[
                (per_sample_df["group"] == grp) &
                (per_sample_df["axis"] == ax_y)
            ].set_index("gsm")["score"]
            common_s = x_samp.index.intersection(y_samp.index)
            if common_s.empty:
                continue

            # Convert raw scores to % using calibration
            def to_pct(scores, axis_key):
                if axis_key in cal.index:
                    c0 = float(cal.loc[axis_key, "ctrl_score"])
                    c1 = float(cal.loc[axis_key, "treat_score"])
                    return 100 * (scores - c0) / (c1 - c0)
                return scores

            xs = to_pct(x_samp[common_s], ax_x)
            ys = to_pct(y_samp[common_s], ax_y)

            ax.scatter(xs, ys,
                       color=COLORS.get(grp, "grey"),
                       s=50, alpha=0.6, edgecolors="white", zorder=3,
                       label=LABELS.get(grp, grp))

        # Draw arrows from control centroid to each biostimulant centroid
        for grp in BIOSTIM_GROUPS:
            if grp not in pct_lookup or "control" not in pct_lookup:
                continue
            tx = pct_lookup[grp][x_key]
            ty = pct_lookup[grp][y_key]
            ax.annotate(
                "", xy=(tx, ty), xytext=(ctrl_x, ctrl_y),
                arrowprops=dict(
                    arrowstyle="-|>",
                    color=COLORS[grp],
                    lw=2.5,
                    mutation_scale=18,
                ),
                zorder=5,
            )
            # Label at tip
            offset_x = 6 if tx >= ctrl_x else -6
            offset_y = 6 if ty >= ctrl_y else -6
            ax.annotate(
                LABELS[grp],
                (tx, ty),
                textcoords="offset points",
                xytext=(offset_x, offset_y),
                fontsize=10, fontweight="bold",
                color=COLORS[grp], zorder=6,
            )

        # Mark control centroid
        ax.scatter([ctrl_x], [ctrl_y], s=150, color=COLORS["control"],
                   edgecolors="black", zorder=4, marker="o")
        ax.annotate("Control", (ctrl_x, ctrl_y),
                    textcoords="offset points", xytext=(0, 8),
                    fontsize=9, ha="center", color="#444444")

        xlab = f"← P-sufficient       {ax_x.replace('_',' ')} (% of ref. range)       P-starved →"
        ylab = f"← low       {ax_y.replace('_',' ')} (% of ref. range)       high →"
        ax.set_xlabel(xlab if ax_x == "P_axis" else
                      f"{ax_x.replace('_',' ')} (% of reference range)", fontsize=9)
        ax.set_ylabel(ylab if ax_y == "ABA_axis" else
                      f"{ax_y.replace('_',' ')} (% of reference range)", fontsize=9)
        ax.set_title(
            f"Biostimulant state-space vectors\n"
            f"{ax_x.replace('_',' ')} vs {ax_y.replace('_',' ')} — GSE297649, Arabidopsis",
            fontsize=12, fontweight="bold"
        )
        # Grid lines at 0 and 100
        ax.axhline(0, color="#cccccc", lw=0.8, ls="--")
        ax.axvline(0, color="#cccccc", lw=0.8, ls="--")
        ax.axhline(100, color="#ffcccc", lw=0.8, ls="--")
        ax.axvline(100, color="#ffcccc", lw=0.8, ls="--")

        handles = [mpatches.Patch(color=COLORS[g], label=LABELS[g])
                   for g in groups_order]
        ax.legend(handles=handles, fontsize=9, loc="upper left")
        plt.tight_layout()
        out = PLOT_DIR / f"GSE297649_vectors_{ax_x}_vs_{ax_y}.png"
        plt.savefig(out, dpi=200)
        plt.close()
        log(f"  Figure 2 saved: {out.name}")

    # ── Figure 3: vector fingerprint — delta bars per biostimulant ───────────────
    # Faceted: one panel per biostimulant, bars show shift per axis
    if has_pct and len(pct_axes) >= 2:
        ctrl_vals = pct_lookup["control"]
        n_bio = len(BIOSTIM_GROUPS)
        fig, axs = plt.subplots(1, n_bio, figsize=(5 * n_bio, 5),
                                sharey=True)
        if n_bio == 1:
            axs = [axs]

        ax_nice = {
            "P_axis":   "P-axis\n(P-starvation)",
            "N_axis":   "N-axis\n(N-supply)",
            "ABA_axis": "ABA-axis\n(stress)",
        }
        x = np.arange(len(pct_axes))
        for i, grp in enumerate(BIOSTIM_GROUPS):
            if grp not in pct_lookup:
                continue
            deltas = [pct_lookup[grp].get(a, 0) - ctrl_vals.get(a, 0)
                      for a in pct_axes]
            bar_colors = [COLORS[grp] if d >= 0 else mcolors.to_rgba(COLORS[grp], 0.5)
                          for d in deltas]
            axs[i].bar(x, deltas, color=bar_colors, edgecolor="white", linewidth=0.5)
            axs[i].axhline(0, color="black", lw=1.0)
            axs[i].set_xticks(x)
            axs[i].set_xticklabels([ax_nice.get(a, a) for a in pct_axes],
                                   fontsize=9)
            axs[i].set_title(LABELS[grp], fontsize=11, fontweight="bold",
                             color=COLORS[grp])
            axs[i].set_ylabel("Shift from control (% of ref. range)" if i == 0 else "",
                              fontsize=9)
            axs[i].grid(axis="y", alpha=0.3)
            # Add value labels on bars
            for j, (bar, d) in enumerate(zip(axs[i].patches, deltas)):
                axs[i].text(bar.get_x() + bar.get_width() / 2,
                            d + (0.3 if d >= 0 else -0.5),
                            f"{d:+.1f}%", ha="center", va="bottom" if d >= 0 else "top",
                            fontsize=8)

        fig.suptitle("Biostimulant physiological fingerprints\n"
                     "Shift from control on absolute calibrated axes",
                     fontsize=12, fontweight="bold")
        plt.tight_layout()
        out = PLOT_DIR / "GSE297649_vector_fingerprints.png"
        plt.savefig(out, dpi=200)
        plt.close()
        log(f"  Figure 3 (fingerprints) saved: {out.name}")

except ImportError:
    log("  matplotlib not available — skipping figures")
except Exception as e:
    import traceback
    log(f"  Figure error: {e}")
    log(traceback.format_exc())


# ══════════════════════════════════════════════════════════════
# Summary
# ══════════════════════════════════════════════════════════════
log("\n" + "=" * 55)
log("PIPELINE COMPLETE")
log("=" * 55)
log(f"Signatures: {[f.name for f in SIG_DIR.glob('*.csv')]}")
log(f"Scores:     {[f.name for f in SCORE_DIR.glob('*.csv')]}")
log(f"Figures:    {[f.name for f in PLOT_DIR.glob('*.png')]}")

log_file = AXES_DIR / f"python_run_log_{datetime.now().strftime('%Y%m%d_%H%M')}.txt"
log_file.write_text("\n".join(LOG), encoding="utf-8")
log(f"Log: {log_file.name}")
