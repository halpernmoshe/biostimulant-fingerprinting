"""
Re-score GSE297649 using log2(CPM+1) instead of raw log2(count+1).
CPM normalization removes library-size depth artifacts.
Compares raw vs normalized results.
"""
import tarfile, gzip
import pandas as pd
import numpy as np
from pathlib import Path

ROOT = Path("C:/Users/moshe/Dropbox/ISF 2025")
TAR  = ROOT / "RNA_seq_data/suppl_files/GSE297649/GSE297649_RAW.tar"
SIG_DIR = ROOT / "RNA_seq_axes/signatures"
CAL_FILE = ROOT / "RNA_seq_axes/calibration.csv"
SCORES_DIR = ROOT / "RNA_seq_axes/scores"

SPECIAL = {"__no_feature", "__ambiguous", "__too_low_aQual",
           "__not_aligned", "__alignment_not_unique"}

groups = {
    "GSM8996165_htseqcount_WTaa_1":       "amino_acid",
    "GSM8996166_htseqcount_WTaa_2":       "amino_acid",
    "GSM8996167_htseqcount_WTaa_3":       "amino_acid",
    "GSM8996168_htseqcount_WTcontrol_1":  "control",
    "GSM8996169_htseqcount_WTcontrol_2":  "control",
    "GSM8996170_htseqcount_WTcontrol_3":  "control",
    "GSM8996171_htseqcount_WTsh_1":       "humic_subst",
    "GSM8996172_htseqcount_WTsh_2":       "humic_subst",
    "GSM8996173_htseqcount_WTsh_3":       "humic_subst",
}

# ── Load count matrix ─────────────────────────────────────────────────────────
counts = {}
with tarfile.open(TAR) as t:
    for member in t.getmembers():
        name = member.name.replace(".txt.gz", "")
        f = t.extractfile(member)
        with gzip.open(f) as gz:
            lines = gz.read().decode("utf-8").strip().split("\n")
        gene2cnt = {}
        for line in lines:
            gene, cnt = line.split("\t")
            if gene not in SPECIAL:
                gene2cnt[gene] = int(cnt)
        counts[name] = gene2cnt

samples = list(counts.keys())
all_genes = sorted(set(g for d in counts.values() for g in d))
df = pd.DataFrame({s: counts[s] for s in samples}, index=all_genes).fillna(0).astype(int)

# CPM normalization
lib_sizes = df.sum(axis=0)
cpm = df.divide(lib_sizes, axis=1) * 1_000_000
log_cpm = np.log2(cpm + 1)  # log2(CPM+1)
log_raw = np.log2(df + 1)   # original raw scoring

print("Library sizes:")
for s in samples:
    print(f"  {s}: {lib_sizes[s]:,}")

# ── Load calibration ──────────────────────────────────────────────────────────
cal = pd.read_csv(CAL_FILE, index_col=0)
print("\nCalibration:", cal)

def calibrate(score, axis):
    lo = cal.loc[axis, "ctrl_score"]
    hi = cal.loc[axis, "treat_score"]
    return (score - lo) / (hi - lo) * 100

# ── Load and score on each axis ───────────────────────────────────────────────
axes = ["P_axis", "N_axis", "ABA_axis"]

results_raw = []
results_cpm = []

for axis in axes:
    sig_file = SIG_DIR / f"{axis}_logFC.csv"
    sig = pd.read_csv(sig_file, index_col=0)
    sig.columns = ["logFC"]

    # Intersect with available genes
    common = sig.index.intersection(df.index)
    sig_c = sig.loc[common, "logFC"]
    print(f"\n{axis}: {len(sig_c)} genes in common")

    for s in samples:
        # Raw score
        raw_score = (log_raw.loc[common, s] * sig_c).sum()
        raw_pct = calibrate(raw_score, axis)
        results_raw.append({"sample": s, "group": groups[s], "axis": axis,
                            "raw_score": raw_score, "raw_pct": raw_pct})

        # CPM score
        cpm_score = (log_cpm.loc[common, s] * sig_c).sum()
        cpm_pct = calibrate(cpm_score, axis)
        results_cpm.append({"sample": s, "group": groups[s], "axis": axis,
                            "cpm_score": cpm_score, "cpm_pct": cpm_pct})

df_raw = pd.DataFrame(results_raw)
df_cpm = pd.DataFrame(results_cpm)

# ── Compare group means: raw vs CPM ──────────────────────────────────────────
print("\n=== Effect vectors: RAW log2(count+1) vs CPM log2(CPM+1) ===")
ctrl_groups = {"P_axis": "control", "N_axis": "control", "ABA_axis": "control"}

for axis in axes:
    print(f"\n--- {axis} ---")
    for grp in ["amino_acid", "control", "humic_subst"]:
        raw_mean = df_raw[(df_raw.axis==axis) & (df_raw.group==grp)]["raw_pct"].mean()
        cpm_mean = df_cpm[(df_cpm.axis==axis) & (df_cpm.group==grp)]["cpm_pct"].mean()
        print(f"  {grp}: raw={raw_mean:.1f}% | CPM={cpm_mean:.1f}%")

    raw_ctrl  = df_raw[(df_raw.axis==axis) & (df_raw.group=="control")]["raw_pct"].mean()
    raw_aa    = df_raw[(df_raw.axis==axis) & (df_raw.group=="amino_acid")]["raw_pct"].mean()
    raw_hs    = df_raw[(df_raw.axis==axis) & (df_raw.group=="humic_subst")]["raw_pct"].mean()
    cpm_ctrl  = df_cpm[(df_cpm.axis==axis) & (df_cpm.group=="control")]["cpm_pct"].mean()
    cpm_aa    = df_cpm[(df_cpm.axis==axis) & (df_cpm.group=="amino_acid")]["cpm_pct"].mean()
    cpm_hs    = df_cpm[(df_cpm.axis==axis) & (df_cpm.group=="humic_subst")]["cpm_pct"].mean()

    print(f"  DELTA amino_acid - control: RAW={raw_aa-raw_ctrl:+.2f}% | CPM={cpm_aa-cpm_ctrl:+.2f}%")
    print(f"  DELTA humic_subst - control: RAW={raw_hs-raw_ctrl:+.2f}% | CPM={cpm_hs-cpm_ctrl:+.2f}%")

# ── Save CPM-normalized scores ────────────────────────────────────────────────
out = SCORES_DIR / "GSE297649_per_sample_scores_CPM.csv"
df_cpm.to_csv(out, index=False)
print(f"\nSaved CPM scores: {out}")

# ── Summary comparison table ──────────────────────────────────────────────────
print("\n=== SUMMARY: Effect vector comparison ===")
print(f"{'Treatment':<20} {'Axis':<12} {'Raw Δ%':>8} {'CPM Δ%':>8} {'Change':>8}")
for axis in axes:
    raw_ctrl = df_raw[(df_raw.axis==axis) & (df_raw.group=="control")]["raw_pct"].mean()
    cpm_ctrl = df_cpm[(df_cpm.axis==axis) & (df_cpm.group=="control")]["cpm_pct"].mean()
    for grp, label in [("amino_acid","Amino acids"), ("humic_subst","Humic subst")]:
        raw_val = df_raw[(df_raw.axis==axis) & (df_raw.group==grp)]["raw_pct"].mean()
        cpm_val = df_cpm[(df_cpm.axis==axis) & (df_cpm.group==grp)]["cpm_pct"].mean()
        raw_d = raw_val - raw_ctrl
        cpm_d = cpm_val - cpm_ctrl
        change = cpm_d - raw_d
        same_dir = "✓" if (raw_d * cpm_d > 0) or (abs(raw_d) < 0.5 and abs(cpm_d) < 0.5) else "⚠️ FLIP"
        print(f"  {label:<18} {axis:<12} {raw_d:>+8.2f} {cpm_d:>+8.2f} {change:>+8.2f}  {same_dir}")
