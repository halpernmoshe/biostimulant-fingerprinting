"""
Score GSE208223 (TiO2 nanoparticles) on the 3 new DESeq2 axes:
  JA_axis, SA_axis, Fe_axis
Then update biostimulant_6axis_summary.csv with TiO2 row.

GSE208223 provides rLog-normalized values directly (not raw counts).
rLog values ARE in log2 space, so use them directly for the dot product.
Columns: E1/E2/E3 = control, F1/F2/F3 = TiO2_treated.
"""
import sys, gzip
sys.stdout.reconfigure(encoding='utf-8')
import pandas as pd
import numpy as np
from pathlib import Path

ROOT     = Path("C:/Users/moshe/Dropbox/ISF 2025")
SIG_DIR  = ROOT / "RNA_seq_axes/signatures"
SCORE_DIR = ROOT / "RNA_seq_axes/scores"
CAL_FILE = ROOT / "RNA_seq_axes/calibration.csv"
GSE208223_FILE = ROOT / "RNA_seq_data/suppl_files/GSE208223/GSE208223_Shoot_norm.txt.gz"

# ── Load data ─────────────────────────────────────────────────────────────────
print("Loading GSE208223 rLog-normalized expression...")
df = pd.read_csv(GSE208223_FILE, sep='\t', index_col=0, compression='gzip')
# Strip isoform suffix from gene IDs: "AT1G01010.1" -> "AT1G01010"
df.index = df.index.str.replace(r'\.\d+$', '', regex=True)
df.index = df.index.str.strip('"')
print(f"  Shape: {df.shape}")
print(f"  Columns: {list(df.columns)}")
print(f"  First 3 gene IDs: {list(df.index[:3])}")
print(f"  Sample values (row 0): {df.iloc[0].to_dict()}")

# Control = E1,E2,E3; TiO2-treated = F1,F2,F3
ctrl_cols  = ['E1', 'E2', 'E3']
treat_cols = ['F1', 'F2', 'F3']

# ── Load calibration ──────────────────────────────────────────────────────────
cal = pd.read_csv(CAL_FILE, index_col=0)
print(f"\nCalibration axes available: {list(cal.index)}")

def calibrate(score, ax_name):
    lo = cal.loc[ax_name, 'ctrl_score']
    hi = cal.loc[ax_name, 'treat_score']
    return (float(score) - lo) / (hi - lo) * 100

# ── Score function ────────────────────────────────────────────────────────────
def score_on_sig(expr_df, sig):
    """Dot-product: sum(expr * logFC) for common genes."""
    common = expr_df.index.intersection(sig.index)
    print(f"  Common genes: {len(common)}")
    return expr_df.loc[common].multiply(sig.loc[common], axis=0).sum(axis=0)

# ── Load and score on each new axis ─────────────────────────────────────────
new_axes = ['JA_axis', 'SA_axis', 'Fe_axis']
results = {}  # axis -> delta_pct

for ax_name in new_axes:
    sig_file = SIG_DIR / f"{ax_name}_logFC.csv"
    sig = pd.read_csv(sig_file, index_col=0).squeeze()
    sig.index = sig.index.str.strip()
    print(f"\n{ax_name}: loaded {len(sig)} genes from signature")

    scores = score_on_sig(df, sig)
    ctrl_score  = scores[ctrl_cols].mean()
    treat_score = scores[treat_cols].mean()
    delta_raw = treat_score - ctrl_score
    delta_pct = calibrate(treat_score, ax_name) - calibrate(ctrl_score, ax_name)

    print(f"  Control score:  {ctrl_score:.2f} ({calibrate(ctrl_score, ax_name):+.1f}%)")
    print(f"  Treated score:  {treat_score:.2f} ({calibrate(treat_score, ax_name):+.1f}%)")
    print(f"  Delta: {delta_raw:+.2f} raw  |  {delta_pct:+.2f}% of calibration range")

    results[ax_name] = delta_pct

print("\n=== TiO2 delta summary (new axes) ===")
for ax, pct in results.items():
    print(f"  {ax}: {pct:+.1f}%")

# ── Update biostimulant_6axis_summary.csv ────────────────────────────────────
summary_file = SCORE_DIR / "biostimulant_6axis_summary.csv"
print(f"\nLoading current 6-axis summary: {summary_file}")
summary = pd.read_csv(summary_file)
print(f"Current shape: {summary.shape}")
print(summary.to_string())

# Check if TiO2 already present
if 'TiO2_treated' in summary['treatment'].values:
    # Update existing JA/SA/Fe rows
    for ax_name, pct in results.items():
        mask = (summary['treatment'] == 'TiO2_treated') & (summary['axis'] == ax_name)
        if mask.any():
            summary.loc[mask, 'delta_pct'] = pct
            print(f"  Updated TiO2_treated {ax_name}: {pct:+.2f}%")
        else:
            new_row = pd.DataFrame([{'treatment': 'TiO2_treated', 'axis': ax_name, 'delta_pct': pct}])
            summary = pd.concat([summary, new_row], ignore_index=True)
            print(f"  Added TiO2_treated {ax_name}: {pct:+.2f}%")
else:
    # Need to get TiO2 P/N/ABA scores from existing data
    # effect_vectors_all_axes.csv stores raw log2-space deltas, must calibrate to %
    ev = pd.read_csv(SCORE_DIR / "effect_vectors_all_axes.csv")
    tio2_existing = ev[ev['treatment'] == 'TiO2_treated'][['axis', 'delta']].copy()

    # Map short axis names to calibration.csv index names
    AX_MAP = {'P': 'P_axis', 'N': 'N_axis', 'ABA': 'ABA_axis'}
    new_rows = []
    for _, row in tio2_existing.iterrows():
        ax_short = row['axis']
        if ax_short not in AX_MAP:
            continue  # skip gene-set axes (already have gene-set scale, not calibrated)
        ax_mapped = AX_MAP[ax_short]
        raw_delta = row['delta']
        # Calibrate: delta_pct = raw_delta / (treat_score - ctrl_score) * 100
        cal_range = cal.loc[ax_mapped, 'treat_score'] - cal.loc[ax_mapped, 'ctrl_score']
        pct = raw_delta / cal_range * 100
        new_rows.append({'treatment': 'TiO2_treated', 'axis': ax_mapped, 'delta_pct': pct})
        print(f"  TiO2_treated {ax_mapped}: raw_delta={raw_delta:+.1f}, cal_range={cal_range:.0f}, pct={pct:+.2f}%")
    for ax_name, pct in results.items():
        new_rows.append({'treatment': 'TiO2_treated', 'axis': ax_name, 'delta_pct': pct})
        print(f"  TiO2_treated {ax_name}: {pct:+.2f}%")
    summary = pd.concat([summary, pd.DataFrame(new_rows)], ignore_index=True)

summary.to_csv(summary_file, index=False)
print(f"\nSaved updated 6-axis summary: {summary_file}")

# ── Print complete 6-axis table ───────────────────────────────────────────────
print("\n=== COMPLETE 6-AXIS BIOSTIMULANT TABLE ===")
axes_order = ['P_axis', 'N_axis', 'ABA_axis', 'JA_axis', 'SA_axis', 'Fe_axis']
summary_pivot = summary.pivot(index='treatment', columns='axis', values='delta_pct')
for ax in axes_order:
    if ax not in summary_pivot.columns:
        summary_pivot[ax] = np.nan

print(f"\n{'Treatment':<25}", end="")
for ax in axes_order:
    print(f" {ax:>10}", end="")
print()
print("-" * 90)
for t in summary_pivot.index:
    print(f"{t:<25}", end="")
    for ax in axes_order:
        val = summary_pivot.loc[t, ax] if ax in summary_pivot.columns else np.nan
        if pd.isna(val):
            print(f" {'NaN':>10}", end="")
        else:
            print(f" {val:>+10.1f}%", end="")
    print()

print("\nDone.")
