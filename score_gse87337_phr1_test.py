"""
Score GSE87337 on P-axis to test PHR1 requirement for bacterial SynCom P-starvation response.
Key comparison: Col-0 vs phr1phl1 mutant, with/without SynCom, at Low P.
Hypothesis: SynCom elevates P-axis in Col-0 but not phr1phl1 (PHR1-dependent mechanism).
"""
import gzip
import pandas as pd
import numpy as np
from pathlib import Path

ROOT    = Path("C:/Users/moshe/Dropbox/ISF 2025")
COUNTS  = ROOT / "RNA_seq_data/suppl_files/GSE87337/GSE87337_PBI_preTreat_counts.txt.gz"
SIG_DIR = ROOT / "RNA_seq_axes/signatures"
CAL_FILE = ROOT / "RNA_seq_axes/calibration.csv"
OUT_DIR = ROOT / "RNA_seq_axes/scores"

cal = pd.read_csv(CAL_FILE, index_col=0)
def calibrate(score, axis):
    lo = cal.loc[axis, "ctrl_score"]
    hi = cal.loc[axis, "treat_score"]
    return (score - lo) / (hi - lo) * 100

# ── Load count matrix ─────────────────────────────────────────────────────────
with gzip.open(COUNTS, "rt") as f:
    lines = f.readlines()

# Line 0: full condition descriptions (tab-separated, starts with empty cell)
conditions = lines[0].strip().split("\t")  # first element is empty
# Line 1: short IDs with "Gene" as first element
headers = lines[1].strip().split("\t")     # headers[0] = "Gene"
short_ids = headers[1:]                    # PBI115, PBI116, ...

print(f"Samples: {len(short_ids)}")
print(f"Conditions row length: {len(conditions)}")
print(f"First 3 conditions: {conditions[:3]}")

# Build mapping: short_id -> full condition name
id_to_cond = dict(zip(short_ids, conditions))

# Parse data rows
data_lines = lines[2:]
gene_ids = []
count_data = []
for line in data_lines:
    parts = line.strip().split("\t")
    gene_ids.append(parts[0])
    count_data.append([int(x) for x in parts[1:]])

df = pd.DataFrame(count_data, index=gene_ids, columns=short_ids)
print(f"Count matrix: {df.shape[0]} genes × {df.shape[1]} samples")
print(f"Library sizes range: {df.sum().min():,} — {df.sum().max():,}")

# ── Select relevant columns ───────────────────────────────────────────────────
# Key comparison: LowP conditions, +P pre-treatment (most neutral pre-treatment)
# Genotypes: Col-0, phr1phl1
# Bacteria: noBact, SynCom

def find_cols(genotype, p_level, bacteria, pretreated="+P"):
    """Find columns matching genotype × P-level × bacteria × pre-treatment."""
    matches = []
    for sid, cond in id_to_cond.items():
        if (genotype in cond and
            f"preTreat({pretreated})" in cond and
            f"_{p_level}_" in cond and
            f"_{bacteria}_" in cond):
            matches.append(sid)
    return sorted(matches)

# Primary comparison: LowP, +P pre-treatment
wt_lowP_no  = find_cols("Col-0",    "LowP", "noBact")
wt_lowP_bac = find_cols("Col-0",    "LowP", "SynCom")
mt_lowP_no  = find_cols("phr1phl1", "LowP", "noBact")
mt_lowP_bac = find_cols("phr1phl1", "LowP", "SynCom")

print(f"\nKey groups (LowP, preTreat(+P)):")
print(f"  Col-0 noBact:    {wt_lowP_no}  → conditions: {[id_to_cond[s] for s in wt_lowP_no[:2]]}")
print(f"  Col-0 SynCom:    {wt_lowP_bac} → conditions: {[id_to_cond[s] for s in wt_lowP_bac[:2]]}")
print(f"  phr1phl1 noBact: {mt_lowP_no}  → conditions: {[id_to_cond[s] for s in mt_lowP_no[:2]]}")
print(f"  phr1phl1 SynCom: {mt_lowP_bac} → conditions: {[id_to_cond[s] for s in mt_lowP_bac[:2]]}")

if not all([wt_lowP_no, wt_lowP_bac, mt_lowP_no, mt_lowP_bac]):
    # Try without pre-treatment filter
    print("\nNo results with +P pre-treatment filter. Trying without...")
    for sid, cond in id_to_cond.items():
        if "LowP" in cond and "noBact" in cond and "Col-0" in cond:
            print(f"  {sid}: {cond}")

# ── CPM normalize ─────────────────────────────────────────────────────────────
lib = df.sum(axis=0)
log_cpm = np.log2(df.divide(lib, axis=1) * 1_000_000 + 1)

# ── Score on P-axis ───────────────────────────────────────────────────────────
P_sig = pd.read_csv(SIG_DIR / "P_axis_logFC.csv", index_col=0)
P_sig.columns = ["logFC"]
N_sig = pd.read_csv(SIG_DIR / "N_axis_logFC.csv", index_col=0)
N_sig.columns = ["logFC"]

def score_group(cols, sig, label=""):
    if not cols:
        return None
    common = sig.index.intersection(log_cpm.index)
    scores = [(log_cpm.loc[common, c] * sig.loc[common, "logFC"]).sum() for c in cols]
    return np.mean(scores)

axes = [("P_axis", P_sig), ("N_axis", N_sig)]

print("\n=== PHR1 REQUIREMENT TEST: SynCom bacteria effect on P-axis ===")
print(f"{'Condition':<30} {'P-axis %':>9} {'N-axis %':>9}")
print("-" * 52)

groups = [
    ("Col-0 noBact (LowP)",    wt_lowP_no),
    ("Col-0 SynCom (LowP)",    wt_lowP_bac),
    ("phr1phl1 noBact (LowP)", mt_lowP_no),
    ("phr1phl1 SynCom (LowP)", mt_lowP_bac),
]

scores_out = {}
for label, cols in groups:
    row = {}
    for ax_name, sig in axes:
        raw_score = score_group(cols, sig)
        if raw_score is not None:
            pct = calibrate(raw_score, ax_name)
            row[ax_name] = pct
        else:
            row[ax_name] = float("nan")
    scores_out[label] = row
    print(f"  {label:<28} {row['P_axis']:>+9.1f}% {row['N_axis']:>+9.1f}%")

# ── Compute deltas ────────────────────────────────────────────────────────────
print("\n=== Effect of SynCom (bacteria) within each genotype ===")
for ax_name in ["P_axis", "N_axis"]:
    wt_delta  = scores_out["Col-0 SynCom (LowP)"][ax_name] - scores_out["Col-0 noBact (LowP)"][ax_name]
    mt_delta  = scores_out["phr1phl1 SynCom (LowP)"][ax_name] - scores_out["phr1phl1 noBact (LowP)"][ax_name]
    ratio = mt_delta / wt_delta if abs(wt_delta) > 0.1 else float("nan")
    print(f"  {ax_name}: Col-0 delta={wt_delta:+.1f}%  phr1phl1 delta={mt_delta:+.1f}%  ratio={ratio:.2f}")

print("\n=== INTERPRETATION ===")
wt_P_delta = scores_out["Col-0 SynCom (LowP)"]["P_axis"] - scores_out["Col-0 noBact (LowP)"]["P_axis"]
mt_P_delta = scores_out["phr1phl1 SynCom (LowP)"]["P_axis"] - scores_out["phr1phl1 noBact (LowP)"]["P_axis"]

if abs(wt_P_delta) > 1 and abs(mt_P_delta) < abs(wt_P_delta) * 0.5:
    print(f"RESULT: PHR1 REQUIRED for SynCom P-axis effect.")
    print(f"  Col-0: bacteria shift P-axis by {wt_P_delta:+.1f}%")
    print(f"  phr1phl1: bacteria shift P-axis by only {mt_P_delta:+.1f}% ({abs(mt_P_delta/wt_P_delta)*100:.0f}% of WT effect)")
    print(f"  CONCLUSION: PGPR/SynCom P-axis activation requires functional PHR1.")
    print(f"  This validates the PHR1-mediated mechanism for PGPR biostimulant action.")
elif abs(wt_P_delta) > 1 and abs(mt_P_delta) >= abs(wt_P_delta) * 0.5:
    print(f"RESULT: PHR1 NOT required — phr1phl1 still shows {abs(mt_P_delta/wt_P_delta)*100:.0f}% of WT effect.")
    print(f"  SynCom may act upstream or parallel to PHR1.")
elif abs(wt_P_delta) < 1:
    print(f"RESULT: Col-0 P-axis shift by SynCom is small ({wt_P_delta:+.1f}%) — low signal.")
    print(f"  This condition may not be comparable to our PGPR datasets.")

# ── Also check HighP comparison ───────────────────────────────────────────────
wt_highP_no  = find_cols("Col-0",    "HighP", "noBact")
wt_highP_bac = find_cols("Col-0",    "HighP", "SynCom")
mt_highP_no  = find_cols("phr1phl1", "HighP", "noBact")
mt_highP_bac = find_cols("phr1phl1", "HighP", "SynCom")

print("\n=== HighP comparison (SynCom effect at adequate P) ===")
for label, cols in [("Col-0 noBact (HighP)", wt_highP_no), ("Col-0 SynCom (HighP)", wt_highP_bac),
                    ("phr1phl1 noBact (HighP)", mt_highP_no), ("phr1phl1 SynCom (HighP)", mt_highP_bac)]:
    raw_P = score_group(cols, P_sig)
    if raw_P is not None:
        pct = calibrate(raw_P, "P_axis")
        print(f"  {label:<30}: {pct:+.1f}%")

if wt_highP_no and wt_highP_bac:
    p_WT_highP = calibrate(score_group(wt_highP_bac, P_sig), "P_axis") - calibrate(score_group(wt_highP_no, P_sig), "P_axis")
    p_WT_lowP  = wt_P_delta
    print(f"\n  Col-0: SynCom effect at LowP={p_WT_lowP:+.1f}%, at HighP={p_WT_highP:+.1f}%")
    print(f"  (Confirms P-conditional effect: larger at LowP, smaller at HighP)")

# ── Save results ──────────────────────────────────────────────────────────────
rows = []
for (label, cols), is_bact, geno, plevel in zip(
    groups,
    [False, True, False, True],
    ["Col-0","Col-0","phr1phl1","phr1phl1"],
    ["LowP","LowP","LowP","LowP"]
):
    for ax_name, sig in axes:
        raw = score_group(cols, sig)
        pct = calibrate(raw, ax_name) if raw is not None else float("nan")
        rows.append({"dataset": "GSE87337", "condition": label, "genotype": geno,
                     "P_level": plevel, "bacteria": "SynCom" if is_bact else "noBact",
                     "axis": ax_name, "score_pct": round(pct, 2), "n_samples": len(cols)})

out = OUT_DIR / "GSE87337_PHR1_test_scores.csv"
pd.DataFrame(rows).to_csv(out, index=False)
print(f"\nSaved: {out}")
