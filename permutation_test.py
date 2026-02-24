"""
Permutation test for axis scores: is the treatment vs control difference
greater than expected by chance?
Uses exact permutation (all C(n_total, n_treat) label assignments).
CPM-normalized counts throughout.
"""
import tarfile, gzip
import pandas as pd
import numpy as np
from itertools import combinations
from pathlib import Path

ROOT    = Path("C:/Users/moshe/Dropbox/ISF 2025")
SIG_DIR = ROOT / "RNA_seq_axes/signatures"
CAL_FILE = ROOT / "RNA_seq_axes/calibration.csv"
SUPPL   = ROOT / "RNA_seq_data/suppl_files"

SPECIAL = {"__no_feature", "__ambiguous", "__too_low_aQual",
           "__not_aligned", "__alignment_not_unique"}

cal = pd.read_csv(CAL_FILE, index_col=0)
def calibrate(score, axis):
    lo = cal.loc[axis, "ctrl_score"]
    hi = cal.loc[axis, "treat_score"]
    return (score - lo) / (hi - lo) * 100

# Pre-load axis signatures
axes = ["P_axis", "N_axis", "ABA_axis"]
sigs = {}
for ax in axes:
    s = pd.read_csv(SIG_DIR / f"{ax}_logFC.csv", index_col=0)
    s.columns = ["logFC"]
    sigs[ax] = s

def score_sample(log_cpm_vec, sig):
    """Dot product of log2(CPM+1) with logFC signature."""
    common = sig.index.intersection(log_cpm_vec.index)
    return (log_cpm_vec.loc[common] * sig.loc[common, "logFC"]).sum()

def permutation_test(log_cpm_df, treat_cols, ctrl_cols, axis, n_perm=None):
    """
    Exact permutation test.
    n_perm=None means use all C(n_total, n_treat) permutations (exact test).
    Returns: observed_delta (%), p_value, n_perms_used
    """
    sig = sigs[axis]
    all_cols = treat_cols + ctrl_cols
    n_treat = len(treat_cols)
    n_total = len(all_cols)

    # Observed scores
    treat_scores = [score_sample(log_cpm_df[c], sig) for c in treat_cols]
    ctrl_scores  = [score_sample(log_cpm_df[c], sig) for c in ctrl_cols]
    obs_delta_raw = np.mean(treat_scores) - np.mean(ctrl_scores)
    obs_delta_pct = calibrate(np.mean(treat_scores), axis) - calibrate(np.mean(ctrl_scores), axis)

    # Null distribution: all possible ways to assign n_treat cols as "treatment"
    from math import comb
    n_possible = comb(n_total, n_treat)

    null_deltas = []
    all_col_scores = {c: score_sample(log_cpm_df[c], sig) for c in all_cols}

    for treat_idx in combinations(range(n_total), n_treat):
        ctrl_idx = [i for i in range(n_total) if i not in treat_idx]
        perm_treat = np.mean([all_col_scores[all_cols[i]] for i in treat_idx])
        perm_ctrl  = np.mean([all_col_scores[all_cols[i]] for i in ctrl_idx])
        null_deltas.append(perm_treat - perm_ctrl)

    null_deltas = np.array(null_deltas)

    # p-value: proportion of permutations where |delta| >= |observed|
    # (two-tailed)
    p_val = np.mean(np.abs(null_deltas) >= np.abs(obs_delta_raw))

    return obs_delta_pct, p_val, n_possible

# ── GSE297649: load from tar ───────────────────────────────────────────────────
print("=== GSE297649: Amino acids + Humic substances ===")
TAR = SUPPL / "GSE297649/GSE297649_RAW.tar"
counts_297 = {}
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
        counts_297[name] = gene2cnt

all_genes_297 = sorted(set(g for d in counts_297.values() for g in d))
df297 = pd.DataFrame(counts_297, index=all_genes_297).fillna(0).astype(int)
lib297 = df297.sum(axis=0)
cpm297 = df297.divide(lib297, axis=1) * 1_000_000
log_cpm297 = np.log2(cpm297 + 1)

# Group assignments
aa_cols  = [c for c in log_cpm297.columns if "WTaa" in c]
hs_cols  = [c for c in log_cpm297.columns if "WTsh" in c]
ctrl_cols = [c for c in log_cpm297.columns if "WTcontrol" in c]
print(f"  amino_acid: {aa_cols}")
print(f"  humic_subst: {hs_cols}")
print(f"  control: {ctrl_cols}")

print(f"\nResults (all C(n,k) exact permutations):")
print(f"{'Treatment':<20} {'Axis':<12} {'Delta%':>8} {'p-value':>8} {'n_perms':>8}")

for label, treat, ctrl in [("amino_acid", aa_cols, ctrl_cols),
                             ("humic_subst", hs_cols, ctrl_cols)]:
    for ax in axes:
        delta, pval, nperms = permutation_test(log_cpm297, treat, ctrl, ax)
        sig_flag = " *" if pval < 0.05 else ("  ~" if pval < 0.10 else "")
        print(f"  {label:<18} {ax:<12} {delta:>+8.2f} {pval:>8.3f} {nperms:>8}{sig_flag}")

# ── GSE138478: load from featureCounts ────────────────────────────────────────
print("\n=== GSE138478: PGPR volatile + diacetyl ===")

def read_fc(gz_path):
    with gzip.open(gz_path, "rt") as f:
        lines = [l for l in f.readlines() if not l.startswith("#")]
    from io import StringIO
    df = pd.read_csv(StringIO("".join(lines)), sep="\t")
    sample_cols = df.columns[6:].tolist()
    return df.set_index("Geneid")[sample_cols]

hmz008 = read_fc(SUPPL / "GSE138478/GSE138478_HMZ008_featureCounts_output.txt.gz")
hmz014 = read_fc(SUPPL / "GSE138478/GSE138478_HMZ014_Dione_featureCounts_output.txt.gz")

def cpm_log(df):
    lib = df.sum(axis=0)
    return np.log2(df.divide(lib, axis=1) * 1_000_000 + 1)

lcpm008 = cpm_log(hmz008)
lcpm014 = cpm_log(hmz014)

# HMZ008 groups: C = control, T = treated
gmv_ctrl  = [c for c in lcpm008.columns if "_C_" in c or "C_5DAT" in c]
gmv_treat = [c for c in lcpm008.columns if "_T_" in c or "T_5DAT" in c]
# HMZ014 groups: CK = control, Dione = treated
dia_ctrl  = [c for c in lcpm014.columns if "CK_" in c]
dia_treat = [c for c in lcpm014.columns if "Dione_" in c]

print(f"  GMV control: {gmv_ctrl}")
print(f"  GMV treated: {gmv_treat}")
print(f"  Diacetyl control: {dia_ctrl}")
print(f"  Diacetyl treated: {dia_treat}")

print(f"\nResults:")
for label, log_cpm, treat, ctrl in [
        ("PGPR GMV", lcpm008, gmv_treat, gmv_ctrl),
        ("PGPR diacetyl", lcpm014, dia_treat, dia_ctrl)]:
    for ax in axes:
        delta, pval, nperms = permutation_test(log_cpm, treat, ctrl, ax)
        sig_flag = " *" if pval < 0.05 else ("  ~" if pval < 0.10 else "")
        print(f"  {label:<18} {ax:<12} {delta:>+8.2f} {pval:>8.3f} {nperms:>8}{sig_flag}")

# ── Summary table ─────────────────────────────────────────────────────────────
print("\n\n=== FULL SUMMARY (for manuscript) ===")
print("Treatment           Axis         DeltaPct   p-value   Significant?")
print("-" * 65)

results_all = []

for label, log_cpm, treat, ctrl in [
        ("AA (GSE297649)",   log_cpm297, aa_cols,   ctrl_cols),
        ("HS (GSE297649)",   log_cpm297, hs_cols,   ctrl_cols),
        ("PGPR GMV",         lcpm008,   gmv_treat, gmv_ctrl),
        ("PGPR diacetyl",    lcpm014,   dia_treat, dia_ctrl)]:
    for ax in axes:
        delta, pval, nperms = permutation_test(log_cpm, treat, ctrl, ax)
        sig = "yes *" if pval < 0.05 else ("marginal" if pval < 0.10 else "no")
        print(f"  {label:<20} {ax:<12} {delta:>+8.2f}   {pval:.3f}    {sig}")
        results_all.append({"treatment": label, "axis": ax,
                            "delta_pct": round(delta, 2),
                            "p_value": round(pval, 3),
                            "n_permutations": nperms,
                            "significant": sig})

# Save
out = ROOT / "RNA_seq_axes/scores/permutation_test_results.csv"
pd.DataFrame(results_all).to_csv(out, index=False)
print(f"\nSaved: {out}")
