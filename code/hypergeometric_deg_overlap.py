"""
Hypergeometric DEG overlap test for amino acid reproducibility.

GSE297649: protein hydrolysate (amino acid class)
GSE298982: glutamate as sole N source

Uses the pre-computed gene-level logFC from amino_acid_gene_level_coherence.csv.
Defines DEGs by |logFC| thresholds (sensitivity analysis at 0.5, 1.0, 1.5).
Reports hypergeometric p-value, observed overlap, expected overlap, and fold enrichment.
"""

import numpy as np
import pandas as pd
from scipy import stats

ROOT = "C:/Users/moshe/Dropbox/ISF 2025"
COHERENCE_CSV = f"{ROOT}/RNA_seq_axes/scores/amino_acid_gene_level_coherence.csv"

print("=" * 70)
print("HYPERGEOMETRIC DEG OVERLAP TEST")
print("GSE297649 (protein hydrolysate) vs GSE298982 (glutamate sole-N)")
print("=" * 70)

df = pd.read_csv(COHERENCE_CSV, index_col=0)
# Drop genes with NaN in either logFC column
df = df[["logFC_hydrolysate", "logFC_glutamate"]].dropna()

M = len(df)  # total expressed genes (population size)
print(f"\nTotal expressed genes (M): {M:,}")

thresholds = [0.5, 1.0, 1.5, 2.0]

results = []
for thresh in thresholds:
    # Define DEG sets
    deg_297 = set(df.index[df["logFC_hydrolysate"].abs() >= thresh])
    deg_298 = set(df.index[df["logFC_glutamate"].abs() >= thresh])

    # Overlap
    overlap = deg_297 & deg_298

    n_297 = len(deg_297)   # "successes in population"
    N_298 = len(deg_298)   # "draws"
    k     = len(overlap)   # observed successes

    # Expected overlap by chance
    expected = n_297 * N_298 / M

    # Fold enrichment
    fold = k / expected if expected > 0 else float("inf")

    # Hypergeometric p-value: P(X >= k)
    # scipy hypergeom: (M, n, N) = (population, success_states, draws)
    p_hyper = stats.hypergeom.sf(k - 1, M, n_297, N_298)

    results.append({
        "threshold": thresh,
        "n_DEG_hydrolysate": n_297,
        "n_DEG_glutamate": N_298,
        "overlap_observed": k,
        "overlap_expected": round(expected, 1),
        "fold_enrichment": round(fold, 2),
        "p_hypergeometric": p_hyper,
    })

    print(f"\n--- |logFC| >= {thresh} ---")
    print(f"  DEGs in GSE297649 (hydrolysate): {n_297:,}")
    print(f"  DEGs in GSE298982 (glutamate):   {N_298:,}")
    print(f"  Overlap observed:                {k:,}")
    print(f"  Overlap expected by chance:      {expected:.1f}")
    print(f"  Fold enrichment:                 {fold:.2f}x")
    print(f"  Hypergeometric p-value:          {p_hyper:.2e}")

# --------------------------------------------------------------------------
# Direction-agreement test (sign test on all expressed genes)
# --------------------------------------------------------------------------
print()
print("=" * 70)
print("SIGN CONSISTENCY TEST (all expressed genes)")
print("=" * 70)

same_sign = ((df["logFC_hydrolysate"] > 0) == (df["logFC_glutamate"] > 0)).sum()
total     = len(df)
frac      = same_sign / total

# Binomial test: expected 50% same sign under null
binom_p = stats.binomtest(same_sign, total, p=0.5, alternative="greater").pvalue

print(f"\n  Same-direction genes: {same_sign:,} / {total:,} ({frac*100:.1f}%)")
print(f"  Expected under null (random): 50.0%")
print(f"  Binomial test p-value:        {binom_p:.2e}")

# --------------------------------------------------------------------------
# Direction test restricted to axis-marker genes
# --------------------------------------------------------------------------
print()
print("=" * 70)
print("SIGN CONSISTENCY TEST (axis marker genes only)")
print("=" * 70)

df_full = pd.read_csv(COHERENCE_CSV, index_col=0).dropna(
    subset=["logFC_hydrolysate", "logFC_glutamate"]
)

for axis in ["P", "N", "JA", "SA", "ABA", "Fe"]:
    col = f"top200_{axis}_axis"
    if col not in df_full.columns:
        continue
    sub = df_full[df_full[col] == 1][["logFC_hydrolysate", "logFC_glutamate"]]
    if len(sub) == 0:
        continue
    same = ((sub["logFC_hydrolysate"] > 0) == (sub["logFC_glutamate"] > 0)).sum()
    n    = len(sub)
    p_b  = stats.binomtest(same, n, p=0.5, alternative="greater").pvalue
    print(f"  {axis}-axis top-200:  {same}/{n} same direction ({same/n*100:.1f}%), "
          f"binomial p = {p_b:.2e}")

# --------------------------------------------------------------------------
# Summary for manuscript
# --------------------------------------------------------------------------
print()
print("=" * 70)
print("MANUSCRIPT-READY NUMBERS")
print("=" * 70)
best = next(r for r in results if r["threshold"] == 1.0)
print(f"""
At |logFC| >= 1.0:
  - {best['overlap_observed']:,} genes were differentially expressed in both
    amino acid experiments (observed overlap)
  - Expected by chance: {best['overlap_expected']:.0f} genes
  - Fold enrichment: {best['fold_enrichment']:.1f}x
  - Hypergeometric p = {best['p_hypergeometric']:.2e}

Suggested sentence for §3.3:
  "Significant overlap was observed between the differentially expressed
  gene sets from the two experiments ({best['overlap_observed']:,} shared
  genes with |log2FC| >= 1, versus {best['overlap_expected']:.0f} expected
  by chance; {best['fold_enrichment']:.1f}-fold enrichment,
  hypergeometric p = {best['p_hypergeometric']:.2e})."
""")

print("Done.")
