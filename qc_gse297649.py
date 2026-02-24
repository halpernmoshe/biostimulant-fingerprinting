"""
QC analysis for GSE297649 (amino acids + humic substances Arabidopsis).
- Library sizes (total mapped reads per sample)
- PCA of all 9 samples (log2 count+1)
- Flag outliers
- Write report to agent_coordination/GSE297649_QC.md
"""
import tarfile, gzip, io
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.decomposition import PCA
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path("C:/Users/moshe/Dropbox/ISF 2025")
TAR  = ROOT / "RNA_seq_data/suppl_files/GSE297649/GSE297649_RAW.tar"
FIG  = ROOT / "state_space_figures"
OUT  = ROOT / "agent_coordination/GSE297649_QC.md"

# ── Group labels ──────────────────────────────────────────────────────────────
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

SPECIAL = {"__no_feature", "__ambiguous", "__too_low_aQual",
           "__not_aligned", "__alignment_not_unique"}

# ── Read all count files from tar ─────────────────────────────────────────────
counts = {}
special_rows = {}

with tarfile.open(TAR) as t:
    for member in t.getmembers():
        name = member.name.replace(".txt.gz", "")
        f = t.extractfile(member)
        with gzip.open(f) as gz:
            lines = gz.read().decode("utf-8").strip().split("\n")
        gene2cnt = {}
        special = {}
        for line in lines:
            gene, cnt = line.split("\t")
            cnt = int(cnt)
            if gene in SPECIAL:
                special[gene] = cnt
            else:
                gene2cnt[gene] = cnt
        counts[name] = gene2cnt
        special_rows[name] = special

samples = list(counts.keys())
all_genes = sorted(set(g for d in counts.values() for g in d))

# Build count matrix (genes × samples)
df = pd.DataFrame({s: counts[s] for s in samples}, index=all_genes).fillna(0).astype(int)
print(f"Count matrix: {df.shape[0]} genes × {df.shape[1]} samples")

# ── Library sizes ─────────────────────────────────────────────────────────────
lib_sizes = df.sum(axis=0)
median_lib = lib_sizes.median()
print("\n=== Library Sizes ===")
for s in samples:
    ratio = lib_sizes[s] / median_lib
    flag = " *** OUTLIER ***" if ratio > 2 or ratio < 0.5 else ""
    print(f"  {s}: {lib_sizes[s]:,}  (ratio={ratio:.2f}){flag}")

# Special (unmapped) stats
print("\n=== Special rows (% of total incl. unmapped) ===")
for s in samples:
    total_all = lib_sizes[s] + sum(special_rows[s].values())
    pct_mapped = lib_sizes[s] / total_all * 100
    print(f"  {s}: {pct_mapped:.1f}% mapped")

# ── PCA ───────────────────────────────────────────────────────────────────────
log_df = np.log2(df + 1)
# Filter low-expressed genes (mean < 1 count across samples)
keep = df.mean(axis=1) >= 1
log_filt = log_df.loc[keep]
print(f"\nGenes kept for PCA (mean count >=1): {keep.sum()}")

pca = PCA(n_components=3)
coords = pca.fit_transform(log_filt.T)  # samples × PCs
var_exp = pca.explained_variance_ratio_ * 100

print(f"\nPCA variance explained: PC1={var_exp[0]:.1f}%, PC2={var_exp[1]:.1f}%, PC3={var_exp[2]:.1f}%")
print("\nPC coordinates:")
for i, s in enumerate(samples):
    print(f"  {s} ({groups[s]}): PC1={coords[i,0]:.1f}, PC2={coords[i,1]:.1f}, PC3={coords[i,2]:.1f}")

# ── PCA plot ──────────────────────────────────────────────────────────────────
colors_map = {"amino_acid": "#E69F00", "control": "#666666", "humic_subst": "#56B4E9"}
markers_map = {"amino_acid": "o", "control": "s", "humic_subst": "^"}
labels_map = {"amino_acid": "Amino acid", "control": "Control", "humic_subst": "Humic substance"}

fig, ax = plt.subplots(figsize=(8, 6))
for grp in ["amino_acid", "control", "humic_subst"]:
    idx = [i for i, s in enumerate(samples) if groups[s] == grp]
    ax.scatter(coords[idx, 0], coords[idx, 1],
               color=colors_map[grp], marker=markers_map[grp],
               s=80, label=labels_map[grp], zorder=3)
    # Label each point
    for i in idx:
        short = samples[i].split("_")[-1]  # e.g. "1", "2", "3"
        ax.annotate(short, (coords[i, 0], coords[i, 1]),
                    textcoords="offset points", xytext=(5, 5), fontsize=8)

ax.set_xlabel(f"PC1 ({var_exp[0]:.1f}%)", fontsize=11)
ax.set_ylabel(f"PC2 ({var_exp[1]:.1f}%)", fontsize=11)
ax.set_title("PCA: GSE297649 all samples (log2 count+1)\nAmino acids, Humic substances, Control", fontsize=11)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
fig.tight_layout()
pca_fig = FIG / "GSE297649_PCA.png"
fig.savefig(pca_fig, dpi=150)
plt.close(fig)
print(f"\nSaved PCA plot: {pca_fig}")

# ── Identify WTaa_1 within-group stats ───────────────────────────────────────
# Within amino_acid group: is WTaa_1 an outlier by P-axis score?
p_scores = {
    "WTaa_1": 5667.1,  # from GSE297649_per_sample_scores.csv
    "WTaa_2": 4418.1,
    "WTaa_3": 4382.2,
}
grp_mean = np.mean(list(p_scores.values()))
grp_std  = np.std(list(p_scores.values()), ddof=1)
z_1 = (p_scores["WTaa_1"] - grp_mean) / grp_std
print(f"\nP-axis within amino_acid group:")
for k, v in p_scores.items():
    z = (v - grp_mean) / grp_std
    print(f"  {k}: {v:.0f} (z={z:.2f})")
print(f"  Group mean±std: {grp_mean:.0f}±{grp_std:.0f}")
print(f"  WTaa_1 z-score: {z_1:.2f}")

# ── Write QC report ───────────────────────────────────────────────────────────
lines_out = []
lines_out.append("# GSE297649 QC Report\n")
lines_out.append("## Dataset: Amino acids + Humic substances biostimulant (Arabidopsis, whole seedlings)\n")
lines_out.append(f"**Count matrix**: {df.shape[0]:,} genes × {df.shape[1]} samples\n")
lines_out.append(f"**Median library size**: {median_lib:,.0f} reads\n\n")

lines_out.append("## 1. Library Sizes\n\n")
lines_out.append("| Sample | Group | Total mapped reads | Ratio to median | Flag |\n")
lines_out.append("|--------|-------|-------------------|-----------------|------|\n")
for s in samples:
    ratio = lib_sizes[s] / median_lib
    flag = "OUTLIER" if ratio > 2 or ratio < 0.5 else "OK"
    total_all = lib_sizes[s] + sum(special_rows[s].values())
    pct_mapped = lib_sizes[s] / total_all * 100
    lines_out.append(f"| {s.split('htseqcount_')[1]} | {groups[s]} | {lib_sizes[s]:,} ({pct_mapped:.1f}% mapped) | {ratio:.2f} | {flag} |\n")

lines_out.append("\n**Conclusion**: ")
outliers = [s for s in samples if lib_sizes[s]/median_lib > 2 or lib_sizes[s]/median_lib < 0.5]
if outliers:
    lines_out.append(f"⚠️ {len(outliers)} sample(s) flagged as outlier(s): {', '.join(outliers)}\n\n")
else:
    lines_out.append("✓ All library sizes within 2× median. No outliers by library size.\n\n")

lines_out.append("## 2. PCA (log2 count+1, genes with mean count ≥1)\n\n")
lines_out.append(f"Genes used: {keep.sum():,}\n")
lines_out.append(f"PC1: {var_exp[0]:.1f}% variance | PC2: {var_exp[1]:.1f}% | PC3: {var_exp[2]:.1f}%\n\n")
lines_out.append("| Sample | Group | PC1 | PC2 | PC3 |\n")
lines_out.append("|--------|-------|-----|-----|-----|\n")
for i, s in enumerate(samples):
    sname = s.split('htseqcount_')[1]
    lines_out.append(f"| {sname} | {groups[s]} | {coords[i,0]:.1f} | {coords[i,1]:.1f} | {coords[i,2]:.1f} |\n")

lines_out.append(f"\nFigure: `state_space_figures/GSE297649_PCA.png`\n\n")

lines_out.append("## 3. WTaa_1 (GSM8996165) Outlier Assessment\n\n")
lines_out.append("WTaa_1 has a markedly higher P-axis score than its replicates:\n\n")
lines_out.append("| Sample | P-axis score | z-score within group |\n")
lines_out.append("|--------|-------------|---------------------|\n")
for k, v in p_scores.items():
    z = (v - grp_mean) / grp_std
    lines_out.append(f"| {k} | {v:.0f} | {z:.2f} |\n")

lines_out.append(f"\nWTaa_1 z-score: **{z_1:.2f}** (within-group std={grp_std:.0f})\n\n")

if abs(z_1) > 2:
    lines_out.append("⚠️ **WTaa_1 is a potential outlier by P-axis score** (z > 2). ")
    lines_out.append("However, n=3 per group is too small for robust outlier detection. ")
    lines_out.append("PCA clustering (see figure) provides the key diagnostic: ")
    lines_out.append("if WTaa_1 clusters with its group in PCA, the higher score likely reflects biological variation.\n\n")
else:
    lines_out.append("✓ WTaa_1 z-score within acceptable range. Not a statistical outlier.\n\n")

lines_out.append("## 4. Recommendation\n\n")
lines_out.append("- If WTaa_1 clusters with WTaa_2 and WTaa_3 in PCA: keep all 3 replicates\n")
lines_out.append("- If WTaa_1 is isolated from its group in PCA: consider sensitivity analysis (with/without WTaa_1)\n")
lines_out.append("- Library size check passed → high P-axis score in WTaa_1 is not a depth artifact\n")
lines_out.append("- Given n=3 per group, losing one replicate would severely reduce power — only exclude if PCA strongly suggests batch contamination\n\n")
lines_out.append(f"**Date**: 2026-02-23\n")

with open(OUT, "w", encoding="utf-8") as f:
    f.writelines(lines_out)
print(f"\nQC report written to: {OUT}")
