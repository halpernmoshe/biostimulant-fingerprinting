"""
PHR1 validation figure: SynCom bacteria P-axis effect in Col-0 vs phr1phl1 mutant.
GSE87337: 144 samples, Arabidopsis + bacterial synthetic community (SynCom), at Low and High P.
"""
import gzip
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

ROOT    = Path("C:/Users/moshe/Dropbox/ISF 2025")
COUNTS  = ROOT / "RNA_seq_data/suppl_files/GSE87337/GSE87337_PBI_preTreat_counts.txt.gz"
SIG_DIR = ROOT / "RNA_seq_axes/signatures"
CAL_FILE = ROOT / "RNA_seq_axes/calibration.csv"
FIG_DIR = ROOT / "state_space_figures"

cal = pd.read_csv(CAL_FILE, index_col=0)
def calibrate(s, ax):
    lo, hi = cal.loc[ax, "ctrl_score"], cal.loc[ax, "treat_score"]
    return (s - lo) / (hi - lo) * 100

with gzip.open(COUNTS, "rt") as f:
    lines = f.readlines()
conditions = lines[0].strip().split("\t")
headers = lines[1].strip().split("\t")
short_ids = headers[1:]
id_to_cond = dict(zip(short_ids, conditions))
data_lines = lines[2:]
gene_ids = [l.split("\t")[0] for l in data_lines]
count_data = [[int(x) for x in l.strip().split("\t")[1:]] for l in data_lines]
df = pd.DataFrame(count_data, index=gene_ids, columns=short_ids)
lib = df.sum(axis=0)
log_cpm = np.log2(df.divide(lib, axis=1) * 1_000_000 + 1)

sig = pd.read_csv(SIG_DIR / "P_axis_logFC.csv", index_col=0)
sig.columns = ["logFC"]
common = sig.index.intersection(log_cpm.index)
sig_c = sig.loc[common, "logFC"]

def get_cols(genotype, p_level, bacteria, pretreated="+P"):
    return sorted([s for s, c in id_to_cond.items()
                   if genotype in c and f"preTreat({pretreated})" in c and
                   f"_{p_level}_" in c and f"_{bacteria}_" in c])

def mean_pct(cols, axis="P_axis"):
    scores = [(log_cpm.loc[common, c] * sig_c).sum() for c in cols]
    return np.mean([calibrate(s, axis) for s in scores])

def per_sample_pct(cols, axis="P_axis"):
    return [calibrate((log_cpm.loc[common, c] * sig_c).sum(), axis) for c in cols]

# ── Build data for figure ─────────────────────────────────────────────────────
# Panel A: LowP — Col-0 and phr1phl1, noBact vs SynCom
groups_A = [
    ("Col-0\nnoBact",    get_cols("Col-0",    "LowP", "noBact")),
    ("Col-0\nSynCom",    get_cols("Col-0",    "LowP", "SynCom")),
    ("phr1phl1\nnoBact", get_cols("phr1phl1", "LowP", "noBact")),
    ("phr1phl1\nSynCom", get_cols("phr1phl1", "LowP", "SynCom")),
]

# Panel B: P-concentration effect in Col-0
groups_B = [
    ("LowP\nnoBact",  get_cols("Col-0", "LowP",  "noBact")),
    ("LowP\nSynCom",  get_cols("Col-0", "LowP",  "SynCom")),
    ("HighP\nnoBact", get_cols("Col-0", "HighP", "noBact")),
    ("HighP\nSynCom", get_cols("Col-0", "HighP", "SynCom")),
]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 6))

# Panel A: PHR1 requirement
colors_A = ["#888888", "#009E73", "#888888", "#D55E00"]
x_A = np.arange(4)
means_A = [mean_pct(cols) for _, cols in groups_A]
samples_A = [per_sample_pct(cols) for _, cols in groups_A]
labels_A = [l for l, _ in groups_A]

bars_A = ax1.bar(x_A, means_A, color=colors_A, width=0.6, alpha=0.85,
                 edgecolor="black", lw=0.7)

# Individual points
for i, pts in enumerate(samples_A):
    x_jitter = np.random.normal(x_A[i], 0.05, len(pts))
    ax1.scatter(x_jitter, pts, color="white", edgecolor="#333333",
                s=30, zorder=5, lw=0.8)

# Delta annotations
wt_delta  = means_A[1] - means_A[0]
mt_delta  = means_A[3] - means_A[2]
ax1.annotate("", xy=(x_A[1], means_A[1]), xytext=(x_A[1], means_A[0]),
             arrowprops=dict(arrowstyle="<->", color="#009E73", lw=2))
ax1.text(x_A[1] + 0.15, (means_A[0] + means_A[1]) / 2, f"+{wt_delta:.0f}%",
         color="#009E73", fontsize=10, fontweight="bold", va="center")
ax1.annotate("", xy=(x_A[3], means_A[3]), xytext=(x_A[3], means_A[2]),
             arrowprops=dict(arrowstyle="<->", color="#D55E00", lw=2))
ax1.text(x_A[3] + 0.15, (means_A[2] + means_A[3]) / 2, f"+{mt_delta:.0f}%",
         color="#D55E00", fontsize=10, fontweight="bold", va="center")

ax1.axvline(1.5, color="#aaaaaa", lw=1, ls="--", alpha=0.7)
ax1.set_xticks(x_A)
ax1.set_xticklabels(labels_A, fontsize=10)
ax1.set_ylabel("P-axis score (% of Pi-starvation range)", fontsize=10)
ax1.set_title("A. PHR1 requirement for SynCom P-axis activation\n(LowP, preTreat +P)",
              fontsize=10, fontweight="bold")
ax1.text(0.75, ax1.get_ylim()[0] + 1, "Col-0", ha="center", fontsize=9,
         color="#009E73", style="italic")
ax1.text(2.5, ax1.get_ylim()[0] + 1, "phr1phl1", ha="center", fontsize=9,
         color="#D55E00", style="italic")
ax1.yaxis.grid(True, alpha=0.3, ls="--")
ax1.set_axisbelow(True)

# Panel B: P-concentration dependence in Col-0
colors_B = ["#888888", "#4477AA", "#aaaaaa", "#88AACC"]
x_B = np.arange(4)
means_B = [mean_pct(cols) for _, cols in groups_B]
samples_B = [per_sample_pct(cols) for _, cols in groups_B]
labels_B = [l for l, _ in groups_B]

bars_B = ax2.bar(x_B, means_B, color=colors_B, width=0.6, alpha=0.85,
                 edgecolor="black", lw=0.7)

for i, pts in enumerate(samples_B):
    x_jitter = np.random.normal(x_B[i], 0.05, len(pts))
    ax2.scatter(x_jitter, pts, color="white", edgecolor="#333333",
                s=30, zorder=5, lw=0.8)

lowP_delta  = means_B[1] - means_B[0]
highP_delta = means_B[3] - means_B[2]
ax2.annotate("", xy=(x_B[1], means_B[1]), xytext=(x_B[1], means_B[0]),
             arrowprops=dict(arrowstyle="<->", color="#4477AA", lw=2))
ax2.text(x_B[1] + 0.15, (means_B[0] + means_B[1]) / 2, f"+{lowP_delta:.0f}%",
         color="#4477AA", fontsize=10, fontweight="bold", va="center")
ax2.annotate("", xy=(x_B[3], means_B[3]), xytext=(x_B[3], means_B[2]),
             arrowprops=dict(arrowstyle="<->", color="#88AACC", lw=2))
ax2.text(x_B[3] + 0.15, (means_B[2] + means_B[3]) / 2, f"+{highP_delta:.0f}%",
         color="#88AACC", fontsize=10, fontweight="bold", va="center")

ax2.axvline(1.5, color="#aaaaaa", lw=1, ls="--", alpha=0.7)
ax2.set_xticks(x_B)
ax2.set_xticklabels(labels_B, fontsize=10)
ax2.set_ylabel("P-axis score (% of Pi-starvation range)", fontsize=10)
ax2.set_title("B. P-concentration dependence of SynCom effect\n(Col-0, preTreat +P)",
              fontsize=10, fontweight="bold")
ax2.text(0.75, ax2.get_ylim()[0] + 1, "Low Pi", ha="center", fontsize=9,
         color="#4477AA", style="italic")
ax2.text(2.5, ax2.get_ylim()[0] + 1, "High Pi", ha="center", fontsize=9,
         color="#88AACC", style="italic")
ax2.yaxis.grid(True, alpha=0.3, ls="--")
ax2.set_axisbelow(True)

fig.suptitle(
    "Bacterial SynCom P-axis activation: PHR1-dependent and Pi-conditional\n"
    "(GSE87337: Arabidopsis whole seedlings, n=4 samples per group)",
    fontsize=11, fontweight="bold"
)
fig.tight_layout()
fname = FIG_DIR / "GSE87337_PHR1_validation.png"
fig.savefig(fname, dpi=150)
print(f"Saved: {fname}")
plt.close(fig)
