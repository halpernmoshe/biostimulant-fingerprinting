"""
Delta-vector bar chart: compare all biostimulant treatments side-by-side.
Shows ΔP-axis, ΔN-axis, ΔABA-axis for each treatment within its experiment.
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

ROOT    = Path("C:/Users/moshe/Dropbox/ISF 2025")
FIG_DIR = ROOT / "state_space_figures"

# ── CPM-corrected delta values (treatment − control, within-experiment) ───────
# Updated 2026-02-23: CPM normalization applied to GSE297649 and GSE138478
# (GSE208223 used rLog-normalized series matrix; GSE199501 used CPM input)
deltas = [
    # label                        ΔP     ΔN     ΔABA  dataset         color
    ("Amino acids\n(GSE297649)",    +0.1,  -2.6,  +2.1, "AA/HA",       "#E69F00"),
    ("Humic substances\n(GSE297649)", -2.3, -1.9, +2.2, "AA/HA",       "#56B4E9"),
    ("PGPR volatile\n(GSE138478)",  +6.5,  -4.3,  +0.3, "PGPR",        "#009E73"),
    ("PGPR diacetyl\n(GSE138478)", +8.9,  +9.7,  -3.9, "PGPR",        "#CC79A7"),
    ("Ti nanoparticles\n(GSE208223)", -1.4, -5.9,  +3.5, "Nanomaterial","#D55E00"),
    ("P. megaterium\n(GSE199501)",  +5.8,  +2.3,  +1.0, "PGPR",        "#8B5CF6"),
]

labels  = [d[0] for d in deltas]
d_p     = [d[1] for d in deltas]
d_n     = [d[2] for d in deltas]
d_aba   = [d[3] for d in deltas]
colors  = [d[5] for d in deltas]

x = np.arange(len(labels))
width = 0.25

fig, ax = plt.subplots(figsize=(13, 6))

bars_p   = ax.bar(x - width,     d_p,   width, label="ΔP-axis",   color="#4477AA", alpha=0.85, edgecolor="black", lw=0.5)
bars_n   = ax.bar(x,             d_n,   width, label="ΔN-axis",   color="#66AADD", alpha=0.85, edgecolor="black", lw=0.5)
bars_aba = ax.bar(x + width,     d_aba, width, label="ΔABA-axis", color="#EE8866", alpha=0.85, edgecolor="black", lw=0.5)

# Value labels
for bar in [bars_p, bars_n, bars_aba]:
    for rect in bar:
        h = rect.get_height()
        ypos = rect.get_y() + h + (0.3 if h >= 0 else -0.8)
        ax.annotate(f"{h:+.1f}",
                    xy=(rect.get_x() + rect.get_width() / 2, ypos),
                    fontsize=7.5, ha='center', va='bottom' if h >= 0 else 'top',
                    color='#333333')

ax.axhline(0, color='black', lw=0.8)
ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=9)
ax.set_ylabel("Δ Axis score (% of biological reference range)", fontsize=10)
ax.set_title(
    "Biostimulant effect vectors: Δ-scores relative to within-experiment control\n"
    "(CPM-normalized; positive = toward reference stress/starvation state)",
    fontsize=11, fontweight='bold'
)
ax.legend(fontsize=10, loc='upper right')
ax.set_ylim(min(min(d_p), min(d_n), min(d_aba)) - 3,
            max(max(d_p), max(d_n), max(d_aba)) + 5)
ax.yaxis.grid(True, alpha=0.4, ls='--')
ax.set_axisbelow(True)

# Separator lines between experiments
ax.axvline(1.5, color='#aaaaaa', lw=1, ls='--', alpha=0.6)
ax.axvline(3.5, color='#aaaaaa', lw=1, ls='--', alpha=0.6)
ax.text(0.5, ax.get_ylim()[1] * 0.95, "GSE297649", ha='center', fontsize=8,
        color='#555555', fontstyle='italic')
ax.text(2.5, ax.get_ylim()[1] * 0.95, "GSE138478", ha='center', fontsize=8,
        color='#555555', fontstyle='italic')
ax.text(4.0, ax.get_ylim()[1] * 0.95, "GSE208223", ha='center', fontsize=8,
        color='#555555', fontstyle='italic')
ax.axvline(4.5, color='#aaaaaa', lw=1, ls='--', alpha=0.6)
ax.text(5.0, ax.get_ylim()[1] * 0.95, "GSE199501", ha='center', fontsize=8,
        color='#555555', fontstyle='italic')

fig.tight_layout()
fname = FIG_DIR / "delta_vectors_all_treatments_CPM.png"
fig.savefig(fname, dpi=150)
print(f"Saved: {fname}")
plt.close(fig)
