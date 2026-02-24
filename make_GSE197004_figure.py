#!/usr/bin/env python3
"""
Create improved figure for GSE197004 Si stress context analysis.
"""
import gzip, numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path

ROOT = Path("C:/Users/moshe/Dropbox/ISF 2025")

# Load data
fpkm_path = ROOT / "RNA_seq_data/suppl_files/GSE197004/GSE197004_All_gene_fpkm.list.txt.gz"
with gzip.open(fpkm_path, 'rt') as f:
    fpkm = pd.read_csv(f, sep='\t', index_col=0)
log2_fpkm = np.log2(fpkm + 1)

contexts = ['Control', 'Drought', 'Salt', 'DroughtSalt']
context_labels = {
    'Control': 'Control',
    'Drought': 'Drought',
    'Salt': 'Salt',
    'DroughtSalt': 'Drought+Salt'
}

sample_map = {}
for tissue in ['A', 'U']:
    sample_map[tissue] = {
        'CK':   [c for c in fpkm.columns if c.startswith('CK-') and c.endswith('-' + tissue)],
        'CKSi': [c for c in fpkm.columns if c.startswith('CKSi-') and c.endswith('-' + tissue)],
        'D':    [c for c in fpkm.columns if c.startswith('D-') and c.endswith('-' + tissue)],
        'DSi':  [c for c in fpkm.columns if c.startswith('DSi-') and c.endswith('-' + tissue)],
        'S':    [c for c in fpkm.columns if c.startswith('S-') and c.endswith('-' + tissue)],
        'SSi':  [c for c in fpkm.columns if c.startswith('SSi-') and c.endswith('-' + tissue)],
        'SD':   [c for c in fpkm.columns if c.startswith('SD-') and c.endswith('-' + tissue)],
        'SDSi': [c for c in fpkm.columns if c.startswith('SDSi-') and c.endswith('-' + tissue)],
    }

means = {}
for tissue in ['A', 'U']:
    means[tissue] = {cond: log2_fpkm[cols].mean(axis=1)
                    for cond, cols in sample_map[tissue].items()}

si_effects = {}
for tissue in ['A', 'U']:
    si_effects[tissue] = {
        'Control':     means[tissue]['CKSi']  - means[tissue]['CK'],
        'Drought':     means[tissue]['DSi']   - means[tissue]['D'],
        'Salt':        means[tissue]['SSi']   - means[tissue]['S'],
        'DroughtSalt': means[tissue]['SDSi']  - means[tissue]['SD'],
    }

def cosine_sim(v1, v2):
    mask = ~(np.isnan(v1) | np.isnan(v2))
    v1m, v2m = v1[mask].values, v2[mask].values
    n1, n2 = np.linalg.norm(v1m), np.linalg.norm(v2m)
    if n1 == 0 or n2 == 0:
        return 0.0
    return float(np.dot(v1m, v2m) / (n1 * n2))

cosine_results = {}
for tissue in ['A', 'U']:
    cos_mat = pd.DataFrame(index=contexts, columns=contexts, dtype=float)
    for c1 in contexts:
        for c2 in contexts:
            cos_mat.loc[c1, c2] = cosine_sim(si_effects[tissue][c1], si_effects[tissue][c2])
    cosine_results[tissue] = cos_mat

si_stats = {}
for tissue in ['A', 'U']:
    si_stats[tissue] = {}
    for cond in contexts:
        vec = si_effects[tissue][cond]
        si_stats[tissue][cond] = {
            'n_up': int((vec > 0.5).sum()),
            'n_down': int((vec < -0.5).sum()),
            'magnitude': float(np.abs(vec).mean()),
        }

# ---- Figure ----
fig = plt.figure(figsize=(16, 9))
fig.patch.set_facecolor('white')
gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.5, wspace=0.4,
                       width_ratios=[1.4, 1, 1])

# Panel a: grouped bar chart of DEG counts
ax_a = fig.add_subplot(gs[:, 0])
x = np.arange(len(contexts))
width = 0.2
bar_sets = [
    ('Aerial up',   [si_stats['A'][c]['n_up'] for c in contexts],    '#1565C0'),
    ('Aerial down', [-si_stats['A'][c]['n_down'] for c in contexts], '#90CAF9'),
    ('Root up',     [si_stats['U'][c]['n_up'] for c in contexts],    '#2E7D32'),
    ('Root down',   [-si_stats['U'][c]['n_down'] for c in contexts], '#A5D6A7'),
]
offsets = [-0.3, -0.1, 0.1, 0.3]

for (label, vals, color), offset in zip(bar_sets, offsets):
    ax_a.bar(x + offset, vals, width, label=label, color=color,
             alpha=0.9, edgecolor='white', linewidth=0.5)

ax_a.axhline(y=0, color='black', linewidth=1.0, zorder=10)
ax_a.set_xticks(x)
ax_a.set_xticklabels([context_labels[c] for c in contexts], fontsize=10)
ax_a.set_ylabel('Number of Si-responsive genes (|log2FC| > 0.5)', fontsize=10)
ax_a.set_title('a. Si-responsive genes per stress context', fontweight='bold', fontsize=11)
ax_a.legend(fontsize=8, loc='upper left')

# Magnitude annotations
for i, cond in enumerate(contexts):
    mag_A = si_stats['A'][cond]['magnitude']
    mag_U = si_stats['U'][cond]['magnitude']
    ymax_val = max(si_stats['A'][cond]['n_up'], si_stats['U'][cond]['n_up'])
    ax_a.text(i, ymax_val + 120,
              'Shoot: %.3f\nRoot: %.3f' % (mag_A, mag_U),
              ha='center', va='bottom', fontsize=7, color='#444444',
              style='italic')

ax_a.set_ylim(-3500, 6200)
ax_a.spines['top'].set_visible(False)
ax_a.spines['right'].set_visible(False)
ax_a.text(0.02, 0.55, 'Up-regulated\nby Si', transform=ax_a.transAxes,
          fontsize=8, color='#555555', style='italic')
ax_a.text(0.02, 0.35, 'Down-regulated\nby Si', transform=ax_a.transAxes,
          fontsize=8, color='#555555', style='italic')

# Panel b1 and b2: cosine similarity matrices
for row, (tissue, tissue_name, panel_label) in enumerate([
    ('A', 'Aerial (shoot)', 'b.'),
    ('U', 'Root', 'c.')
]):
    ax_cos = fig.add_subplot(gs[row, 1:])
    cos_mat = cosine_results[tissue].astype(float)
    im = ax_cos.imshow(cos_mat.values, cmap='RdYlGn', vmin=-0.5, vmax=1.0, aspect='auto')
    ax_cos.set_xticks(range(len(contexts)))
    ax_cos.set_yticks(range(len(contexts)))
    short_labels = ['Control', 'Drought', 'Salt', 'D+S']
    ax_cos.set_xticklabels(short_labels, fontsize=9, rotation=30, ha='right')
    ax_cos.set_yticklabels(short_labels, fontsize=9)

    for i in range(len(contexts)):
        for j in range(len(contexts)):
            val = cos_mat.iloc[i, j]
            color = 'white' if abs(val) > 0.6 else 'black'
            bold = 'bold' if i == j else 'normal'
            ax_cos.text(j, i, '%.3f' % val,
                       ha='center', va='center', fontsize=10,
                       color=color, fontweight=bold)

    cbar = plt.colorbar(im, ax=ax_cos, shrink=0.85, pad=0.03)
    cbar.set_label('Cosine similarity', fontsize=8)
    ax_cos.set_title('%s Cosine similarity of Si-effect vectors: %s' % (panel_label, tissue_name),
                    fontweight='bold', fontsize=10)

    # Add interpretation note
    ctrl_drought = float(cos_mat.loc['Control', 'Drought'])
    ctrl_ds = float(cos_mat.loc['Control', 'DroughtSalt'])
    note = ('Si effect mostly stable\nunder single stress' if ctrl_drought > 0.1
            else 'Si effect changes\nwith stress context')
    ax_cos.text(1.25, 0.5, note, transform=ax_cos.transAxes,
               fontsize=7, ha='left', va='center', color='#333333',
               bbox=dict(boxstyle='round,pad=0.3', facecolor='#F5F5F5', alpha=0.8))

plt.suptitle(
    'GSE197004: Silicon biostimulant transcriptomics in G. uralensis\n'
    'under drought, salt, and combined stress (Ma et al. 2022 Front Plant Sci)',
    fontsize=12, fontweight='bold', y=1.01
)

fig_path = ROOT / "state_space_figures/GSE197004_Si_stress_context.png"
plt.savefig(fig_path, dpi=150, bbox_inches='tight', facecolor='white')
plt.close()
print("Figure saved: %s" % fig_path)

# Print final numbers
print("\nKey results:")
for tissue, tissue_name in [('A', 'Aerial'), ('U', 'Root')]:
    print("\n  %s:" % tissue_name)
    for cond in contexts:
        s = si_stats[tissue][cond]
        print("    %s: %d up, %d down, magnitude=%.4f" % (
            cond, s['n_up'], s['n_down'], s['magnitude']))
    print("  Cosine similarities (vs Control):")
    for c2 in ['Drought', 'Salt', 'DroughtSalt']:
        sim = float(cosine_results[tissue].loc['Control', c2])
        print("    Control vs %s: %.4f" % (c2, sim))
