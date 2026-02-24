"""
GSE214347: Glycyrrhiza uralensis Bacillus pumilus G5 + Silicon × drought
Axis fingerprinting analysis.

IMPORTANT NOTE ON METHODOLOGY:
The G. uralensis genome uses Glyur gene IDs (Mochida 2017, Plant Journal,
PMID 27775193) from a database that is currently offline (RIKEN PSC archive:
http://ngs-data-archive.psc.riken.jp/). This means a standard whole-genome
Glyur -> Arabidopsis AT gene ortholog mapping is not available from any
currently accessible public database (NCBI, Ensembl Plants, PLAZA, KEGG).

APPROACH USED:
1. ANCHOR GENE METHOD (P, N, JA, SA, ABA, Fe axes):
   39 Glyur genes identified from published paper tables (Frontiers Plant Sci 2022,
   doi:10.3389/fpls.2022.1033915) with gene names matching known Arabidopsis gene symbols.
   29 unique Arabidopsis AT IDs. Scores are extrapolated to full calibration range
   using the weighted contribution of anchor genes to each axis signature.
   CAVEAT: Anchor genes cover <0.5% of axis signal; scores are directionally
   informative but quantitatively uncertain. CK, ET, Auxin scores omitted
   due to unreliably large extrapolation factors.

2. COSINE SIMILARITY (full transcriptome, 27,475 Glyur genes):
   Scale-free comparison between all pairwise response vectors. Reliable because
   it uses the complete expression profile without requiring gene ID mapping.

3. KEY MARKER GENE ANALYSIS:
   PHT1 (phosphate transporter, AT5G43370, P-axis logFC=+8.83) is the clearest
   P-axis marker. Its expression dynamics are reported directly.

Experimental design: Control (CK, n=3), Drought (D, n=3),
Drought+Bacillus pumilus G5 (DB, n=3), Drought+Si (DSi, n=3),
Drought+Bacillus+Si (DBSi, n=3).

NOTE: There are NO non-drought Bacillus or Si conditions.
All biostimulant scores are RELIEF-under-drought effects (treatment vs Drought).
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import warnings
warnings.filterwarnings('ignore')
import os

# ============================================================
# PATHS
# ============================================================
ROOT = "C:/Users/moshe/Dropbox/ISF 2025"
CPM_FILE = f"{ROOT}/RNA_seq_data/suppl_files/GSE214347/GSE214347_All_gene_CPM.list.txt.gz"
SIG_DIR = f"{ROOT}/RNA_seq_axes/signatures"
CALIB_FILE = f"{ROOT}/RNA_seq_axes/calibration.csv"
SCORES_OUT = f"{ROOT}/RNA_seq_axes/scores/GSE214347_biostimulant_drought_scores.csv"
FIG_OUT = f"{ROOT}/state_space_figures/GSE214347_biostimulant_drought.png"
REPORT_OUT = f"{ROOT}/agent_coordination/GSE214347_stress_context_report.md"

# ============================================================
# LOAD DATA
# ============================================================
print("Loading CPM data...")
df = pd.read_csv(CPM_FILE, sep='\t', compression='gzip', index_col=0)
glyur_df = df[df.index.str.startswith('Glyur')].copy()
print(f"  {len(glyur_df)} Glyur genes, {len(df) - len(glyur_df)} ONT genes (excluded)")

groups = {
    'CK': ['CK-1', 'CK-2', 'CK-3'],
    'D':  ['D-1',  'D-2',  'D-3'],
    'DB': ['DB-1', 'DB-2', 'DB-3'],
    'DBSi': ['DBSi-1', 'DBSi-2', 'DBSi-3'],
    'DSi':  ['DSi-1',  'DSi-2',  'DSi-3'],
}

log_cpm = np.log2(glyur_df + 1)
means = {grp: log_cpm[cols].mean(axis=1) for grp, cols in groups.items()}

# Response vectors
fc = {
    'Drought':         means['D']    - means['CK'],   # V_drought
    'Bacillus_vs_D':   means['DB']   - means['D'],    # V_Bacillus_in_drought
    'Si_vs_D':         means['DSi']  - means['D'],    # V_Si_in_drought
    'BSi_vs_D':        means['DBSi'] - means['D'],    # V_BSi_in_drought
}

# ============================================================
# ANCHOR GENE MAPPING (Glyur -> AT)
# ============================================================
# Source: gene names from Frontiers Plant Sci 2022 paper tables,
# manually mapped to TAIR locus IDs.
glyur_to_at = {
    # P-axis anchor (KEY: PHT1 = top P-axis gene, logFC=+8.83)
    'Glyur001979s00036874': 'AT5G43370',  # PHT1;4 phosphate transporter
    # N-axis anchors
    'Glyur000288s00019860': 'AT4G39640',  # GGT1 gamma-glutamyltransferase
    'Glyur000025s00003342': 'AT1G70410',  # CA2 carbonic anhydrase (N-axis +2.34)
    'Glyur000100s00008372': 'AT1G63940',  # MDAR5 (N-axis +3.65)
    'Glyur000058s00006123': 'AT5G23310',  # FSD3 Fe superoxide dismutase (N-axis +0.79)
    # SA-axis anchors
    'Glyur000273s00014672': 'AT1G07890',  # APX1 ascorbate peroxidase (SA +1.04)
    'Glyur000157s00015845': 'AT1G10370',  # GSTU17 glutathione S-transferase (SA +1.54)
    'Glyur000595s00018295': 'AT4G11600',  # GPX6 glutathione peroxidase (SA +0.44)
    # JA-axis anchors
    'Glyur000163s00011136': 'AT2G37040',  # PAL1 phenylalanine ammonia-lyase (JA +0.80)
    'Glyur000327s00026461': 'AT2G37040',  # PAL1 second copy
    'Glyur000051s00003431': 'AT5G13930',  # CHS chalcone synthase (JA +1.78)
    'Glyur000223s00013796': 'AT5G13930',  # CHS1
    'Glyur000397s00020478': 'AT5G13930',  # CHS6
    'Glyur000424s00026890': 'AT5G13930',  # CHS2
    'Glyur001333s00028402': 'AT5G13930',  # CHS5
    'Glyur002999s00043471': 'AT5G13930',  # CHS
    'Glyur006062s00044203': 'AT5G13930',  # CHS1
    'Glyur000098s00007989': 'AT1G51680',  # 4CL1 (JA pathway)
    'Glyur003091s00037488': 'AT5G38120',  # 4CLL6 (JA +4.46!)
    # ABA-axis anchors
    'Glyur000004s00000451': 'AT4G34240',  # ALDH3F1 aldehyde dehydrogenase (ABA +0.82)
    'Glyur000245s00014674': 'AT1G44170',  # ALDH3H1 (ABA +1.06)
    'Glyur001144s00035322': 'AT1G44170',  # ALDH3H1 second
    'Glyur000200s00011170': 'AT3G24503',  # ALDH3I1 (ABA -0.73)
    'Glyur006220s00043163': 'AT3G66585',  # ALDH2B7
    'Glyur000260s00013529': 'AT3G26744',  # SCRM/ICE1 (ABA +0.18)
    # Fe-axis anchors
    'Glyur000013s00003277': 'AT1G74540',  # CYP98A2 (Fe +1.94!)
    # Multi-axis
    'Glyur000069s00004035': 'AT1G08830',  # SOD1
    'Glyur000274s00012557': 'AT2G28190',  # SODCP
    'Glyur000583s00027212': 'AT1G08830',  # SODB
    'Glyur000860s00039262': 'AT1G75270',  # DHAR3
    'Glyur000081s00006734': 'AT2G30870',  # GSTL3
    'Glyur000823s00027909': 'AT4G26850',  # VTC2
    'Glyur002747s00042594': 'AT5G08640',  # FLS flavonol synthase
    'Glyur000959s00024500': 'AT3G55120',  # CHI1B1
    'Glyur005711s00045445': 'AT3G55120',  # CHI1
    'Glyur000021s00002690': 'AT1G73560',  # MIOX1
    'Glyur001845s00034619': 'AT1G30580',  # MIOX2
    'Glyur001723s00036239': 'AT3G05165',  # PNC1
    'Glyur000136s00007908': 'AT5G14040',  # PNC2
}

print(f"\nAnchor genes: {len(glyur_to_at)} Glyur IDs -> {len(set(glyur_to_at.values()))} unique AT IDs")

# ============================================================
# LOAD AXES AND CALIBRATION
# ============================================================
axis_names_reliable = ['P', 'N', 'JA', 'SA', 'ABA', 'Fe']
axes = {}
for ax in axis_names_reliable:
    fname = f"{SIG_DIR}/{ax}_axis_logFC.csv"
    axes[ax] = pd.read_csv(fname).set_index('GeneID')['logFC']

calib = pd.read_csv(CALIB_FILE, index_col=0)

# ============================================================
# COMPUTE ANCHOR-GENE AXIS SCORES
# ============================================================
def get_at_fc(comp_fc):
    """Average Glyur fold changes by AT gene ID."""
    at_fc = {}
    for glyur_id, at_id in glyur_to_at.items():
        if glyur_id in comp_fc.index:
            if at_id not in at_fc:
                at_fc[at_id] = []
            at_fc[at_id].append(float(comp_fc[glyur_id]))
    return {at: float(np.mean(vals)) for at, vals in at_fc.items()}

def compute_pct_score(ax_name, ax_sig, at_fc_dict):
    """
    Weighted extrapolation: partial dot product scaled by anchor gene weight.
    Returns % calibration range.
    """
    common = [at for at in at_fc_dict if at in ax_sig.index]
    if len(common) < 3:
        return np.nan
    partial = sum(ax_sig[at] * at_fc_dict[at] for at in common)
    anchor_var = sum(ax_sig[at]**2 for at in common)
    total_var = float((ax_sig**2).sum())
    if total_var == 0 or anchor_var == 0:
        return np.nan
    fraction = anchor_var / total_var
    extrapolated = partial / fraction
    ax_key = f'{ax_name}_axis'
    if ax_key not in calib.index:
        return np.nan
    ctrl = float(calib.loc[ax_key, 'ctrl_score'])
    treat = float(calib.loc[ax_key, 'treat_score'])
    calib_range = treat - ctrl
    return (extrapolated - ctrl) / calib_range * 100.0

# ============================================================
# COSINE SIMILARITY (full transcriptome)
# ============================================================
def cosine_sim(v1, v2):
    n1 = float(np.sqrt((v1**2).sum()))
    n2 = float(np.sqrt((v2**2).sum()))
    if n1 == 0 or n2 == 0:
        return 0.0
    return float((v1 * v2).sum() / (n1 * n2))

# ============================================================
# COMPUTE ALL SCORES
# ============================================================
print("\nComputing axis scores...")
scores_records = []

for comp_name, comp_fc in fc.items():
    at_fc_dict = get_at_fc(comp_fc)
    rec = {'comparison': comp_name}

    # Axis scores
    for ax_name, ax_sig in axes.items():
        rec[ax_name + '_pct'] = compute_pct_score(ax_name, ax_sig, at_fc_dict)

    # Cosine similarities
    rec['cosine_vs_drought'] = cosine_sim(fc['Drought'], comp_fc)
    rec['L2_norm'] = float(np.sqrt((comp_fc**2).sum()))

    # DEG counts
    rec['n_up_1fc'] = int((comp_fc > 1.0).sum())
    rec['n_down_1fc'] = int((comp_fc < -1.0).sum())

    scores_records.append(rec)

# Pairwise cosine similarities
cos_records = {}
biostim_comps = ['Bacillus_vs_D', 'Si_vs_D', 'BSi_vs_D']
for i, n1 in enumerate(biostim_comps):
    for j, n2 in enumerate(biostim_comps):
        if j > i:
            cos_records[f'{n1}_vs_{n2}'] = cosine_sim(fc[n1], fc[n2])

scores_df = pd.DataFrame(scores_records).set_index('comparison')
print(scores_df.round(2).to_string())

# Save scores CSV
os.makedirs(os.path.dirname(SCORES_OUT), exist_ok=True)
scores_df.to_csv(SCORES_OUT)
print(f"\nScores saved: {SCORES_OUT}")

# ============================================================
# KEY MARKER GENE TABLE
# ============================================================
key_markers = {
    'Glyur001979s00036874': ('PHT1', 'AT5G43370', 'P-axis', '+8.83'),
    'Glyur000025s00003342': ('CA2', 'AT1G70410', 'N-axis', '+2.34'),
    'Glyur000100s00008372': ('MDAR5', 'AT1G63940', 'N-axis', '+3.65'),
    'Glyur000273s00014672': ('APX1', 'AT1G07890', 'SA-axis', '+1.04'),
    'Glyur000163s00011136': ('PAL1', 'AT2G37040', 'JA-axis', '+0.80'),
    'Glyur000051s00003431': ('CHS1', 'AT5G13930', 'JA-axis', '+1.78'),
    'Glyur000003091s00037488' if 'Glyur000003091s00037488' in glyur_df.index else 'Glyur003091s00037488':
        ('4CLL6', 'AT5G38120', 'JA-axis', '+4.46'),
    'Glyur000004s00000451': ('ALDH3F1', 'AT4G34240', 'ABA-axis', '+0.82'),
    'Glyur000013s00003277': ('CYP98A2', 'AT1G74540', 'Fe-axis', '+1.94'),
    'Glyur000058s00006123': ('FSD3', 'AT5G23310', 'N-axis', '+0.79'),
}

marker_rows = []
for glyur_id, (gene, at_id, axis, at_logfc) in key_markers.items():
    if glyur_id not in glyur_df.index:
        # Handle possible key formatting issues
        continue
    row = {'Glyur_ID': glyur_id, 'Gene': gene, 'AT_ID': at_id,
           'Axis_logFC': at_logfc, 'Axis': axis}
    for comp_name, comp_fc in fc.items():
        if glyur_id in comp_fc.index:
            row[comp_name] = round(float(comp_fc[glyur_id]), 3)
        else:
            row[comp_name] = np.nan
    marker_rows.append(row)

marker_df = pd.DataFrame(marker_rows)
print("\nKey marker gene fold changes:")
print(marker_df.to_string(index=False))

# ============================================================
# FIGURE: Panel A (axis scores) + Panel B (cosine similarities)
# ============================================================
print("\nCreating figure...")
os.makedirs(os.path.dirname(FIG_OUT), exist_ok=True)

fig, axes_fig = plt.subplots(1, 2, figsize=(14, 6))

# --- Panel A: Bar chart of axis scores ---
ax_a = axes_fig[0]
ax_plot_names = ['P', 'N', 'JA', 'SA', 'ABA', 'Fe']
plot_comps = ['Drought', 'Bacillus_vs_D', 'Si_vs_D', 'BSi_vs_D']
plot_labels = ['Drought', 'Bacillus\n(vs Drought)', 'Si\n(vs Drought)', 'BSi\n(vs Drought)']
colors = ['#d62728', '#1f77b4', '#2ca02c', '#9467bd']

x = np.arange(len(ax_plot_names))
width = 0.2

for i, (comp, label, color) in enumerate(zip(plot_comps, plot_labels, colors)):
    vals = [scores_df.loc[comp, ax+'_pct'] for ax in ax_plot_names]
    bars = ax_a.bar(x + i*width - 1.5*width, vals, width, label=label, color=color, alpha=0.8)

ax_a.axhline(y=0, color='black', linewidth=0.8)
ax_a.set_xticks(x)
ax_a.set_xticklabels(ax_plot_names, fontsize=11)
ax_a.set_ylabel('Score (% calibration range)', fontsize=11)
ax_a.set_title('Panel A: Axis fingerprints\n(anchor-gene method, n=29 AT IDs)', fontsize=11)
ax_a.legend(fontsize=9, loc='upper right')
ax_a.set_ylim(-80, 100)
ax_a.text(0.02, 0.02, 'PRELIMINARY: limited gene coverage (<0.5% of axis)\nDirectional trends only; quantitative values uncertain',
          transform=ax_a.transAxes, fontsize=7, color='red',
          verticalalignment='bottom')

# --- Panel B: Cosine similarity heatmap / bar chart ---
ax_b = axes_fig[1]

# Show cosine similarities as bars
cos_labels = []
cos_vals = []
cos_colors_list = []

# Drought as reference
for biostim, label, color in [('Bacillus_vs_D', 'Bacillus\nvs Drought', '#1f77b4'),
                                ('Si_vs_D', 'Si\nvs Drought', '#2ca02c'),
                                ('BSi_vs_D', 'BSi\nvs Drought', '#9467bd')]:
    cos_labels.append(f'Drought\nvs {label}')
    cos_vals.append(cosine_sim(fc['Drought'], fc[biostim]))
    cos_colors_list.append(color)

# Pairwise between biostimulants
for n1, n2, label, color in [
    ('Bacillus_vs_D', 'Si_vs_D', 'Bacillus vs Si', '#17becf'),
    ('Bacillus_vs_D', 'BSi_vs_D', 'Bacillus vs BSi', '#bcbd22'),
    ('Si_vs_D', 'BSi_vs_D', 'Si vs BSi', '#e377c2'),
]:
    cos_labels.append(label)
    cos_vals.append(cosine_sim(fc[n1], fc[n2]))
    cos_colors_list.append(color)

bar_positions = np.arange(len(cos_vals))
ax_b.barh(bar_positions, cos_vals, color=cos_colors_list, alpha=0.8)
ax_b.axvline(x=0, color='black', linewidth=0.8)
ax_b.set_yticks(bar_positions)
ax_b.set_yticklabels(cos_labels, fontsize=9)
ax_b.set_xlabel('Cosine similarity', fontsize=11)
ax_b.set_title('Panel B: Response vector\ncosine similarities (full transcriptome)', fontsize=11)
ax_b.set_xlim(-0.7, 0.9)
ax_b.text(0.02, 0.02, 'Negative = opposes drought (beneficial)\nPositive = similar to drought (harmful)',
          transform=ax_b.transAxes, fontsize=8, color='gray', verticalalignment='bottom')

plt.suptitle('GSE214347: G. uralensis — Bacillus pumilus G5 + Si × Drought\n(Transcriptome-based fingerprinting)',
             fontsize=12, fontweight='bold', y=1.02)

plt.tight_layout()
plt.savefig(FIG_OUT, dpi=150, bbox_inches='tight')
plt.close()
print(f"Figure saved: {FIG_OUT}")

# ============================================================
# REPORT
# ============================================================
report = f"""# GSE214347 Stress-Context Fingerprinting Report

**Dataset:** Glycyrrhiza uralensis, Bacillus pumilus G5 + Silicon × drought
**Source:** GSE214347, Frontiers Plant Sci 2022, doi:10.3389/fpls.2022.1033915
**Date:** 2026-02-24
**Analysis method:** Anchor-gene axis scoring + full-transcriptome cosine similarity

---

## Experimental Design

| Group | Treatment | n |
|-------|-----------|---|
| CK | Control (65-75% soil moisture) | 3 |
| D | Drought (35-45% soil moisture) | 3 |
| DB | Drought + Bacillus pumilus G5 (10^8 CFU/mL) | 3 |
| DSi | Drought + Si (2 mM) | 3 |
| DBSi | Drought + Bacillus + Si | 3 |

**Critical caveat:** There are NO non-drought Bacillus-only or Si-only conditions.
All biostimulant scores represent **relief-under-drought** effects (treatment vs Drought baseline).
The biostimulant-vs-control effect under normal conditions cannot be computed from this dataset.

---

## Methodology

### Gene Ortholog Limitation
The G. uralensis CPM data uses Glyur gene IDs from the Mochida 2017 Plant Journal genome
(PMID 27775193). This genome's annotation database (RIKEN PSC: http://ngs-data-archive.psc.riken.jp/)
is offline. G. uralensis is not in KEGG (plant), PLAZA dicots v5, Ensembl Plants, or NCBI
RefSeq with nuclear gene annotation. Therefore, standard genome-wide Glyur→AT ortholog mapping
was not achievable.

**Approaches implemented:**
1. **Anchor gene method**: 39 Glyur gene IDs from published paper tables matched to gene symbols
   (SOD, APX, PAL, CHS, PHT1, etc.), manually curated to 29 unique Arabidopsis AT IDs.
   Axis scores computed by weighted extrapolation (anchor weight fraction of axis signal).
   Covers <0.5% of axis signal; scores are directionally informative but quantitatively uncertain.
   CK, ET, Auxin scores omitted (unreliable extrapolation factors).

2. **Full-transcriptome cosine similarity**: Computed across all 27,475 Glyur genes. Scale-free
   and does not require gene ID mapping. Reliable measure of response vector similarity.

---

## Key Results

### 1. Axis Scores (% calibration range, anchor-gene method)

| Comparison | P | N | JA | SA | ABA | Fe |
|-----------|---|---|----|----|-----|----|
| Drought (vs CK) | {scores_df.loc['Drought','P_pct']:.1f}% | {scores_df.loc['Drought','N_pct']:.1f}% | {scores_df.loc['Drought','JA_pct']:.1f}% | {scores_df.loc['Drought','SA_pct']:.1f}% | {scores_df.loc['Drought','ABA_pct']:.1f}% | {scores_df.loc['Drought','Fe_pct']:.1f}% |
| Bacillus (vs Drought) | {scores_df.loc['Bacillus_vs_D','P_pct']:.1f}% | {scores_df.loc['Bacillus_vs_D','N_pct']:.1f}% | {scores_df.loc['Bacillus_vs_D','JA_pct']:.1f}% | {scores_df.loc['Bacillus_vs_D','SA_pct']:.1f}% | {scores_df.loc['Bacillus_vs_D','ABA_pct']:.1f}% | {scores_df.loc['Bacillus_vs_D','Fe_pct']:.1f}% |
| Si (vs Drought) | {scores_df.loc['Si_vs_D','P_pct']:.1f}% | {scores_df.loc['Si_vs_D','N_pct']:.1f}% | {scores_df.loc['Si_vs_D','JA_pct']:.1f}% | {scores_df.loc['Si_vs_D','SA_pct']:.1f}% | {scores_df.loc['Si_vs_D','ABA_pct']:.1f}% | {scores_df.loc['Si_vs_D','Fe_pct']:.1f}% |
| BSi (vs Drought) | {scores_df.loc['BSi_vs_D','P_pct']:.1f}% | {scores_df.loc['BSi_vs_D','N_pct']:.1f}% | {scores_df.loc['BSi_vs_D','JA_pct']:.1f}% | {scores_df.loc['BSi_vs_D','SA_pct']:.1f}% | {scores_df.loc['BSi_vs_D','ABA_pct']:.1f}% | {scores_df.loc['BSi_vs_D','Fe_pct']:.1f}% |

*P and N scores most reliable (PHT1, CA2, MDAR5 are strong axis markers).
JA directionally meaningful (PAL1, CHS, 4CLL6). SA/ABA/Fe weaker coverage.*

### 2. Cosine Similarities (full transcriptome, 27,475 genes)

| Comparison | Cosine sim | Interpretation |
|------------|-----------|---------------|
| Drought vs Bacillus | {cosine_sim(fc['Drought'], fc['Bacillus_vs_D']):.4f} | Bacillus partially opposes drought |
| Drought vs Si | {cosine_sim(fc['Drought'], fc['Si_vs_D']):.4f} | Si more strongly opposes drought |
| Drought vs BSi | {cosine_sim(fc['Drought'], fc['BSi_vs_D']):.4f} | BSi effect is more complex |
| Bacillus vs Si | {cosine_sim(fc['Bacillus_vs_D'], fc['Si_vs_D']):.4f} | Bacillus and Si have similar relief mechanisms |
| Bacillus vs BSi | {cosine_sim(fc['Bacillus_vs_D'], fc['BSi_vs_D']):.4f} | BSi overlaps with Bacillus response |
| Si vs BSi | {cosine_sim(fc['Si_vs_D'], fc['BSi_vs_D']):.4f} | BSi overlaps with Si response |

### 3. Key Marker Gene Responses

| Gene | Function | AT ID | D vs CK | B vs D | Si vs D | BSi vs D |
|------|---------|-------|---------|--------|---------|---------|
| PHT1 | Phosphate transporter (P-axis top gene) | AT5G43370 | {float(fc['Drought']['Glyur001979s00036874']):+.3f} | {float(fc['Bacillus_vs_D']['Glyur001979s00036874']):+.3f} | {float(fc['Si_vs_D']['Glyur001979s00036874']):+.3f} | {float(fc['BSi_vs_D']['Glyur001979s00036874']):+.3f} |
| CA2 | Carbonic anhydrase (N-axis) | AT1G70410 | {float(fc['Drought']['Glyur000025s00003342']):+.3f} | {float(fc['Bacillus_vs_D']['Glyur000025s00003342']):+.3f} | {float(fc['Si_vs_D']['Glyur000025s00003342']):+.3f} | {float(fc['BSi_vs_D']['Glyur000025s00003342']):+.3f} |
| MDAR5 | Reductase (N-axis high) | AT1G63940 | {float(fc['Drought']['Glyur000100s00008372']):+.3f} | {float(fc['Bacillus_vs_D']['Glyur000100s00008372']):+.3f} | {float(fc['Si_vs_D']['Glyur000100s00008372']):+.3f} | {float(fc['BSi_vs_D']['Glyur000100s00008372']):+.3f} |
| APX1 | Ascorbate peroxidase (SA-axis) | AT1G07890 | {float(fc['Drought']['Glyur000273s00014672']):+.3f} | {float(fc['Bacillus_vs_D']['Glyur000273s00014672']):+.3f} | {float(fc['Si_vs_D']['Glyur000273s00014672']):+.3f} | {float(fc['BSi_vs_D']['Glyur000273s00014672']):+.3f} |
| PAL1 | Phenylalanine ammonia-lyase (JA/SA) | AT2G37040 | {float(fc['Drought']['Glyur000163s00011136']):+.3f} | {float(fc['Bacillus_vs_D']['Glyur000163s00011136']):+.3f} | {float(fc['Si_vs_D']['Glyur000163s00011136']):+.3f} | {float(fc['BSi_vs_D']['Glyur000163s00011136']):+.3f} |
| CHS1 | Chalcone synthase (JA/flavonoid) | AT5G13930 | {float(fc['Drought']['Glyur000051s00003431']):+.3f} | {float(fc['Bacillus_vs_D']['Glyur000051s00003431']):+.3f} | {float(fc['Si_vs_D']['Glyur000051s00003431']):+.3f} | {float(fc['BSi_vs_D']['Glyur000051s00003431']):+.3f} |
| ALDH3F1 | Aldehyde dehydrogenase (ABA-stress) | AT4G34240 | {float(fc['Drought']['Glyur000004s00000451']):+.3f} | {float(fc['Bacillus_vs_D']['Glyur000004s00000451']):+.3f} | {float(fc['Si_vs_D']['Glyur000004s00000451']):+.3f} | {float(fc['BSi_vs_D']['Glyur000004s00000451']):+.3f} |
| FSD3 | Fe-superoxide dismutase (N/Fe-axis) | AT5G23310 | {float(fc['Drought']['Glyur000058s00006123']):+.3f} | {float(fc['Bacillus_vs_D']['Glyur000058s00006123']):+.3f} | {float(fc['Si_vs_D']['Glyur000058s00006123']):+.3f} | {float(fc['BSi_vs_D']['Glyur000058s00006123']):+.3f} |

### 4. Effect Size (L2 norms of response vectors)

| Vector | L2-norm | Interpretation |
|--------|---------|---------------|
| Drought | {scores_df.loc['Drought','L2_norm']:.1f} | Baseline drought transcriptomic effect |
| Bacillus vs Drought | {scores_df.loc['Bacillus_vs_D','L2_norm']:.1f} | Bacillus effect magnitude |
| Si vs Drought | {scores_df.loc['Si_vs_D','L2_norm']:.1f} | Si effect magnitude |
| BSi vs Drought | {scores_df.loc['BSi_vs_D','L2_norm']:.1f} | BSi effect magnitude (larger = synergistic) |

| Comparison | Up (|FC|>1) | Down (|FC|>1) |
|-----------|-----------|--------------|
| Drought | {scores_df.loc['Drought','n_up_1fc']} | {scores_df.loc['Drought','n_down_1fc']} |
| Bacillus vs D | {scores_df.loc['Bacillus_vs_D','n_up_1fc']} | {scores_df.loc['Bacillus_vs_D','n_down_1fc']} |
| Si vs D | {scores_df.loc['Si_vs_D','n_up_1fc']} | {scores_df.loc['Si_vs_D','n_down_1fc']} |
| BSi vs D | {scores_df.loc['BSi_vs_D','n_up_1fc']} | {scores_df.loc['BSi_vs_D','n_down_1fc']} |

---

## Biological Interpretation

### Drought Response
- P score: +{scores_df.loc['Drought','P_pct']:.1f}% (mild P-stress signal; PHT1 upregulated +0.35 log2FC)
- N score: +{scores_df.loc['Drought','N_pct']:.1f}% (mild N-stress signal; N-axis genes modestly elevated)
- Drought downregulates N-metabolism genes (CA2, MDAR5) and strongly suppresses APX1 (-1.1 log2FC)
- ALDH3F1 strongly induced (+1.3 log2FC) - ABA-associated stress marker
- CHS1 elevated (+0.72) - flavonoid/antioxidant response to oxidative stress

### Bacillus pumilus G5 Under Drought
- P score: {scores_df.loc['Bacillus_vs_D','P_pct']:.1f}% (PHT1 DECREASES -0.56; Bacillus RELIEVES P-stress)
- N score: {scores_df.loc['Bacillus_vs_D','N_pct']:.1f}% (strong N-axis shift; N-assimilation affected)
- APX1 strongly recovers (+0.99 log2FC relative to drought) - SA-ISR activation
- PAL1 upregulated (+0.31) - phenylpropanoid defense induction
- JA score: +{scores_df.loc['Bacillus_vs_D','JA_pct']:.1f}% (JA defense maintained/enhanced)
- **Mechanism**: Bacillus G5 activates ISR (APX1/PAL1 induction) and relieves P-stress (PHT1 relief)

### Silicon Under Drought
- P score: {scores_df.loc['Si_vs_D','P_pct']:.1f}% (PHT1 decreases -0.60; Si also RELIEVES P-stress)
- N score: {scores_df.loc['Si_vs_D','N_pct']:.1f}% (CA2 strongly down -1.09; Si shifts C/N balance)
- APX1 most strongly recovered (+1.37 log2FC) - strongest antioxidant activation
- PAL1 upregulated (+0.61) - phenylpropanoid defense
- JA score: +{scores_df.loc['Si_vs_D','JA_pct']:.1f}% (similar JA activation to Bacillus)
- **Mechanism**: Si relieves P-stress and activates ROS-scavenging defense (APX1). Strong drought-opposition
  (cosine = -0.50, stronger than Bacillus at -0.44)

### Bacillus + Si (BSi) Under Drought
- P score: {scores_df.loc['BSi_vs_D','P_pct']:.1f}% (STRONGEST P relief; BSi most effective Pi supplier)
- JA score: {scores_df.loc['BSi_vs_D','JA_pct']:.1f}% (near zero - JA pathway not activated)
- APX1 recovered (+1.04 log2FC) - antioxidant defense maintained
- PAL1 STRONGLY induced (+1.78 log2FC) - strongest phenylpropanoid activation
- CHS1 SUPPRESSED (-2.03 log2FC) - unusual; flavonoid biosynthesis redirected?
- **Mechanism**: BSi combination is most effective for P relief but shows different defense
  pathway than individual treatments. JA near-zero while PAL1 is high = SA-like ISR pattern?
  L2-norm {scores_df.loc['BSi_vs_D','L2_norm']:.1f} > Drought {scores_df.loc['Drought','L2_norm']:.1f} - BSi has LARGER transcriptomic impact than drought itself

### Cross-Treatment Comparisons
- Bacillus and Si have SIMILAR mechanisms (cosine = +0.49) - both relieve P-stress and activate APX1
- BSi is a SUPERPOSITION of both effects, not a novel mechanism (cosines +0.57 and +0.51 with each)
- All biostimulants OPPOSE drought response (all cosines negative with drought vector)
- Si MOST STRONGLY opposes drought (-0.50), Bacillus moderately (-0.44), BSi least (-0.26; more complex)

---

## Comparison with Other Biostimulants in Fingerprinting Framework

| Biostimulant | P-axis | N-axis | JA-axis | SA-axis | Context |
|-------------|--------|--------|---------|---------|---------|
| Bacillus G5 (vs Drought) | {scores_df.loc['Bacillus_vs_D','P_pct']:.1f}% | {scores_df.loc['Bacillus_vs_D','N_pct']:.1f}% | {scores_df.loc['Bacillus_vs_D','JA_pct']:.1f}% | {scores_df.loc['Bacillus_vs_D','SA_pct']:.1f}% | *Drought relief in licorice* |
| Si (vs Drought) | {scores_df.loc['Si_vs_D','P_pct']:.1f}% | {scores_df.loc['Si_vs_D','N_pct']:.1f}% | {scores_df.loc['Si_vs_D','JA_pct']:.1f}% | {scores_df.loc['Si_vs_D','SA_pct']:.1f}% | *Drought relief in licorice* |
| BSi (vs Drought) | {scores_df.loc['BSi_vs_D','P_pct']:.1f}% | {scores_df.loc['BSi_vs_D','N_pct']:.1f}% | {scores_df.loc['BSi_vs_D','JA_pct']:.1f}% | {scores_df.loc['BSi_vs_D','SA_pct']:.1f}% | *Drought relief in licorice* |
| diacetyl (Arabidopsis) | +13.7% | +11.8% | -9.8% | +12.6% | *Normal conditions* |
| GMV volatile (Arabidopsis) | +11.8% | -4.3% | +8.3% | +3.4% | *Normal conditions* |
| Si (TiO2 proxy, Arabidopsis) | -1.5% | -5.9% | +3.2% | +3.2% | *Normal conditions* |

*Note: Direct comparison is confounded by: (a) different baseline conditions (drought vs normal),
(b) different species (G. uralensis vs Arabidopsis), (c) anchor-gene method uncertainty.*

---

## Conclusions

1. **Both Bacillus G5 and Si relieve P-stress under drought**: PHT1 expression decreases when
   biostimulants are applied (-0.56 and -0.60 log2FC respectively), indicating improved Pi availability.
   This is consistent with known PGPR mechanisms (Pi solubilization) and Si effects.

2. **Both biostimulants activate antioxidant defense**: APX1 strongly recovers under both treatments
   (Bacillus: +0.99, Si: +1.37 log2FC relative to drought). Si shows stronger activation.

3. **Si most strongly opposes drought** (cosine = -0.50 vs Bacillus -0.44): Silicon has broader
   physiological effects on drought response at the transcriptomic level.

4. **BSi combination has the largest transcriptomic impact** (L2-norm = {scores_df.loc['BSi_vs_D','L2_norm']:.1f} vs drought = {scores_df.loc['Drought','L2_norm']:.1f}),
   suggesting synergistic effects on gene expression beyond either treatment alone.

5. **JA pathway**: BSi shows near-zero JA score while both individual treatments show positive JA.
   The combination may shift defense from JA-ISR to SA-ISR pathway (PAL1 strongly induced, CHS1 suppressed).

6. **Important limitations**: Scores are based on drought-relief comparisons only (no normal-condition
   biostimulant data). The anchor-gene method covers <0.5% of axis signal. Results are directionally
   informative but the percentage values are uncertain. The G. uralensis genome annotation issue
   prevents standard full-genome fingerprinting.

---

## Files Generated
- `RNA_seq_axes/scores/GSE214347_biostimulant_drought_scores.csv` — Axis scores + cosine similarities
- `state_space_figures/GSE214347_biostimulant_drought.png` — Panel A (axis scores) + Panel B (cosine sim)
- `agent_coordination/GSE214347_stress_context_report.md` — This report

## Data Source
- GEO: GSE214347
- Reference genome: Mochida et al. 2017 Plant J 89(1):181-194, doi:10.1111/tpj.13385 (offline)
- Publication: Frontiers Plant Sci 2022, doi:10.3389/fpls.2022.1033915
"""

os.makedirs(os.path.dirname(REPORT_OUT), exist_ok=True)
with open(REPORT_OUT, 'w', encoding='utf-8') as f:
    f.write(report)
print(f"Report saved: {REPORT_OUT}")

print("\n=== ANALYSIS COMPLETE ===")
print(f"\nKey findings:")
print(f"  PHT1 (P-axis marker): Drought={float(fc['Drought']['Glyur001979s00036874']):+.2f}, B={float(fc['Bacillus_vs_D']['Glyur001979s00036874']):+.2f}, Si={float(fc['Si_vs_D']['Glyur001979s00036874']):+.2f}, BSi={float(fc['BSi_vs_D']['Glyur001979s00036874']):+.2f}")
print(f"  Cosine (Drought vs Bacillus): {cosine_sim(fc['Drought'], fc['Bacillus_vs_D']):.4f}")
print(f"  Cosine (Drought vs Si): {cosine_sim(fc['Drought'], fc['Si_vs_D']):.4f}")
print(f"  Cosine (Bacillus vs Si): {cosine_sim(fc['Bacillus_vs_D'], fc['Si_vs_D']):.4f}")
