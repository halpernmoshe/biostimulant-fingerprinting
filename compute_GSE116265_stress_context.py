"""
GSE116265 Stress-context cosine similarity analysis
Cucumber silicon × salt (NaCl) dataset

4 conditions: Control (CT), Si-alone, NaCl-alone (Na), Si+NaCl (NaSi)
  CT:   3 replicates
  Na:   3 replicates
  Si:   2 replicates  ← NOTE: only 2 replicates (caveat flagged in report)
  NaSi: 3 replicates

Fingerprint vectors:
  V_Si         = Si_alone   - Control     (pure Si effect)
  V_NaCl       = NaCl_alone - Control     (pure salt effect)
  V_Si_salt    = Si+NaCl    - NaCl_alone  (Si effect within salt stress)
  V_interaction= Si+NaCl    - NaCl_alone - Si_alone + Control  (interaction term)

Axes used: SA, JA, P, N, ABA, Auxin, Fe, ET (Arabidopsis Arabidopsis-axis signatures,
           projected through Arabidopsis-cucumber orthologs already embedded in
           Cucumber_Si_detailed_scores.csv columns *_sig)

Cosine similarity: cos(V_Si, V_Si_salt)
Rotation angle   = arccos(cosine similarity) in degrees
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.stats import pearsonr
import warnings
warnings.filterwarnings('ignore')

# ── paths ──────────────────────────────────────────────────────────────────
BASE     = "C:/Users/moshe/Dropbox/ISF 2025"
SCORES   = f"{BASE}/RNA_seq_axes/scores"
SIGS_DIR = f"{BASE}/RNA_seq_axes/signatures"
FIGS_DIR = f"{BASE}/state_space_figures"
AGT_DIR  = f"{BASE}/agent_coordination"

INPUT_CSV = f"{SCORES}/Cucumber_Si_detailed_scores.csv"
CAL_CSV   = f"{BASE}/RNA_seq_axes/calibration.csv"
OUT_CSV   = f"{SCORES}/GSE116265_Si_salt_stress_context.csv"
OUT_FIG   = f"{FIGS_DIR}/GSE116265_Si_salt_stress_context.png"
OUT_MD    = f"{AGT_DIR}/GSE116265_stress_context_report.md"

# ── calibration ranges ──────────────────────────────────────────────────────
cal = pd.read_csv(CAL_CSV)
cal_range = {}
for _, row in cal.iterrows():
    cal_range[row['axis']] = abs(row['treat_score'] - row['ctrl_score'])

# Axis metadata: (column_in_detailed_scores, calibration_key, display_label)
AXES = [
    ('p_sig',   'P_axis',    'P'),
    ('n_sig',   'N_axis',    'N'),
    ('ja_sig',  'JA_axis',   'JA'),
    ('sa_sig',  'SA_axis',   'SA'),
    ('aba_sig', 'ABA_axis',  'ABA'),
    ('aux_sig', 'Auxin_axis','Auxin'),
    ('fe_sig',  'Fe_axis',   'Fe'),
    ('et_sig',  'ET_axis',   'ET'),
]

# ── load data ──────────────────────────────────────────────────────────────
print("Loading Cucumber_Si_detailed_scores.csv ...")
df = pd.read_csv(INPUT_CSV)
print(f"  Loaded {len(df):,} genes, columns: {df.columns.tolist()}")

# Check which axis signature columns already exist
present_sigs = [c for c in df.columns if c.endswith('_sig')]
print(f"  Present sig columns: {present_sigs}")

# ── add Fe and ET axis signatures if not already present ──────────────────
# Use csa1g_to_tair_final.csv which uses 'Csa1G...' format matching gene_ids in data
ortho = pd.read_csv(f"{BASE}/RNA_seq_axes/csa1g_to_tair_final.csv")
print(f"  Ortholog file columns: {ortho.columns.tolist()}")
ortho = ortho.rename(columns={'Csa1G': 'csa_id', 'TAIR_ID': 'tair_id',
                               'conf': 'confidence', 'pct_id': 'pct_id'})
ortho = ortho[['csa_id','tair_id','confidence','pct_id']].dropna(subset=['tair_id'])
# One-to-one: best ortholog per cucumber gene
ortho_best = (ortho.sort_values(['confidence','pct_id'], ascending=False)
                   .drop_duplicates('csa_id', keep='first')
                   [['csa_id','tair_id']])
print(f"  Best-ortholog pairs: {len(ortho_best):,}")

def add_axis_sig(df, col_name, sig_file, ortho_best):
    """Merge axis logFC signature into df via ortholog mapping, return df."""
    if col_name in df.columns:
        return df
    sig = pd.read_csv(sig_file)
    # Standardise GeneID column
    gcol = [c for c in sig.columns if c.lower() in ('geneid','gene_id','gene')][0]
    sig = sig.rename(columns={gcol: 'tair_id', 'logFC': 'logFC'})
    merged = ortho_best.merge(sig[['tair_id','logFC']], on='tair_id', how='inner')
    # Aggregate if multiple cucumber genes map to same TAIR (take mean)
    merged = merged.groupby('csa_id')['logFC'].mean().reset_index()
    merged = merged.rename(columns={'logFC': col_name, 'csa_id': 'gene_id'})
    df = df.merge(merged, on='gene_id', how='left')
    n = df[col_name].notna().sum()
    print(f"  Added {col_name}: {n:,} genes with signature value")
    return df

if 'fe_sig' not in df.columns:
    df = add_axis_sig(df, 'fe_sig', f"{SIGS_DIR}/Fe_axis_logFC.csv", ortho_best)
if 'et_sig' not in df.columns:
    df = add_axis_sig(df, 'et_sig', f"{SIGS_DIR}/ET_axis_logFC.csv", ortho_best)

print(f"  Available sig columns now: {[c for c in df.columns if c.endswith('_sig')]}")

# ── RPKM → log2(RPKM + 1) transformation ─────────────────────────────────
# Use pre-computed mean RPKM columns from the detailed scores file
# (already computed from raw xlsx files)
eps = 1.0
df['logRPKM_CT']   = np.log2(df['mean_CT']   + eps)
df['logRPKM_Si']   = np.log2(df['mean_Si']   + eps)
df['logRPKM_Na']   = np.log2(df['mean_Na']   + eps)
df['logRPKM_NaSi'] = np.log2(df['mean_NaSi'] + eps)

# ── compute per-gene effect vectors ───────────────────────────────────────
# All log2FC comparisons; use log-space difference = log2FC
df['lfc_Si']       = df['logRPKM_Si']   - df['logRPKM_CT']    # Si vs Control
df['lfc_NaCl']     = df['logRPKM_Na']   - df['logRPKM_CT']    # NaCl vs Control
df['lfc_Si_salt']  = df['logRPKM_NaSi'] - df['logRPKM_Na']    # Si+NaCl vs NaCl  ← key
df['lfc_interact'] = df['logRPKM_NaSi'] - df['logRPKM_Na'] - df['logRPKM_Si'] + df['logRPKM_CT']

# ── project onto each axis ────────────────────────────────────────────────
# score for a given vector = sum(per_gene_lfc * axis_logFC_signature)
results = {}

for (sig_col, cal_key, label) in AXES:
    if sig_col not in df.columns:
        print(f"  WARNING: {sig_col} not in df — skipping {label}")
        continue
    mask = df[sig_col].notna() & df['lfc_Si'].notna()
    sub  = df[mask].copy()
    n    = len(sub)

    # Raw dot-product scores
    s_Si       = (sub['lfc_Si']      * sub[sig_col]).sum()
    s_NaCl     = (sub['lfc_NaCl']    * sub[sig_col]).sum()
    s_Si_salt  = (sub['lfc_Si_salt'] * sub[sig_col]).sum()
    s_interact = (sub['lfc_interact']* sub[sig_col]).sum()

    # Calibration: express as % of calibration range
    cr = cal_range.get(cal_key, None)
    if cr is None or cr == 0:
        print(f"  WARNING: no calibration range for {cal_key}")
        cr = 1.0

    results[label] = {
        'axis'          : label,
        'n_genes'       : n,
        'raw_Si'        : s_Si,
        'raw_NaCl'      : s_NaCl,
        'raw_Si_salt'   : s_Si_salt,
        'raw_interact'  : s_interact,
        'pct_Si'        : 100 * s_Si        / cr,
        'pct_NaCl'      : 100 * s_NaCl      / cr,
        'pct_Si_salt'   : 100 * s_Si_salt   / cr,
        'pct_interact'  : 100 * s_interact  / cr,
        'calibration_range': cr,
    }
    print(f"  {label:6s}: V_Si={100*s_Si/cr:+7.2f}%  V_Si_salt={100*s_Si_salt/cr:+7.2f}%  "
          f"V_NaCl={100*s_NaCl/cr:+7.2f}%  V_interact={100*s_interact/cr:+7.2f}%  "
          f"(n={n:,})")

# ── cosine similarity between V_Si and V_Si_salt ─────────────────────────
# Use the raw dot-product scores vector (one value per axis)
axis_order = [label for (_, _, label) in AXES if label in results]
v_Si       = np.array([results[a]['pct_Si']      for a in axis_order])
v_NaCl     = np.array([results[a]['pct_NaCl']    for a in axis_order])
v_Si_salt  = np.array([results[a]['pct_Si_salt'] for a in axis_order])
v_interact = np.array([results[a]['pct_interact'] for a in axis_order])

def cosine_similarity(a, b):
    norm_a = np.linalg.norm(a)
    norm_b = np.linalg.norm(b)
    if norm_a == 0 or norm_b == 0:
        return 0.0
    return np.dot(a, b) / (norm_a * norm_b)

cos_Si_vs_salt = cosine_similarity(v_Si, v_Si_salt)
cos_Si_vs_NaCl = cosine_similarity(v_Si, v_NaCl)
angle_deg      = np.degrees(np.arccos(np.clip(cos_Si_vs_salt, -1, 1)))

# Pearson r between V_Si and V_Si_salt axis scores
pearson_r, pearson_p = pearsonr(v_Si, v_Si_salt)

print(f"\nCosine similarity (V_Si vs V_Si_salt): {cos_Si_vs_salt:.4f}")
print(f"Rotation angle:                         {angle_deg:.1f} degrees")
print(f"Cosine similarity (V_Si vs V_NaCl):    {cos_Si_vs_NaCl:.4f}")
print(f"Pearson r (V_Si vs V_Si_salt):          {pearson_r:.4f}  p={pearson_p:.4f}")

# ── build output dataframe ────────────────────────────────────────────────
out_rows = []
for a in axis_order:
    r = results[a]
    out_rows.append({
        'axis'             : a,
        'n_genes'          : r['n_genes'],
        'pct_Si'           : round(r['pct_Si'],       3),
        'pct_NaCl'         : round(r['pct_NaCl'],     3),
        'pct_Si_salt'      : round(r['pct_Si_salt'],  3),
        'pct_interact'     : round(r['pct_interact'], 3),
        'raw_Si'           : round(r['raw_Si'],       4),
        'raw_NaCl'         : round(r['raw_NaCl'],     4),
        'raw_Si_salt'      : round(r['raw_Si_salt'],  4),
        'raw_interact'     : round(r['raw_interact'], 4),
        'calibration_range': round(r['calibration_range'], 4),
    })

out_df = pd.DataFrame(out_rows)
# Append summary row
summary = pd.DataFrame([{
    'axis'        : 'SUMMARY',
    'n_genes'     : '',
    'pct_Si'      : '',
    'pct_NaCl'    : '',
    'pct_Si_salt' : '',
    'pct_interact': '',
    'raw_Si'      : '',
    'raw_NaCl'    : '',
    'raw_Si_salt' : '',
    'raw_interact': '',
    'calibration_range': '',
}])
out_df.to_csv(OUT_CSV, index=False)

# Append cosine/angle metadata at bottom as comments
with open(OUT_CSV, 'a') as f:
    f.write(f"\n# cosine_similarity_Si_vs_Si_salt,{cos_Si_vs_salt:.6f}\n")
    f.write(f"# rotation_angle_degrees,{angle_deg:.2f}\n")
    f.write(f"# cosine_similarity_Si_vs_NaCl,{cos_Si_vs_NaCl:.6f}\n")
    f.write(f"# pearson_r_Si_vs_Si_salt,{pearson_r:.6f}\n")
    f.write(f"# pearson_p_Si_vs_Si_salt,{pearson_p:.6f}\n")
    f.write(f"# CAVEAT: Si-alone group has only 2 replicates\n")
print(f"Saved: {OUT_CSV}")

# ── FIGURE ────────────────────────────────────────────────────────────────
COLORS = {
    'V_Si'       : '#2196F3',  # blue
    'V_Si_salt'  : '#FF5722',  # deep orange
    'V_NaCl'     : '#9C27B0',  # purple
    'V_interact' : '#4CAF50',  # green
}

fig = plt.figure(figsize=(14, 6))
gs  = gridspec.GridSpec(1, 2, width_ratios=[1.2, 1], wspace=0.35,
                        left=0.08, right=0.97, top=0.88, bottom=0.15)

# ── Panel A: grouped bar chart ────────────────────────────────────────────
ax_a = fig.add_subplot(gs[0])
n_axes = len(axis_order)
x      = np.arange(n_axes)
w      = 0.20

bars_Si       = ax_a.bar(x - 1.5*w, v_Si,       w, label='V_Si (Si − Control)',
                          color=COLORS['V_Si'],       alpha=0.85, zorder=3)
bars_NaCl     = ax_a.bar(x - 0.5*w, v_NaCl,     w, label='V_NaCl (NaCl − Control)',
                          color=COLORS['V_NaCl'],     alpha=0.85, zorder=3)
bars_Si_salt  = ax_a.bar(x + 0.5*w, v_Si_salt,  w, label='V_Si|NaCl (Si+NaCl − NaCl)',
                          color=COLORS['V_Si_salt'],  alpha=0.85, zorder=3)
bars_interact = ax_a.bar(x + 1.5*w, v_interact, w, label='V_interact',
                          color=COLORS['V_interact'], alpha=0.85, zorder=3)

ax_a.axhline(0, color='black', linewidth=0.8, zorder=2)
ax_a.set_xticks(x)
ax_a.set_xticklabels(axis_order, fontsize=10)
ax_a.set_ylabel('Axis score (% calibration range)', fontsize=10)
ax_a.set_xlabel('Axis', fontsize=10)
ax_a.set_title('Panel A: Fingerprint vectors across axes\nGSE116265 Cucumber Si × NaCl', fontsize=10, pad=8)
ax_a.legend(fontsize=7.5, loc='upper right', framealpha=0.9)
ax_a.grid(axis='y', linewidth=0.4, alpha=0.5, zorder=1)
ax_a.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(5))

# value labels on bars (only for bars whose abs > 0.5%)
for bars_set in [bars_Si, bars_NaCl, bars_Si_salt, bars_interact]:
    for bar in bars_set:
        h = bar.get_height()
        if abs(h) > 0.5:
            ax_a.text(bar.get_x() + bar.get_width()/2., h,
                      f'{h:+.1f}', ha='center',
                      va='bottom' if h >= 0 else 'top',
                      fontsize=6.5, rotation=90)

# ── Panel B: scatter V_Si vs V_Si_salt ───────────────────────────────────
ax_b = fig.add_subplot(gs[1])

scatter = ax_b.scatter(v_Si, v_Si_salt,
                       s=80, c=range(n_axes),
                       cmap='tab10', zorder=4, edgecolors='white', linewidths=0.5)

# Label each point with axis name
for i, a in enumerate(axis_order):
    ax_b.annotate(a, (v_Si[i], v_Si_salt[i]),
                  textcoords='offset points', xytext=(5, 3), fontsize=9)

# Diagonal y=x line (perfect conservation)
all_vals = np.concatenate([v_Si, v_Si_salt])
lim_min  = min(all_vals) * 1.3
lim_max  = max(all_vals) * 1.3
diag_lo  = min(lim_min, -0.5)
diag_hi  = max(lim_max,  0.5)
ax_b.plot([diag_lo, diag_hi], [diag_lo, diag_hi],
          'k--', linewidth=0.9, alpha=0.5, label='y = x (perfect conservation)')

ax_b.axhline(0, color='grey', linewidth=0.5, alpha=0.4)
ax_b.axvline(0, color='grey', linewidth=0.5, alpha=0.4)

ax_b.set_xlabel('V_Si score (% calibration range)', fontsize=10)
ax_b.set_ylabel('V_Si|NaCl score (% calibration range)', fontsize=10)
ax_b.set_title('Panel B: Si fingerprint conservation\nunder NaCl stress', fontsize=10, pad=8)

# Annotation box with metrics
textstr = (f'cos(V_Si, V_Si|NaCl) = {cos_Si_vs_salt:.3f}\n'
           f'Rotation angle = {angle_deg:.1f}°\n'
           f'Pearson r = {pearson_r:.3f}  p = {pearson_p:.3f}')
props = dict(boxstyle='round', facecolor='lightyellow', alpha=0.85)
ax_b.text(0.04, 0.97, textstr, transform=ax_b.transAxes,
          fontsize=8.5, verticalalignment='top', bbox=props)

ax_b.legend(fontsize=8, loc='lower right')
ax_b.grid(linewidth=0.3, alpha=0.4)

# footnote caveat
fig.text(0.01, 0.01,
         '* Si-alone: n=2 replicates (Control, NaCl, Si+NaCl: n=3). '
         'Fingerprint vectors from log2(RPKM+1) differences. '
         'Axes = Arabidopsis signatures projected via Arabidopsis-cucumber orthologs.',
         fontsize=7, color='grey')

plt.savefig(OUT_FIG, dpi=150, bbox_inches='tight')
print(f"Saved: {OUT_FIG}")
plt.close()

# ── MARKDOWN REPORT ────────────────────────────────────────────────────────
report_lines = [
    "# GSE116265 Stress-Context Cosine Similarity Report",
    "",
    "**Dataset**: GSE116265 — Cucumber (*Cucumis sativus*) silicon × salt stress",
    "**Analysis date**: 2026-02-24",
    "**Author**: automated via compute_GSE116265_stress_context.py",
    "",
    "## Experimental design",
    "",
    "| Condition | Replicates | Abbreviation |",
    "|-----------|------------|--------------|",
    "| Control   | 3          | CT           |",
    "| Si alone  | **2** (CAVEAT) | Si        |",
    "| NaCl alone| 3          | Na           |",
    "| Si + NaCl | 3          | NaSi         |",
    "",
    "> **CAVEAT**: The Si-alone group has only **n=2 replicates**, compared to n=3 for",
    "> all other conditions. This reduces statistical power for V_Si and V_Si_under_NaCl",
    "> computations. Axis scores should be interpreted with caution; magnitude estimates",
    "> may be inflated by sampling noise.",
    "",
    "## Fingerprint vectors (% calibration range)",
    "",
    "- **V_Si** = log2(Si) − log2(Control)  ← pure Si effect in normal conditions",
    "- **V_NaCl** = log2(NaCl) − log2(Control)  ← pure salt effect",
    "- **V_Si|NaCl** = log2(Si+NaCl) − log2(NaCl)  ← Si effect *within* salt stress",
    "- **V_interaction** = log2(Si+NaCl) − log2(NaCl) − log2(Si) + log2(Control)  ← synergy/antagonism",
    "",
    "All values expressed as % of Arabidopsis calibration range for each axis.",
    "Expression data: log2(RPKM + 1), using precomputed condition means.",
    "Axis signatures: Arabidopsis logFC vectors, mapped to cucumber via one-to-one orthologs",
    "(Ensembl Plants Compara, confidence ≥ 0; best-match by confidence then %identity).",
    "",
    "### Axis-by-axis scores",
    "",
    "| Axis | n genes | V_Si (%) | V_NaCl (%) | V_Si\\|NaCl (%) | V_interact (%) |",
    "|------|---------|----------|------------|----------------|----------------|",
]

for a in axis_order:
    r = results[a]
    report_lines.append(
        f"| {a:6s} | {r['n_genes']:,} | "
        f"{r['pct_Si']:+7.2f} | "
        f"{r['pct_NaCl']:+7.2f} | "
        f"{r['pct_Si_salt']:+7.2f} | "
        f"{r['pct_interact']:+7.2f} |"
    )

report_lines += [
    "",
    "## Key metric: cosine similarity",
    "",
    f"| Metric | Value |",
    f"|--------|-------|",
    f"| cos(V_Si, V_Si\\|NaCl) | **{cos_Si_vs_salt:.4f}** |",
    f"| Rotation angle | **{angle_deg:.1f}°** |",
    f"| cos(V_Si, V_NaCl) | {cos_Si_vs_NaCl:.4f} |",
    f"| Pearson r (V_Si vs V_Si\\|NaCl) | {pearson_r:.4f} (p = {pearson_p:.4f}) |",
    "",
    "## Interpretation",
    "",
    f"The cosine similarity between the Si-alone fingerprint (V_Si) and the Si",
    f"fingerprint under salt stress (V_Si|NaCl) is **{cos_Si_vs_salt:.3f}**,",
    f"corresponding to a rotation angle of **{angle_deg:.1f}°**.",
    "",
]

# Qualitative interpretation based on angle
if angle_deg < 20:
    interp = (
        "This is a **small rotation** (< 20°), indicating that the Si transcriptomic "
        "fingerprint is **largely conserved** under NaCl stress. The plant's Si response "
        "engages similar biological pathways regardless of whether salt stress is present."
    )
elif angle_deg < 45:
    interp = (
        f"This is a **moderate rotation** ({angle_deg:.0f}°), indicating **partial conservation** "
        "of the Si fingerprint under NaCl stress. Some axes shift in direction or magnitude, "
        "suggesting context-dependent modulation of Si's transcriptomic mode of action."
    )
elif angle_deg < 90:
    interp = (
        f"This is a **substantial rotation** ({angle_deg:.0f}°), indicating that the Si "
        "transcriptomic fingerprint is **substantially altered** under NaCl stress. Si engages "
        "partially different biological axes when the plant is simultaneously experiencing "
        "salt stress, suggesting stress-context dependence of the biostimulant response."
    )
else:
    interp = (
        f"This is an **extreme rotation** ({angle_deg:.0f}° ≥ 90°), indicating that the Si "
        "fingerprint under NaCl stress is **orthogonal or reversed** relative to the normal-condition "
        "Si fingerprint. Si's transcriptomic mode of action is fundamentally different under salt "
        "stress, pointing to strong stress-context dependence."
    )

report_lines.append(interp)
report_lines += [
    "",
    "### V_NaCl vs V_Si comparison",
    f"cos(V_Si, V_NaCl) = {cos_Si_vs_NaCl:.3f}: indicates whether Si-alone response",
    "resembles the salt-stress response. Values near 0 indicate orthogonality (distinct modes);",
    "values near 1 indicate overlap (Si may prime for salt tolerance via shared pathways).",
    "",
    "### Interaction vector V_interact",
    "V_interact captures non-additive effects. Axes with large |V_interact| values signal",
    "synergistic or antagonistic interaction between Si and NaCl beyond their independent effects.",
    "",
    "## Caveats",
    "",
    "1. **n=2 replicates for Si-alone**: This is the most critical caveat. With only 2 replicates,",
    "   there is no estimate of within-group variance for the Si condition. The log2FC vectors",
    "   (V_Si, V_Si|NaCl) cannot be filtered by statistical significance per gene. The fingerprint",
    "   should be treated as a directional estimate rather than a precise quantity.",
    "2. **RPKM-based analysis**: Data are RPKM (not count-based). Projection uses log2(RPKM+1)",
    "   differences, which approximate log fold-change but are not TMM- or DESeq2-normalized.",
    "3. **Cross-species ortholog mapping**: Arabidopsis axis signatures projected onto cucumber",
    "   via Ensembl Compara orthologs. Some axes may lose coverage in cucumber.",
    "4. **No statistical testing**: With only 2 Si replicates, formal differential expression",
    "   testing was not feasible; all per-gene values are arithmetic means.",
    "",
    "## Output files",
    "",
    "| File | Description |",
    "|------|-------------|",
    f"| `RNA_seq_axes/scores/GSE116265_Si_salt_stress_context.csv` | Axis scores for all 4 vectors |",
    f"| `state_space_figures/GSE116265_Si_salt_stress_context.png` | Dual-panel figure |",
    f"| `agent_coordination/GSE116265_stress_context_report.md` | This report |",
    "",
    "---",
    "*Generated automatically by compute_GSE116265_stress_context.py*",
]

with open(OUT_MD, 'w', encoding='utf-8') as f:
    f.write('\n'.join(report_lines) + '\n')
print(f"Saved: {OUT_MD}")
print("\nAll done.")
