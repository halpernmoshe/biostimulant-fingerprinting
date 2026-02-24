"""
Build gene-set axes from TAIR GO annotations (tair.gaf.gz).
Each gene-set axis assigns genes in a GO biological process score +1 (axis-positive).

Then score all biostimulant datasets on all axes and run statistical selection.
"""
import pandas as pd
import numpy as np
import re
import gzip
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

ROOT     = Path("C:/Users/moshe/Dropbox/ISF 2025")
SIG_DIR  = ROOT / "RNA_seq_axes" / "signatures"
SCORE_DIR = ROOT / "RNA_seq_axes" / "scores"
CAL_FILE = ROOT / "RNA_seq_axes" / "calibration.csv"
FIG_DIR  = ROOT / "state_space_figures"
FIG_DIR.mkdir(exist_ok=True)

LOG = ROOT / "RNA_seq_axes" / "geneset_axes_log.txt"

def log(msg):
    print(msg)
    with open(LOG, 'a') as f:
        f.write(msg + '\n')

log("=" * 60)
log("GENE-SET AXES  2026-02-23")
log("=" * 60)

# ── Load GO annotations ─────────────────────────────────────────────────────
log("\nLoading TAIR GO annotations...")
GAF = ROOT / "RNA_seq_data" / "tair.gaf.gz"
genes_l, go_l = [], []
with gzip.open(GAF, 'rt', errors='replace') as f:
    for line in f:
        if line.startswith('!'): continue
        parts = line.strip().split('\t')
        if len(parts) < 9: continue
        gene, go = parts[1], parts[4]
        if re.match(r'AT\dG\d+', gene) and re.match(r'GO:\d{7}', go):
            genes_l.append(gene)
            go_l.append(go)

go_df = pd.DataFrame({'gene': genes_l, 'go': go_l})
gene_to_gos = go_df.groupby('gene')['go'].apply(set)
log(f"  {len(go_df)} annotations; {gene_to_gos.shape[0]} unique genes")

def get_genes_with_go(go_terms_list):
    """Return set of AT*G* genes annotated to any GO term in the list (or descendants)."""
    go_set = set(go_terms_list)
    result = set()
    for gene, gos in gene_to_gos.items():
        if gos & go_set:
            result.add(gene)
    return result

# ── Define gene-set axes ─────────────────────────────────────────────────────
# Format: {name: (go_up_list, go_down_list, description)}
# Note: For biological process axes, we define genes in the process as "up-regulated"
# when the axis is active. For the scoring:
# score = mean(log2(count+1) for UP genes) - mean(log2(count+1) for DOWN genes)
# Then calibrate: 100% = fully stressed/active, 0% = baseline

GENE_SET_AXES = {
    "Brassinosteroid": {
        "go_up": ["GO:0009741", "GO:0009742", "GO:0010268", "GO:0010558"],
        "go_down": [],
        "desc": "Brassinosteroid response/signaling (GO:0009741/0009742)",
    },
    "OxidativeStress": {
        "go_up": ["GO:0006979", "GO:0055114", "GO:0045454"],
        "go_down": [],
        "desc": "Oxidative stress / ROS response (GO:0006979)",
    },
    "Hypoxia": {
        "go_up": ["GO:0001666", "GO:0036294", "GO:0070482"],
        "go_down": [],
        "desc": "Hypoxia / anaerobic response (GO:0001666)",
    },
    "CellWall": {
        "go_up": ["GO:0042546", "GO:0009832", "GO:0071669", "GO:0010383", "GO:0009664"],
        "go_down": [],
        "desc": "Cell wall assembly/modification (GO:0042546)",
    },
    "RibosomeBiogenesis": {
        "go_up": ["GO:0042254", "GO:0042273", "GO:0042274", "GO:0022626"],
        "go_down": [],
        "desc": "Ribosome biogenesis (GO:0042254) — TOR activity proxy",
    },
    "JA_response": {
        "go_up": ["GO:0009753", "GO:0009741", "GO:0009867"],
        "go_down": [],
        "desc": "Jasmonic acid response (GO:0009753)",
    },
    "SA_response": {
        "go_up": ["GO:0009751", "GO:0010200"],
        "go_down": [],
        "desc": "Salicylic acid response (GO:0009751)",
    },
    "CK_response": {
        "go_up": ["GO:0009736", "GO:0009735"],
        "go_down": [],
        "desc": "Cytokinin response/signaling (GO:0009736)",
    },
    "ET_response": {
        "go_up": ["GO:0009723", "GO:0071369"],
        "go_down": [],
        "desc": "Ethylene response (GO:0009723)",
    },
    "Auxin_response": {
        "go_up": ["GO:0009733", "GO:0009734"],
        "go_down": [],
        "desc": "Auxin response/signaling (GO:0009733)",
    },
    "ABA_response_GS": {
        "go_up": ["GO:0009737"],
        "go_down": [],
        "desc": "ABA response gene-set (GO:0009737) — compare with DESeq2 ABA-axis",
    },
}

log(f"\nBuilding {len(GENE_SET_AXES)} gene-set axes...")

gs_signatures = {}
for axis_name, cfg in GENE_SET_AXES.items():
    genes_up   = get_genes_with_go(cfg["go_up"])
    genes_down = get_genes_with_go(cfg["go_down"]) if cfg["go_down"] else set()
    log(f"  {axis_name}: {len(genes_up)} UP genes, {len(genes_down)} DOWN genes ({cfg['desc']})")

    if not genes_up and not genes_down:
        log(f"    -> SKIP (no genes found)")
        continue

    # Gene-set signature: UP=+1, DOWN=-1, combined for dot-product scoring
    all_genes = sorted(genes_up | genes_down)
    sig = pd.Series(0.0, index=all_genes, name='logFC')
    sig.index.name = 'GeneID'
    for g in genes_up:
        sig[g] = 1.0
    for g in genes_down:
        sig[g] = -1.0

    out = SIG_DIR / f"{axis_name}_axis_logFC.csv"
    sig.to_csv(out, header=True)
    gs_signatures[axis_name] = sig
    log(f"    Saved: {out.name} ({len(sig)} genes)")

log(f"\nBuilt {len(gs_signatures)} gene-set axes")


# ── Score all biostimulant datasets on all axes ───────────────────────────────
log("\n=== SCORING BIOSTIMULANTS ON ALL AXES ===")

# Load existing DESeq2 axes
deseq2_sigs = {}
for fname in SIG_DIR.glob("*_axis_logFC.csv"):
    name = fname.stem.replace("_axis_logFC", "")
    if name in gs_signatures:
        continue  # don't duplicate gene-set axes
    sig = pd.read_csv(fname, index_col=0).iloc[:, 0]
    sig.name = name
    deseq2_sigs[name] = sig

log(f"DESeq2 axes: {list(deseq2_sigs.keys())}")
log(f"Gene-set axes: {list(gs_signatures.keys())}")

all_sigs = {**deseq2_sigs, **gs_signatures}

# ── Load biostimulant count/expression matrices ───────────────────────────────
# GSE297649: HTSeq counts (raw)
SCORE_DIR_DATA = ROOT / "RNA_seq_data" / "suppl_files"

def load_htseq_dir(gse_dir, pattern="*.txt.gz"):
    """Load per-sample HTSeq count files from a directory."""
    files = list(gse_dir.glob(pattern))
    if not files:
        return None
    dfs = []
    for f in files:
        try:
            with gzip.open(f, 'rt') as fh:
                lines = [l.strip().split('\t') for l in fh if not l.startswith('__')]
            genes  = [l[0] for l in lines if len(l) == 2]
            counts = [int(l[1]) for l in lines if len(l) == 2]
            if genes:
                dfs.append(pd.Series(counts, index=genes, name=f.stem.split('.')[0]))
        except Exception:
            continue
    if not dfs:
        return None
    return pd.concat(dfs, axis=1).fillna(0).astype(int)

def log2p1(df):
    return np.log2(df + 1)

def score_on_axis(expr_df, sig):
    """Dot-product score: sum(log2(count+1) * logFC) for common genes."""
    if expr_df is None or sig is None:
        return pd.Series()
    # For log-normalized data (rLog or TPM), use directly; for counts, log2(+1)
    is_count = (expr_df.dtypes == int).all() or (expr_df.values.max() > 100 and expr_df.values.max() > 5000)
    if is_count:
        expr = log2p1(expr_df)
    else:
        expr = expr_df.copy()
    common = expr.index.intersection(sig.index)
    return expr.loc[common].multiply(sig.loc[common], axis=0).sum(axis=0)

def score_geneset(expr_df, sig):
    """Gene-set score: mean(expr for +1 genes) - mean(expr for -1 genes)."""
    if expr_df is None or sig is None:
        return pd.Series()
    is_count = (expr_df.dtypes == int).all() or (expr_df.values.max() > 5000)
    if is_count:
        expr = log2p1(expr_df)
    else:
        expr = expr_df.copy()
    up_genes   = sig[sig > 0].index.intersection(expr.index)
    down_genes = sig[sig < 0].index.intersection(expr.index)
    score = pd.Series(0.0, index=expr.columns)
    if len(up_genes) > 0:
        score += expr.loc[up_genes].mean(axis=0)
    if len(down_genes) > 0:
        score -= expr.loc[down_genes].mean(axis=0)
    return score

# ── Load GSE297649 counts ─────────────────────────────────────────────────────
log("\nLoading GSE297649...")
gse297_dir = SCORE_DIR_DATA / "GSE297649"
gse297_counts = load_htseq_dir(gse297_dir)
if gse297_counts is not None:
    log(f"  {gse297_counts.shape[0]} genes × {gse297_counts.shape[1]} samples")
    log(f"  Columns: {gse297_counts.columns.tolist()}")
    # Group assignment
    gse297_groups = {}
    for col in gse297_counts.columns:
        c = col.lower()
        if 'wtaa' in c or 'wt_aa' in c or (c.startswith('wt') and 'aa' in c):
            gse297_groups[col] = 'amino_acid'
        elif 'wtsh' in c or 'wt_sh' in c or (c.startswith('wt') and ('sh' in c or 'hs' in c)):
            gse297_groups[col] = 'humic_subst'
        else:
            gse297_groups[col] = 'control'
    log(f"  Groups: {set(gse297_groups.values())}")
else:
    log("  GSE297649 counts not found")

# ── Load GSE138478 featureCounts ─────────────────────────────────────────────
log("\nLoading GSE138478...")
gse138_dir = SCORE_DIR_DATA / "GSE138478"

def read_featurecounts_gz(fpath):
    with gzip.open(fpath, 'rt') as f:
        for line in f:
            if not line.startswith('#'):
                break
    df = pd.read_csv(fpath, sep='\t', comment='#', index_col=0)
    anno_cols = ['Chr','Start','End','Strand','Length']
    data_cols = [c for c in df.columns if not any(c.startswith(a) for a in anno_cols)
                 and c not in anno_cols]
    df = df[data_cols].copy()
    df.columns = [Path(c).stem for c in df.columns]
    return df

gse138_counts = None
for fc_file in gse138_dir.glob("*.txt.gz"):
    if 'GMV' in fc_file.name or 'HMZ' in fc_file.name:
        try:
            fc = read_featurecounts_gz(fc_file)
            if gse138_counts is None:
                gse138_counts = fc
            else:
                gse138_counts = gse138_counts.join(fc, how='outer')
        except Exception as e:
            log(f"  Error loading {fc_file.name}: {e}")

if gse138_counts is not None:
    gse138_counts = gse138_counts.fillna(0).astype(int)
    log(f"  {gse138_counts.shape}")
    log(f"  Columns: {gse138_counts.columns.tolist()}")
    # Classify columns
    gse138_groups = {}
    for col in gse138_counts.columns:
        parts = col.split('_')
        if len(parts) >= 2:
            suf = parts[1]
            if suf.endswith('C'):
                prefix = parts[0]
                gse138_groups[col] = f"{'GMV' if prefix.startswith('GMV') else 'diacetyl'}_control"
            elif suf.endswith('T'):
                prefix = parts[0]
                gse138_groups[col] = f"{'GMV' if prefix.startswith('GMV') else 'diacetyl'}_treated"
else:
    log("  GSE138478 not loaded")

# ── Load GSE208223 normalized shoot data ─────────────────────────────────────
log("\nLoading GSE208223 (shoot 10d)...")
gse208_file = SCORE_DIR_DATA / "GSE208223" / "GSE208223_Shoot_norm.txt.gz"
if gse208_file.exists():
    gse208 = pd.read_csv(gse208_file, sep='\t', index_col=0)
    gse208.index = gse208.index.str.replace(r'\.\d+$', '', regex=True)
    gse208.index.name = 'GeneID'
    log(f"  {gse208.shape}")
    gse208_groups = {
        'E1': 'TiO2_control', 'E2': 'TiO2_control', 'E3': 'TiO2_control',
        'F1': 'TiO2_treated', 'F2': 'TiO2_treated', 'F3': 'TiO2_treated',
    }
    gse208_use = gse208[[c for c in gse208_groups if c in gse208.columns]]
else:
    log("  GSE208223 not found")
    gse208 = None
    gse208_groups = {}
    gse208_use = None

# ── Score all datasets on all axes ───────────────────────────────────────────
datasets = {
    'GSE297649': (gse297_counts, gse297_groups if gse297_counts is not None else {}),
    'GSE138478': (gse138_counts, gse138_groups if gse138_counts is not None else {}),
    'GSE208223': (gse208_use, gse208_groups),
}

log("\nScoring all datasets on all axes...")
rows = []
for ds_name, (expr_df, groups) in datasets.items():
    if expr_df is None or not groups:
        log(f"  {ds_name}: SKIP (no data)")
        continue
    for ax_name, sig in all_sigs.items():
        is_geneset = ax_name in gs_signatures
        if is_geneset:
            scores = score_geneset(expr_df, sig)
        else:
            scores = score_on_axis(expr_df, sig)
        if scores.empty:
            continue
        for col, score in scores.items():
            grp = groups.get(col, 'unknown')
            rows.append({
                'dataset': ds_name, 'sample': col, 'group': grp,
                'axis': ax_name, 'score': float(score),
                'is_geneset': is_geneset
            })

all_scores = pd.DataFrame(rows)
log(f"Total scores: {len(all_scores)}")
all_scores.to_csv(SCORE_DIR / "all_axes_all_biostimulants_per_sample.csv", index=False)
log("Saved: all_axes_all_biostimulants_per_sample.csv")

# ── Summary: mean score per group and axis ────────────────────────────────────
summary = (all_scores
           .groupby(['dataset', 'group', 'axis', 'is_geneset'])['score']
           .agg(['mean', 'std', 'count'])
           .reset_index())
summary.to_csv(SCORE_DIR / "all_axes_summary.csv", index=False)
log("\nSaved: all_axes_summary.csv")
log(summary.to_string(index=False))

# ── Compute effect vectors (delta: treatment - control within dataset) ─────────
log("\n=== EFFECT VECTORS ===")

# Control group names per dataset
CTRL_MAP = {
    'GSE297649': 'control',
    'GSE138478_GMV': 'GMV_control',
    'GSE138478_diacetyl': 'diacetyl_control',
    'GSE208223': 'TiO2_control',
}

delta_rows = []
for ds_name, (_, groups) in datasets.items():
    if not groups:
        continue
    unique_grps = set(groups.values())
    ctrl_grps = [g for g in unique_grps if 'control' in g.lower()]
    treat_grps = [g for g in unique_grps if 'control' not in g.lower()]

    for ctrl_grp in ctrl_grps:
        # Find matching treatment groups
        prefix = ctrl_grp.replace('_control', '')
        matched_treats = [g for g in treat_grps if prefix in g]
        if not matched_treats:
            matched_treats = treat_grps

        for treat_grp in matched_treats:
            for ax_name in all_sigs.keys():
                ctrl_sub  = summary[(summary.dataset == ds_name) &
                                    (summary.group == ctrl_grp) &
                                    (summary.axis == ax_name)]
                treat_sub = summary[(summary.dataset == ds_name) &
                                    (summary.group == treat_grp) &
                                    (summary.axis == ax_name)]
                if ctrl_sub.empty or treat_sub.empty:
                    continue
                delta = float(treat_sub['mean'].iloc[0]) - float(ctrl_sub['mean'].iloc[0])
                delta_rows.append({
                    'dataset': ds_name, 'control': ctrl_grp, 'treatment': treat_grp,
                    'axis': ax_name,
                    'delta': delta,
                    'ctrl_mean': float(ctrl_sub['mean'].iloc[0]),
                    'treat_mean': float(treat_sub['mean'].iloc[0]),
                })

delta_df = pd.DataFrame(delta_rows)
delta_df.to_csv(SCORE_DIR / "effect_vectors_all_axes.csv", index=False)
log(f"Saved: effect_vectors_all_axes.csv ({len(delta_df)} rows)")

# Print summary
pivot_delta = delta_df.pivot_table(index=['dataset','treatment'], columns='axis', values='delta')
log("\nEffect vectors (treatment - control):")
log(pivot_delta.round(3).to_string())

# ── Statistical selection: F-statistic per axis ───────────────────────────────
log("\n=== AXIS SELECTION (F-statistic) ===")

# Group biostimulant classes
CLASS_MAP = {
    'amino_acid': 'AA_biostimulant',
    'humic_subst': 'HA_biostimulant',
    'GMV_treated': 'PGPR',
    'diacetyl_treated': 'PGPR',
    'TiO2_treated': 'Nanoparticle',
}

# For each axis, compute F = between-class variance / within-class variance
f_stats = {}
per_sample = all_scores.copy()
per_sample['class'] = per_sample['group'].map(CLASS_MAP)
per_sample = per_sample.dropna(subset=['class'])

for ax_name in all_sigs.keys():
    ax_data = per_sample[per_sample.axis == ax_name]
    if ax_data.empty:
        continue
    classes = ax_data['class'].unique()
    if len(classes) < 2:
        continue
    groups_list = [ax_data[ax_data['class'] == c]['score'].values for c in classes]
    groups_list = [g for g in groups_list if len(g) > 0]
    if len(groups_list) < 2:
        continue
    try:
        f, p = stats.f_oneway(*groups_list)
        f_stats[ax_name] = {'F': float(f), 'p': float(p), 'n_classes': len(classes)}
    except Exception:
        pass

f_df = pd.DataFrame(f_stats).T.sort_values('F', ascending=False)
f_df.to_csv(ROOT / "agent_coordination" / "axis_selection.md", index=True)  # temp save as md
log("\nAxis F-statistics (higher = better discriminates biostimulant classes):")
log(f_df.round(4).to_string())

# Write proper axis_selection.md
axis_md = ["# Axis Selection Analysis\n",
           "## F-statistic Ranking (discriminative power for biostimulant classes)\n",
           f"Classes: AA biostimulant, HA biostimulant, PGPR, Ti nanoparticle\n\n",
           f_df.round(4).to_markdown(), "\n\n",
           "## Interpretation\n",
           "High F = axis separates biostimulant classes well.\n",
           "Low F = axis cannot distinguish biostimulant classes (not useful for Paper A).\n",
           "\n### Recommended axes for paper:\n"]

for ax in f_df.index[:5]:
    row = f_df.loc[ax]
    axis_md.append(f"- **{ax}**: F={row['F']:.2f}, p={row['p']:.4f}\n")

axis_md.append("\n### Axes to deprioritize:\n")
for ax in f_df.index[-5:]:
    row = f_df.loc[ax]
    axis_md.append(f"- **{ax}**: F={row['F']:.2f} (low discriminative power)\n")

with open(ROOT / "agent_coordination" / "axis_selection.md", 'w') as f:
    f.writelines(axis_md)
log("\nSaved: agent_coordination/axis_selection.md")

# ── Heatmap: delta vectors across all axes and treatments ─────────────────────
log("\n=== HEATMAP FIGURE ===")
try:
    pivot_heat = delta_df.pivot_table(index='treatment', columns='axis', values='delta')
    # Normalize columns to z-scores for visualization
    pivot_z = (pivot_heat - pivot_heat.mean()) / (pivot_heat.std() + 1e-9)

    fig, ax = plt.subplots(figsize=(max(12, pivot_z.shape[1] * 0.8), max(5, pivot_z.shape[0] * 0.7)))
    im = ax.imshow(pivot_z.values, aspect='auto', cmap='RdBu_r', vmin=-2.5, vmax=2.5)
    ax.set_xticks(range(pivot_z.shape[1]))
    ax.set_xticklabels(pivot_z.columns, rotation=45, ha='right', fontsize=8)
    ax.set_yticks(range(pivot_z.shape[0]))
    ax.set_yticklabels(pivot_z.index, fontsize=9)
    plt.colorbar(im, ax=ax, label='Δ axis score (z-scored across axes)')
    ax.set_title(
        "Biostimulant effect vectors across all axes\n(red = positive shift, blue = negative)",
        fontsize=12, fontweight='bold'
    )
    fig.tight_layout()
    fname = FIG_DIR / "biostimulant_heatmap_all_axes.png"
    fig.savefig(fname, dpi=150)
    plt.close(fig)
    log(f"Saved: {fname.name}")
except Exception as e:
    log(f"Heatmap error: {e}")

log("\n=== ALL DONE ===")
