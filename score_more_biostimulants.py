#!/usr/bin/env python3
"""
score_more_biostimulants.py
===========================
Downloads and scores additional biostimulant datasets on existing axes.

Datasets:
  GSE138478  PGPR volatiles (GMV) + diacetyl, Arabidopsis, featureCounts
  GSE208223  TiO2/Si nanoparticles, Arabidopsis, RNA-seq
  GSE217631  Protein hydrolysate PSI-475, tomato — gene overlap check only

Outputs combined multi-biostimulant state-space figure.
"""
import sys, gzip, urllib.request, warnings
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime

ROOT      = Path(__file__).parent
SUPPL     = ROOT / "RNA_seq_data" / "suppl_files"
MATRIX    = ROOT / "RNA_seq_data" / "series_matrices"
SIG_DIR   = ROOT / "RNA_seq_axes" / "signatures"
SCORE_DIR = ROOT / "RNA_seq_axes" / "scores"
PLOT_DIR  = ROOT / "state_space_figures"
CAL_FILE  = ROOT / "RNA_seq_axes" / "calibration.csv"
for d in [SUPPL, SCORE_DIR, PLOT_DIR]:
    d.mkdir(parents=True, exist_ok=True)

LOG = []
def log(msg):
    ts = datetime.now().strftime("%H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line)
    LOG.append(line)


# ── Utilities ────────────────────────────────────────────────

def download_file(url, dest, min_bytes=500):
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > min_bytes:
        log(f"  Cached: {dest.name}")
        return True
    log(f"  Downloading: {dest.name} ...")
    try:
        urllib.request.urlretrieve(url, dest)
        log(f"  OK ({dest.stat().st_size // 1024} KB)")
        return True
    except Exception as e:
        log(f"  FAIL: {e}")
        return False


def load_signatures():
    sigs = {}
    for f in SIG_DIR.glob("*.csv"):
        ax = f.stem.replace("_logFC", "")   # e.g. "P_axis"
        df = pd.read_csv(f, index_col=0)
        if "logFC" in df.columns:
            sigs[ax] = df["logFC"].dropna()
    log(f"  Loaded signatures: {list(sigs.keys())}")
    return sigs


def load_calibration():
    if CAL_FILE.exists():
        df = pd.read_csv(CAL_FILE)
        if "axis" in df.columns:
            df = df.set_index("axis")
        return df
    return pd.DataFrame()


def to_pct(scores, ax_name, cal):
    if ax_name in cal.index:
        c0 = float(cal.loc[ax_name, "ctrl_score"])
        c1 = float(cal.loc[ax_name, "treat_score"])
        rng = c1 - c0
        if abs(rng) > 0:
            return 100 * (scores - c0) / rng
    return scores   # raw if no calibration


def score_on_axis(counts_df, sig):
    common = sig.index.intersection(counts_df.index)
    if len(common) < 100:
        log(f"    Only {len(common)} genes in common — skipping")
        return pd.Series(dtype=float)
    s = sig[common]
    e = np.log2(counts_df.loc[common].astype(float) + 1)
    return e.multiply(s, axis=0).sum(axis=0)


def parse_series_matrix(path):
    with gzip.open(path, 'rb') as f:
        raw = f.read()
    text = raw.decode('latin-1')
    meta = {}
    for line in text.split('\n'):
        if line.startswith('!Sample_'):
            parts = line.strip().split('\t')
            key = parts[0].lstrip('!')
            vals = [v.strip('"') for v in parts[1:]]
            if key not in meta:
                meta[key] = vals
    return meta


# ── Load existing axes and calibration ──────────────────────
log("=" * 55)
log("Loading axes and calibration")
log("=" * 55)
sigs = load_signatures()
cal  = load_calibration()
log(f"  Calibration:\n{cal.to_string() if not cal.empty else '  (none)'}")

all_scored_rows = []


# ══════════════════════════════════════════════════════════════
# GSE297649: already scored — reload scores
# ══════════════════════════════════════════════════════════════
score_f = SCORE_DIR / "GSE297649_per_sample_scores.csv"
if score_f.exists():
    gse297 = pd.read_csv(score_f)
    all_scored_rows.append(gse297)
    log(f"  GSE297649: loaded {len(gse297)} rows")


# ══════════════════════════════════════════════════════════════
# GSE138478: PGPR volatiles (GMV) + diacetyl — featureCounts
# ══════════════════════════════════════════════════════════════
log("\n" + "=" * 55)
log("GSE138478: PGPR volatiles + diacetyl (Arabidopsis)")
log("=" * 55)

BASE138 = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138478/suppl"
d138 = SUPPL / "GSE138478"
d138.mkdir(exist_ok=True)

f_hmz008 = d138 / "GSE138478_HMZ008_featureCounts_output.txt.gz"
f_hmz014 = d138 / "GSE138478_HMZ014_Dione_featureCounts_output.txt.gz"

ok008 = download_file(f"{BASE138}/GSE138478_HMZ008_featureCounts_output.txt.gz", f_hmz008)
ok014 = download_file(f"{BASE138}/GSE138478_HMZ014_Dione_featureCounts_output.txt.gz", f_hmz014)

def read_featurecounts_gz(path):
    """featureCounts: header line starts with #, then col names, then data.
    Columns: Geneid, Chr, Start, End, Strand, Length, sample1, sample2, ...
    """
    with gzip.open(path, 'rb') as f:
        raw = f.read()
    text = raw.decode('latin-1')
    lines = [l for l in text.split('\n') if l and not l.startswith('#')]
    if not lines:
        return pd.DataFrame()
    header = lines[0].split('\t')
    rows = []
    for line in lines[1:]:
        parts = line.split('\t')
        if parts:
            rows.append(parts)
    df = pd.DataFrame(rows, columns=header[:len(rows[0])] if rows else header)
    df = df.set_index('Geneid')
    # Drop annotation columns, keep count columns
    count_cols = [c for c in df.columns if c not in ['Chr','Start','End','Strand','Length']]
    df = df[count_cols].apply(pd.to_numeric, errors='coerce').fillna(0).astype(int)
    # Simplify column names: take filename stem
    df.columns = [Path(c).stem for c in df.columns]
    return df


rows_138 = []

if ok008 and f_hmz008.exists():
    try:
        fc008 = read_featurecounts_gz(f_hmz008)
        log(f"  HMZ008 matrix: {fc008.shape[0]} genes x {fc008.shape[1]} samples")
        log(f"  Columns: {list(fc008.columns)}")
        # Classify by 2nd field after split: '0.05C' or '0.5C' = control; '0.05T' or '0.5T' = GMV
        # e.g. HMZ008_0.05C_5DAT_rep1 -> split('_')[1] = '0.05C' -> endswith('C')
        ctrl_cols  = [c for c in fc008.columns if c.split('_')[1].endswith('C')]
        treat_cols = [c for c in fc008.columns if c.split('_')[1].endswith('T')]
        log(f"  Control: {ctrl_cols}")
        log(f"  GMV treated: {treat_cols}")

        if ctrl_cols and treat_cols:
            for ax_name, sig in sigs.items():
                for grp, grp_cols, cat in [
                    ('GMV_control', ctrl_cols, 'Control (PGPR)'),
                    ('GMV_treated', treat_cols, 'PGPR volatile (GMV)')
                ]:
                    sc = score_on_axis(fc008[grp_cols], sig)
                    for samp, val in sc.items():
                        rows_138.append({
                            'gse': 'GSE138478', 'gsm': samp,
                            'group': grp, 'category': cat,
                            'axis': ax_name, 'score': float(val)
                        })
    except Exception as e:
        log(f"  HMZ008 error: {e}")
        import traceback; log(traceback.format_exc())

if ok014 and f_hmz014.exists():
    try:
        fc014 = read_featurecounts_gz(f_hmz014)
        log(f"  HMZ014 matrix: {fc014.shape[0]} genes x {fc014.shape[1]} samples")
        log(f"  Columns: {list(fc014.columns)}")
        ctrl_cols = [c for c in fc014.columns if 'CK' in c or 'mock' in c.lower() or '_C_' in c]
        treat_cols = [c for c in fc014.columns if 'Dione' in c or 'diacetyl' in c.lower() or '_T_' in c]
        log(f"  Control: {ctrl_cols}")
        log(f"  Diacetyl treated: {treat_cols}")

        if ctrl_cols and treat_cols:
            for ax_name, sig in sigs.items():
                for grp, grp_cols, cat in [
                    ('diacetyl_control', ctrl_cols, 'Control (diacetyl)'),
                    ('diacetyl_treated', treat_cols, 'PGPR diacetyl')
                ]:
                    sc = score_on_axis(fc014[grp_cols], sig)
                    for samp, val in sc.items():
                        rows_138.append({
                            'gse': 'GSE138478', 'gsm': samp,
                            'group': grp, 'category': cat,
                            'axis': ax_name, 'score': float(val)
                        })
    except Exception as e:
        log(f"  HMZ014 error: {e}")
        import traceback; log(traceback.format_exc())

if rows_138:
    df138 = pd.DataFrame(rows_138)
    df138.to_csv(SCORE_DIR / "GSE138478_per_sample_scores.csv", index=False)
    all_scored_rows.append(df138)
    log(f"  GSE138478: {len(df138)} rows saved")
    summ = df138.groupby(['gse','category','group','axis'])['score'].agg(['mean','std','count'])
    log(f"\n{summ.to_string()}")


# ══════════════════════════════════════════════════════════════
# GSE208223: TiO2/Si nanoparticles — Arabidopsis RNA-seq
# ══════════════════════════════════════════════════════════════
log("\n" + "=" * 55)
log("GSE208223: TiO2 nanoparticles (Arabidopsis)")
log("=" * 55)

sm208 = MATRIX / "GSE208223_series_matrix.txt.gz"
if not sm208.exists():
    download_file(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE208nnn/GSE208223/matrix/GSE208223_series_matrix.txt.gz",
        sm208)

rows_208 = []
if sm208.exists():
    meta208 = parse_series_matrix(sm208)
    titles208 = meta208.get('Sample_title', [])
    gsms208   = meta208.get('Sample_geo_accession', [])
    chars208  = [v for k, v in meta208.items() if 'characteristics' in k.lower()]
    log(f"  {len(gsms208)} samples: {titles208[:6]}")

    # Get supplementary file URLs
    suppl_urls = meta208.get('Sample_supplementary_file_1', [])
    ftp_urls   = [u for u in suppl_urls if u.startswith('ftp') and u.endswith('.gz')]
    log(f"  Suppl URLs (first 3): {ftp_urls[:3]}")

    # Convert ftp:// to https://
    https_urls = [u.replace('ftp://', 'https://') for u in ftp_urls]

    d208 = SUPPL / "GSE208223"
    d208.mkdir(exist_ok=True)

    # Download count files
    downloaded_208 = []
    for url in https_urls[:20]:  # cap at 20 to avoid huge datasets
        fname = url.split('/')[-1]
        dest  = d208 / fname
        if download_file(url, dest):
            downloaded_208.append(dest)

    if downloaded_208:
        # Read count files (HTSeq format)
        dfs = []
        for f in downloaded_208:
            try:
                with gzip.open(f, 'rt') as fh:
                    lines = [l.strip().split('\t') for l in fh if not l.startswith('__')]
                genes  = [l[0] for l in lines if len(l) == 2]
                counts = [int(l[1]) for l in lines if len(l) == 2]
                if genes:
                    dfs.append(pd.Series(counts, index=genes, name=f.stem.split('.')[0]))
            except Exception:
                continue
        if dfs:
            counts208 = pd.concat(dfs, axis=1).fillna(0).astype(int)
            log(f"  Loaded {counts208.shape[0]} genes x {counts208.shape[1]} samples")

            # Classify samples from title
            gsm_to_title = dict(zip(gsms208, titles208))
            ctrl_cols = [c for c in counts208.columns
                         if any(gsm in c for gsm, t in gsm_to_title.items()
                                if 'control' in t.lower() or 'mock' in t.lower()
                                or 'ctrl' in t.lower() or 'untreated' in t.lower())]
            treat_cols = [c for c in counts208.columns if c not in ctrl_cols]
            log(f"  Control: {ctrl_cols}")
            log(f"  Treated: {treat_cols}")

            if ctrl_cols and treat_cols:
                for ax_name, sig in sigs.items():
                    for grp, grp_cols, cat in [
                        ('TiO2_control', ctrl_cols, 'Control (TiO2)'),
                        ('TiO2_treated', treat_cols, 'Si/TiO2 nanoparticles')
                    ]:
                        if grp_cols:
                            sc = score_on_axis(counts208[grp_cols], sig)
                            for samp, val in sc.items():
                                rows_208.append({
                                    'gse': 'GSE208223', 'gsm': samp,
                                    'group': grp, 'category': cat,
                                    'axis': ax_name, 'score': float(val)
                                })

if rows_208:
    df208 = pd.DataFrame(rows_208)
    df208.to_csv(SCORE_DIR / "GSE208223_per_sample_scores.csv", index=False)
    all_scored_rows.append(df208)
    log(f"  GSE208223: {len(df208)} rows saved")
    summ = df208.groupby(['gse','category','group','axis'])['score'].agg(['mean','std','count'])
    log(f"\n{summ.to_string()}")


# ══════════════════════════════════════════════════════════════
# GSE217631: Protein hydrolysate PSI-475 — Tomato (FPKM)
# Gene overlap check with Arabidopsis axes
# ══════════════════════════════════════════════════════════════
log("\n" + "=" * 55)
log("GSE217631: Protein hydrolysate PSI-475 (tomato)")
log("=" * 55)

BASE217 = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE217nnn/GSE217631/suppl"
d217 = SUPPL / "GSE217631"
d217.mkdir(exist_ok=True)
f_fpkm = d217 / "GSE217631_fpkm_sample.xls.gz"

ok217 = download_file(f"{BASE217}/GSE217631_fpkm_sample.xls.gz", f_fpkm)

if ok217:
    try:
        # Try reading as tab-separated (despite .xls extension, often TSV)
        with gzip.open(f_fpkm, 'rb') as f:
            raw = f.read()
        text = raw.decode('latin-1')
        lines = text.strip().split('\n')
        header = lines[0].split('\t')
        log(f"  Columns (first 8): {header[:8]}")
        log(f"  Rows: {len(lines) - 1}")
        # Gene IDs
        first_ids = [lines[i].split('\t')[0] for i in range(1, min(6, len(lines)))]
        log(f"  First gene IDs: {first_ids}")

        # Check overlap with any axis
        all_tair = set(sigs.get('P_axis', pd.Series()).index[:5000])
        tomato_ids = set(lines[i].split('\t')[0] for i in range(1, len(lines)))
        overlap = len(tomato_ids & all_tair)
        log(f"  Gene ID overlap with P-axis (TAIR): {overlap} / {len(tomato_ids)}")
        if overlap < 100:
            log("  -> Gene IDs incompatible (tomato ITAG vs TAIR). Need orthology mapping.")
            log("  -> Skipping scoring for now. Note: use DIOPT or OrthoFinder for ITAG->TAIR mapping.")
        else:
            log("  -> Sufficient overlap for direct scoring.")
    except Exception as e:
        log(f"  Error reading GSE217631: {e}")


# ══════════════════════════════════════════════════════════════
# Combined state-space: all scored datasets
# ══════════════════════════════════════════════════════════════
log("\n" + "=" * 55)
log("COMBINED MULTI-BIOSTIMULANT STATE-SPACE FIGURE")
log("=" * 55)

if not all_scored_rows:
    log("  No data to combine")
else:
    combined = pd.concat(all_scored_rows, ignore_index=True)
    log(f"  Combined: {len(combined)} rows, {combined['gse'].nunique()} datasets")

    # Compute pct for all rows
    def apply_pct(df, cal):
        rows = []
        for _, r in df.iterrows():
            ax = r['axis']
            sc = r['score']
            pct = float(to_pct(pd.Series([sc]), ax, cal).iloc[0]) if not cal.empty else sc
            rows.append({**r.to_dict(), 'pct': pct})
        return pd.DataFrame(rows)

    combined_pct = apply_pct(combined, cal)

    # Pivot: one row per (gse, group), one col per axis pct
    pivot_all = (combined_pct
                 .groupby(['gse', 'category', 'group', 'axis'])['pct']
                 .mean()
                 .unstack('axis')
                 .reset_index())

    pivot_all.to_csv(ROOT / "all_biostimulants_scores.csv", index=False)
    log(f"  Saved: all_biostimulants_scores.csv")
    log(f"\n{pivot_all.to_string(index=False)}")

    # ── Combined state-space vector figure ──
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches

        # Color/style palette by dataset + treatment
        PALETTE = {
            # GSE297649
            'amino_acid':        ('#E69F00', 'o', 'Amino acids (AA)'),
            'humic_subst':       ('#56B4E9', 'o', 'Humic substances'),
            # GSE138478
            'GMV_treated':       ('#009E73', 's', 'PGPR volatile (GMV)'),
            'diacetyl_treated':  ('#CC79A7', 's', 'PGPR diacetyl'),
            # GSE208223
            'TiO2_treated':      ('#D55E00', '^', 'Si/TiO2 nanoparticles'),
            # Controls (grey, all merged)
            'control':           ('#999999', 'x', 'Control'),
            'GMV_control':       ('#bbbbbb', 'x', 'Control (PGPR)'),
            'diacetyl_control':  ('#cccccc', 'x', 'Control (diacetyl)'),
            'TiO2_control':      ('#dddddd', 'x', 'Control (TiO2)'),
        }

        axis_pairs = [('P_axis', 'N_axis'),
                      ('P_axis', 'ABA_axis'),
                      ('N_axis', 'ABA_axis')]

        for ax_x, ax_y in axis_pairs:
            if ax_x not in pivot_all.columns or ax_y not in pivot_all.columns:
                continue

            fig, ax = plt.subplots(figsize=(10, 8))

            # Reference landmarks
            ax.scatter([0], [0], s=100, color='#eeeeee', edgecolors='#aaaaaa',
                       marker='s', zorder=1)
            ax.scatter([100], [100], s=100, color='#ffeeee', edgecolors='#cc8888',
                       marker='s', zorder=1)

            # Identify control rows for each dataset
            ctrl_groups = ['control', 'GMV_control', 'diacetyl_control', 'TiO2_control']

            # Plot each group
            for _, row in pivot_all.iterrows():
                grp = row['group']
                if pd.isna(row.get(ax_x)) or pd.isna(row.get(ax_y)):
                    continue
                color, marker, label = PALETTE.get(grp, ('#aaaaaa', 'o', grp))
                ax.scatter([row[ax_x]], [row[ax_y]],
                           color=color, marker=marker, s=80,
                           alpha=0.8, edgecolors='black', linewidths=0.5,
                           zorder=3)

            # Draw arrows: from matching control to each treatment
            ctrl_map = {}  # gse -> (ctrl_x, ctrl_y)
            for _, row in pivot_all.iterrows():
                if row['group'] in ctrl_groups and pd.notna(row.get(ax_x)) and pd.notna(row.get(ax_y)):
                    gse = row['gse']
                    if gse not in ctrl_map:
                        ctrl_map[gse] = (row[ax_x], row[ax_y])

            for _, row in pivot_all.iterrows():
                grp = row['group']
                if grp in ctrl_groups:
                    continue
                if pd.isna(row.get(ax_x)) or pd.isna(row.get(ax_y)):
                    continue
                gse = row['gse']
                if gse not in ctrl_map:
                    continue
                cx, cy = ctrl_map[gse]
                tx, ty = row[ax_x], row[ax_y]
                color, marker, label = PALETTE.get(grp, ('#aaaaaa', 'o', grp))

                ax.annotate("", xy=(tx, ty), xytext=(cx, cy),
                            arrowprops=dict(arrowstyle="-|>", color=color,
                                           lw=2, mutation_scale=15),
                            zorder=4)
                # Label at tip
                off_x = 5 if tx >= cx else -5
                off_y = 5 if ty >= cy else -5
                ax.annotate(label, (tx, ty),
                            textcoords='offset points', xytext=(off_x, off_y),
                            fontsize=8, color=color, fontweight='bold', zorder=5)

            # Mark control centroids
            for gse, (cx, cy) in ctrl_map.items():
                ax.scatter([cx], [cy], s=120, color='#666666',
                           edgecolors='black', marker='o', zorder=4)
                ax.annotate(f"Ctrl\n({gse})", (cx, cy),
                            textcoords='offset points', xytext=(3, 6),
                            fontsize=7, color='#444444')

            ax.set_xlabel(f"{ax_x.replace('_', ' ')} (% of reference range)", fontsize=11)
            ax.set_ylabel(f"{ax_y.replace('_', ' ')} (% of reference range)", fontsize=11)
            ax.set_title(f"Multi-biostimulant physiological state-space\n"
                         f"{ax_x.replace('_', ' ')} vs {ax_y.replace('_', ' ')}",
                         fontsize=12, fontweight='bold')

            # Reference lines
            ax.axhline(0, color='#dddddd', lw=0.8, ls='--')
            ax.axvline(0, color='#dddddd', lw=0.8, ls='--')
            ax.axhline(100, color='#ffcccc', lw=0.8, ls='--')
            ax.axvline(100, color='#ffcccc', lw=0.8, ls='--')

            # Legend
            legend_handles = [
                mpatches.Patch(color=c, label=lbl)
                for grp, (c, m, lbl) in PALETTE.items()
                if any(pivot_all['group'] == grp)
            ]
            ax.legend(handles=legend_handles, fontsize=8, loc='best',
                      framealpha=0.8)
            plt.tight_layout()
            out = PLOT_DIR / f"combined_state_space_{ax_x}_vs_{ax_y}.png"
            plt.savefig(out, dpi=200)
            plt.close()
            log(f"  Saved: {out.name}")

    except Exception as e:
        import traceback
        log(f"  Figure error: {e}")
        log(traceback.format_exc())


# ══════════════════════════════════════════════════════════════
# Summary
# ══════════════════════════════════════════════════════════════
log("\n" + "=" * 55)
log("DONE")
log("=" * 55)
log(f"Scores: {[f.name for f in SCORE_DIR.glob('*.csv')]}")
log(f"Figures: {[f.name for f in PLOT_DIR.glob('combined_*.png')]}")

log_path = ROOT / "RNA_seq_axes" / f"more_biostim_log_{datetime.now().strftime('%Y%m%d_%H%M')}.txt"
log_path.write_text('\n'.join(LOG), encoding='utf-8')
log(f"Log: {log_path.name}")
