"""
Build additional axes for exploratory axis expansion:
1. Drought-axis from GSE226105 (rLog-normalized counts)
2. Gene-set axes using TAIR GO annotations:
   - Brassinosteroid (GO:0009741)
   - Oxidative stress (GO:0006979)
   - Hypoxia (GO:0001666)
   - Cell wall remodeling (GO:0042546)
   - Ribosome biogenesis (GO:0042254) — TOR activity proxy

Then score all biostimulant datasets on all axes and run statistical selection.
"""
import pandas as pd
import numpy as np
from pathlib import Path
import urllib.request
import gzip
import re

ROOT    = Path("C:/Users/moshe/Dropbox/ISF 2025")
SIG_DIR = ROOT / "RNA_seq_axes" / "signatures"
SCORE_DIR = ROOT / "RNA_seq_axes" / "scores"
SUPPL   = ROOT / "RNA_seq_data" / "suppl_files"
CAL_FILE = ROOT / "RNA_seq_axes" / "calibration.csv"
SIG_DIR.mkdir(parents=True, exist_ok=True)

LOG = ROOT / "RNA_seq_axes" / "build_more_axes_log.txt"

def log(msg):
    print(msg)
    with open(LOG, 'a') as f:
        f.write(msg + '\n')

log("=" * 60)
log("AXIS EXPANSION PIPELINE  2026-02-23")
log("=" * 60)

# ── Load biostimulant expression matrices for scoring ─────────────────────────
# We'll build these after we know which scoring approach to use

# ── 1. DROUGHT AXIS from GSE226105 rLog-normalized ───────────────────────────
log("\n=== DROUGHT AXIS (GSE226105) ===")

drought_file = SUPPL / "GSE226105" / "GSE226105_rLog-Normalized_counts.csv.gz"
if drought_file.exists():
    df_d = pd.read_csv(drought_file, index_col=0)
    log(f"  Loaded: {df_d.shape[0]} genes × {df_d.shape[1]} samples")
    log(f"  Columns: {df_d.columns.tolist()[:10]}")

    # Identify control vs drought columns
    # From series matrix: "CONT" = control, "endpoint_DROUGHT" = drought endpoint
    ctrl_cols    = [c for c in df_d.columns if c.startswith('CONT')]
    drought_cols = [c for c in df_d.columns if 'DROUGHT' in c and 'endpoint' in c.lower()]
    log(f"  Control cols: {ctrl_cols}")
    log(f"  Drought cols: {drought_cols}")

    if ctrl_cols and drought_cols:
        # Compute fold-change: log2((mean_drought + ε) / (mean_ctrl + ε))
        # But data is already rLog-normalized → simply use mean_drought - mean_ctrl
        # (log2FC from rLog values = rLog_treated - rLog_control)
        mean_ctrl    = df_d[ctrl_cols].mean(axis=1)
        mean_drought = df_d[drought_cols].mean(axis=1)
        logfc = mean_drought - mean_ctrl  # difference of rLog values ≈ log2FC

        # Remove low-variance/zero genes
        sig_drought = logfc[logfc.abs() > 0.01].copy()
        sig_drought.name = 'logFC'
        sig_drought.index.name = 'GeneID'

        out = SIG_DIR / "Drought_axis_logFC.csv"
        sig_drought.to_csv(out, header=True)
        log(f"  Saved: {out.name} ({len(sig_drought)} genes)")
        log(f"  logFC range: {sig_drought.min():.2f} to {sig_drought.max():.2f}")
        log(f"  Top 10 genes UP (drought): {sig_drought.nlargest(10).to_dict()}")
        log(f"  Top 10 genes DOWN: {sig_drought.nsmallest(10).to_dict()}")
    else:
        log("  ERROR: Could not identify control/drought columns")
else:
    log(f"  File not found: {drought_file}")


# ── 2. GENE-SET AXES from TAIR GO annotations ─────────────────────────────────
log("\n=== GENE-SET AXES (TAIR GO annotations) ===")

# Download TAIR GO annotation (GOSLIM)
GO_FILE = ROOT / "RNA_seq_data" / "ATH_GO_GOSLIM.txt"

if not GO_FILE.exists():
    log("  Downloading TAIR GO SLIM annotations...")
    try:
        url = "https://www.arabidopsis.org/download_files/GO_and_PO_Annotations/Gene_Ontology_Annotations/ATH_GO_GOSLIM.txt"
        urllib.request.urlretrieve(url, GO_FILE)
        log(f"  Downloaded: {GO_FILE.stat().st_size // 1024} KB")
    except Exception as e:
        log(f"  GOSLIM download failed: {e}")
        # Try alternative URL
        try:
            url2 = "https://www.arabidopsis.org/download_files/GO_and_PO_Annotations/Gene_Ontology_Annotations/gene_association.tair.gz"
            urllib.request.urlretrieve(url2, str(GO_FILE) + ".gz")
            with gzip.open(str(GO_FILE) + ".gz", 'rb') as f_in:
                with open(GO_FILE, 'wb') as f_out:
                    f_out.write(f_in.read())
            log(f"  Downloaded via alternative URL")
        except Exception as e2:
            log(f"  Alternative URL also failed: {e2}")

# Try full GO annotation if GOSLIM not available
FULL_GO_FILE = ROOT / "RNA_seq_data" / "gene_association.tair"
FULL_GO_GZ   = ROOT / "RNA_seq_data" / "gene_association.tair.gz"

if not GO_FILE.exists() and not FULL_GO_GZ.exists():
    log("  Trying full TAIR GO annotation download...")
    try:
        url3 = "https://www.arabidopsis.org/download_files/GO_and_PO_Annotations/Gene_Ontology_Annotations/gene_association.tair.gz"
        urllib.request.urlretrieve(url3, FULL_GO_GZ)
        log(f"  Downloaded: {FULL_GO_GZ.stat().st_size // 1024} KB")
    except Exception as e3:
        log(f"  Full GO download failed: {e3}")

# Load whichever GO file is available
go_df = None
if GO_FILE.exists():
    try:
        # GOSLIM format: tab-separated with columns including gene ID and GO terms
        go_df = pd.read_csv(GO_FILE, sep='\t', comment='!', header=None,
                            low_memory=False, on_bad_lines='skip')
        log(f"  GO SLIM loaded: {go_df.shape}")
        log(f"  Columns sample: {go_df.iloc[0].tolist()[:5]}")
    except Exception as e:
        log(f"  Error loading GOSLIM: {e}")

elif FULL_GO_GZ.exists():
    try:
        go_df = pd.read_csv(FULL_GO_GZ, sep='\t', comment='!', header=None,
                            low_memory=False, on_bad_lines='skip')
        log(f"  Full GO loaded: {go_df.shape}")
    except Exception as e:
        log(f"  Error loading full GO: {e}")

if go_df is None:
    log("  No GO annotation available — skipping gene-set axes")
else:
    # Gene associations typically: col[1]=gene_ID (AT1G...), col[3]=GO_id
    # or col[2]=gene_name, col[4]=GO_id depending on format
    # Try to find which column has AT*G* gene IDs and which has GO:* terms
    gene_col = None
    go_col   = None
    for i in range(min(10, go_df.shape[1])):
        vals = go_df.iloc[:, i].astype(str)
        if vals.str.contains(r'AT\dG', na=False).mean() > 0.3:
            gene_col = i
        if vals.str.contains(r'GO:\d{7}', na=False).mean() > 0.3:
            go_col = i

    log(f"  Detected gene column: {gene_col}, GO column: {go_col}")

    if gene_col is not None and go_col is not None:
        go_df.columns = range(go_df.shape[1])
        gene_to_gos = go_df.groupby(gene_col)[go_col].apply(set)
        log(f"  Unique genes with GO annotations: {len(gene_to_gos)}")

        # Define GO term → axis name mapping
        # For each axis: (name, GO_terms_UP, GO_terms_DOWN, description)
        GO_AXES = {
            "Brassinosteroid": (
                ["GO:0009741", "GO:0009742", "GO:0010268"],  # response to BR, signal transduction
                [],
                "Brassinosteroid response (GO:0009741)"
            ),
            "OxidativeStress": (
                ["GO:0006979", "GO:0055114"],  # response to oxidative stress, redox
                [],
                "Oxidative stress response (GO:0006979)"
            ),
            "Hypoxia": (
                ["GO:0001666", "GO:0036294"],  # response to hypoxia, anaerobic
                [],
                "Hypoxia / anaerobic response (GO:0001666)"
            ),
            "CellWall": (
                ["GO:0042546", "GO:0009832", "GO:0071669"],  # cell wall assembly/deposition/modification
                [],
                "Cell wall assembly (GO:0042546)"
            ),
            "RibosomeBiogenesis": (
                ["GO:0042254", "GO:0042273"],  # ribosome biogenesis, ribosomal large subunit
                [],
                "Ribosome biogenesis (GO:0042254) — TOR activity proxy"
            ),
        }

        for axis_name, (go_up, go_down, desc) in GO_AXES.items():
            go_up_set  = set(go_up)
            go_down_set = set(go_down)

            genes_up   = set()
            genes_down = set()

            for gene, gos in gene_to_gos.items():
                gene_str = str(gene).strip()
                if not re.match(r'AT\dG\d+', gene_str):
                    continue
                if gos & go_up_set:
                    genes_up.add(gene_str)
                if gos & go_down_set:
                    genes_down.add(gene_str)

            log(f"\n  {axis_name}: {len(genes_up)} UP genes, {len(genes_down)} DOWN genes")
            log(f"    ({desc})")

            if not genes_up and not genes_down:
                log(f"    -> No genes found — GO terms may not be in GOSLIM")
                continue

            # Create gene-set signature: UP genes = +1, DOWN genes = -1
            all_genes = genes_up | genes_down
            sig = pd.Series(0.0, index=sorted(all_genes))
            for g in genes_up:
                sig[g] = 1.0
            for g in genes_down:
                sig[g] = -1.0

            sig.name = 'logFC'
            sig.index.name = 'GeneID'
            out = SIG_DIR / f"{axis_name}_axis_logFC.csv"
            sig.to_csv(out, header=True)
            log(f"    Saved: {out.name}")
    else:
        log("  Could not identify gene ID and GO columns — manual check needed")

log("\n=== DONE ===")
