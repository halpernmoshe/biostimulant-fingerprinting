#!/usr/bin/env python3
import numpy as np, pandas as pd, gzip
from pathlib import Path

ROOT      = Path("C:/Users/moshe/Dropbox/ISF 2025")
AXES_DIR  = ROOT / "RNA_seq_axes"
SIG_DIR   = AXES_DIR / "signatures"
SCORE_DIR = AXES_DIR / "scores"
SCORE_DIR.mkdir(exist_ok=True)

print("=" * 60)
print("Cross-species scoring: GSE217631 tomato -> Arabidopsis axes")
print("=" * 60)

# 1. Load tomato FPKM data
print("\n1. Loading GSE217631 FPKM data...")
fpkm_path = ROOT / "RNA_seq_data/suppl_files/GSE217631/GSE217631_fpkm_sample.xls.gz"
with gzip.open(fpkm_path, "rt") as fh:
    fpkm = pd.read_csv(fh, sep="\t", index_col=0)
print(f"   Shape: {fpkm.shape}")
print(f"   Columns: {list(fpkm.columns)}")
print(f"   First gene IDs: {list(fpkm.index[:5])}")

fpkm.index = fpkm.index.str.replace(r"\.\d+$", "", regex=True)
print(f"   After stripping versions: e.g. {list(fpkm.index[:3])}")

leaf_ctrl_cols = [c for c in fpkm.columns if "LeafCtr" in c or "Leaf_Ctrl" in c or ("leaf" in c.lower() and "ctr" in c.lower())]
leaf_trt_cols  = [c for c in fpkm.columns if "LeafTrt" in c or "Leaf_Trt" in c or ("leaf" in c.lower() and "trt" in c.lower())]
root_ctrl_cols = [c for c in fpkm.columns if "RootCtr" in c or "Root_Ctrl" in c or ("root" in c.lower() and "ctr" in c.lower())]
root_trt_cols  = [c for c in fpkm.columns if "RootTrt" in c or "Root_Trt" in c or ("root" in c.lower() and "trt" in c.lower())]

if not leaf_ctrl_cols:
    cols = list(fpkm.columns)
    leaf_ctrl_cols = cols[0:3]
    leaf_trt_cols  = cols[3:6]
    root_ctrl_cols = cols[6:9]
    root_trt_cols  = cols[9:12]
    print("   Used positional column assignment")

print(f"   Leaf ctrl: {leaf_ctrl_cols}")
print(f"   Leaf trt:  {leaf_trt_cols}")
print(f"   Root ctrl: {root_ctrl_cols}")
print(f"   Root trt:  {root_trt_cols}")

# 2. Load ortholog tables
print("\n2. Loading ortholog tables...")
solyc_orth = pd.read_csv(AXES_DIR / "tair_to_solyc_orthologs.tsv", sep="\t")
solyc_orth.columns = ["tair_id", "solyc_id", "solyc_name", "orthology_type"]
solyc_one2one = solyc_orth[solyc_orth["orthology_type"] == "ortholog_one2one"].copy()
solyc_one2one = solyc_one2one.dropna(subset=["solyc_id"])
solyc_one2one["solyc_id_base"] = solyc_one2one["solyc_id"].str.replace(r"\.\d+$", "", regex=True)
solyc_one2one = solyc_one2one.drop_duplicates(subset=["tair_id"])
solyc_one2one = solyc_one2one.drop_duplicates(subset=["solyc_id_base"])
print(f"   Tomato one2one orthologs: {len(solyc_one2one)} pairs")

rice_orth = pd.read_csv(AXES_DIR / "tair_to_rice_orthologs.tsv", sep="\t")
rice_orth.columns = ["tair_id", "rice_id", "orthology_type"]
rice_one2one = rice_orth[rice_orth["orthology_type"] == "ortholog_one2one"].dropna(subset=["rice_id"])
rice_one2one = rice_one2one.drop_duplicates(subset=["tair_id"])
print(f"   Rice one2one orthologs:   {len(rice_one2one)} pairs")

maize_orth = pd.read_csv(AXES_DIR / "tair_to_maize_orthologs.tsv", sep="\t")
maize_orth.columns = ["tair_id", "maize_id", "orthology_type"]
maize_one2one = maize_orth[maize_orth["orthology_type"] == "ortholog_one2one"].dropna(subset=["maize_id"])
maize_one2one = maize_one2one.drop_duplicates(subset=["tair_id"])
print(f"   Maize one2one orthologs:  {len(maize_one2one)} pairs")

# 3. Load axis signatures
print("\n3. Loading axis signatures...")
P_sig   = pd.read_csv(SIG_DIR / "P_axis_logFC.csv", index_col=0)["logFC"]
N_sig   = pd.read_csv(SIG_DIR / "N_axis_logFC.csv", index_col=0)["logFC"]
ABA_sig = pd.read_csv(SIG_DIR / "ABA_axis_logFC.csv", index_col=0)["logFC"]
print(f"   P-axis:   {len(P_sig)} genes")
print(f"   N-axis:   {len(N_sig)} genes")
print(f"   ABA-axis: {len(ABA_sig)} genes")

# 4. Score function (FPKM -> log2 -> dot product)
def score_axis_via_orthologs(fpkm_df, axis_sig, ortho_df, species="tomato"):
    if "solyc_id_base" in ortho_df.columns:
        id_col = "solyc_id_base"
    elif "rice_id" in ortho_df.columns:
        id_col = "rice_id"
    else:
        id_col = ortho_df.columns[1]
    mapping     = ortho_df.set_index("tair_id")[id_col].to_dict()
    axis_tair   = axis_sig.index.intersection(list(mapping.keys()))
    species_ids = [mapping[g] for g in axis_tair]
    avail       = [sid for sid in species_ids if sid in fpkm_df.index]
    tair_avail  = [axis_tair[i] for i, sid in enumerate(species_ids) if sid in fpkm_df.index]
    n_overlap = len(avail)
    n_orth    = len(axis_tair)
    n_total   = len(axis_sig)
    print(f"   [{species}] axis genes: {n_total} total -> {n_orth} with orthologs -> {n_overlap} in dataset ({n_overlap/n_total*100:.1f}%)")
    if n_overlap < 50:
        print(f"   WARNING: too few genes ({n_overlap}) -- skipping")
        return None, n_total, n_orth, n_overlap
    sig_subset  = axis_sig[tair_avail].values
    expr_subset = np.log2(fpkm_df.loc[avail].astype(float) + 1)
    scores      = expr_subset.multiply(sig_subset, axis=0).sum(axis=0)
    return scores, n_total, n_orth, n_overlap

# 5. Load calibration
calib = pd.read_csv(AXES_DIR / "calibration.csv", index_col=0)
print("\n4. Calibration:")
print(calib)

def to_pct(score, axis_name):
    row = calib.loc[axis_name]
    return (score - row["ctrl_score"]) / (row["treat_score"] - row["ctrl_score"]) * 100

# 6. Score each group
print("\n5. Scoring GSE217631 groups...")
groups = {
    "LeafCtrl": leaf_ctrl_cols,
    "LeafTrt":  leaf_trt_cols,
    "RootCtrl": root_ctrl_cols,
    "RootTrt":  root_trt_cols,
}
results = []
overlap_stats = {}

for axis_name, axis_sig in [("P_axis", P_sig), ("N_axis", N_sig), ("ABA_axis", ABA_sig)]:
    print(f"\n  --- {axis_name} ---")
    scores, n_tot, n_orth, n_ovl = score_axis_via_orthologs(fpkm, axis_sig, solyc_one2one, "tomato")
    overlap_stats[axis_name] = {"total": n_tot, "with_ortholog": n_orth, "in_dataset": n_ovl,
                                "pct_covered": round(n_ovl / n_tot * 100, 1)}
    if scores is None:
        continue
    for group_name, cols in groups.items():
        cols_avail = [c for c in cols if c in scores.index]
        if not cols_avail:
            cols_avail = [c for c in scores.index if any(g.lower() in c.lower() for g in cols)]
        if not cols_avail:
            print(f"   WARNING: no columns found for {group_name}")
            continue
        grp_mean = scores[cols_avail].mean()
        grp_pct  = to_pct(grp_mean, axis_name)
        results.append({
            "dataset":      "GSE217631",
            "species":      "tomato",
            "biostimulant": "Protein_hydrolysate_seaweed",
            "axis":         axis_name,
            "group":        group_name,
            "tissue":       "Leaf" if "Leaf" in group_name else "Root",
            "treatment":    "PSI-475_treated" if "Trt" in group_name else "Control",
            "raw_score":    grp_mean,
            "pct_of_range": grp_pct,
        })
        print(f"   {group_name}: raw={grp_mean:.1f}  pct={grp_pct:.2f}%")

results_df = pd.DataFrame(results)
print("\n6. Per-group scores:")
print(results_df.to_string())

# 7. Compute within-experiment deltas
print("\n7. Treatment effects (delta = Trt - Ctrl):")
delta_rows = []
for tissue in ["Leaf", "Root"]:
    ctrl_rows = results_df[(results_df["tissue"] == tissue) & (results_df["treatment"] == "Control")]
    trt_rows  = results_df[(results_df["tissue"] == tissue) & (results_df["treatment"] == "PSI-475_treated")]
    for axis_name in ["P_axis", "N_axis", "ABA_axis"]:
        ctrl_pct = ctrl_rows[ctrl_rows["axis"] == axis_name]["pct_of_range"].values
        trt_pct  = trt_rows[trt_rows["axis"] == axis_name]["pct_of_range"].values
        if len(ctrl_pct) and len(trt_pct):
            delta = trt_pct[0] - ctrl_pct[0]
            print(f"   {tissue} {axis_name}: Ctrl={ctrl_pct[0]:.2f}%  Trt={trt_pct[0]:.2f}%  Delta={delta:+.2f}%")
            delta_rows.append({"tissue": tissue, "axis": axis_name,
                               "ctrl_pct": ctrl_pct[0], "trt_pct": trt_pct[0], "delta_pct": delta})

delta_df = pd.DataFrame(delta_rows)
print("\nDelta summary:")
print(delta_df.to_string())

# 8. Save scores
out_path = SCORE_DIR / "GSE217631_projected_scores.csv"
results_df.to_csv(out_path, index=False)
print(f"\n8. Saved: {out_path}")

# 9. Conserved gene module (Task F)
print("\n" + "=" * 60)
print("Task F: Conserved gene module analysis")
print("=" * 60)

tomato_tair   = set(solyc_one2one["tair_id"])
rice_tair     = set(rice_one2one["tair_id"])
maize_tair    = set(maize_one2one["tair_id"])
all_3_species = tomato_tair & rice_tair & maize_tair
print(f"\nTAIR genes with one2one ortholog in:")
print(f"  Tomato: {len(tomato_tair)}")
print(f"  Rice:   {len(rice_tair)}")
print(f"  Maize:  {len(maize_tair)}")
print(f"  ALL 3:  {len(all_3_species)}")

conserved_stats = {}
for axis_name, sig in [("P_axis", P_sig), ("N_axis", N_sig), ("ABA_axis", ABA_sig)]:
    strong    = sig[sig.abs() > 1].index
    conserved = strong.intersection(all_3_species)
    pct = len(conserved) / len(strong) * 100 if len(strong) > 0 else 0
    conserved_stats[axis_name] = {
        "strong_responders":  len(strong),
        "conserved_in_all_3": len(conserved),
        "pct_conserved":      round(pct, 1),
    }
    print(f"\n  {axis_name}: {len(strong)} strong responders -> {len(conserved)} conserved in tomato+rice+maize ({pct:.1f}%)")

conserved_df  = pd.DataFrame(conserved_stats).T
conserved_csv = AXES_DIR / "conserved_core_modules.csv"
conserved_df.to_csv(conserved_csv)
print(f"\nSaved conserved modules: {conserved_csv}")

print("\nAxis coverage via ortholog mapping:")
for axis, stats in overlap_stats.items():
    print(f"  {axis}: {stats['in_dataset']}/{stats['total']} genes covered ({stats['pct_covered']}%)")

print("\n=== DONE ===")
