import pandas as pd, re

# The refmap file maps TraesCS transcript IDs -> TGACv1 transcript IDs
# Columns: ref_id (TraesCS transcript), ccode, tid (TGACv1 transcript), gid (TGACv1 gene), ..., ref_gene (TraesCS gene)
print("Reading refmap...")
df = pd.read_csv(
    "C:/Users/moshe/Dropbox/ISF 2025/RNA_seq_axes/iwgsc_refseqv1.0_vs_TGACv1.refmap",
    sep="\t", low_memory=False
)
print(f"Total rows: {len(df)}")
print(f"Columns: {list(df.columns)}")
print(f"\nFirst 3 rows:")
print(df.head(3).to_string())

# Extract gene-level: ref_gene (TraesCS*) <-> gid (TRIAE_CS42_*) 
# Keep only rows where both are not NA
gene_map = df[["ref_gene", "gid"]].dropna()
gene_map.columns = ["TraesCS_gene", "TGACv1_gene"]
# Strip transcript suffix from TGACv1 gene IDs (they already show gene IDs in gid)
gene_map = gene_map[gene_map["TGACv1_gene"] != "NA"]
gene_map = gene_map.drop_duplicates()
print(f"\nGene-level pairs (non-NA): {len(gene_map)}")

# Check format
print(f"\nTraesCS sample: {list(gene_map['TraesCS_gene'][:3])}")
print(f"TGACv1 sample:  {list(gene_map['TGACv1_gene'][:3])}")

# Save
out = "C:/Users/moshe/Dropbox/ISF 2025/RNA_seq_axes/traescs_to_tgacv1_gene_xref.csv"
gene_map.to_csv(out, index=False)
print(f"\nSaved: {out}")
print(f"Unique TraesCS genes: {gene_map['TraesCS_gene'].nunique()}")
print(f"Unique TGACv1 genes:  {gene_map['TGACv1_gene'].nunique()}")
