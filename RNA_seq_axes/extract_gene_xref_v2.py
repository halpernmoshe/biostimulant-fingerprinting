import pandas as pd, re

print("Reading refmap...")
df = pd.read_csv(
    "C:/Users/moshe/Dropbox/ISF 2025/RNA_seq_axes/iwgsc_refseqv1.0_vs_TGACv1.refmap",
    sep="\t", low_memory=False, usecols=["ref_gene", "gid"]
)
print(f"Total rows: {len(df)}")

# Drop NA (ccode = NA rows)
df = df.dropna()
df = df[df["gid"] != "NA"]

# The gid column contains transcript-level IDs like:
#   TRIAE_CS42_1AS_TGACv1_021658_AA0082030.1.mrna1
# Gene IDs in TGACv1 look like: TRIAE_CS42_1AS_TGACv1_021658_AA0082030
# Strip everything from .1.mrna onwards OR from the last dot-number pattern
def to_gene_id(tid):
    # Pattern: everything up to the last .N.mrnaN portion
    # e.g. TRIAE_CS42_1AS_TGACv1_021658_AA0082030.1.mrna1 -> TRIAE_CS42_1AS_TGACv1_021658_AA0082030
    return re.sub(r'\.\d+\.mrna\d+$', '', str(tid))

df["TGACv1_gene"] = df["gid"].apply(to_gene_id)
df = df.rename(columns={"ref_gene": "TraesCS_gene"})
df = df[["TraesCS_gene", "TGACv1_gene"]].drop_duplicates()

print(f"\nGene-level pairs after cleaning: {len(df)}")
print(f"\nTraesCS sample: {list(df['TraesCS_gene'][:5])}")
print(f"TGACv1 sample:  {list(df['TGACv1_gene'][:5])}")
print(f"\nUnique TraesCS genes: {df['TraesCS_gene'].nunique()}")
print(f"Unique TGACv1 genes:  {df['TGACv1_gene'].nunique()}")

# Check one-to-one stats
traes_counts = df.groupby("TraesCS_gene")["TGACv1_gene"].nunique()
tgac_counts  = df.groupby("TGACv1_gene")["TraesCS_gene"].nunique()
print(f"\nTraesCS genes with exactly 1 TGACv1 match: {(traes_counts == 1).sum()}")
print(f"TraesCS genes with >1 TGACv1 match:         {(traes_counts > 1).sum()}")
print(f"TGACv1 genes with exactly 1 TraesCS match:  {(tgac_counts == 1).sum()}")

out = "C:/Users/moshe/Dropbox/ISF 2025/RNA_seq_axes/traescs_to_tgacv1_gene_xref.csv"
df.to_csv(out, index=False)
print(f"\nSaved: {out}")
