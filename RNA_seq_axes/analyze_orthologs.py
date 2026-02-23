import pandas as pd
df = pd.read_csv("C:/Users/moshe/Dropbox/ISF 2025/RNA_seq_axes/tair_to_wheat_orthologs.tsv", sep="\t")
df.columns = ["tair_id", "wheat_id", "orthology_type"]
wheat_ids = df["wheat_id"].dropna()
print(f"Total rows: {len(df)}")
print(f"Wheat IDs with values: {len(wheat_ids)}")
print("First 10 wheat IDs:")
for i in list(wheat_ids[:10]):
    print(f"  {repr(i)}")
one2one = df[df["orthology_type"] == "ortholog_one2one"].dropna(subset=["wheat_id"])
print(f"\none2one orthologs: {len(one2one)}")
traes_cs = wheat_ids.str.startswith("TraesCS").sum()
traes_old = wheat_ids.str.startswith("Traes_").sum()
print(f"TraesCS format (RefSeq v1.0): {traes_cs}")
print(f"Traes_ format (TGACv1/v2.2):  {traes_old}")
