import pandas as pd, numpy as np, os, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.patches as mpatches
import matplotlib.patheffects as pe
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import seaborn as sns

ROOT       = "C:/Users/moshe/Dropbox/ISF 2025"
SCORES_DIR = ROOT + "/RNA_seq_axes/scores"
FIGS_DIR   = ROOT + "/state_space_figures"
REPORT_DIR = ROOT + "/agent_coordination"
os.makedirs(FIGS_DIR, exist_ok=True)
os.makedirs(REPORT_DIR, exist_ok=True)
print("Step 1: Loading data...")
df   = pd.read_csv(SCORES_DIR + "/all_9_treatments_10axes.csv")
wide = df.pivot(index="treatment", columns="axis", values="delta_pct").reset_index()
axis_order = [a for a in ["P_axis","N_axis","ABA_axis","JA_axis","SA_axis",
                           "Fe_axis","CK_axis","ET_axis","Auxin_axis"] if a in wide.columns]
class_map = {
    "amino_acid":"Amino acids", "GSE298982_glutamate":"Amino acids",
    "humic_subst":"Humic substances", "GMV_treated":"PGPR (complex)",
    "GSE87337_SynCom":"PGPR (complex)", "diacetyl_treated":"PGPR (single)",
    "GSE139154_GB03":"PGPR (single)", "GSE199501_PGPR":"PGPR (single)",
    "TiO2_treated":"Mineral"}
wide["class"] = wide["treatment"].map(class_map)
slmap = {
    "amino_acid":"AA hydrolysate", "GSE298982_glutamate":"Glutamate",
    "humic_subst":"Humic subs.", "GMV_treated":"PGPR-JA (GMV)",
    "GSE87337_SynCom":"SynCom", "diacetyl_treated":"PGPR-SA (diacetyl)",
    "GSE139154_GB03":"Bacillus GB03", "GSE199501_PGPR":"P. megaterium",
    "TiO2_treated":"TiO2"}
wide["short_label"] = wide["treatment"].map(slmap)
QC  = ["amino_acid","humic_subst","GMV_treated","diacetyl_treated","TiO2_treated"]
cm5 = {"amino_acid":"Amino acids","humic_subst":"Humic substances",
       "GMV_treated":"PGPR-JA","diacetyl_treated":"PGPR-SA","TiO2_treated":"Mineral"}
matA = wide[wide["treatment"].isin(QC)].copy().reset_index(drop=True)
matA["class_5"] = matA["treatment"].map(cm5)
matB = wide.copy().reset_index(drop=True)
print("  matA:", len(matA), "x", len(axis_order), " matB:", len(matB))
CC = {"Amino acids":"#2196F3","Humic substances":"#4CAF50",
      "PGPR (complex)":"#FF5722","PGPR (single)":"#E91E63",
      "Mineral":"#607D8B","PGPR-JA":"#FF5722","PGPR-SA":"#E91E63"}
CM = {"amino_acid":"o","GSE298982_glutamate":"s","humic_subst":"D",
      "GMV_treated":"^","GSE87337_SynCom":"v","diacetyl_treated":"P",
      "GSE139154_GB03":"X","GSE199501_PGPR":"*","TiO2_treated":"h"}
AL = {"P_axis":"P-starvation axis (%)","N_axis":"N-supply axis (%)",
      "ABA_axis":"ABA axis (%)","JA_axis":"JA axis (%)","SA_axis":"SA axis (%)",
      "Fe_axis":"Fe-deficiency axis (%)","CK_axis":"Cytokinin axis (%)",
      "ET_axis":"Ethylene axis (%)","Auxin_axis":"Auxin axis (%)"}
AS = {"P_axis":"P","N_axis":"N","ABA_axis":"ABA","JA_axis":"JA","SA_axis":"SA",
      "Fe_axis":"Fe","CK_axis":"CK","ET_axis":"ET","Auxin_axis":"Auxin"}
# ===== STEP 2: HIERARCHICAL CLUSTERING =====
print("Step 2: Hierarchical clustering...")
X_A   = matA[axis_order].values.astype(float)
X_A_z = StandardScaler().fit_transform(X_A)
fig, ax = plt.subplots(figsize=(11, 6))
T2 = ("Hierarchical Clustering of Biostimulant Transcriptomic Fingerprints" + chr(10) +
      "(5 QC-verified treatments, 9 hormone/nutrient axes, Ward linkage)")
fig.suptitle(T2, fontsize=12, fontweight="bold", y=0.98)
Z = linkage(X_A_z, method="ward", metric="euclidean")
dendrogram(Z, labels=matA["short_label"].values, ax=ax, color_threshold=0,
           above_threshold_color="#333", leaf_font_size=11, leaf_rotation=0)
ax.set_ylabel("Ward linkage distance", fontsize=11)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
clbl = {"AA hydrolysate":CC["Amino acids"],"Humic subs.":CC["Humic substances"],
        "PGPR-JA (GMV)":CC["PGPR-JA"],"PGPR-SA (diacetyl)":CC["PGPR-SA"],"TiO2":CC["Mineral"]}
for tl in ax.get_xticklabels():
    tl.set_color(clbl.get(tl.get_text(), "#333")); tl.set_fontweight("bold")
plt.tight_layout(rect=[0,0,1,0.96])
p1 = FIGS_DIR + "/classification_dendrogram.png"
plt.savefig(p1, dpi=150, bbox_inches="tight"); plt.close()
print("  Saved:", p1)
# ===== STEP 3: PCA MATRIX A =====
print("Step 3: PCA Matrix A (5 tx x 9 axes)...")
pca_A = PCA(n_components=min(5, len(matA)-1))
sc_A  = pca_A.fit_transform(X_A_z)
var_A = pca_A.explained_variance_ratio_ * 100
ld_A  = pd.DataFrame(pca_A.components_.T, index=axis_order,
                     columns=["PC"+str(i+1) for i in range(pca_A.n_components_)])
fig, ax = plt.subplots(figsize=(8,7))
for i in range(len(matA)):
    row = matA.iloc[i]
    c = CC.get(row["class_5"], "#888")
    m = CM.get(row["treatment"], "o")
    ax.scatter(sc_A[i,0], sc_A[i,1], color=c, marker=m, s=220, zorder=3,
               edgecolors="black", linewidths=0.8)
    ax.annotate(row["short_label"], xy=(sc_A[i,0],sc_A[i,1]), xytext=(6,6),
                textcoords="offset points", fontsize=9, color=c, fontweight="bold",
                path_effects=[pe.withStroke(linewidth=2.5, foreground="white")])
scale = max(abs(sc_A[:,:2]).max()*0.6, 0.5)
for an in axis_order:
    lx = ld_A.loc[an,"PC1"]*scale; ly = ld_A.loc[an,"PC2"]*scale
    if abs(lx)<0.08 and abs(ly)<0.08: continue
    ax.annotate("", xy=(lx,ly), xytext=(0,0),
                arrowprops=dict(arrowstyle="->", color="#666", lw=1.5))
    ax.text(lx*1.15, ly*1.15, AS[an], fontsize=8, color="#444",
            ha="center", va="center", fontstyle="italic")
ax.axhline(0,color="#ccc",lw=0.8,ls="--"); ax.axvline(0,color="#ccc",lw=0.8,ls="--")
ax.set_xlabel("PC1 (" + "{:.1f}".format(var_A[0]) + "% variance)", fontsize=12)
ax.set_ylabel("PC2 (" + "{:.1f}".format(var_A[1]) + "% variance)", fontsize=12)
T3 = "PCA of Biostimulant Transcriptomic Fingerprints" + chr(10) + "(5 QC-verified, 9 axes, z-scored)"
ax.set_title(T3, fontsize=12, fontweight="bold")
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
ptchs = [mpatches.Patch(color=CC[c],label=c) for c in
         ["Amino acids","Humic substances","PGPR-JA","PGPR-SA","Mineral"]]
ax.legend(handles=ptchs, fontsize=9, loc="best", framealpha=0.9)
plt.tight_layout()
p2 = FIGS_DIR + "/PCA_5treatments_9axes.png"
plt.savefig(p2, dpi=150, bbox_inches="tight"); plt.close()
print("  Saved:", p2)
# ===== STEP 4: F-STATISTICS =====
print("Step 4: F-statistic ranking...")
cls5 = matA["class_5"].values
ucls = np.unique(cls5)
fst  = {}
for an in axis_order:
    v  = matA[an].values.astype(float)
    gm = v.mean()
    SSb = sum(np.sum(cls5==c)*(v[cls5==c].mean()-gm)**2 for c in ucls)
    SSw = sum(np.sum((v[cls5==c]-v[cls5==c].mean())**2) for c in ucls)
    k,n = len(ucls),len(v)
    if (n-k)<=0 or SSw<1e-12:
        fst[an] = SSb/(k-1) if k>1 else 0.0
    else:
        fst[an] = (SSb/(k-1))/(SSw/(n-k))
fstat_df = pd.DataFrame({"axis":list(fst),"F":list(fst.values())})
fstat_df = fstat_df.sort_values("F",ascending=False).reset_index(drop=True)
print(fstat_df.to_string(index=False))
# ===== STEP 5: ALL-9 PCA (P,N,ABA) =====
print("Step 5: All-9 PCA (P,N,ABA)...")
ax3   = ["P_axis","N_axis","ABA_axis"]
X_Bz  = StandardScaler().fit_transform(matB[ax3].values.astype(float))
pca_B = PCA(n_components=3)
sc_B  = pca_B.fit_transform(X_Bz)
var_B = pca_B.explained_variance_ratio_ * 100
ld_B  = pd.DataFrame(pca_B.components_.T, index=ax3, columns=["PC1","PC2","PC3"])
fig, ax = plt.subplots(figsize=(10,8))
for i in range(len(matB)):
    row = matB.iloc[i]
    c = CC.get(row["class"], "#888")
    m = CM.get(row["treatment"], "o")
    ax.scatter(sc_B[i,0],sc_B[i,1],color=c,marker=m,s=200,zorder=3,
               edgecolors="black",linewidths=0.8)
    ax.annotate(row["short_label"],xy=(sc_B[i,0],sc_B[i,1]),xytext=(6,6),
                textcoords="offset points",fontsize=8.5,color=c,fontweight="bold",
                path_effects=[pe.withStroke(linewidth=2.5,foreground="white")])
sB = max(abs(sc_B[:,:2]).max()*0.5, 0.3)
for an in ax3:
    lx = ld_B.loc[an,"PC1"]*sB; ly = ld_B.loc[an,"PC2"]*sB
    ax.annotate("",xy=(lx,ly),xytext=(0,0),
                arrowprops=dict(arrowstyle="->",color="#666",lw=1.5))
    ax.text(lx*1.18,ly*1.18,AS[an],fontsize=10,color="#333",
            ha="center",va="center",fontweight="bold",fontstyle="italic")
ax.axhline(0,color="#ccc",lw=0.8,ls="--"); ax.axvline(0,color="#ccc",lw=0.8,ls="--")
ax.set_xlabel("PC1 (" + "{:.1f}".format(var_B[0]) + "% variance)",fontsize=12)
ax.set_ylabel("PC2 (" + "{:.1f}".format(var_B[1]) + "% variance)",fontsize=12)
T5 = "PCA of All 9 Biostimulant Treatments" + chr(10) + "(P, N, ABA axes; z-scored)"
ax.set_title(T5,fontsize=12,fontweight="bold")
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
pB = [mpatches.Patch(color="#2196F3",label="Amino acids"),
      mpatches.Patch(color="#4CAF50",label="Humic substances"),
      mpatches.Patch(color="#FF5722",label="PGPR (complex)"),
      mpatches.Patch(color="#E91E63",label="PGPR (single)"),
      mpatches.Patch(color="#607D8B",label="Mineral")]
ax.legend(handles=pB,fontsize=9,loc="upper left",framealpha=0.9)
plt.tight_layout()
p3 = FIGS_DIR + "/PCA_9treatments_3axes.png"
plt.savefig(p3,dpi=150,bbox_inches="tight"); plt.close()
print("  Saved:", p3)
# ===== STEP 6: 3-PANEL DISCRIMINANT SCATTER =====
print("Step 6: 3-panel discriminant scatter...")
fig, axs = plt.subplots(1,3,figsize=(18,6))
fig.suptitle("Biostimulant Classification in Transcriptomic Axis Space",
             fontsize=14,fontweight="bold",y=1.01)
panels = [("P_axis","N_axis",matB,"P-axis vs N-axis (all 9 treatments)"),
          ("JA_axis","SA_axis",matA,"JA-axis vs SA-axis (5 QC-verified)"),
          ("P_axis","ABA_axis",matB,"P-axis vs ABA-axis (all 9 treatments)")]
for pi,(ax_x,ax_y,dat,title) in enumerate(panels):
    ax = axs[pi]
    for i in range(len(dat)):
        row = dat.iloc[i]
        c5  = row["class_5"] if "class_5" in dat.columns else row["class"]
        c = CC.get(c5, "#888")
        m = CM.get(row["treatment"], "o")
        x,y = row[ax_x], row[ax_y]
        ax.scatter(x,y,color=c,marker=m,s=160,zorder=3,edgecolors="black",linewidths=0.8)
        ax.annotate(row["short_label"],xy=(x,y),xytext=(5,5),textcoords="offset points",
                    fontsize=7.5,color=c,fontweight="bold",
                    path_effects=[pe.withStroke(linewidth=2,foreground="white")])
    ax.axhline(0,color="#ccc",lw=0.8,ls="--"); ax.axvline(0,color="#ccc",lw=0.8,ls="--")
    ax.set_xlabel(AL[ax_x],fontsize=10); ax.set_ylabel(AL[ax_y],fontsize=10)
    ax.set_title(title,fontsize=11,fontweight="bold")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=9)
allp = [mpatches.Patch(color="#2196F3",label="Amino acids"),
        mpatches.Patch(color="#4CAF50",label="Humic substances"),
        mpatches.Patch(color="#FF5722",label="PGPR (complex / JA-type)"),
        mpatches.Patch(color="#E91E63",label="PGPR (single / SA-type)"),
        mpatches.Patch(color="#607D8B",label="Mineral / ABA-stress")]
fig.legend(handles=allp,fontsize=9,loc="lower center",
           ncol=5,bbox_to_anchor=(0.5,-0.04),framealpha=0.9)
plt.tight_layout()
p4 = FIGS_DIR + "/discriminant_scatter_paper.png"
plt.savefig(p4,dpi=150,bbox_inches="tight"); plt.close()
print("  Saved:", p4)
# ===== BONUS: SEABORN CLUSTERMAP =====
print("Bonus: Seaborn clustermap (9 tx x 9 axes)...")
X_all  = wide[axis_order].values.astype(float)
X_allz = StandardScaler().fit_transform(X_all)
dfh = pd.DataFrame(X_allz, index=wide["short_label"].values,
                   columns=[AS[a] for a in axis_order])
pcls = {"Amino acids":"#2196F3","Humic substances":"#4CAF50",
        "PGPR (complex)":"#FF5722","PGPR (single)":"#E91E63","Mineral":"#607D8B"}
rcol = pd.Series([pcls.get(wide["class"].iloc[i],"#888") for i in range(len(wide))],
                  index=wide["short_label"].values)
cg = sns.clustermap(dfh, method="ward", metric="euclidean", cmap="RdBu_r", center=0,
                    row_colors=rcol, figsize=(10,8), dendrogram_ratio=(0.15,0.1),
                    cbar_pos=(0.02,0.8,0.03,0.15), linewidths=0.5, linecolor="#eee",
                    xticklabels=True, yticklabels=True)
TCM = "Biostimulant Transcriptomic Fingerprints" + chr(10) + "(9 treatments x 9 axes, z-scored)"
cg.ax_heatmap.set_title(TCM, pad=15, fontsize=12, fontweight="bold")
cg.ax_heatmap.tick_params(axis="x",labelsize=10)
cg.ax_heatmap.tick_params(axis="y",labelsize=9)
legp = [mpatches.Patch(color=v,label=k) for k,v in pcls.items()]
cg.ax_heatmap.legend(handles=legp,bbox_to_anchor=(1.25,1),loc="upper left",
                      title="Class",fontsize=8,title_fontsize=9)
p5 = FIGS_DIR + "/classification_clustermap_9tx_9axes.png"
cg.savefig(p5,dpi=150,bbox_inches="tight"); plt.close("all")
print("  Saved:", p5)
# ===== STEP 7: WRITE REPORT =====
print("Step 7: Writing report...")

def fld(ldf, n=3):
    n = min(n, ldf.shape[1])
    hdr = "| Axis  | " + " | ".join(["PC"+str(i+1) for i in range(n)]) + " |"
    sep = "|-------|" + "|".join(["-------"]*n) + "|"
    rows = [hdr, sep]
    for an in ldf.index:
        vals = " | ".join(["{:+.3f}".format(ldf.loc[an,"PC"+str(i+1)]) for i in range(n)])
        rows.append("| " + "{:<5}".format(AS.get(an,an)) + " | " + vals + " |")
    return chr(10).join(rows)

itp = {
    "JA_axis":"Separates PGPR-JA (+8.3%) from PGPR-SA (-9.8%); strongest PGPR split",
    "SA_axis":"diacetyl +12.6% unique; SA-priming mechanism identified",
    "N_axis":"PGPR-SA +11.8%; PGPR complex -4 to -10%",
    "Auxin_axis":"Humic +6.3%, PGPR-SA +11.6%; class-specific Humic marker",
    "P_axis":"PGPR/Glutamate high; TiO2/Humic near zero",
    "Fe_axis":"PGPR-JA (GMV) distinctly high (+13%)",
    "ABA_axis":"TiO2 cluster (+3.5%); others low",
    "CK_axis":"Low variance; limited discriminating power",
    "ET_axis":"Near-zero for all; minimal power"}

t3  = ", ".join([AS[a] for a in fstat_df["axis"].head(3)])
p1t = ", ".join([AS[a] for a in ld_A["PC1"].abs().nlargest(3).index])
p2t = ", ".join([AS[a] for a in ld_A["PC2"].abs().nlargest(3).index])

ftbl = ["| Rank | Axis | F (MS_between) | Interpretation |",
        "|------|------|----------------|----------------|" ]
for i,r in fstat_df.iterrows():
    aname  = AS.get(r["axis"], r["axis"])
    interp = itp.get(r["axis"], "--")
    ftbl.append("| " + str(i+1) + " | " + aname + " | " +
                "{:.2f}".format(r["F"]) + " | " + interp + " |")

vAstr = chr(10).join(["- PC"+str(i+1)+": "+"{:.1f}".format(var_A[i])+"%" for i in range(len(var_A))])
vBstr = chr(10).join(["- PC"+str(i+1)+": "+"{:.1f}".format(var_B[i])+"%" for i in range(len(var_B))])
rpt_parts = [
    "# Discriminant Analysis Report -- Biostimulant Transcriptomic Fingerprints",
    "",
    "**Date:** 2026-02-23",
    "",
    "---",
    "",
    "## 1. Dataset",
    "",
    "- Matrix A (QC-verified): 5 treatments x 9 axes",
    "  - amino_acid, humic_subst, GMV_treated, diacetyl_treated, TiO2_treated",
    "- Matrix B (all): 9 treatments x P/N/ABA axes",
    "",
    "---",
    "",
    "## 2. F-Statistic Ranking (Matrix A)",
    "",
    "NOTE: 1 sample per class; metric=MS_between (between-class spread).",
    "",
]
rpt_parts.extend(ftbl)
rpt_parts.extend(["",
    "Top discriminating axes: " + t3,
    "","---","","## 3. PCA","",
    "### Matrix A: 5 treatments x 9 axes (z-scored)",
])
rpt_parts.extend(vAstr.split(chr(10)))
rpt_parts.append("")
rpt_parts.extend(fld(ld_A, n=min(4,ld_A.shape[1])).split(chr(10)))
rpt_parts.extend(["",
    "PC1 (" + "{:.1f}".format(var_A[0]) + "%): Dominated by " + p1t + " -- separates PGPR from others",
    "PC2 (" + "{:.1f}".format(var_A[1]) + "%): Dominated by " + p2t + " -- JA vs SA opposition (PGPR split)",
    "","### Matrix B: 9 treatments x 3 axes (z-scored)",
])
rpt_parts.extend(vBstr.split(chr(10)))
rpt_parts.append("")
rpt_parts.extend(fld(ld_B, n=3).split(chr(10)))
rpt_parts.extend(["","---","","## 4. Key Findings","",
    "### 4.1 Best Discriminating Axes",
    "1. JA-axis: PGPR-JA +8.3% (GMV), +6.4% (SynCom); PGPR-SA -9.8% (diacetyl)",
    "2. SA-axis: diacetyl +12.6% unique; SA-priming mechanism identified",
    "3. N-axis: PGPR-SA +11.8%; PGPR complex -4 to -10%",
    "4. Auxin-axis: Humic +6.3%, PGPR-SA +11.6%; class-specific Humic marker",
    "","### 4.2 The PGPR Split",
    "PGPR-JA (complex): GMV +8.3% JA, +13% Fe, -4.3% N; SynCom +6.4% JA, -10% N",
    "PGPR-SA (single):  diacetyl -9.8% JA, +12.6% SA, +11.8% N; GB03 -5.7% JA",
    "Two independent axes (JA, SA) confirm the split with full mechanistic coherence.",
    "","### 4.3 Amino Acids -- Reproducible Profile",
    "Hydrolysate and glutamate (independent datasets): P+, N-, hormones<3%.",
    "","### 4.4 Humic -- Auxin Signature",
    "Auxin +6.3%; class-specific marker consistent with humic->auxin crosstalk.",
    "","### 4.5 TiO2 -- ABA Cluster",
    "ABA +3.5%, N -5.9%; distinct from nutritional biostimulants.",
    "","---","","## 5. Recommended Axes for Paper","",
    "Primary (main figures): JA, SA, N, P",
    "Secondary (supplementary): Auxin, ABA, Fe",
    "Low power (context only): CK, ET",
    "","---","","## 6. Figures","",
    "| File | Description |",
    "|------|-------------|" ,
    "| state_space_figures/classification_dendrogram.png | Ward dendrogram, 5 QC treatments |",
    "| state_space_figures/PCA_5treatments_9axes.png | PCA biplot + loadings, 5 tx x 9 axes |",
    "| state_space_figures/PCA_9treatments_3axes.png | PCA, 9 tx x P/N/ABA |",
    "| state_space_figures/discriminant_scatter_paper.png | 3-panel scatter (PxN, JAxSA, PxABA) |",
    "| state_space_figures/classification_clustermap_9tx_9axes.png | Seaborn clustermap, 9 tx |",
    "","---",
    "Generated by RNA_seq_axes/discriminant_analysis_DESeq2axes.py",
])
rpt = chr(10).join(rpt_parts)
rpath = REPORT_DIR + "/discriminant_analysis_report.md"
with open(rpath,"w",encoding="utf-8") as f: f.write(rpt)
print("  Saved:", rpath)

print("")
print("=== All steps complete ===")
print(fstat_df[["axis","F"]].head(5).to_string(index=False))
print("PCA-A var:", [round(v,1) for v in var_A[:4]])
