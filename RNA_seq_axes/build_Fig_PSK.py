import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec

BLUE   = "#1f77b4"
RED    = "#d62728"
GREEN  = "#2ca02c"
PURPLE = "#7b2d8b"
GRAY   = "#888888"
ORANGE = "#ff7f0e"

fig = plt.figure(figsize=(12, 5.2))
gs = gridspec.GridSpec(1, 2, figure=fig, wspace=0.42, left=0.07, right=0.97, top=0.86, bottom=0.14)

axes_radar = ["P", "N", "ABA", "JA", "SA", "Auxin"]
N = len(axes_radar)
amino = dict(P=0.07, N=-2.0, ABA=2.2, JA=0.4, SA=0.8, Auxin=2.2)
psk_d = dict(P=-5.69, N=-3.86, ABA=14.37, JA=4.07, SA=0.29, Auxin=1.78)
humic = dict(P=-1.37, N=-1.36, ABA=2.65, JA=2.71, SA=0.55, Auxin=6.33)
r_min, r_max = -20, 20
def norm(v):
    return [(xi-r_min)/(r_max-r_min) for xi in v]
angles = np.linspace(0, 2*np.pi, N, endpoint=False).tolist()
ac = angles + [angles[0]]
def mc(d):
    v = [d[a] for a in axes_radar]; return v + [v[0]]
a_n = norm(mc(amino))
p_n = norm(mc(psk_d))
h_n = norm(mc(humic))

axA = fig.add_subplot(gs[0], polar=True)
axA.set_ylim(0, 1)
for ref in [r_min, -10, 0, 10, r_max]:
    nv = (ref-r_min)/(r_max-r_min)
    axA.plot(ac, [nv]*(N+1), color="#dddddd", lw=0.7)
n0 = (0-r_min)/(r_max-r_min)
axA.plot(ac, [n0]*(N+1), color="#888888", lw=1.0, ls="--")
for ang in angles: axA.plot([ang,ang],[0,1],color="#cccccc",lw=0.7)
axA.plot(ac, a_n, color=BLUE,   lw=2.2, ls="-",  label="Amino acids (hydrolysate)")
axA.fill(ac, a_n, color=BLUE,   alpha=0.10)
axA.plot(ac, p_n, color=RED,    lw=2.5, ls="-",  label="PSK (mean)")
axA.fill(ac, p_n, color=RED,    alpha=0.12)
axA.plot(ac, h_n, color=GREEN,  lw=2.0, ls="--", label="Humic substances")
axA.fill(ac, h_n, color=GREEN,  alpha=0.08)
axA.set_xticks(angles)
axA.set_xticklabels(axes_radar, fontsize=10, fontweight="bold")
axA.set_yticks([])
for ref,lbl in [(r_min,str(r_min)+"%"),(0,"0%"),(r_max,str(r_max)+"%")]:
    nv=(ref-r_min)/(r_max-r_min)
    axA.text(angles[0]+0.12, nv, lbl, fontsize=5.8, ha="left", va="center", color="#666666")
axA.legend(fontsize=7.5, loc="upper right", bbox_to_anchor=(1.5, 0.22), framealpha=0.9, edgecolor="#cccccc")
axA.set_title("A  PSK clusters with ABA/stress, not amino acids", fontweight="bold", fontsize=8.5, pad=22, loc="left", x=-0.05)

axB = fig.add_subplot(gs[1])
pts = [
    ("Amino acids",         0.07,   2.2,  BLUE,   "o", 80),
    ("Glutamate",          29.51,   8.18, BLUE,   "^", 80),
    ("Humic substances",   -1.37,   2.65, GREEN,  "o", 80),
    ("PGPR GMV",           11.78,  -0.50, ORANGE, "o", 80),
    ("PGPR diacetyl",      13.72,  -2.86, ORANGE, "D", 80),
    ("TiO2",               -1.46,   3.46, GRAY,   "o", 80),
    ("PSK (mean)",          -5.69,  14.37, PURPLE, "*",200),
]
axB.axhline(0, color="#cccccc", lw=0.8)
axB.axvline(0, color="#cccccc", lw=0.8)
r1=patches.FancyBboxPatch((-13,1.2),16,14.5,boxstyle="round,pad=0.5",lw=1.2,edgecolor=RED,facecolor=RED,alpha=0.07)
axB.add_patch(r1)
axB.text(-12, 15.3, "ABA/stress cluster", fontsize=7.5, color=RED, fontweight="bold")
r2=patches.FancyBboxPatch((9,-4.5),24,14.0,boxstyle="round,pad=0.5",lw=1.2,edgecolor=ORANGE,facecolor=ORANGE,alpha=0.07)
axB.add_patch(r2)
axB.text(10, 10.3, "P-acquisition cluster", fontsize=7.5, color=ORANGE, fontweight="bold")
loff = {
    "Amino acids":        (0.8, 0.5),
    "Glutamate":          (0.8, 0.5),
    "Humic substances":   (-0.5, -1.8),
    "PGPR GMV":           (0.5, -1.8),
    "PGPR diacetyl":      (0.5, -1.8),
    "TiO2":               (0.5,  0.5),
    "PSK (mean)":         (0.8,  0.5),
}
for (lbl,P,ABA,col,mk,sz) in pts:
    axB.scatter([P],[ABA],color=col,marker=mk,s=sz,zorder=5,edgecolors="white",linewidths=0.5)
    ox,oy = loff.get(lbl,(0.5,0.5))
    axB.text(P+ox, ABA+oy, lbl, fontsize=7.0, color=col, fontweight="bold")
axB.annotate("Origin: amino acid derivative", xy=(-5.69,14.37), xytext=(-13,8.5),
    fontsize=7.0, color=PURPLE, fontweight="bold",
    arrowprops=dict(arrowstyle="->", color=PURPLE, lw=1.2))
axB.set_xlabel("P-axis score (%)", fontsize=10)
axB.set_ylabel("ABA-axis score (%)", fontsize=10)
axB.set_xlim(-15, 37)
axB.set_ylim(-8, 20)
axB.spines["top"].set_visible(False)
axB.spines["right"].set_visible(False)
axB.grid(True, alpha=0.20)
axB.set_title("B  PSK fingerprint: ABA/stress cluster (not amino acid cluster)", fontweight="bold", fontsize=8.5, pad=4, loc="left")

fig.suptitle("PSK (phytosulfokine): origin-based classification fails transcriptomic reality", fontsize=10.5, fontweight="bold", y=0.985)
out = "C:/Users/moshe/Dropbox/ISF 2025/state_space_figures/Fig_PSK_origin_classification_failure.png"
fig.savefig(out, dpi=150, facecolor="white")
print("Saved:", out)