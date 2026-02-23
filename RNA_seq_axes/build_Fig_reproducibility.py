import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec

axes_9 = ["P", "N", "ABA", "JA", "SA", "Fe", "CK", "ET", "Auxin"]
hydro = {"P": 0.07, "N": -2.0, "ABA": 2.2, "JA": 0.4, "SA": 0.8, "Fe": -0.0, "CK": 1.4, "ET": -0.2, "Auxin": 2.2}
glut  = {"P": 29.51, "N": -17.74, "ABA": 8.18, "JA": 19.88, "SA": 38.25, "Fe": 24.71, "CK": -1.46, "ET": 2.35, "Auxin": 9.16}
hydro_vals = [hydro[a] for a in axes_9]
glut_vals  = [glut[a]  for a in axes_9]
TEAL = "#1a9e8f"
GREEN = "#2ca02c"
RED = "#d62728"
ORANGE = "#ff7f0e"

fig = plt.figure(figsize=(14, 5.2))
gs = gridspec.GridSpec(1, 3, figure=fig, wspace=0.42, left=0.07, right=0.97, top=0.86, bottom=0.19)

# Panel A
ax_a = fig.add_subplot(gs[0])
x = np.arange(len(axes_9))
w = 0.35
ax_a.bar(x - w/2, hydro_vals, w, color=TEAL, label="Protein hydrolysate (GSE297649)", edgecolor="white", lw=0.5, zorder=3)
ax_a.bar(x + w/2, glut_vals,  w, color=TEAL, label="Glutamate (GSE298982)", hatch="////", edgecolor="white", lw=0.5, alpha=0.65, zorder=3)
ax_a.axhline(0, color="black", lw=0.8)
ax_a.set_xticks(x)
ax_a.set_xticklabels(axes_9, fontsize=9)
ax_a.set_ylabel("Score (% calibration range)", fontsize=9.5)
ax_a.set_ylim(-25, 55)
ax_a.set_title("A", fontweight="bold", loc="left", fontsize=13, pad=4)
ax_a.set_xlabel("Axis", fontsize=9.5)
for i, (h, g) in enumerate(zip(hydro_vals, glut_vals)):
    same = (h >= 0) == (g >= 0)
    if same and not (abs(h) < 0.15 and abs(g) < 0.15):
        ax_a.text(x[i], max(h, g) + 2.0, "*", ha="center", fontsize=11, color="#333333")
ax_a.legend(fontsize=7.5, loc="upper right", framealpha=0.9, edgecolor="#cccccc")
ax_a.spines["top"].set_visible(False)
ax_a.spines["right"].set_visible(False)
ax_a.yaxis.grid(True, alpha=0.3)
ax_a.set_axisbelow(True)
ax_a.annotate("*P-axis corrected for library-size artifact", xy=(0.01, 0.01), xycoords="axes fraction", fontsize=5.8, color="#888888", va="bottom")

# Panel B
ax_b = fig.add_subplot(gs[1])
vectors = [
    ("Protein hydrolysate",  0.07,  -2.0,   TEAL,   "-"),
    ("Glutamate",           29.51, -17.74,  TEAL,   "--"),
    ("Humic substances",    -1.37,  -1.36,  GREEN,  "-"),
    ("PGPR GMV",            11.78,  -4.32,  RED,    "-"),
    ("PGPR diacetyl",       13.72,  11.77,  ORANGE, "-"),
]
label_offsets = {
    "Protein hydrolysate": (0.5, -2.8),
    "Glutamate":           (0.8, -2.0),
    "Humic substances":    (-7.2,  0.8),
    "PGPR GMV":            (0.5,   1.8),
    "PGPR diacetyl":       (0.5,   1.8),
}
ax_b.axhline(0, color="#bbbbbb", lw=0.8, zorder=1)
ax_b.axvline(0, color="#bbbbbb", lw=0.8, zorder=1)
for (lbl, px, ny, col, ls) in vectors:
    ax_b.annotate("", xy=(px, ny), xytext=(0, 0),
                  arrowprops=dict(arrowstyle="->", color=col, lw=2.0, linestyle=ls, mutation_scale=14))
    ox, oy = label_offsets.get(lbl, (0.5, 1.5))
    ax_b.text(px + ox, ny + oy, lbl, color=col, fontsize=6.8, fontweight="bold")
ax_b.scatter([0], [0], color="black", s=40, zorder=5)
ax_b.set_xlabel("P-axis score (%)", fontsize=9.5)
ax_b.set_ylabel("N-axis score (%)", fontsize=9.5)
ax_b.set_title("B", fontweight="bold", loc="left", fontsize=13, pad=4)
ax_b.spines["top"].set_visible(False)
ax_b.spines["right"].set_visible(False)
ax_b.grid(True, alpha=0.25, zorder=0)
leg_b = [mpatches.Patch(color=c, label=l) for (l, px, ny, c, ls) in vectors]
ax_b.legend(handles=leg_b, fontsize=7.0, loc="lower right", framealpha=0.9, edgecolor="#cccccc")

# Panel C: Radar
axes_radar = ["P", "N", "ABA", "JA", "SA", "Auxin"]
hydro_r = [hydro[a] for a in axes_radar]
glut_r  = [glut[a]  for a in axes_radar]
N = len(axes_radar)
angles = np.linspace(0, 2 * np.pi, N, endpoint=False).tolist()
angles_cl = angles + [angles[0]]
h_plot = hydro_r + [hydro_r[0]]
g_plot = glut_r  + [glut_r[0]]
r_min, r_max = -15, 55
def norm(v):
    return [(xi - r_min) / (r_max - r_min) for xi in v]
h_n = norm(h_plot)
g_n = norm(g_plot)
ax_c = fig.add_subplot(gs[2], polar=True)
ax_c.set_ylim(0, 1)
for ref in [r_min, 0, 15, 30, r_max]:
    nv = (ref - r_min) / (r_max - r_min)
    ax_c.plot(angles_cl, [nv] * (N + 1), color="#dddddd", lw=0.7)
n0 = (0 - r_min) / (r_max - r_min)
ax_c.plot(angles_cl, [n0] * (N + 1), color="#888888", lw=1.0, ls="--")
for ang in angles:
    ax_c.plot([ang, ang], [0, 1], color="#cccccc", lw=0.7)
ax_c.plot(angles_cl, h_n, color=TEAL, lw=2.2, ls="-", label="Protein hydrolysate (GSE297649)")
ax_c.fill(angles_cl, h_n, color=TEAL, alpha=0.12)
ax_c.plot(angles_cl, g_n, color=TEAL, lw=2.0, ls="--", label="Glutamate (GSE298982)")
ax_c.fill(angles_cl, g_n, color=TEAL, alpha=0.06)
ax_c.set_xticks(angles)
ax_c.set_xticklabels(axes_radar, fontsize=9.5, fontweight="bold")
ax_c.set_yticks([])
ax_c.set_title("C", fontweight="bold", loc="left", fontsize=13, pad=12, x=-0.08)
ax_c.legend(fontsize=7.0, loc="upper right", bbox_to_anchor=(1.38, 0.18), framealpha=0.9, edgecolor="#cccccc")
for ref, lbl in [(r_min, str(r_min)+"%"), (0, "0%"), (30, "30%")]:
    nv = (ref - r_min) / (r_max - r_min)
    ax_c.text(angles[0] + 0.12, nv, lbl, fontsize=5.8, ha="left", va="center", color="#666666")

fig.suptitle("Amino acid biostimulants: fingerprint reproducibility across independent experiments", fontsize=10.5, fontweight="bold", y=0.985)

out = "C:/Users/moshe/Dropbox/ISF 2025/state_space_figures/Fig_reproducibility_amino_acids.png"
fig.savefig(out, dpi=150, facecolor="white")
print("Saved:", out)