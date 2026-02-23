import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, Ellipse
import warnings
warnings.filterwarnings("ignore")

C_BOX    = "#F0F0F0"
C_BOX_ED = "#CCCCCC"
C_ARROW  = "#555555"
C_BLUE   = "#2C7BB6"
C_ORANGE = "#D7191C"
C_AA     = "#4393C3"
C_HS     = "#4DAF4A"
C_PGPR   = "#E41A1C"
C_TIO2   = "#969696"
FONT     = "DejaVu Sans"

def rounded_box(ax, x, y, w, h,
                color=None, edgecolor=None, lw=1.2,
                radius=0.04, zorder=2):
    if color is None: color = C_BOX
    if edgecolor is None: edgecolor = C_BOX_ED
    box = FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0.0,rounding_size={radius}",
        facecolor=color, edgecolor=edgecolor, linewidth=lw, zorder=zorder,
        transform=ax.transAxes, clip_on=False
    )
    ax.add_patch(box)
    return box

def arrow_ax(ax, x0, y0, x1, y1, color=None, lw=1.5, hw=0.012, hl=0.018):
    if color is None: color = C_ARROW
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
        xycoords="axes fraction", textcoords="axes fraction",
        arrowprops=dict(arrowstyle=f"->,head_width={hw},head_length={hl}",
            color=color, lw=lw), zorder=5)

def box_title(ax, x, y, w, text, fontsize=7.5):
    ax.text(x + w/2, y, text, ha="center", va="top",
            fontsize=fontsize, fontweight="bold", fontfamily=FONT,
            transform=ax.transAxes, zorder=6)

def draw_panel_a(ax):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    BOX_H = 0.72; BOX_Y = 0.14; BOX_TOP = BOX_Y + BOX_H
    GAP = 0.025; ARR_W = 0.035; n_boxes = 4; n_gaps = 3
    total_arrow = n_gaps * (GAP*2 + ARR_W)
    BW = (0.98 - total_arrow) / n_boxes
    xs = []; cursor = 0.01
    for i in range(n_boxes):
        xs.append(cursor); cursor += BW
        if i < n_gaps: cursor += GAP*2 + ARR_W
    for x in xs: rounded_box(ax, x, BOX_Y, BW, BOX_H)
    for i in range(n_gaps):
        x0 = xs[i]+BW+GAP; x1 = x0+ARR_W
        arrow_ax(ax, x0, BOX_Y+BOX_H/2, x1, BOX_Y+BOX_H/2)
    mg = 0.015
    # Step 1
    s1x = xs[0]
    ax1 = ax.inset_axes([s1x+mg, BOX_Y+0.18, BW-2*mg, BOX_H-0.32], transform=ax.transAxes)
    ax1.bar([0,1],[3.2,1.0],color=[C_ORANGE,C_BLUE],width=0.55,edgecolor="white",linewidth=0.8)
    ax1.set_xlim(-0.6,1.6); ax1.set_ylim(0,4.0)
    ax1.set_xticks([0,1]); ax1.set_xticklabels(["-Pi","+Pi"],fontsize=6.5,fontfamily=FONT)
    ax1.set_yticks([])
    ax1.spines[["top","right","left"]].set_visible(False)
    ax1.spines["bottom"].set_linewidth(0.7)
    ax1.tick_params(bottom=False,length=0)
    ax1.text(0.5,0.97,"logFC signature",ha="center",va="top",fontsize=5.8,
             fontstyle="italic",fontfamily=FONT,transform=ax1.transAxes,color="#444444")
    box_title(ax,s1x,BOX_TOP-0.01,BW,"1. Reference experiment",fontsize=7.2)
    ax.text(s1x+BW/2,BOX_Y+0.04,"DESeq2 on public RNA-seq",ha="center",va="bottom",
            fontsize=5.8,fontfamily=FONT,color="#666666",transform=ax.transAxes)
    # Step 2
    s2x = xs[1]
    ax2 = ax.inset_axes([s2x+mg, BOX_Y+0.28, BW-2*mg, 0.18], transform=ax.transAxes)
    ax2.set_xlim(-0.1,1.1); ax2.set_ylim(-0.5,0.5)
    ax2.axhline(0,color="#888888",lw=1.2,zorder=1)
    ax2.scatter([0,1],[0,0],s=60,c=[C_ORANGE,C_BLUE],zorder=3,edgecolors="white",linewidths=0.8)
    ax2.set_xticks([0,1]); ax2.set_xticklabels(["0%","100%"],fontsize=6.5,fontfamily=FONT)
    ax2.set_yticks([])
    ax2.spines[["top","right","left","bottom"]].set_visible(False)
    ax2.tick_params(bottom=False,length=0)
    ax2.text(0.0,-0.48,"-Pi state",ha="center",va="top",fontsize=5.5,fontfamily=FONT,color=C_ORANGE)
    ax2.text(1.0,-0.48,"+Pi state",ha="center",va="top",fontsize=5.5,fontfamily=FONT,color=C_BLUE)
    box_title(ax,s2x,BOX_TOP-0.01,BW,"2. Axis calibration",fontsize=7.2)
    ax.text(s2x+BW/2,BOX_Y+0.04,"% calibration range",ha="center",va="bottom",
            fontsize=5.8,fontfamily=FONT,color="#666666",transform=ax.transAxes)
    # Step 3
    s3x = xs[2]
    ax3 = ax.inset_axes([s3x+mg, BOX_Y+0.16, BW-2*mg, BOX_H-0.28], transform=ax.transAxes)
    ax3.bar([0,1,2],[2.4,-1.4,8.3],color=[C_AA,C_HS,C_PGPR],width=0.55,edgecolor="white",linewidth=0.8)
    ax3.axhline(0,color="#888888",lw=0.8,linestyle="--")
    ax3.set_xlim(-0.6,2.6); ax3.set_ylim(-4,12)
    ax3.set_xticks([0,1,2]); ax3.set_xticklabels(["AA","HS","PGPR"],fontsize=6.5,fontfamily=FONT)
    ax3.set_yticks([0,5,10]); ax3.set_yticklabels(["0%","5%","10%"],fontsize=5.5,fontfamily=FONT)
    ax3.spines[["top","right"]].set_visible(False)
    ax3.spines[["left","bottom"]].set_linewidth(0.7)
    ax3.tick_params(length=0)
    ax3.text(0.5,0.97,"P-axis score",ha="center",va="top",fontsize=5.8,
             fontstyle="italic",fontfamily=FONT,transform=ax3.transAxes,color="#444444")
    box_title(ax,s3x,BOX_TOP-0.01,BW,"3. Score biostimulants",fontsize=7.2)
    # Step 4 radar
    s4x = xs[3]
    rcx = s4x + BW/2; rcy = BOX_Y + BOX_H/2 - 0.01
    rr = min(BW,BOX_H)*0.31
    n_ax = 5; lbs = ["P","N","JA","SA","ABA"]
    angs = np.linspace(0,2*np.pi,n_ax,endpoint=False)
    fp1 = np.array([0.12,0.04,0.83,0.34,0.05])
    fp2 = np.array([0.14,0.12,0.02,1.00,0.29])
    for frac in [0.33,0.67,1.0]:
        cpts = np.linspace(0,2*np.pi,100)
        ax.plot(rcx+rr*frac*np.sin(cpts),rcy+rr*frac*np.cos(cpts),
                color="#DDDDDD",lw=0.6,transform=ax.transAxes,zorder=3)
    for a in angs:
        ax.plot([rcx,rcx+rr*np.sin(a)],[rcy,rcy+rr*np.cos(a)],
                color="#BBBBBB",lw=0.6,transform=ax.transAxes,zorder=3)
    for fp,col,alp in [(fp1,C_PGPR,0.25),(fp2,"#FF7F00",0.25)]:
        px=[rcx+rr*fp[i]*np.sin(angs[i]) for i in range(n_ax)]
        py=[rcy+rr*fp[i]*np.cos(angs[i]) for i in range(n_ax)]
        px.append(px[0]); py.append(py[0])
        ax.fill(px,py,color=col,alpha=alp,transform=ax.transAxes,zorder=4)
        ax.plot(px,py,color=col,lw=1.0,transform=ax.transAxes,zorder=5)
    for a,lb in zip(angs,lbs):
        ax.text(rcx+rr*1.18*np.sin(a),rcy+rr*1.18*np.cos(a),lb,
                ha="center",va="center",fontsize=6.0,fontweight="bold",
                fontfamily=FONT,transform=ax.transAxes,zorder=6)
    box_title(ax,s4x,BOX_TOP-0.01,BW,"4. Fingerprint",fontsize=7.2)
    lx=s4x+BW*0.05; ly=BOX_Y+0.06
    for col,lab in [(C_PGPR,"JA-ISR"),("#FF7F00","SA-ISR")]:
        ax.plot([lx,lx+0.035],[ly,ly],color=col,lw=2,transform=ax.transAxes,zorder=6)
        ax.text(lx+0.042,ly,lab,va="center",fontsize=5.5,fontfamily=FONT,
                transform=ax.transAxes,zorder=6)
        ly+=0.055
    ax.text(0.01,0.98,"(a)",ha="left",va="top",fontsize=11,fontweight="bold",
            fontfamily=FONT,transform=ax.transAxes)

def draw_panel_b(ax, standalone=False):
    data = [
        ("GMV",     8.3,  3.4, C_PGPR, "o", 9),
        ("Diacetyl",-9.8,12.6, C_PGPR, "D", 8),
        ("AA",      0.4,  0.8, C_AA,   "o", 9),
        ("HS",      2.7,  0.6, C_HS,   "o", 9),
        ("TiO2",    3.2,  3.2, C_TIO2, "o", 8),
    ]
    ax.set_xlim(-15,15); ax.set_ylim(-5,20)
    ax.set_xlabel("JA-axis score (% calibration range)",fontsize=9,fontfamily=FONT,labelpad=4)
    ax.set_ylabel("SA-axis score (% calibration range)",fontsize=9,fontfamily=FONT,labelpad=4)
    ax.axhline(0,color="#CCCCCC",lw=0.8,zorder=1)
    ax.axvline(0,color="#CCCCCC",lw=0.8,zorder=1)
    ax.set_facecolor("white")
    for sp in ax.spines.values(): sp.set_linewidth(0.8); sp.set_color("#AAAAAA")
    ax.tick_params(labelsize=8,length=3,color="#AAAAAA")
    ja_ell=Ellipse((7.5,3.8),width=9,height=6,angle=20,facecolor="#FFDDDD",
                  edgecolor=C_PGPR,linestyle="--",linewidth=1.0,alpha=0.55,zorder=2)
    ax.add_patch(ja_ell)
    ax.text(12.0,0.8,"JA-ISR cluster",ha="center",va="center",fontsize=7.5,
            color=C_PGPR,fontfamily=FONT,fontstyle="italic",zorder=6)
    sa_ell=Ellipse((-9.0,12.8),width=9,height=6,angle=-15,facecolor="#FFE8CC",
                  edgecolor="#FF7F00",linestyle="--",linewidth=1.0,alpha=0.55,zorder=2)
    ax.add_patch(sa_ell)
    ax.text(-12.5,16.5,"SA-ISR cluster",ha="center",va="center",fontsize=7.5,
            color="#D06000",fontfamily=FONT,fontstyle="italic",zorder=6)
    dx=np.array([-14,14]); dy=-0.7*dx+2.5
    ax.plot(dx,dy,color="#BBBBBB",lw=1.0,linestyle=(0,(6,4)),zorder=2)
    ax.text(-1.5,4.0,"SA-JA antagonism",rotation=-28,ha="center",va="center",
            fontsize=7,color="#999999",fontfamily=FONT,fontstyle="italic",zorder=5)
    offs = {"GMV":(1.4,0.8),"Diacetyl":(-1.4,1.0),"AA":(1.4,0.6),"HS":(-0.5,-1.5),"TiO2":(1.4,-1.1)}
    for nm,ja,sa,col,mk,ms in data:
        ax.scatter(ja,sa,s=ms**2*0.9,c=col,marker=mk,edgecolors="white",linewidths=0.8,zorder=7)
        ddx,ddy=offs.get(nm,(1.0,0.5))
        ax.annotate(nm,xy=(ja,sa),xytext=(ja+ddx,sa+ddy),fontsize=8,fontfamily=FONT,
                    color=col,fontweight="bold",
                    arrowprops=dict(arrowstyle="-",color=col,lw=0.6,shrinkA=0,shrinkB=4),zorder=8)
    legend_items=[
        mpatches.Patch(facecolor=C_PGPR,edgecolor="white",label="PGPR"),
        mpatches.Patch(facecolor=C_AA,edgecolor="white",label="Amino acids"),
        mpatches.Patch(facecolor=C_HS,edgecolor="white",label="Humic subst."),
        mpatches.Patch(facecolor=C_TIO2,edgecolor="white",label="TiO2 (control)"),
    ]
    ax.legend(handles=legend_items,loc="lower right",fontsize=7.5,framealpha=0.9,
              edgecolor="#CCCCCC",handlelength=1.2,handleheight=0.9,borderpad=0.6,
              prop={"family":FONT,"size":7.5})
    ax.set_title("PGPR biostimulants reveal two distinct ISR modes",
                 fontsize=10,fontfamily=FONT,fontweight="bold",pad=8)
    ax.text(-0.09,1.03,"(b)",ha="left",va="bottom",fontsize=11,fontweight="bold",
            fontfamily=FONT,transform=ax.transAxes)
    ax.text(0.3,-4.5,"JA = 0",fontsize=6.5,color="#AAAAAA",fontfamily=FONT,ha="center")
    ax.text(-14.5,0.4,"SA = 0",fontsize=6.5,color="#AAAAAA",fontfamily=FONT,ha="left")

OUT="C:/Users/moshe/Dropbox/ISF 2025/state_space_figures"
fig,axes=plt.subplots(1,2,figsize=(14,6),gridspec_kw={"width_ratios":[1.15,1.0]})
fig.patch.set_facecolor("white")
draw_panel_a(axes[0])
draw_panel_b(axes[1])
plt.tight_layout(pad=1.2,w_pad=2.0)
p1=f"{OUT}/Fig1_framework_schematic.png"
fig.savefig(p1,dpi=300,bbox_inches="tight",facecolor="white",edgecolor="none")
plt.close(fig); print("Saved:",p1)
fig2,ax2=plt.subplots(figsize=(7,6))
fig2.patch.set_facecolor("white")
draw_panel_b(ax2,standalone=True)
plt.tight_layout(pad=1.2)
p2=f"{OUT}/Fig1b_PGPR_split_scatter.png"
fig2.savefig(p2,dpi=300,bbox_inches="tight",facecolor="white",edgecolor="none")
plt.close(fig2); print("Saved:",p2)
print("Done.")
