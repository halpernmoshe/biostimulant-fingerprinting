import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import math

SCORES = 'C:/Users/moshe/Dropbox/ISF 2025/RNA_seq_axes/scores/'
FIGS   = 'C:/Users/moshe/Dropbox/ISF 2025/state_space_figures/'

df = pd.read_csv(SCORES + 'biostimulant_6axis_summary.csv')
M  = df.pivot(index='treatment', columns='axis', values='delta_pct')
print('Loaded:', M.shape)

TI = {
    'amino_acid':       {'label':'Amino acids (hydrolysate)',     'color':'#2196F3', 'marker':'o'},
    'humic_subst':      {'label':'Humic substances',              'color':'#9C27B0', 'marker':'s'},
    'GMV_treated':      {'label':'PGPR GMV (complex mix)',        'color':'#FF5722', 'marker':'D'},
    'diacetyl_treated': {'label':'Pure diacetyl (single strain)', 'color':'#E91E63', 'marker':'^'},
    'TiO2_treated':     {'label':'TiO2 nanoparticles',            'color':'#607D8B', 'marker':'P'},
}
CC = '#FF5722'
SC = '#E91E63'
ORD = ['amino_acid','humic_subst','GMV_treated','diacetyl_treated','TiO2_treated']

def qbox(ax, txt, x, y, ha, va, col, fs=9):
    ax.text(x, y, txt, transform=ax.transAxes, ha=ha, va=va,
            fontsize=fs, color=col, style='italic', fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.25', fc='white', ec=col, alpha=0.75))

LO_JA = {'amino_acid':(8,5),'humic_subst':(8,5),'GMV_treated':(8,-16),
         'diacetyl_treated':(-130,8),'TiO2_treated':(8,8)}
LO_PN = {'amino_acid':(8,5),'humic_subst':(8,5),'GMV_treated':(8,-16),
         'diacetyl_treated':(8,8),'TiO2_treated':(8,5)}

# ===== FIGURE 1: 4-panel =====
fig, axs = plt.subplots(2, 2, figsize=(13, 11))
fig.suptitle('PGPR biostimulants split into two mechanistically distinct sub-classes',
             fontsize=14, fontweight='bold', y=0.98)

# Panel A: JA vs SA scatter
ax = axs[0,0]
for t in ORD:
    if t not in M.index: continue
    info = TI[t]
    ja = M.loc[t,'JA_axis']; sa = M.loc[t,'SA_axis']
    ax.scatter(ja, sa, c=info['color'], marker=info['marker'], s=220,
               zorder=5, edgecolors='black', linewidths=0.8)
    ax.annotate(info['label'], (ja, sa),
                textcoords='offset points', xytext=LO_JA.get(t,(8,5)), fontsize=8.5)
ax.axhline(0, color='gray', ls='--', alpha=0.5, lw=1)
ax.axvline(0, color='gray', ls='--', alpha=0.5, lw=1)
ax.set_xlabel('JA-axis score (delta-% of MeJA range)', fontsize=11)
ax.set_ylabel('SA-axis score (delta-% of BTH range)',   fontsize=11)
ax.set_title('A.  JA vs SA: PGPR sub-class separation', fontweight='bold', fontsize=11)
ax.set_xlim(-14, 13); ax.set_ylim(-1, 17)
ax.fill_betweenx([0,17], -14, 0, alpha=0.06, color=SC)
ax.fill_betweenx([-1,0],  0, 13, alpha=0.06, color=CC)
qbox(ax, 'SA-ISR' + chr(10) + 'SA up, JA down', 0.02, 0.97, 'left',  'top',    SC)
qbox(ax, 'JA-ISR' + chr(10) + 'JA up, SA mod.', 0.97, 0.04, 'right', 'bottom', CC)

# Panel B: P vs N scatter
ax = axs[0,1]
for t in ORD:
    if t not in M.index: continue
    info = TI[t]
    p = M.loc[t,'P_axis']; n = M.loc[t,'N_axis']
    ax.scatter(p, n, c=info['color'], marker=info['marker'], s=220,
               zorder=5, edgecolors='black', linewidths=0.8)
    ax.annotate(info['label'], (p, n),
                textcoords='offset points', xytext=LO_PN.get(t,(8,5)), fontsize=8.5)
ax.axhline(0, color='gray', ls='--', alpha=0.5, lw=1)
ax.axvline(0, color='gray', ls='--', alpha=0.5, lw=1)
ax.set_xlabel('P-axis score (delta-% of Pi-starvation range)', fontsize=11)
ax.set_ylabel('N-axis score (delta-% of N-supply range)',       fontsize=11)
ax.set_title('B.  P vs N: Nutrient fingerprint split', fontweight='bold', fontsize=11)
ax.set_xlim(-4, 18); ax.set_ylim(-9, 15)
ax.fill_betweenx([0,15],  0, 18, alpha=0.06, color=SC)
ax.fill_betweenx([-9,0],  0, 18, alpha=0.06, color=CC)
qbox(ax, 'P up, N up' + chr(10) + 'SA-ISR',   0.97, 0.97, 'right', 'top',    SC)
qbox(ax, 'P up, N down' + chr(10) + 'JA-ISR', 0.97, 0.04, 'right', 'bottom', CC)

# Panel C: JA bar chart
ax = axs[1,0]
TRE = [t for t in ORD if t in M.index]
JV  = [M.loc[t,'JA_axis'] for t in TRE]
SV  = [M.loc[t,'SA_axis'] for t in TRE]
COL = [TI[t]['color'] for t in TRE]
x   = np.arange(len(TRE))
bars = ax.bar(x, JV, color=COL, edgecolor='black', linewidth=0.6, width=0.6)
ax.axhline(0, color='black', lw=1)
ax.set_xticks(x)
ax.set_xticklabels([TI[t]['label'] for t in TRE], fontsize=9)
ax.set_ylabel('JA-axis score (delta-% of MeJA range)', fontsize=10)
ax.set_title('C.  JA-axis: activated by GMV (JA-ISR)' + chr(10) + 'suppressed by pure diacetyl (SA-JA antagonism)',
             fontweight='bold', fontsize=10)
ymn = min(JV)*1.35-1; ymx = max(JV)*1.35+1
ax.set_ylim(ymn, ymx)
for bar,v in zip(bars,JV):
    o = 0.3 if v>=0 else -1.0
    ax.text(bar.get_x()+bar.get_width()/2, v+o, f'{v:.1f}%', ha='center', fontsize=8.5, fontweight='bold')
gi=TRE.index('GMV_treated'); di=TRE.index('diacetyl_treated')
ya=ymx*0.82
ax.annotate('', xy=(di,ya), xytext=(gi,ya),
            arrowprops=dict(arrowstyle='<->', color='black', lw=1.5))
ax.text((gi+di)/2, ya+0.35, 'PGPR split', ha='center', fontsize=8, style='italic')

# Panel D: SA bar chart
ax = axs[1,1]
bars = ax.bar(x, SV, color=COL, edgecolor='black', linewidth=0.6, width=0.6)
ax.axhline(0, color='black', lw=1)
ax.set_xticks(x)
ax.set_xticklabels([TI[t]['label'] for t in TRE], fontsize=9)
ax.set_ylabel('SA-axis score (delta-% of BTH range)', fontsize=10)
ax.set_title('D.  SA-axis: strongly activated by pure diacetyl' + chr(10) + '(SA-ISR confirmed)',
             fontweight='bold', fontsize=10)
ax.set_ylim(-0.5, max(SV)*1.35)
for bar,v in zip(bars,SV):
    o = 0.2 if v>=0 else -0.8
    ax.text(bar.get_x()+bar.get_width()/2, v+o, f'{v:.1f}%', ha='center', fontsize=8.5, fontweight='bold')

handles = [mpatches.Patch(color=TI[t]['color'], label=TI[t]['label']) for t in ORD]
fig.legend(handles=handles, loc='lower center', ncol=5,
           fontsize=9, bbox_to_anchor=(0.5,-0.01), framealpha=0.9, edgecolor='gray')
plt.tight_layout(rect=[0,0.05,1,0.97])
plt.savefig(FIGS+'PGPR_split_JA_SA_4panel.png', dpi=150, bbox_inches='tight')
print('Saved: PGPR_split_JA_SA_4panel.png')
plt.close()

# ===== FIGURE 2: 2-panel =====
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
fig.suptitle('PGPR transcriptomic fingerprints reveal two mechanistic sub-classes',
             fontsize=13, fontweight='bold')
for t in ORD:
    if t not in M.index: continue
    info = TI[t]
    ja = M.loc[t,'JA_axis']; sa = M.loc[t,'SA_axis']
    ax1.scatter(ja, sa, c=info['color'], marker=info['marker'], s=220,
                zorder=5, edgecolors='black', linewidths=0.8, label=info['label'])
    ax1.annotate(info['label'], (ja, sa),
                 textcoords='offset points', xytext=LO_JA.get(t,(8,5)), fontsize=9)
ax1.axhline(0, color='gray', ls='--', alpha=0.5, lw=1)
ax1.axvline(0, color='gray', ls='--', alpha=0.5, lw=1)
ax1.set_xlabel('JA-axis score (delta-% of MeJA range)', fontsize=11)
ax1.set_ylabel('SA-axis score (delta-% of BTH range)',   fontsize=11)
ax1.set_title('JA vs SA: PGPR sub-class separation', fontweight='bold')
ax1.set_xlim(-14, 13); ax1.set_ylim(-1, 17)
qbox(ax1, 'SA-ISR' + chr(10) + 'SA up, JA down', 0.02, 0.97, 'left',  'top',    SC)
qbox(ax1, 'JA-ISR' + chr(10) + 'JA up, SA mod.', 0.97, 0.04, 'right', 'bottom', CC)
for t in ORD:
    if t not in M.index: continue
    info = TI[t]
    p = M.loc[t,'P_axis']; n = M.loc[t,'N_axis']
    ax2.scatter(p, n, c=info['color'], marker=info['marker'], s=220,
                zorder=5, edgecolors='black', linewidths=0.8)
    ax2.annotate(info['label'], (p, n),
                 textcoords='offset points', xytext=LO_PN.get(t,(8,5)), fontsize=9)
ax2.axhline(0, color='gray', ls='--', alpha=0.5, lw=1)
ax2.axvline(0, color='gray', ls='--', alpha=0.5, lw=1)
ax2.set_xlabel('P-axis score (delta-% of Pi-starvation range)', fontsize=11)
ax2.set_ylabel('N-axis score (delta-% of N-supply range)',       fontsize=11)
ax2.set_title('P vs N: Nutrient fingerprint split', fontweight='bold')
ax2.set_xlim(-4, 18); ax2.set_ylim(-9, 15)
qbox(ax2, 'P up, N up' + chr(10) + 'SA-ISR',   0.97, 0.97, 'right', 'top',    SC)
qbox(ax2, 'P up, N down' + chr(10) + 'JA-ISR', 0.97, 0.04, 'right', 'bottom', CC)
handles2 = [mpatches.Patch(color=TI[t]['color'], label=TI[t]['label']) for t in ORD]
fig.legend(handles=handles2, loc='lower center', ncol=5,
           fontsize=9.5, bbox_to_anchor=(0.5,-0.08), framealpha=0.9, edgecolor='gray')
plt.tight_layout(rect=[0,0.08,1,0.96])
plt.savefig(FIGS+'PGPR_split_2panel.png', dpi=150, bbox_inches='tight')
print('Saved: PGPR_split_2panel.png')
plt.close()

# ===== FIGURE 3: Radar chart =====
AL   = ['P_axis','N_axis','ABA_axis','JA_axis','SA_axis','Fe_axis']
ALAB = ['P-starvation','N-supply','ABA stress','JA (defense)','SA (defense)','Fe-starvation']
AV   = [(a,l) for a,l in zip(AL,ALAB) if a in M.columns]
AL2  = [a for a,l in AV]; ALAB2 = [l for a,l in AV]
NA   = len(AL2)
ang  = [n/float(NA)*2*math.pi for n in range(NA)] + [0]
fig, ax = plt.subplots(figsize=(8,8), subplot_kw=dict(polar=True))
CR = {'amino_acid':'#2196F3','humic_subst':'#9C27B0',
      'GMV_treated':'#FF5722','diacetyl_treated':'#E91E63','TiO2_treated':'#607D8B'}
LR = {'amino_acid':'Amino acids (hydrolysate)','humic_subst':'Humic substances',
      'GMV_treated':'PGPR GMV (complex mix)','diacetyl_treated':'Pure diacetyl (single strain)',
      'TiO2_treated':'TiO2 nanoparticles'}
for t in ORD:
    if t not in M.index: continue
    vals = [M.loc[t,a] if a in M.columns else 0 for a in AL2]
    vals = vals + [vals[0]]
    ax.plot(ang, vals, 'o-', linewidth=2, color=CR[t], label=LR.get(t,t), alpha=0.85)
    ax.fill(ang, vals, alpha=0.08, color=CR[t])
ax.set_xticks(ang[:-1])
ax.set_xticklabels(ALAB2, fontsize=11)
ax.set_title('Biostimulant physiological fingerprints (6-axis radar chart)',
             fontsize=12, fontweight='bold', pad=25)
ax.legend(loc='upper right', bbox_to_anchor=(1.38,1.15), fontsize=10)
ax.yaxis.set_tick_params(labelsize=8)
plt.tight_layout()
plt.savefig(FIGS+'fingerprint_radar_5treatments.png', dpi=150, bbox_inches='tight')
print('Saved: fingerprint_radar_5treatments.png')
plt.close()
print('All 3 figures done.')
