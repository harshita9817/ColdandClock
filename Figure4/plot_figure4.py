#!/usr/bin/env python3
"""Figure 4: raffinose pathway induction and raffinose accumulation.

A: diurnal-adjusted expression change 27 h after dawn- and dusk-onset chilling
   for four raffinose-pathway enzymes (data/generatable/
   diurnal_adjusted_response_categories.csv.gz). Values plotted are written to
   Figure4_pathway_expression.csv.
B: raffinose change 27 h after each onset, from Figure4_raffinose_model.csv
   (run fit_figure4_raffinose.R first). Stars: * p<0.05, ** p<0.01,
   *** p<0.001; chilled vs control two-sided, dusk vs dawn one-sided.
"""
from pathlib import Path
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, pandas as pd
from matplotlib import font_manager as fm
from matplotlib.patches import Patch, Polygon
FIGDIR=Path(__file__).resolve().parent; ROOT=FIGDIR.parent.parent
for f in ['NimbusSans-Regular.otf','NimbusSans-Bold.otf']:
    p=Path('/usr/share/fonts/opentype/urw-base35')/f
    if p.exists(): fm.fontManager.addfont(str(p))
_inst={f.name for f in fm.fontManager.ttflist}
FONT=next((f for f in ['Nimbus Sans','Helvetica','Arial'] if f in _inst),'DejaVu Sans')
plt.rcParams.update({'font.family':FONT,'font.size':7,'axes.linewidth':0.8,'xtick.major.width':0.8,'ytick.major.width':0.8,'xtick.major.size':2.5,'ytick.major.size':2.5,'xtick.major.pad':1.5,'ytick.major.pad':1.5,'axes.labelpad':2})
DAWN,DUSK,SHARED,FAINT,INK='#D9822B','#7B3294','#1F8A70','#B5DDD1','#000000'
FS=7;FW,FH=6.5,6.5/3
c=pd.read_csv(ROOT/'data/generatable/diurnal_adjusted_response_categories.csv.gz',low_memory=False);c=c[c.hours_after_onset==27]
pct=lambda g,col:(2**c[c.gene_id==g][col].iloc[0]-1)*100
lab=lambda v:f'{round(v):+,d}%'.replace('-','−')
fig=plt.figure(figsize=(FW,FH))
ax=fig.add_axes([0,0,4.333/FW,1]);ax.set_xlim(0,4.333);ax.set_ylim(0,FH);ax.axis('off')
Y=1.67
for x,t in [(0.31,'Glucose-\n6-phosphate'),(1.31,'myo-Inositol-\n3-phosphate'),(2.31,'myo-Inositol'),(3.25,'Galactinol'),(4.09,'Raffinose')]:
    ax.text(x,Y,t,ha='center',va='center',fontsize=FS,color=INK,fontweight='bold',linespacing=1.0)
E=[(0.63,0.97,'myo-inositol-3-\nphosphate synthase',DUSK,'Sobic.001G472800',0.86,False),(1.64,1.98,'inositol\nmonophosphatase',DUSK,'Sobic.009G048900',1.84,False),(2.63,2.97,'galactinol\nsynthase',SHARED,'Sobic.002G423600',2.80,True),(3.51,3.85,'raffinose\nsynthase',FAINT,'Sobic.003G052300',3.68,False)]
base=0.19;k=0.0009;BR_LO=650;BR_HI=1500;gap=0.055
ht=lambda v,brk:v*k if (v<=BR_LO or not brk) else BR_LO*k+gap+(v-BR_HI)*k
rows=[]
for a,b,nm,col,g,m,brk in E:
    ax.annotate('',xy=(b,Y),xytext=(a,Y),arrowprops=dict(arrowstyle='-|>,head_width=0.22,head_length=0.32',lw=2.4,color=col,shrinkA=0,shrinkB=0))
    tx=(a+b)/2
    ax.text(tx,Y-0.13,nm,ha='center',va='top',fontsize=FS,color=INK,linespacing=1.0)
    ax.plot([m-.41,m+.41],[base,base],color=INK,lw=.8,zorder=5,solid_capstyle='butt')
    for off,cc,cl in [(-.2,'mean_adjusted_dawn_response_log2',DAWN),(.2,'mean_adjusted_dusk_response_log2',DUSK)]:
        v=pct(g,cc);H=ht(v,brk);ax.bar(m+off,H,bottom=base,width=.34,color=cl,lw=0,zorder=2)
        if brk:
            y0=base+BR_LO*k;x0=m+off-.19;x1=m+off+.19
            ax.add_patch(Polygon([(x0,y0-.012),(x1,y0+.012),(x1,y0+gap+.012),(x0,y0+gap-.012)],closed=True,fc='white',ec='none',zorder=3))
            ax.plot([x0,x1],[y0-.012,y0+.012],color=INK,lw=.8,zorder=4);ax.plot([x0,x1],[y0+gap-.012,y0+gap+.012],color=INK,lw=.8,zorder=4)
        rows.append({'gene_id':g,'enzyme':nm.replace(chr(10),' ').replace('- ','-'),'onset':cc.split('_')[2],'percent_change':v})
        ax.text(m+off,base+H+(.025 if v>=0 else -.025),lab(v),ha='center',va='bottom' if v>=0 else 'top',fontsize=FS,color=INK)
# y-axis on the leftmost expression plot only
xa=0.43
ax.plot([xa,xa],[base-50*k,base+600*k],color=INK,lw=.8,solid_capstyle='butt')
for v in [0,200,400,600]:
    ax.plot([xa-.035,xa],[base+v*k]*2,color=INK,lw=.8);ax.text(xa-.05,base+v*k,f'{v}',ha='right',va='center',fontsize=FS,color=INK)
ax.text(0.115,base+275*k,'Expression\nchange (%)',rotation=90,ha='center',va='center',fontsize=FS,color=INK,linespacing=1.0)
for x,t in [(2.80,'UDP-galactose'),(3.68,'Sucrose')]:
    ax.text(x,2.03,t,ha='center',va='center',fontsize=FS,color=INK,fontweight='bold')
    ax.annotate('',xy=(x,Y+.05),xytext=(x,1.96),arrowprops=dict(arrowstyle='-|>,head_width=0.14,head_length=0.22',lw=.9,color=INK,shrinkA=0,shrinkB=0))
ax.legend(handles=[Patch(fc=DAWN,label='Dawn onset'),Patch(fc=DUSK,label='Dusk onset')],loc='center left',bbox_to_anchor=(0.27/4.333,2.05/FH),ncol=2,frameon=False,fontsize=FS,handlelength=1.0,handleheight=.8,handletextpad=.4,columnspacing=.9,borderpad=0,borderaxespad=0,labelcolor=INK)
fig.text(0.012/FW,(FH-0.02)/FH,'A',ha='left',va='top',fontsize=9,fontweight='bold')
fig.text(4.42/FW,(FH-0.02)/FH,'B',ha='left',va='top',fontsize=9,fontweight='bold')
# Panel B: robust model estimates at 27 hours
pd.DataFrame(rows).to_csv(FIGDIR/'Figure4_pathway_expression.csv',index=False)
mod=pd.read_csv(FIGDIR/'Figure4_raffinose_model.csv').set_index('contrast')
star=lambda p:'***' if p<0.001 else '**' if p<0.01 else '*' if p<0.05 else 'n.s.'
v=[mod.loc['dawn_onset','percent_change'],mod.loc['dusk_onset','percent_change']]
s=[star(mod.loc['dawn_onset','p_two_sided']),star(mod.loc['dusk_onset','p_two_sided'])]
sdiff=star(mod.loc['dusk_minus_dawn','p_one_sided_greater'])
b=fig.add_axes([4.92/FW,0.19/FH,1.5/FW,1.72/FH])
b.bar([0,1],v,width=.72,color=[DAWN,DUSK],lw=0)
for x,val,st in zip([0,1],v,s):
    b.text(x,val+4,lab(val),ha='center',va='bottom',fontsize=FS,color=INK)
    b.text(x,val+21,st,ha='center',va='bottom',fontsize=FS,color=INK)
b.plot([0,0,1,1],[222,228,228,222],color=INK,lw=.8);b.text(.5,229,sdiff,ha='center',va='bottom',fontsize=FS,color=INK)
b.set_xlim(-.6,1.6);b.set_ylim(0,250);b.set_xticks([]);b.set_yticks([0,50,100,150,200]);b.set_ylabel('Raffinose change (%)',fontsize=FS,color=INK);b.tick_params(labelsize=FS)
for sp in ['top','right']:b.spines[sp].set_visible(False)
fig.savefig(FIGDIR/'Figure4.png',dpi=300,facecolor='white');fig.savefig(FIGDIR/'Figure4.svg',facecolor='white')
print('Figure 4 written to',FIGDIR)
