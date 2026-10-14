# Figure 6 (v1.9): post hoc donor-paired cases GSE236156 (MDM) and GSE237879 (pDC). Reads saved per-donor outputs only.
from pathlib import Path
import argparse
import pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'pdf.fonttype':42,'font.size':7,'axes.spines.top':False,'axes.spines.right':False})
D=pd.read_csv(Path(__file__).resolve().parent.parent/'analyses/v1.10_replay/outputs/per_donor_primary.csv'); a=D[D.dataset=='GSE236156'].copy(); a=a[a.stim.isin(['IFNb','IFNg'])].replace({'IFNb':'IFN-β','IFNg':'IFN-γ'})
b=D[D.dataset=='GSE237879'].copy(); b=b[b.stim.isin(['beta','gamma'])].replace({'beta':'IFN-β','gamma':'IFN-γ'})
W=170/25.4; f,ax=plt.subplots(1,2,figsize=(W,62/25.4),sharey=False)
col={'IFN-β':'#0072B2','IFN-γ':'#D55E00'}; xs={'CPM':0,'CP10k':1}
for k,(df,ttl) in enumerate([(a,'a  GSE236156: monocyte-derived macrophages, 5 donors'),(b,'b  GSE237879: plasmacytoid DCs, 4 donors')]):
    A=ax[k]; A.axhline(0,color='k',lw=0.6)
    for j,stim in enumerate(['IFN-β','IFN-γ']):
        for d in sorted(df.donor.unique()):
            s=df[(df.stim==stim)&(df.donor==d)].set_index('scale').IFNG_minus_IFN1
            x=[xs['CPM']+3*j,xs['CP10k']+3*j]; A.plot(x,[s['CPM'],s['CP10k']],'-o',color=col[stim],ms=3,lw=0.8,alpha=0.85)
        for sc in ['CPM','CP10k']:
            v=df[(df.stim==stim)&(df.scale==sc)].IFNG_minus_IFN1; ok=int(((v>0) if stim=='IFN-γ' else (v<0)).sum())
            A.text(xs[sc]+3*j,A.get_ylim()[1] if False else max(df.IFNG_minus_IFN1)+3,f'{ok}/{len(v)}',ha='center',fontsize=6,color=col[stim])
    A.set_xticks([0,1,3,4]); A.set_xticklabels(['CPM*','CP10k']*2,fontsize=5.5)
    A.text(0.5,-0.22,'IFN-β stimulation',transform=A.get_xaxis_transform(),ha='center',color=col['IFN-β'])
    A.text(3.5,-0.22,'IFN-γ stimulation',transform=A.get_xaxis_transform(),ha='center',color=col['IFN-γ'])
    A.set_xlim(-0.6,4.6); A.set_ylim(min(df.IFNG_minus_IFN1)-3,max(df.IFNG_minus_IFN1)+6)
    A.set_ylabel('CytoSig z(IFNG) − z(IFN1)\n(>0 favours IFN-γ; <0 favours type I)'); A.set_title(ttl,loc='left',fontsize=7)
f.tight_layout(rect=(0,0,1,1)); p=argparse.ArgumentParser();p.add_argument('--out',default=str(Path(__file__).resolve().parent));z=p.parse_args();Path(z.out).mkdir(parents=True,exist_ok=True);O=str(Path(z.out).resolve()/'Fig6_posthoc_donor_paired_cases')
f.savefig(O+'.pdf'); f.savefig(O+'.png',dpi=300)
