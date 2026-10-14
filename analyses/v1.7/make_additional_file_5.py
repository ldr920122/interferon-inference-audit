# v1.6: regenerate Additional file 5 (recall; reproduced vs official signature, 1,306 vs 1,307 columns)
import pandas as pd,numpy as np,json,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
off=pd.read_csv('/workspace/cytosig_repro/CytoSig/CytoSig/signature.centroid',sep='\t',index_col=0)
def diffs(f):
    c=pd.read_csv(f,sep='\t',index_col=0); g=off.index.intersection(c.index)
    return pd.Series({k:float((c.loc[g,k]-off.loc[g,k]).abs().max()) for k in off.columns})
a=diffs('diff.centroid_1306'); b=diffs('diff.centroid_1307')
rec=[('Project_known_IFN_sources',10,10),('Aybey_IFN_stim_sets',9,9),('TxCyto_SI2_IFN',27,29),('CytoSig_IFNB_or_IFNG',116,124),('CytoSig_IFN_any',160,174)]
with PdfPages('../../manuscript/additional_files/Additional_file_5_recall_and_CytoSig_reproduction.pdf') as pdf:
    f,ax=plt.subplots(figsize=(8,3.5))
    for i,(n,k,m) in enumerate(rec[::-1]): ax.barh(i,100*k/m,color='#4C72B0'); ax.text(100*k/m+0.3,i,f'{k}/{m} ({100*k/m:.1f}%)',va='center',fontsize=8)
    ax.set_yticks(range(5)); ax.set_yticklabels([r[0] for r in rec[::-1]],fontsize=8); ax.set_xlim(80,108); ax.set_xlabel('Recall (%)')
    ax.set_title('Additional file 5a | Search recall against five reference sets',fontsize=10); f.tight_layout(); pdf.savefig(f); plt.close(f)
    f,ax=plt.subplots(figsize=(8,9)); o=b.sort_values().index; y=np.arange(len(o)); fl=1e-17
    ax.scatter(a[o].clip(lower=fl),y,marker='o',facecolors='none',edgecolors='grey',label='1,306 reconstructed columns (v1.4)')
    ax.scatter(b[o].clip(lower=fl),y,marker='x',color='#C44E52',label='1,307 resolved columns (+GSE13887 NO)')
    ax.set_xscale('log'); ax.axvline(1e-6,ls='--',c='k',lw=0.7); ax.set_yticks(y); ax.set_yticklabels(o,fontsize=7)
    ax.set_xlabel('Max |official − reproduced| over 4,881 genes (zeros plotted at 1e-17)')
    ax.set_title(f'Additional file 5b | Reproduced vs official CytoSig signature\n1,306 columns: {int((a<1e-6).sum())}/43 identical within 1e-6; 1,307 columns: {int((b<1e-6).sum())}/43 (max {b.max():.1e})',fontsize=9)
    ax.legend(fontsize=7,loc='lower right'); f.tight_layout(); pdf.savefig(f); plt.close(f)
pd.DataFrame({'max_abs_diff_1306':a,'max_abs_diff_1307':b}).to_csv('centroid_1306_vs_1307_vs_official.csv')
print((a<1e-6).sum(),(b<1e-6).sum(),a.idxmax(),a.max())
