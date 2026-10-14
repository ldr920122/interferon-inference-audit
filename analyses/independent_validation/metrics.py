import pandas as pd, numpy as np
from math import comb
class _B:
    def __init__(s,p): s.pvalue=p
def binomtest(k,n,alternative='greater'): return _B(sum(comb(n,i) for i in range(k,n+1))/2**n)
def holm(p):
    p=np.asarray(p,float); o=np.argsort(p); a=np.empty_like(p); run=0
    for i,j in enumerate(o): run=max(run,min(1,(len(p)-i)*p[j])); a[j]=run
    return a
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
r=pd.read_csv('tables/per_contrast_z.csv'); zc=[c for c in r if c.startswith('z_')]; cyt=[c[2:] for c in zc]
rng=np.random.default_rng(1)
def auroc(z,t):
    pos=z[t]; neg=z.drop(t); return ((neg<pos).sum()+0.5*(neg==pos).sum())/len(neg)
rows=[]
for _,x in r.iterrows():
    z=pd.Series(x[zc].values.astype(float),index=cyt)
    d=dict(x[['signature','dataset','stim','target','donor','indep']]); d.update(z_IFN1=z['IFN1'],z_IFNG=z['IFNG'],z_IFNL=z['IFNL'],top1=z.idxmax(),top3=','.join(z.nlargest(3).index))
    if x.target in ('IFN1','IFNG'):
        d.update(rank=int((z>z[x.target]).sum()+1),auroc=auroc(z,x.target),
                 specific=(z['IFN1']>z['IFNG']) if x.target=='IFN1' else (z['IFNG']>z['IFN1']))
        d['perm_auroc']=np.mean([auroc(z,rng.choice(cyt)) for _ in range(200)])
    rows.append(d)
m=pd.DataFrame(rows); m.to_csv('tables/per_contrast_metrics.csv',index=False)
S=[]
for (s,ds,t,ind),g in m[m.target.isin(['IFN1','IFNG'])].groupby(['signature','dataset','target','indep']):
    dz=g.groupby('donor')[f'z_{t}'].median()
    p=binomtest(int((dz>0).sum()),len(dz),alternative='greater').pvalue
    S.append(dict(signature=s,dataset=ds,target=t,indep=ind,n_contrasts=len(g),n_donor_units=len(dz),top1_rate=(g['rank']==1).mean(),top3_rate=(g['rank']<=3).mean(),
      median_auroc=g.auroc.median(),median_z_target=g[f'z_{t}'].median(),specificity_rate=(g.specific==True).mean(),median_margin=(g['z_'+t]-g['z_'+('IFNG' if t=='IFN1' else 'IFN1')]).median(),perm_auroc_median=g.perm_auroc.median(),sign_test_p=p))
S=pd.DataFrame(S)
# pooled per cytokine
P=[]
for (s,t),g in m[m.target.isin(['IFN1','IFNG'])].groupby(['signature','target']):
    for lab,gg in [('all',g),('full_indep_only',g[~g.indep.str.startswith('partial')])]:
        dz=gg.groupby(['dataset','donor'])[f'z_{t}'].median()
        P.append(dict(signature=s,target=t,set=lab,n_contrasts=len(gg),n_units=len(dz),top3_rate=(gg['rank']<=3).mean(),median_auroc=gg.auroc.median(),median_z=gg[f'z_{t}'].median(),specificity=(gg.specific==True).mean(),median_margin=(gg['z_'+t]-gg['z_'+('IFNG' if t=='IFN1' else 'IFN1')]).median(),
          sign_p=binomtest(int((dz>0).sum()),len(dz),alternative='greater').pvalue))
P=pd.DataFrame(P)
P['sign_p_holm']=np.nan
for s in P.signature.unique():
    i=P.signature==s; P.loc[i,'sign_p_holm']=holm(P.loc[i,'sign_p'])
P['P1_top3>=0.75']=P.top3_rate>=0.75; P['P2_AUROC>=0.80']=P.median_auroc>=0.8; P['P3_z>2&p<.05']=(P.median_z>2)&(P.sign_p_holm<0.05); P['P4_spec>=0.80']=P.specificity>=0.8
S.to_csv('tables/per_dataset_summary.csv',index=False); P.to_csv('tables/pooled_preregistered_endpoints.csv',index=False)
N=m[~m.target.isin(['IFN1','IFNG'])]
neg=N.groupby(['signature','dataset','stim','target']).agg(n=('z_IFN1','size'),frac_IFN1_z_gt2=('z_IFN1',lambda v:(v>2).mean()),frac_IFNG_z_gt2=('z_IFNG',lambda v:(v>2).mean()),median_z_IFN1=('z_IFN1','median'),median_z_IFNG=('z_IFNG','median')).reset_index()
neg.to_csv('tables/negative_controls.csv',index=False)
pd.set_option('display.width',250)
print(S.round(3).to_string()); print(P.round(4).to_string()); print(neg.round(2).to_string())
# figures
o=m[m.signature=='official']
fig,ax=plt.subplots(1,2,figsize=(12,5))
for a,t in zip(ax,['IFN1','IFNG']):
    g=o[o.target==t]; dss=sorted(g.dataset.unique())
    for i,ds in enumerate(dss):
        gg=g[g.dataset==ds]; a.scatter(gg.z_IFN1,gg.z_IFNG,label=f'{ds} (n={len(gg)})',s=25)
    lim=[min(o.z_IFN1.min(),o.z_IFNG.min()),max(o.z_IFN1.max(),o.z_IFNG.max())]; a.plot(lim,lim,'k--',lw=.7)
    a.set_xlabel('z IFN1'); a.set_ylabel('z IFNG'); a.set_title(f'{t}-stimulated contrasts (official signature)'); a.legend(fontsize=7)
plt.tight_layout(); plt.savefig('figures/fig1_IFN1_vs_IFNG_z.png',dpi=150)
fig,ax=plt.subplots(figsize=(9,4.5))
for i,(t,c) in enumerate([('IFN1','tab:blue'),('IFNG','tab:red')]):
    for j,s in enumerate(['official','reproduced']):
        g=m[(m.target==t)&(m.signature==s)]; ax.scatter(np.full(len(g),i*2+j)+rng.uniform(-.2,.2,len(g)),g['rank'],c=c,alpha=.5,s=12)
ax.set_xticks(range(4)); ax.set_xticklabels(['IFN1 official','IFN1 reproduced','IFNG official','IFNG reproduced']); ax.invert_yaxis(); ax.axhline(3.5,ls='--',c='k',lw=.7); ax.set_ylabel('rank of target among 43 cytokines')
plt.tight_layout(); plt.savefig('figures/fig2_target_rank.png',dpi=150)
fig,ax=plt.subplots(figsize=(10,4.5)); nn=neg[neg.signature=='official']
lab=nn.dataset+'|'+nn.stim; x=np.arange(len(nn)); ax.bar(x-.2,nn.median_z_IFN1,.4,label='IFN1'); ax.bar(x+.2,nn.median_z_IFNG,.4,label='IFNG'); ax.axhline(2,ls='--',c='k',lw=.7)
ax.set_xticks(x); ax.set_xticklabels(lab,rotation=75,fontsize=6); ax.set_ylabel('median z'); ax.legend(); ax.set_title('Negative / induced controls (official)')
plt.tight_layout(); plt.savefig('figures/fig3_negative_controls.png',dpi=150)
