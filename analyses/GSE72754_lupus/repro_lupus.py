"""Reproduce Jiang 2021 Nat Methods Fig.4b (GSE72754, IFN-K vaccine in SLE). See 00_prestated_tests.md.
CytoSig ridge (alpha=1e4, flag_normalize as in CytoSig Util.py) re-implemented in numpy, with 1000 gene-label
permutations for z-scores (the C/GSL extension is not available on this host)."""
import gzip,re,json,numpy as np,pandas as pd
from scipy import stats
rng=np.random.default_rng(0)
# ---- metadata
L=[l.rstrip('\n').split('\t') for l in gzip.open('GSE72754_series_matrix.txt.gz','rt')]
gsm=[x.strip('"') for x in next(l for l in L if l[0]=='!Sample_geo_accession')[1:]]
title=[x.strip('"') for x in next(l for l in L if l[0]=='!Sample_title')[1:]]
meta=pd.DataFrame({'gsm':gsm,'title':title})
for l in L:
    if l[0]=='!Sample_characteristics_ch1':
        for i,v in enumerate(l[1:]):
            v=v.strip('"')
            if ':' in v: k,val=v.split(':',1); meta.loc[i,k.strip()]=val.strip()
meta.to_csv('sample_metadata.csv',index=False)
tcol=[c for c in meta.columns if c.startswith('neutralizing')]
print('titer columns:',tcol)
meta['month']=meta['time'].str.extract(r'(\d+)').astype(float)
meta.loc[meta['time'].fillna('').str.contains('baseline'),'month']=0
# ---- expression
i0=next(i for i,l in enumerate(L) if l[0]=='!series_matrix_table_begin')
i1=next(i for i,l in enumerate(L) if l[0]=='!series_matrix_table_end')
E_=pd.DataFrame([l[1:] for l in L[i0+2:i1]],index=[l[0].strip('"') for l in L[i0+2:i1]],columns=gsm).apply(pd.to_numeric,errors='coerce')
print('expr',E_.shape,'range',float(E_.min().min()),float(E_.max().max()))
if E_.max().max()>50: E_=np.log2(E_.clip(lower=1))
ann=pd.read_csv('GPL570.annot.gz',sep='\t',skiprows=27,usecols=['ID','Gene symbol'],dtype=str).dropna()
ann=ann[~ann['Gene symbol'].str.contains('///')].set_index('ID')['Gene symbol']
G=E_.loc[E_.index.intersection(ann.index)].groupby(ann).mean()
print('genes',G.shape)
# ---- differential profiles: post minus same patient's baseline (IFN-K arm only)
pat=meta.dropna(subset=['patient'])
base=pat[pat.month==0].set_index('patient')['gsm']
post=pat[(pat.month>0)&(pat['agent']=='IFN-K')].copy()
post['baseline_gsm']=post.patient.map(base); post=post.dropna(subset=['baseline_gsm'])
Y=pd.DataFrame({r.gsm:G[r.gsm]-G[r.baseline_gsm] for r in post.itertuples()})
print('n post samples',Y.shape[1],'patients',post.patient.nunique())
def cytosig(S,Y,alpha=1e4,nrand=1000):
    c=S.index.intersection(Y.index); X=S.loc[c]; Yc=Y.loc[c]
    X=(X-X.mean())/X.std(); Yc=(Yc-Yc.mean())/Yc.std()
    Xa=X.values; P=np.linalg.solve(Xa.T@Xa+alpha*np.eye(Xa.shape[1]),Xa.T)
    b=P@Yc.values; R=np.zeros((nrand,)+b.shape)
    for k in range(nrand): R[k]=P@Yc.values[rng.permutation(len(c))]
    z=(b-R.mean(0))/R.std(0)
    return pd.DataFrame(z,index=X.columns,columns=Y.columns),len(c)
off=pd.read_csv('/workspace/cytosig_repro/CytoSig/CytoSig/signature.centroid',sep='\t',index_col=0)
rep=pd.read_csv('/workspace/cytosig_repro/output/diff.centroid',sep='\t',index_col=0)
res={}
for name,S in (('official',off),('reproduced',rep)):
    z,n=cytosig(S,Y); res[name]=z; z.T.to_csv(f'activity_z_{name}.csv'); print(name,'genes used',n)
post['titer']=pd.to_numeric(post[tcol[0]],errors='coerce')
for name in res: post[f'IFN1_{name}']=post.gsm.map(res[name].loc['IFN1']); post[f'IFNG_{name}']=post.gsm.map(res[name].loc['IFNG'])
# merge original source data
o=pd.read_excel('orig_source_data_fig4.xlsx',sheet_name='b')
o.columns=['id','titer_orig','IFN1_orig']
o['patient']=o.id.str.extract(r'blood&(\d+)')[0].astype(str); o['month']=o.id.str.extract(r'_(\d+) mon')[0].astype(float)
post['patient']=post.patient.astype(str).str.lstrip('0')
m=post.merge(o[['patient','month','titer_orig','IFN1_orig']],on=['patient','month'],how='left')
m.to_csv('per_sample_results.csv',index=False)
def ci_fisher(r,n):
    zr=np.arctanh(r); se=1/np.sqrt(n-3); return float(np.tanh(zr-1.96*se)),float(np.tanh(zr+1.96*se))
def boot(df,x,y,B=2000):
    ps=df.patient.unique(); out=[]
    for _ in range(B):
        s=pd.concat([df[df.patient==p] for p in rng.choice(ps,len(ps))]); 
        if s[x].nunique()>1: out.append(stats.spearmanr(s[x],s[y])[0])
    return [float(np.nanpercentile(out,2.5)),float(np.nanpercentile(out,97.5))]
S={'n_samples':int(len(m)),'n_patients':int(m.patient.nunique())}
r,p=stats.spearmanr(m.titer_orig,m.IFN1_orig); S['original_sourcedata_spearman']=dict(rho=float(r),p=float(p),n=int(m.IFN1_orig.notna().sum()))
for name in res:
    r,p=stats.spearmanr(m.titer,m[f'IFN1_{name}'])
    S[f'T1T2_{name}']=dict(rho=float(r),p_two_sided=float(p),ci95_fisher=ci_fisher(r,len(m)),ci95_patient_bootstrap=boot(m,'titer',f'IFN1_{name}'))
    g=m.groupby('patient')[['titer',f'IFN1_{name}',f'IFNG_{name}']].median()
    r2,p2=stats.spearmanr(g.titer,g[f'IFN1_{name}']); S[f'T3_patient_level_{name}']=dict(rho=float(r2),p=float(p2),n=len(g))
    d=m[f'IFN1_{name}']-m[f'IFNG_{name}']; w=stats.wilcoxon(m[f'IFN1_{name}'],m[f'IFNG_{name}'])
    S[f'T4_IFN1_vs_IFNG_{name}']=dict(median_IFN1=float(m[f'IFN1_{name}'].median()),median_IFNG=float(m[f'IFNG_{name}'].median()),
        n_samples_IFN1_below_IFNG=int((d<0).sum()),wilcoxon_p=float(w.pvalue),
        n_patients_median_IFN1_below_IFNG=int((g[f'IFN1_{name}']<g[f'IFNG_{name}']).sum()))
    r3,p3=stats.spearmanr(m[f'IFN1_{name}'],m.IFN1_orig,nan_policy='omit'); S[f'agreement_with_original_values_{name}']=dict(spearman=float(r3),p=float(p3))
    rI,pI=stats.spearmanr(m.titer,m[f'IFNG_{name}']); S[f'IFNG_vs_titer_{name}']=dict(rho=float(rI),p=float(pI))
    hi=m[m.titer>100]; lo=m[m.titer<=100]
    S[f'T5_group_{name}']=dict(mean_IFN1_titer_detected=float(hi[f'IFN1_{name}'].mean()),n_detected=len(hi),
        mean_IFN1_titer_undetected=float(lo[f'IFN1_{name}'].mean()),n_undetected=len(lo),mannwhitney_p=float(stats.mannwhitneyu(hi[f'IFN1_{name}'],lo[f'IFN1_{name}']).pvalue))
S['titer_match_with_original']=int((m.titer==m.titer_orig).sum())
json.dump(S,open('summary.json','w'),indent=1); print(json.dumps(S,indent=1))
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig,ax=plt.subplots(1,3,figsize=(15,4.5))
for a,(col,lab) in zip(ax,(('IFN1_orig','Original source data'),('IFN1_official','Reproduction, official signature'),('IFN1_reproduced','Reproduction, reproduced signature'))):
    for p_,g in m.groupby('patient'): a.scatter(g.titer,g[col],label=p_,s=25)
    a.set_xscale('log'); a.set_xlabel('neutralizing anti-IFNα titer (100 = undetected)'); a.set_ylabel('IFN1 activity (post − baseline)')
    r,p=stats.spearmanr(m.titer,m[col],nan_policy='omit'); a.set_title(f'{lab}\nSpearman ρ={r:.2f}, p={p:.2g}, n={m[col].notna().sum()} samples')
ax[0].legend(title='patient',fontsize=7); plt.tight_layout(); plt.savefig('fig_lupus_repro.png',dpi=150)
