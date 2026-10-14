"""Leakage pilot: CytoSig-style inference on IFN-beta/gamma training profiles, with the full (leaky)
reproduced signature vs a signature rebuilt without the profile's own accession (leave-one-accession-out).
Ridge closed form (alpha=1e4), signature and input z-scored across genes (approximation of
CytoSig flag_normalize); no permutation z-scores (C extension unavailable). No new model is fitted:
signatures are medians of the reproduced final columns exactly as in CytoSig_data_filter."""
import pandas as pd, numpy as np, json, re
sig=pd.read_csv('/workspace/cytosig_repro/CytoSig_data_filter/data/diff.merge.gz',sep='\t',index_col=0)
sig.columns=[c.strip('"') for c in sig.columns]
final=[l.strip('"') for l in open('/workspace/ici/deliverables/routeB/cytosig_repro/final_columns_in_families_ge5.txt').read().split('\n') if l]
off=pd.read_csv('/workspace/cytosig_repro/CytoSig/CytoSig/signature.centroid',sep='\t',index_col=0)
rn={'IFNA':'IFN1','IFNB':'IFN1','IL36A':'IL36','IL36B':'IL36','IL36G':'IL36'}
lig=lambda c:re.split(r'[&@]',c)[0]; fam=lambda c:rn.get(lig(c),lig(c)); acc=lambda c:c.split('@')[-1].split('.')[0]
D=sig[final]; F=pd.Series([fam(c) for c in final],index=final); A=pd.Series([acc(c) for c in final],index=final)
genes=off.index
def centroid(cols):
    out={}
    for f,cs in F[cols].groupby(F[cols]).groups.items():
        if len(cs)<5: continue
        x=D[list(cs)]; x=x.loc[x.notnull().sum(axis=1)>=5]; out[f]=x.median(axis=1)
    return pd.DataFrame(out).reindex(genes)
full=centroid(final); print('full vs official max diff (non-NO):',float((full[off.columns.drop('NO')]-off[off.columns.drop('NO')]).abs().max().max()))
def infer(S,y):
    m=S.notnull().all(axis=1)&y.notnull(); X=S[m].values; Y=y[m].values
    X=(X-X.mean(0))/X.std(0); Y=(Y-Y.mean())/Y.std()
    b=np.linalg.solve(X.T@X+1e4*np.eye(X.shape[1]),X.T@Y); return pd.Series(b,index=S.columns)
tests=[c for c in final if lig(c) in ('IFNB','IFNG')]
rows=[]
for a in sorted(set(A[tests])):
    keep=[c for c in final if A[c]!=a]; S=centroid(keep).reindex(columns=full.columns)
    for c in [t for t in tests if A[t]==a]:
        y=D[c].reindex(genes); tf=fam(c); r={}
        for k,Sx in (('leaky',full),('loao',S)):
            b=infer(Sx.dropna(axis=1,how='all'),y); rk=b.rank(ascending=False)
            r[k+'_rank_true']=int(rk.get(tf,np.nan)) if tf in rk else None; r[k+'_top']=b.idxmax()
            r[k+'_IFNG_minus_IFN1']=float(b['IFNG']-b['IFN1'])
        rows.append(dict(accession=a,column=c,true_family=tf,**r))
df=pd.DataFrame(rows); df.to_csv('loao_results.csv',index=False)
s={}
for tf,g in df.groupby('true_family'):
    s[tf]=dict(n_profiles=len(g),n_accessions=g.accession.nunique(),
      top1_leaky=int((g.leaky_rank_true==1).sum()),top1_loao=int((g.loao_rank_true==1).sum()),
      median_rank_leaky=float(g.leaky_rank_true.median()),median_rank_loao=float(g.loao_rank_true.median()),
      correct_sign_IFNG_minus_IFN1_leaky=int(((g.leaky_IFNG_minus_IFN1>0)==(tf=='IFNG')).sum()),
      correct_sign_IFNG_minus_IFN1_loao=int(((g.loao_IFNG_minus_IFN1>0)==(tf=='IFNG')).sum()),
      rank_worse_after_removal=int((g.loao_rank_true>g.leaky_rank_true).sum()))
json.dump(s,open('loao_summary.json','w'),indent=1); print(json.dumps(s,indent=1))
