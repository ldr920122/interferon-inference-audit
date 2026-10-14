"""GSE78193 case: AMG 811 (anti-IFN-g) SLE arms, IFNG inference with full vs leave-series-out reproduced CytoSig signature.
Reuses the pilot's signature build and ridge approximation (no permutation z)."""
import pandas as pd, numpy as np, gzip, re, json, io
from scipy.stats import spearmanr, mannwhitneyu
exec(open('/workspace/ici/deliverables/leakage_pilot/loao_pilot.py').read().split('tests=[')[0])  # sig, final, centroid, infer, full
lso=centroid([c for c in final if A[c]!='GSE78193']).reindex(columns=full.columns)
removed=[c for c in final if A[c]=='GSE78193']
# expression + metadata
lines=gzip.open('GSE78193-GPL6480_series_matrix.txt.gz','rt').read().split('\n')
meta={}
for l in lines:
    if l.startswith('!Sample_characteristics_ch1') or l.startswith('!Sample_geo_accession') or l.startswith('!Sample_title'):
        v=[x.strip('"') for x in l.split('\t')[1:]]
        k=v[0].split(':')[0] if ':' in v[0] else l.split('\t')[0]
        meta[k]=[x.split(': ',1)[-1] for x in v]
M=pd.DataFrame(meta); M.columns=['title','gsm']+list(M.columns[2:])
i0=lines.index('!series_matrix_table_begin'); i1=lines.index('!series_matrix_table_end')
X=pd.read_csv(io.StringIO('\n'.join(lines[i0+1:i1])),sep='\t',index_col=0)
ann=pd.read_csv('GPL6480.annot.gz',sep='\t',skiprows=28,usecols=['ID','Gene symbol'],dtype=str).dropna()
ann=ann[~ann['Gene symbol'].str.contains('///')].set_index('ID')['Gene symbol']
G=X.loc[X.index.intersection(ann.index)].groupby(ann).mean()
M=M.set_index('gsm'); M['day']=M['treatment day']; M['dose']=M['dose']
rows=[]
for sid,g in M.groupby('subject id'):
    base=g[g.day=='day 1']
    if len(base)!=1: continue
    b=G[base.index[0]]
    for gsm,r in g[g.day!='day 1'].iterrows():
        y=(G[gsm]-b).reindex(genes)
        out=dict(subject=sid,gsm=gsm,day=r.day,dose=r.dose)
        for k,S in (('full',full),('lso',lso)):
            be=infer(S.dropna(axis=1,how='all'),y); out[k+'_IFNG']=float(be['IFNG']); out[k+'_IFNG_rank_low']=int(be.rank()['IFNG'])
        rows.append(out)
R=pd.DataFrame(rows); R['arm']=np.where(R.dose=='0mg','placebo','AMG811')
R.to_csv('per_sample_results.csv',index=False)
s={'removed_training_columns':removed,'n_subjects':int(R.subject.nunique()),'n_post_samples':len(R),
   'n_genes_used':int(G.index.intersection(genes).size)}
s['T1_spearman_full_vs_lso']=float(spearmanr(R.full_IFNG,R.lso_IFNG)[0])
s['T1_medians']={f'{a}|{d}|{k}':float(v) for (a,d),g in R.groupby(['arm','day']) for k,v in (('full',g.full_IFNG.median()),('lso',g.lso_IFNG.median()))}
subj=R.groupby(['subject','arm'])[['full_IFNG','lso_IFNG']].mean().reset_index()
for k in ('full','lso'):
    a=subj[subj.arm=='AMG811'][k+'_IFNG']; p=subj[subj.arm=='placebo'][k+'_IFNG']
    s[f'T2_{k}']=dict(n_amg=len(a),n_pbo=len(p),median_amg=float(a.median()),median_pbo=float(p.median()),p=float(mannwhitneyu(a,p).pvalue))
for k in ('full','lso'):
    S=full if k=='full' else lso
    for d in ('day 15','day 56','end of study'):
        gs=R[(R.arm=='AMG811')&(R.day==d)]
        if not len(gs): continue
        prof=pd.concat([(G[x]-G[M[(M['subject id']==sb)&(M.day=='day 1')].index[0]]) for x,sb in zip(gs.gsm,gs.subject)],axis=1).mean(axis=1).reindex(genes)
        be=infer(S.dropna(axis=1,how='all'),prof); s[f'T3_{k}_{d}']=dict(IFNG_rank_most_negative=int(be.rank()['IFNG']),n_cytokines=len(be),most_negative=be.idxmin())
json.dump(s,open('summary.json','w'),indent=1); print(json.dumps(s,indent=1))
