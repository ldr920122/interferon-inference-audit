# v1.6: rebuild the reproduced signature from the 1,307 resolved final columns
# (1,306 reconstructed columns + NO&Dose:600 uM_Duration:24h@Condition:T cell@GSE13887; see missing_column_search.csv)
# and compare every family with the official CytoSig signature.centroid. Same family merge and >=5 rule as repro.py.
import pandas as pd,numpy as np,hashlib
D='/workspace/cytosig_repro'
cols=[l.strip().strip('"') for l in open('final_columns_resolved_1307.txt') if l.strip()]
sig=pd.read_csv(f'{D}/CytoSig_data_filter/data/diff.merge.gz',sep='\t',index_col=0)
assert all(c in sig.columns for c in cols), 'missing column'
rn={'IFNA':'IFN1','IFNB':'IFN1','IL36A':'IL36','IL36B':'IL36','IL36G':'IL36'}
d=sig[cols]; fam=[rn.get(v.split('@')[0].split('&')[0],v.split('@')[0].split('&')[0]) for v in d.columns]
cents=[]
for f,cc in pd.Series(d.columns,index=fam).groupby(level=0):
    x=d[list(cc)]
    if x.shape[1]<5: continue
    x=x.loc[(~x.isnull()).sum(axis=1)>=5]; a=x.median(axis=1); a.name=f; cents.append(a)
cen=pd.concat(cents,axis=1,join='outer').dropna(); cen.to_csv('diff.centroid_1307',sep='\t')
off=pd.read_csv(f'{D}/CytoSig/CytoSig/signature.centroid',sep='\t',index_col=0)
rows=[]
for f in off.columns:
    g=off.index.intersection(cen.index); m=float((cen.loc[g,f]-off.loc[g,f]).abs().max()) if f in cen else np.nan
    rows.append(dict(family=f,genes_compared=len(g),max_abs_diff=m,identical=bool(m<1e-12)))
r=pd.DataFrame(rows); r.to_csv('centroid_1307_vs_official.csv',index=False)
print('columns',len(cols),'centroid',cen.shape,'official',off.shape,'genes equal',set(cen.index)==set(off.index),
      'families identical',int(r.identical.sum()),'/',len(r),'max diff',r.max_abs_diff.max())
for fn in ['final_columns_resolved_1307.txt','diff.centroid_1307','centroid_1307_vs_official.csv']:
    print(fn,'sha256',hashlib.sha256(open(fn,'rb').read()).hexdigest())
