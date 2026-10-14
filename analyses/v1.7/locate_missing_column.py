# M7: try to locate the CytoSig final column missing from the reconstruction (NO family)
import pandas as pd,numpy as np,sys
D='/workspace/cytosig_repro'
sig=pd.read_csv(sys.argv[1],sep='\t',index_col=0)
off=pd.read_csv(f'{D}/CytoSig/CytoSig/signature.centroid',sep='\t',index_col=0)
inc=set(pd.read_csv(f'{D}/output/included_columns_after_correlation.txt',header=None)[0])
fam=lambda v:v.split('@')[0].split('&')[0]
nocols=[c for c in sig.columns if fam(c)=='NO']
ours=[c for c in nocols if c in inc]; others=[c for c in nocols if c not in inc]
print('NO columns in diff.merge',len(nocols),'in reconstruction',len(ours),'excluded',len(others))
genes=off.index.intersection(sig.index)
def err(cols):
    x=sig.loc[:,cols]; x=x.loc[(~x.isnull()).sum(axis=1)>=5]; m=x.median(axis=1).reindex(genes).dropna()
    return float((m-off.loc[m.index,'NO']).abs().max()), int(len(m))
print('current',err(ours))
res=[]
for c in others: res.append(('add',c)+err(ours+[c]))
for c in ours: res.append(('drop',c)+err([o for o in ours if o!=c]))
r=pd.DataFrame(res,columns=['op','column','max_abs_diff','genes']).sort_values('max_abs_diff')
print(r.head(8).to_string()); r.to_csv('/workspace/ici/deliverables/v15_analyses/missing_column_search.csv',index=False)
