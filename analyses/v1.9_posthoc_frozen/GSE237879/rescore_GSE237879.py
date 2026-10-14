"""GSE237879 post-hoc per-donor rescoring (v1.9 re-audit). Same contract as GSE236156 rescore.py:
log2(1+x/libsize*scale) stimulated minus same-donor medium; all gene rows in denominator; official CytoSig
signature, ridge alpha=10000, 1000 gene-row permutations, seed 20261014. Inputs: GEO supplementary counts (Salmon) + metadata."""
import pandas as pd, numpy as np, sys, json, hashlib
import os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); sys.path.insert(1,'/workspace/ici/v19/gse236156')
c=pd.read_csv('GSE237879_counts.txt.gz',sep='\t',index_col=0); c.index=c.index.str.split('.').str[0]
m=pd.read_csv('GSE237879_metadata.txt.gz',sep='\t').set_index('ID_library'); assert set(m.index)==set(c.columns)
h=pd.read_csv('/workspace/ici/v19/gse236156/data/hgnc_complete_set.txt',sep='\t',usecols=['symbol','ensembl_gene_id'],dtype=str).dropna()
k=h.ensembl_gene_id.value_counts(); h=h[h.ensembl_gene_id.map(k)==1]; e2s=dict(zip(h.ensembl_gene_id,h.symbol))
SIG=pd.read_csv(__import__('os').environ.get('CYTOSIG_SIG','/workspace/cytosig_repro/CytoSig/CytoSig/signature.centroid'),sep='\t',index_col=0)
from rescore_lib import ridge
rows=[]
for scale,lab in [(1e4,'CP10k'),(1e6,'CPM')]:
    lg=np.log2(c/c.sum()*scale+1); lg=lg[lg.index.isin(e2s)]; lg.index=lg.index.map(e2s); lg=lg.groupby(level=0).mean()
    cols={}
    for d in sorted(m.Donor.unique()):
        md=m[m.Donor==d]; ctrl=md.index[md.Stim=='medium'][0]
        for s in ['beta','gamma','alpha','lambda']: cols[(d,s)]=lg[md.index[md.Stim==s][0]]-lg[ctrl]
    Y=pd.DataFrame(cols); g=SIG.index.intersection(Y.index); Z=ridge(SIG.loc[g],Y.loc[g],seed=20261014)
    for (d,s) in Y.columns:
        z=Z[(d,s)]; o=list(z.sort_values(ascending=False).index)
        rows.append(dict(scale=lab,donor=d,stim=s,n_genes=len(g),z_IFN1=z['IFN1'],z_IFNG=z['IFNG'],IFNG_minus_IFN1=z['IFNG']-z['IFN1'],rank_IFN1=o.index('IFN1')+1,rank_IFNG=o.index('IFNG')+1,top1=o[0]))
R=pd.DataFrame(rows); R.to_csv('per_donor_cytosig_GSE237879.csv',index=False)
print(R[R.stim.isin(['beta','gamma'])].round(2).to_string())
for lab in ['CP10k','CPM']:
    for s,w in [('gamma',1),('beta',-1)]:
        x=R[(R.scale==lab)&(R.stim==s)].IFNG_minus_IFN1; print(lab,s,'correct',int((x*w>0).sum()),'/',len(x))
json.dump({f:hashlib.sha256(open(f,'rb').read()).hexdigest() for f in ['GSE237879_counts.txt.gz','GSE237879_metadata.txt.gz']},open('input_sha256.json','w'),indent=1)
