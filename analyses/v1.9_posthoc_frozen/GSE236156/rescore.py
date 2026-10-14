"""GSE236156 post-hoc per-donor rescoring (v1.9). Inputs: author featureCounts + metadata
(hawn-lab/AM_MDM_TB_public @0076bde4, MDM_IFN/data_raw), HGNC complete set, official CytoSig signature.centroid.
Response = log2(1+count/libsize*scale) treatment minus same-donor untreated; denominator = all 60,664 gene rows.
CytoSig: ridge alpha=10000, 1000 gene-row permutations, fixed seed. No tuning."""
import pandas as pd, numpy as np, json, hashlib
D='data/'
fc=pd.read_csv(D+'MDM-IFN.featurecounts.tsv',sep='\t',comment='#',index_col=0).iloc[:,5:]
fc.columns=[c.split('/')[-1].replace('_filter.bam','').replace('-','.') for c in fc.columns]
meta=pd.read_csv(D+'MDM-IFN.metadata.csv').set_index('libID')
assert set(meta.index)==set(fc.columns) and fc.shape[0]==60664
h=pd.read_csv(D+'hgnc_complete_set.txt',sep='\t',usecols=['symbol','ensembl_gene_id'],dtype=str).dropna()
cnt=h.ensembl_gene_id.value_counts(); h=h[h.ensembl_gene_id.map(cnt)==1]; e2s=dict(zip(h.ensembl_gene_id,h.symbol))
SIG=pd.read_csv(__import__('os').environ.get('CYTOSIG_SIG','/workspace/cytosig_repro/CytoSig/CytoSig/signature.centroid'),sep='\t',index_col=0)
def ridge(X,Y,seed,alpha=10000,nrand=1000):
    rng=np.random.default_rng(seed)
    X=(X-X.mean())/X.std(); Y=(Y-Y.mean())/Y.std(); Xa=X.values
    P=np.linalg.solve(Xa.T@Xa+alpha*np.eye(Xa.shape[1]),Xa.T); b=P@Y.values
    s1=np.zeros_like(b); s2=np.zeros_like(b); Yv=Y.values
    for _ in range(nrand):
        br=P@Yv[rng.permutation(len(Yv))]; s1+=br; s2+=br**2
    mu=s1/nrand; sd=np.sqrt(s2/nrand-mu**2)
    return pd.DataFrame((b-mu)/sd,index=X.columns,columns=Y.columns)
stims=['IFNb','IFNg','IFNA8','IFNe']; donors=sorted(meta.ptID.unique())
rows=[]; allz=[]
for scale,lab in [(1e4,'CP10k'),(1e6,'CPM')]:
    lg=np.log2(fc/fc.sum()*scale+1)
    lg=lg[lg.index.isin(e2s)]; lg.index=lg.index.map(e2s); lg=lg.groupby(level=0).mean()
    cols={}; 
    for d in donors:
        m=meta[meta.ptID==d]; ctrl=m.index[m.IFN=='untreated'][0]
        for s in stims: cols[(d,s)]=lg[m.index[m.IFN==s][0]]-lg[ctrl]
    Y=pd.DataFrame(cols); g=SIG.index.intersection(Y.index)
    Z=ridge(SIG.loc[g],Y.loc[g],seed=20261014)
    for (d,s) in Y.columns:
        z=Z[(d,s)]; ai=sorted(z.index,key=lambda k:-z[k])
        rows.append(dict(scale=lab,donor=d,stim=s,n_genes=len(g),z_IFN1=z['IFN1'],z_IFNG=z['IFNG'],
            IFNG_minus_IFN1=z['IFNG']-z['IFN1'],rank_IFN1=ai.index('IFN1')+1,rank_IFNG=ai.index('IFNG')+1,top1=ai[0]))
        allz.append(pd.Series(z,name=f'{lab}|{d}|{s}'))
R=pd.DataFrame(rows); R.to_csv('per_donor_cytosig.csv',index=False)
pd.concat(allz,axis=1).T.to_csv('cytosig_all43_per_contrast.csv')
print(R[R.stim.isin(['IFNb','IFNg'])].round(2).to_string())
for lab in ['CP10k','CPM']:
    for s,want in [('IFNg',1),('IFNb',-1)]:
        x=R[(R.scale==lab)&(R.stim==s)].IFNG_minus_IFN1
        print(lab,s,'correct',int(((x*want)>0).sum()),'/5','mean',round(x.mean(),2))
json.dump({f:hashlib.sha256(open(D+f,'rb').read()).hexdigest() for f in ['MDM-IFN.featurecounts.tsv','MDM-IFN.metadata.csv','hgnc_complete_set.txt']},open('input_sha256.json','w'),indent=1)
