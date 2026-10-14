import pandas as pd, numpy as np, glob, gzip, os, json
rng=np.random.default_rng(0)
D='data/'
# ---- signatures
SIG={'reproduced':pd.read_csv('../cytosig_repro/diff.centroid',sep='\t',index_col=0),
     'official':pd.read_csv('/workspace/cytosig_repro/CytoSig/CytoSig/signature.centroid',sep='\t',index_col=0)}
# ensembl->symbol map
m=pd.read_csv(D+'GSE125352_Raw_counts_IFNB_all_donors_2h.csv.gz',index_col=0)
e2s=m.iloc[:,0].to_dict()
def sym(df):
    df=df.copy(); df.index=[e2s.get(i,None) for i in df.index]; df=df[df.index.notnull()]
    return df.groupby(level=0).sum()
def logcpm(c): return np.log2(c/c.sum()*1e6+1)
def ridge(X,Y,alpha=10000,nrand=1000):
    X=(X-X.mean())/X.std(); Y=(Y-Y.mean())/Y.std()
    Xa=X.values; P=np.linalg.solve(Xa.T@Xa+alpha*np.eye(Xa.shape[1]),Xa.T)
    b=P@Y.values; cnt=np.zeros_like(b); s1=np.zeros_like(b); s2=np.zeros_like(b)
    Yv=Y.values
    for _ in range(nrand):
        br=P@Yv[rng.permutation(len(Yv))]; s1+=br; s2+=br**2
    mu=s1/nrand; sd=np.sqrt(s2/nrand-mu**2)
    return pd.DataFrame((b-mu)/sd,index=X.columns,columns=Y.columns)
contrasts=[]  # (dataset, label, stim_class, stim_series, ctrl_series, donor, note)
def add(ds,stim,cls,s,c,donor,indep='full'):
    contrasts.append(dict(dataset=ds,stim=stim,target=cls,s=s,c=c,donor=str(donor),indep=indep))
# GSE125352 IFN-b MDM; controls media from GSE82227 (normalized; in CytoSig pool) -> partial
for t in ['2h','6h','24h']:
    r=pd.read_csv(D+f'GSE125352_Raw_counts_IFNB_all_donors_{t}.csv.gz',index_col=0).iloc[:,1:]
    r.columns=[c.split()[1] for c in r.columns]; r=logcpm(sym(r))
    n=pd.read_csv(D+'GSE82227_counts_DESeq.normalized.csv.gz',index_col=0)
    idx={'2h':'01','6h':'05','24h':'09'}[t]
    for d in r.columns:
        cc=[c for c in n.columns if c.startswith(f'SR{int(d):03d}_{idx}')]
        if not cc: continue
        cv=logcpm(sym(n[cc]))[cc[0]]
        add('GSE125352',f'IFNB_{t}','IFN1',r[d],cv,d,'partial(ctrl in CytoSig pool; normalized ctrl)')
# GSE330780 THP1 WT
tp=pd.read_csv(D+'GSE330780_TPMs_allSamples.csv.gz',index_col=0); tp=tp.groupby(level=0).mean(); tp=np.log2(tp+1)
ctrl=tp[['moWT_R1','moWT_R2','moWT_R3','moWT_R4']].mean(1)
for c in ['moWT_IFNB1_R1','moWT_IFNB1_R2','moWT_IFNB1_24h','moWT_IFNB1_1x_R1','moWT_IFNB1_100x_R1']: add('GSE330780','IFNB_THP1','IFN1',tp[c],ctrl,c,'cell line; TPM')
for c in ['moWT_LPS_R1','moWT_LPS_R2']: add('GSE330780','LPS','induced',tp[c],ctrl,c,'cell line')
# IFNAR2 KO + IFNB: expected loss of response
ko=tp[['moIFNAR2KO_R1','moIFNAR2KO_R2']].mean(1); add('GSE330780','IFNB_in_IFNAR2KO','negKO',tp['moIFNAR2KO_IFNB1_R1'],ko,'KO','cell line')
# GSE174380 Jurkat normalized (assume log-scale already? check)
j=pd.read_csv(D+'GSE174380_HIV_Latency_Normalized_Count_Matrix_GEO.csv.gz',index_col=0)
if j.max().max()>50: j=np.log2(j+1)
for line in ['Jurkat','JLat92']:
    cm=j[[c for c in j if line in c and 'Mock' in c]].mean(1)
    for c in [c for c in j if line in c and 'IFNB' in c]: add('GSE174380','IFNB_Jurkat','IFN1',j[c],cm,c,'cell line; normalized')
# GSE327707 PBMC IFNa etc 48h, paired donors
r=pd.read_csv(D+'GSE327707_raw_counts.csv.gz',index_col=0)
import re
sm=gzip.open(D+'GSE327707_series_matrix.txt.gz','rt').read().splitlines()
tit=[l for l in sm if l.startswith('!Sample_title')][0].split('\t')[1:]
lib=[l for l in sm if l.startswith('!Sample_description')][0].split('\t')[1:]
meta=pd.DataFrame({'t':[x.strip('"') for x in tit],'lib':[x.strip('"').replace('Library name: ','') for x in lib]})
meta[['donor','cell','trt']]=meta.t.str.extract(r'Donor (\d+), (.+), (.+)')
r=logcpm(sym(r[meta.lib]))
cmap={'IFNa':'IFN1','LPS':'induced','Cytostim':'induced','HEK BCS (CD79B)':'nonIFN'}
for _,x in meta[meta.trt!='Baseline'].iterrows():
    b=meta[(meta.donor==x.donor)&(meta.cell==x.cell)&(meta.trt=='Baseline')].lib.iloc[0]
    add('GSE327707',f'{x.trt}|{x.cell}',cmap[x.trt],r[x.lib],r[b],x.donor)
# GSE218535 neutrophil 1h, paired donors
r=pd.read_csv(D+'GSE218535_rsubread_featureCounts.csv.gz',index_col=0).set_index('symbol').drop(columns='entrezid')
r=logcpm(r.groupby(level=0).sum())
for c in r.columns:
    if c.startswith('unstim'): continue
    d=c.split('_')[-1]; u=[x for x in r if x.startswith('unstim') and x.endswith(d)][0]
    k=c.split('_')[0]; k='IL13' if c.startswith('IL_13') else ('IL4' if c.startswith('IL_4') else 'IFNG')
    add('GSE218535',k,'IFNG' if k=='IFNG' else 'nonIFN',r[c],r[u],d)
# GSE288312 neutrophil 1h, unpaired
r=pd.read_csv(D+'GSE288312_Gene_counts_Wright_2025.csv.gz',index_col=0).set_index('GeneID'); r=logcpm(r.groupby(level=0).sum())
ut=r[[c for c in r if c.startswith('UT')]].mean(1)
for c in r:
    if c.startswith('UT'): continue
    k=c.rsplit('_',1)[0]; add('GSE288312',k,{'IFNA':'IFN1','IFNG':'IFNG'}.get(k,'nonIFN'),r[c],ut,c,'unpaired ctrl mean')
# GSE268354 hMDM 18h paired
fs=sorted(glob.glob(D+'GSM829*.tab.gz'))
r=pd.concat({os.path.basename(f).split('_',1)[1].split('.')[0]:pd.read_csv(f,sep='\t',index_col=0)['count'] for f in fs},axis=1)
r=logcpm(sym(r))
for d in ['D37','D91','D92']:
    for c in [c for c in r if c.startswith(d) and 'ctrl' not in c]:
        add('GSE268354','IFNG' if c.endswith('IFNg') else 'IFNG+TNF','IFNG',r[c],r[f'{d}_NTC_ctrl'],d)
# ---- run
out=[]
for sname,S in SIG.items():
    for grp in pd.Series([c['dataset'] for c in contrasts]).unique():
        cs=[c for c in contrasts if c['dataset']==grp]
        Y=pd.concat([c['s']-c['c'] for c in cs],axis=1).dropna(); Y.columns=range(len(cs))
        g=S.index.intersection(Y.index); Z=ridge(S.loc[g],Y.loc[g])
        for i,c in enumerate(cs):
            z=Z[i]; row={k:c[k] for k in ['dataset','stim','target','donor','indep']}
            row.update(signature=sname,n_genes=len(g)); row.update({f'z_{k}':v for k,v in z.items()}); out.append(row)
res=pd.DataFrame(out); res.to_csv('tables/per_contrast_z.csv',index=False)
print(res.groupby(['signature','dataset']).size())
