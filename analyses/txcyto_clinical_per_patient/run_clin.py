import sys,gzip,json,re,numpy as np,pandas as pd
from scipy import stats
sys.path.insert(0,'/workspace/cytosig_repro/data_significance'); sys.path.insert(0,'/workspace/cytosig_repro/CytoSig/CytoSig')
import Util as U
S=pd.read_csv('/workspace/cytosig_repro/CytoSig/CytoSig/signature.centroid',sep='\t',index_col=0)
def load(g,gpl):
    meta={};rows=[]
    with gzip.open(f'{g}_series_matrix.txt.gz','rt',errors='ignore') as f:
        for l in f:
            if l.startswith('!Sample_geo_accession'): meta['gsm']=[x.strip('"') for x in l.rstrip('\n').split('\t')[1:]]
            elif l.startswith('!Sample_title'): meta['title']=[x.strip('"') for x in l.rstrip('\n').split('\t')[1:]]
            elif l.startswith('!Sample_characteristics_ch1'):
                v=[x.strip('"') for x in l.rstrip('\n').split('\t')[1:]]
                k=v[0].split(':')[0].strip().lower(); meta[k]=[x.split(':',1)[1].strip() if ':' in x else x for x in v]
            elif l.startswith('!series_matrix_table_begin'): break
        E=pd.read_csv(f,sep='\t',index_col=0,comment='!')
    E.columns=[c.strip('"') for c in E.columns]; E.index=E.index.astype(str).str.strip('"')
    M=pd.DataFrame(meta)
    a=pd.read_csv(f'{gpl}.annot.gz',sep='\t',skiprows=27,usecols=['ID','Gene symbol'],dtype=str).dropna()
    a=a[~a['Gene symbol'].str.contains('///')]; mp=dict(zip(a.ID,a['Gene symbol']))
    E=E.apply(pd.to_numeric,errors='coerce'); E['sym']=E.index.map(mp); E=E.dropna(subset=['sym']).groupby('sym').mean()
    if E.max().max()>50: E=np.log2(E.clip(lower=1))
    return E,M
def infer(E,pairs):
    Y=pd.DataFrame({p:E[p]-E[b] for p,b in pairs.items()})
    c=S.index.intersection(Y.index)
    return U.ridge_significance_test(S.loc[c],Y.loc[c],1e4,nrand=1000,verbose=False)[2].T,len(c)
def tests(d,targets,arm_col,treated,comp):
    out=[]
    for t in targets:
        for lab,sub in [('treated',d[d[arm_col]==treated])]:
            ps=stats.wilcoxon(sub[t]).pvalue if len(sub)>1 else np.nan
            pm=sub.groupby('patient')[t].median(); pp=stats.wilcoxon(pm).pvalue if len(pm)>1 else np.nan
            out.append(dict(target=t,test='A treated change vs 0',n_samples=len(sub),n_patients=len(pm),median_sample=sub[t].median(),p_per_sample=ps,p_per_patient=pp))
        if comp is not None:
            a=d[d[arm_col]==treated]; b=d[d[arm_col]==comp]
            if len(b)>1:
                ps=stats.mannwhitneyu(a[t],b[t]).pvalue
                am=a.groupby('patient')[t].median(); bm=b.groupby('patient')[t].median()
                pp=stats.mannwhitneyu(am,bm).pvalue
                out.append(dict(target=t,test=f'B {treated} vs {comp}',n_samples=f'{len(a)} vs {len(b)}',n_patients=f'{len(am)} vs {len(bm)}',median_sample=a[t].median()-b[t].median(),p_per_sample=ps,p_per_patient=pp))
    return out
R=[];log={}
def run(g,gpl,build,targets,arm_col,treated,comp):
    try:
        E,M=load(g,gpl); M,pairs=build(M)
        Z,n=infer(E,pairs); d=M.set_index('gsm').loc[list(pairs)].join(Z[targets])
        d.to_csv(f'{g}_per_sample.csv')
        for r in tests(d,targets,arm_col,treated,comp): r.update(series=g,n_genes=n); R.append(r)
        log[g]='ok'
    except Exception as e: log[g]=f'FAILED: {type(e).__name__}: {e}'
# GSE100093: lesional, baseline D-1; drug = AMG811 on Day1 sequence (D15,D57 within first period)
def b100093(M):
    M=M[M['skin lesion status']=='Lesional'].copy(); M['patient']=M['sbj']
    M['arm']=np.where(M['sequence'].str.startswith('AMG811'),'AMG811','Placebo')
    base=M[M['visit name']=='D-1'].set_index('patient').gsm.to_dict()
    post=M[M['visit name']!='D-1']; post=post[post.patient.isin(base)]
    return post,{r.gsm:base[r.patient] for r in post.itertuples()}
run('GSE100093','GPL13158',b100093,['IFNG'],'arm','AMG811','Placebo')
def b130588(M):
    M=M[(M.tissue=='LS')].copy(); M['patient']=M['subject id']
    base=M[M.time=='Week0'].set_index('patient').gsm.to_dict()
    post=M[(M.time!='Week0')&M.patient.isin(base)]
    return post,{r.gsm:base[r.patient] for r in post.itertuples()}
run('GSE130588','GPL570',b130588,['IL4','IL13'],'treatment','Dupilumab','Placebo')
def b59294(M):
    M=M[M.tissue=='LS skin'].copy(); M['arm']=np.where(M.treatment.str.contains('Placebo'),'Placebo','Dupilumab')
    base=M[M.time=='Pre'].set_index('patient').gsm.to_dict()
    post=M[(M.time=='Post')&M.patient.isin(base)]
    return post,{r.gsm:base[r.patient] for r in post.itertuples()}
run('GSE59294','GPL570',b59294,['IL4','IL13'],'arm','Dupilumab','Placebo')
def b11903(M):
    M=M[M.condition=='lesional'].copy(); M['patient']=M['patient']; M['arm']=M['group']
    base=M[M.time=='0'].set_index('patient').gsm.to_dict()
    post=M[(M.time!='0')&M.patient.isin(base)]
    return post,{r.gsm:base[r.patient] for r in post.itertuples()}
run('GSE11903','GPL571',b11903,['TNFA','IL17A'],'arm','responder','non-responder')
def b80060(M):
    M=M[M['disease status']=='SJIA'].copy(); M['patient']=M['subject id']; M['arm']=M.treatment
    base=M[M.visit=='Day1_BL'].set_index('patient').gsm.to_dict()
    post=M[(M.visit=='Day3')&M.patient.isin(base)]
    return post,{r.gsm:base[r.patient] for r in post.itertuples()}
run('GSE80060','GPL570',b80060,['IL1B'],'arm','Canakinumab','Placebo')
def b99802(M):
    M=M[M.tissue=='LS'].copy(); M['patient']=M['patient id']; M['arm']=M.treatment
    base=M[M.week=='Wk 0'].set_index('patient').gsm.to_dict()
    post=M[(M.week!='Wk 0')&M.patient.isin(base)]
    return post,{r.gsm:base[r.patient] for r in post.itertuples()}
run('GSE99802','GPL570',b99802,['IL22'],'arm','Drug','Placebo')
log['GSE42048']='excluded in advance: ACHN cell-line xenografts in nude mice, no patients'
pd.DataFrame(R).to_csv('results_per_sample_vs_per_patient.csv',index=False)
json.dump(log,open('run_log.json','w'),indent=1)
print(json.dumps(log,indent=1)); print(pd.DataFrame(R).to_string())
