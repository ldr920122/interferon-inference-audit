# Faithful numpy port of CytoSig_data_filter@5e4cf05b src/correlation.py.
# Only change: CytoSig.ridge_significance_test(alpha=0,nrand=0,flag_const=True,flag_normalize=False) replaced by
# closed-form OLS t-statistic of the signature coefficient (data_significance C/GSL extension could not be built: no make/GSL/sudo).
import re,os,sys,time,numpy as np,pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
D='/workspace/cytosig_repro/CytoSig_data_filter/data'; OUT='/workspace/cytosig_repro/output'; os.makedirs(OUT,exist_ok=True)
def load_gene_map():
  sp=re.compile('[,+]'); inc=set(); gm=[]
  for title in ['Cytokine','Chemokine','Growth_Factor','Inhibitory']:
    xls=pd.ExcelFile(f'{D}/annotation/{title}.xlsx',engine='openpyxl')
    for fam in xls.sheet_names:
      if fam=='Extra': continue
      d=xls.parse(fam,index_col=0)[['Gene','Receptor']]; d=d.loc[~d.index.isnull()]
      for v in d:
        f=~d[v].isnull(); d.loc[f,v]=d.loc[f,v].apply(lambda x:set(sp.split(x.replace(' ','').replace('*',''))))
      d=d.loc[d.index.difference(inc)]; gm.append(d); inc.update(d.index)
  return pd.concat(gm)
def zscores(sig,expr):
  common=expr.index.intersection(sig.index); S=sig.loc[common].values; Y=expr.loc[common].values.astype(np.float64)
  M=(~np.isnan(S)).astype(np.float64); S0=np.nan_to_num(S)
  n=M.sum(0); Ss=S0.sum(0); Sss=(S0*S0).sum(0)
  Sy=M.T@Y; Syy=M.T@(Y*Y); Ssy=S0.T@Y
  ss=(Sss-Ss**2/n)[:,None]; sy=Ssy-Ss[:,None]*Sy/n[:,None]; yy=Syy-Sy**2/n[:,None]
  beta=sy/ss; rss=yy-beta*sy; se=np.sqrt(rss/(n[:,None]-2)/ss)
  return pd.DataFrame((beta/se).T,index=expr.columns,columns=sig.columns)
def run_cohort(cohort,f,sig,gm):
  out=f'{OUT}/{cohort}.{f}'
  if os.path.exists(out): return
  expr=pd.read_csv(f'{D}/{cohort}/{f}.gz',sep='\t',index_col=0)
  Z=zscores(sig,expr); rows={}
  for title in sig.columns:
    gid=title.split('@')[0].split('&')[0]
    if gid not in gm.index: continue
    piv=gm.loc[gid].dropna(); r={}
    for cat,s in piv.items():
      g=expr.index.intersection(list(s))
      if len(g)==0: continue
      p=expr.loc[g].mean()
      if p.std()==0: continue
      r['pearson.'+cat]=Z[title].corr(p); r['spearman.'+cat]=Z[title].corr(p,method='spearman')
    rows[title]=r
  pd.DataFrame(rows).T.reindex(columns=['pearson.Gene','pearson.Receptor','spearman.Gene','spearman.Receptor']).to_csv(out,sep='\t')
def merge(run):
  for piv in ['Gene','Receptor']:
    cm={}
    for c,f in run:
      if 'ALL' in f or 'AML' in f: continue
      d=pd.read_csv(f'{OUT}/{c}.{f}',sep='\t',index_col=0)['pearson.'+piv].dropna(); d.name=f; cm.setdefault(c,[]).append(d)
    for c,l in cm.items():
      d=pd.concat(l,axis=1,join='outer',sort=False); d=d.loc[d.isnull().mean(axis=1)<0.1]
      d=d.apply(lambda v:pd.Series([v.dropna().median(),stats.wilcoxon(v.dropna(),alternative='greater')[1]],index=['med','p']),axis=1)
      d['FDR']=multipletests(d['p'],method='fdr_bh')[1]; d.to_csv(f'{OUT}/merge.{piv}.{c}.stat',sep='\t')
def filt(sig,cnt=5,q=0.05):
  rn={'IFNA':'IFN1','IFNB':'IFN1','IL36A':'IL36','IL36B':'IL36','IL36G':'IL36'}; inc=None
  for c in ['TCGA','GTEx']:
    g=pd.read_csv(f'{OUT}/merge.Gene.{c}.stat',sep='\t',index_col=0); r=pd.read_csv(f'{OUT}/merge.Receptor.{c}.stat',sep='\t',index_col=0)
    s=g.index[g['FDR']<q].union(r.index[r['FDR']<q]); inc=s if inc is None else inc.intersection(s)
  d=sig[inc]; fam=[rn.get(v.split('@')[0].split('&')[0],v.split('@')[0].split('&')[0]) for v in d.columns]
  famcount=pd.Series(fam).value_counts(); kept=[c for c,f in zip(d.columns,fam) if famcount[f]>=cnt]
  cents=[]
  for f,cols in pd.Series(d.columns,index=fam).groupby(level=0):
    x=d[list(cols)]
    if x.shape[1]<cnt: continue
    x=x.loc[(~x.isnull()).sum(axis=1)>=cnt]; a=x.median(axis=1); a.name=f; cents.append(a)
  cen=pd.concat(cents,axis=1,join='outer').dropna()
  pd.Series(list(inc)).to_csv(f'{OUT}/included_columns_after_correlation.txt',index=False,header=False)
  pd.Series(kept).to_csv(f'{OUT}/final_columns_in_families_ge5.txt',index=False,header=False)
  cen.to_csv(f'{OUT}/diff.centroid',sep='\t'); print('included',len(inc),'kept',len(kept),'centroid',cen.shape)
if __name__=='__main__':
  gm=load_gene_map(); run=[]
  for c in ['TCGA','GTEx']:
    for f in pd.read_csv(f'{D}/{c}/file.list',header=None,sep='\t').iloc[:,0]: run.append((c,f))
  sig=pd.read_csv(f'{D}/diff.merge.gz',sep='\t',index_col=0); print('sig',sig.shape,flush=True)
  t0=time.time()
  for i,(c,f) in enumerate(run):
    run_cohort(c,f,sig,gm); print(i,c,f,round(time.time()-t0),flush=True)
  merge(run); filt(sig); print('total s',round(time.time()-t0))
