# v1.6 (two-sided sign-flip fixed; v1.5 computed one-sided): per-dataset, per-donor stratified summary of loose independent validation (official signature)
import pandas as pd,numpy as np,itertools
d=pd.read_csv('../routeB/validation/tables/per_contrast_metrics.csv')
meta={'GSE327707':('IFN-α','primary','paired','independent'),'GSE330780':('IFN-β','cell line','unpaired (control mean)','independent'),
'GSE174380':('IFN-β','cell line','technical replicates','independent'),'GSE125352':('IFN-β','primary','external control (GSE82227)','partial overlap'),
'GSE218535':('IFN-γ','primary','paired','independent'),'GSE268354':('IFN-γ','primary','paired','independent'),'GSE288312':(None,'primary','unpaired (control mean)','independent')}
def signflip(x,two=True):
    x=np.asarray(x); n=len(x); obs=x.sum(); c=0;t=0
    for s in itertools.product([1,-1],repeat=n): t+=1; v=np.dot(s,np.abs(x)); c+= (abs(v)>=abs(obs)-1e-12) if two else (v>=obs-1e-12)
    return c/t
rows=[]
for sig in ['official','reproduced']:
  e=d[(d.signature==sig)&d.target.isin(['IFN1','IFNG'])].copy()
  e['other']=np.where(e.target=='IFN1',e.z_IFNG,e.z_IFN1); e['ztar']=np.where(e.target=='IFN1',e.z_IFN1,e.z_IFNG)
  e['margin']=e.ztar-e.other
  for (ds,tar),g in e.groupby(['dataset','target']):
    sub,prim,pair,ov=meta[ds]; sub=sub or ('IFN-α' if tar=='IFN1' else 'IFN-γ')
    dn=g.groupby('donor').agg(ztar=('ztar','median'),margin=('margin','median'),top3=('top3',lambda s: all(tar in x for x in s)))
    rows.append(dict(signature=sig,dataset=ds,stimulus=sub,cells=prim,control_design=pair,overlap=ov,contrasts=len(g),units=len(dn),
      units_target_top3_all=int(dn.top3.sum()),median_unit_z=round(dn.ztar.median(),1),
      unit_z_range=f'{dn.ztar.min():.1f} to {dn.ztar.max():.1f}',median_unit_margin_vs_other_IFN=round(dn.margin.median(),1),
      unit_margin_range=f'{dn.margin.min():.1f} to {dn.margin.max():.1f}',units_margin_positive=int((dn.margin>0).sum()),
      exact_signflip_P_two_sided=signflip(dn.margin.values) if len(dn)<=18 else None, exact_signflip_P_one_sided=signflip(dn.margin.values,False) if len(dn)<=18 else None))
r=pd.DataFrame(rows); r.to_csv('stratified_validation_per_dataset_per_unit.csv',index=False)
print(r[r.signature=='official'].drop(columns='signature').to_string())
