"""Fixed-contract replay of the v1.9 post hoc cases; no fitting/tuning of new signatures.
Requires Python, pandas and numpy. Input names are listed in input_sha256.json.
Output summaries remain donor paired, separately by scale, stimulus, source and seed.
"""
from pathlib import Path
import argparse,json,hashlib,sys,time
import numpy as np,pandas as pd

def infer(X,Y,seed,alpha=10000,nrand=1000):
 X=(X-X.mean())/X.std(ddof=1);Y=(Y-Y.mean())/Y.std(ddof=1)
 x=X.to_numpy();y=Y.to_numpy();p=np.linalg.solve(x.T@x+alpha*np.eye(x.shape[1]),x.T);b=p@y
 rng=np.random.default_rng(seed);s1=np.zeros_like(b);s2=np.zeros_like(b)
 for _ in range(nrand):
  z=p@y[rng.permutation(len(y))];s1+=z;s2+=z*z
 mu=s1/nrand;sd=np.sqrt(s2/nrand-mu*mu)
 return pd.DataFrame((b-mu)/sd,index=X.columns,columns=Y.columns),pd.DataFrame(b,index=X.columns,columns=Y.columns)

def run(inp,out,seeds):
 out.mkdir(parents=True,exist_ok=True)
 h=pd.read_csv(inp/'hgnc_complete_set.txt',sep='\t',usecols=['symbol','ensembl_gene_id'],dtype=str).dropna()
 counts=h.ensembl_gene_id.value_counts();h=h[h.ensembl_gene_id.map(counts)==1];e2s=dict(zip(h.ensembl_gene_id,h.symbol))
 h.to_csv(out/'unique_ensembl_hgnc_map.tsv',sep='\t',index=False)
 sig=pd.read_csv(inp/'signature.centroid',sep='\t',index_col=0)
 fc=pd.read_csv(inp/'MDM-IFN.featurecounts.tsv',sep='\t',comment='#',index_col=0).iloc[:,5:]
 fc.columns=[c.split('/')[-1].replace('_filter.bam','').replace('-','.') for c in fc.columns]
 md=pd.read_csv(inp/'MDM-IFN.metadata.csv').set_index('libID');assert set(md.index)==set(fc.columns) and fc.shape==(60664,30)
 pc=pd.read_csv(inp/'GSE237879_counts.txt.gz',sep='\t',index_col=0);pc.index=pc.index.str.split('.').str[0]
 pm=pd.read_csv(inp/'GSE237879_metadata.txt.gz',sep='\t').set_index('ID_library');assert set(pm.index)==set(pc.columns) and pc.shape==(57084,28)
 all_rows=[];channels=[];checks=[];inputs=[]
 for p in sorted(inp.iterdir()):
  if p.is_file():inputs.append({'file':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size})
 for ds,c,m,dk,sk,ctrl,stims in [
  ('GSE236156',fc,md,'ptID','IFN','untreated',['IFNb','IFNg','IFNA8','IFNe']),
  ('GSE237879',pc,pm,'Donor','Stim','medium',['beta','gamma','alpha','lambda'])]:
  yy={}
  for scale,label in [(1e4,'CP10k'),(1e6,'CPM')]:
   lg=np.log2(c/c.sum()*scale+1);lg=lg[lg.index.isin(e2s)];lg.index=lg.index.map(e2s);lg=lg.groupby(level=0).mean()
   contrasts={};metadata=[]
   for donor in sorted(m[dk].unique()):
    dm=m[m[dk]==donor];controls=dm.index[dm[sk]==ctrl];assert len(controls)==1
    for stim in stims:
     ss=dm.index[dm[sk]==stim];assert len(ss)==1
     contrasts[(donor,stim)]=lg[ss[0]]-lg[controls[0]]
     metadata.append({'dataset':ds,'donor':donor,'stim':stim,'stim_library':ss[0],'control_library':controls[0]})
   y=pd.DataFrame(contrasts);genes=sig.index.intersection(y.index);y=y.loc[genes];yy[label]=y
   y.to_csv(out/f'{ds}_{label}_processed_contrasts.csv.gz',compression='gzip')
   pd.DataFrame(metadata).to_csv(out/f'{ds}_library_pairs.csv',index=False)
   print(ds,label,'genes',len(genes),'contrasts',len(y.columns),flush=True)
   for seed in seeds:
    z,b=infer(sig.loc[genes],y,seed)
    for donor,stim in y.columns:
     v=z[(donor,stim)];order=list(v.sort_values(ascending=False).index)
     row=dict(dataset=ds,scale=label,seed=seed,donor=donor,stim=stim,n_genes=len(genes),z_IFN1=v['IFN1'],z_IFNG=v['IFNG'],IFNG_minus_IFN1=v['IFNG']-v['IFN1'],rank_IFN1=order.index('IFN1')+1,rank_IFNG=order.index('IFNG')+1,top1=order[0])
     all_rows.append(row)
     for k in sig.columns:channels.append({'dataset':ds,'scale':label,'seed':seed,'donor':donor,'stim':stim,'channel':k,'z':v[k],'ridge_coefficient':b.loc[k,(donor,stim)]})
    print(ds,label,'seed',seed,'done',flush=True)
  # A unit change with matched pseudocount is algebraically invariant after subtraction.
  a=np.log2(c/c.sum()*1e6+1);u=np.log2(c/c.sum()*1e4+0.01)
  residual=float(((a-u)-np.log2(100)).abs().to_numpy().max())
  checks.append({'dataset':ds,'check':'matched_pseudocount_unit_conversion_identity','max_abs_residual':residual,'tolerance':1e-10,'pass':residual<1e-10})
 rows=pd.DataFrame(all_rows);rows.to_csv(out/'per_donor_all_seeds.csv',index=False)
 pd.DataFrame(channels).to_csv(out/'all43_all_seeds.csv.gz',index=False,compression='gzip')
 primary=rows[rows.seed==20261014];primary.to_csv(out/'per_donor_primary.csv',index=False)
 pd.DataFrame(channels).query('seed==20261014').to_csv(out/'all43_primary.csv',index=False)
 dirs=[]
 for (ds,sc,se,st),g in rows[rows.stim.isin(['IFNb','IFNg','beta','gamma'])].groupby(['dataset','scale','seed','stim']):
  wanted=1 if st in ['IFNg','gamma'] else -1
  dirs.append({'dataset':ds,'scale':sc,'seed':int(se),'stim':st,'n_donors':len(g),'correct_relative_order':int((g.IFNG_minus_IFN1*wanted>0).sum()),'target_rank1':int(((g.rank_IFNG if wanted==1 else g.rank_IFN1)==1).sum()),'target_top3':int(((g.rank_IFNG if wanted==1 else g.rank_IFN1)<=3).sum())})
 pd.DataFrame(dirs).to_csv(out/'direction_counts.csv',index=False)
 json.dump({'inputs':inputs,'seed_primary':20261014,'alternative_seeds':[s for s in seeds if s!=20261014],'alpha':10000,'permutations':1000,'signature_genes':len(sig),'signature_channels':len(sig.columns),'software':{'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__},'checks':checks},open(out/'replay_contract_and_checks.json','w'),indent=2)
 print('REPLAY COMPLETE',flush=True)

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--inputs',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--seeds',type=int,nargs='+',default=[20261014,1,2,3,4,5]);a=ap.parse_args();run(a.inputs,a.out,a.seeds)
