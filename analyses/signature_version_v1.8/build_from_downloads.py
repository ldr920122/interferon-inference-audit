from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,gzip,re,time
b=Path(__file__).parent.resolve()
import argparse
parser=argparse.ArgumentParser();parser.add_argument('--input-dir',required=True);parser.add_argument('--resource-dir',required=True);parser.add_argument('--output-dir',required=True);args=parser.parse_args()
d=Path(args.input_dir).resolve();out=Path(args.output_dir).resolve();out.mkdir(parents=True,exist_ok=True)
r=Path(args.resource_dir).resolve();old=pd.read_csv(r/'cytosig_reproduction/diff.centroid_1306',sep='\t',index_col=0);new=pd.read_csv(r/'cytosig_reproduction/diff.centroid_1307',sep='\t',index_col=0);off=pd.read_csv(d/'signature.centroid',sep='\t',index_col=0)
assert old.index.equals(new.index) and new.index.equals(off.index) and new.columns.equals(off.columns)
summary={'official_commit':'ce53cc7f','seed':20261013,'permutations':1000,'new_vs_official_max_centroid_abs':float((new-off).abs().max().max()),'n_genes':len(new),'n_channels':len(new.columns),'datasets':[]};assert summary['new_vs_official_max_centroid_abs']<1e-12
# Build exactly the pre-existing seven-dataset contrasts, without running the old top-level scoring code.
prefix=(r/'analyses/independent_validation/run_validation.py').read_text().split('# ---- run')[0]
prefix=re.sub(r"SIG=\{.*?\}\n# ensembl",'# ensembl',prefix,flags=re.S)
prefix=prefix.replace("D='data/'",f'D={str(d)+chr(47)!r}')
ns={};exec(compile(prefix,'legacy_contrast_builder','exec'),ns);contrasts=ns['contrasts']
# Add the already-used lupus case; no source-data spreadsheet is required for the signature comparison.
L=[l.rstrip('\n').split('\t') for l in gzip.open(d/'GSE72754_series_matrix.txt.gz','rt')];gsm=[x.strip('"') for x in next(l for l in L if l[0]=='!Sample_geo_accession')[1:]];meta=pd.DataFrame({'gsm':gsm})
for l in L:
 if l[0]=='!Sample_characteristics_ch1':
  for i,v in enumerate(l[1:]):
   v=v.strip('"')
   if ':' in v:k,val=v.split(':',1);meta.loc[i,k.strip()]=val.strip()
meta['month']=meta['time'].str.extract(r'(\d+)').astype(float);meta.loc[meta['time'].fillna('').str.contains('baseline'),'month']=0
i0=next(i for i,l in enumerate(L) if l[0]=='!series_matrix_table_begin');i1=next(i for i,l in enumerate(L) if l[0]=='!series_matrix_table_end')
e=pd.DataFrame([l[1:] for l in L[i0+2:i1]],index=[l[0].strip('"') for l in L[i0+2:i1]],columns=gsm).apply(pd.to_numeric,errors='coerce')
if e.max().max()>50:e=np.log2(e.clip(lower=1))
ann=pd.read_csv(d/'GPL570.annot.gz',sep='\t',skiprows=27,usecols=['ID','Gene symbol'],dtype=str).dropna();ann=ann[~ann['Gene symbol'].str.contains('///')].set_index('ID')['Gene symbol'];G=e.loc[e.index.intersection(ann.index)].groupby(ann).mean();pat=meta.dropna(subset=['patient']);base=pat[pat.month==0].set_index('patient').gsm;post=pat[(pat.month>0)&(pat.agent=='IFN-K')].copy();post['baseline_gsm']=post.patient.map(base);post=post.dropna(subset=['baseline_gsm']);tc=[c for c in meta if c.startswith('neutralizing')][0]
for row in post.itertuples():contrasts.append(dict(dataset='GSE72754',stim='IFN-K',target='clinical',donor=str(row.patient),indep='clinical',s=G[row.gsm],c=G[row.baseline_gsm],gsm=row.gsm,titer=float(meta.loc[meta.gsm==row.gsm,tc].iloc[0])))
allrows=[];comparisons=[];channels=list(new.columns);started=time.time()
for gi,ds in enumerate(dict.fromkeys(c['dataset'] for c in contrasts)):
 cs=[c for c in contrasts if c['dataset']==ds];Y=pd.concat([c['s']-c['c'] for c in cs],axis=1).dropna();g=new.index.intersection(Y.index);Y=Y.loc[g];yv=((Y-Y.mean())/Y.std()).to_numpy();ops=[]
 for s in [old,new,off]:
  X=s.loc[g];x=((X-X.mean())/X.std()).to_numpy();ops.append(np.linalg.solve(x.T@x+10000*np.eye(x.shape[1]),x.T))
 P=np.vstack(ops);coef=P@yv;sum1=np.zeros_like(coef);sum2=np.zeros_like(coef);rng=np.random.default_rng(20261013+gi)
 for it in range(1000):
  br=P@yv[rng.permutation(len(g))];sum1+=br;sum2+=br**2
 mu=sum1/1000;sd=np.sqrt(sum2/1000-mu**2);Z=(coef-mu)/sd
 zold,znew,zoff=np.split(Z,3);summary['datasets'].append(dict(dataset=ds,n_contrasts=len(cs),n_genes=len(g),new_vs_official_max_z_abs=float(abs(znew-zoff).max()),old_vs_new_max_z_abs=float(abs(znew-zold).max()),elapsed_s=round(time.time()-started,1)))
 Y.to_csv(out/f'{ds}_expression_contrasts.tsv.gz',sep='\t',index_label='gene',compression='gzip')
 for j,c in enumerate(cs):
  for name,z in [('legacy1306',zold),('resolved1307',znew),('official',zoff)]:
   v=z[:,j];o=np.argsort(-v,kind='stable');rec={k:c[k] for k in ['dataset','stim','target','donor','indep']};rec.update(signature=name,gsm=c.get('gsm',''),titer=c.get('titer',None),n_genes=len(g));rec.update({f'z_{ch}':float(v[k]) for k,ch in enumerate(channels)});rec.update(top1=channels[o[0]],top3=','.join(channels[k] for k in o[:3]));allrows.append(rec)
  a={k:c[k] for k in ['dataset','stim','target','donor']};a.update(max_abs_z_delta=float(abs(znew[:,j]-zold[:,j]).max()),IFN1_delta=float(znew[channels.index('IFN1'),j]-zold[channels.index('IFN1'),j]),IFNG_delta=float(znew[channels.index('IFNG'),j]-zold[channels.index('IFNG'),j]),full_top1_changed=bool(np.argmax(znew[:,j])!=np.argmax(zold[:,j])),IFN_threshold_changed=bool(((znew[channels.index('IFN1'),j]>2)|(znew[channels.index('IFNG'),j]>2))!=((zold[channels.index('IFN1'),j]>2)|(zold[channels.index('IFNG'),j]>2))),margin_sign_changed=bool((znew[channels.index('IFN1'),j]>znew[channels.index('IFNG'),j])!=(zold[channels.index('IFN1'),j]>zold[channels.index('IFNG'),j])));comparisons.append(a)
 print(ds,len(cs),summary['datasets'][-1],flush=True)
pd.DataFrame(allrows).to_csv(out/'fixed_permutation_all_channels.csv',index=False);pd.DataFrame(comparisons).to_csv(out/'legacy1306_vs_resolved1307_per_contrast.csv',index=False)
summary['n_total_contrasts']=len(comparisons);summary['sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [r/'cytosig_reproduction/diff.centroid_1306',r/'cytosig_reproduction/diff.centroid_1307',d/'signature.centroid']};(out/'summary.json').write_text(json.dumps(summary,indent=2));print('COMPLETE',summary,flush=True)
