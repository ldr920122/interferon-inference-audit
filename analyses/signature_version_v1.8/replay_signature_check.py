"""Replay the signature-version comparison from the archived processed contrasts.
This independently written checker is MIT licensed; matrices retain upstream terms.
Run: python replay_signature_check.py --out replay_outputs
"""
from pathlib import Path
import argparse,json,hashlib
import numpy as np
import pandas as pd
p=argparse.ArgumentParser();p.add_argument('--out',default='replay_outputs');a=p.parse_args();root=Path(__file__).resolve().parent;out=Path(a.out).resolve();out.mkdir(parents=True,exist_ok=True)
sigs={k:pd.read_csv(root/'matrices'/f,sep='\t',index_col=0) for k,f in [('legacy1306','diff.centroid_1306'),('resolved1307','diff.centroid_1307'),('official','signature.centroid')]}
meta=pd.read_csv(root/'contrast_metadata.csv',dtype={'donor':str,'gsm':str});reference=pd.read_csv(root/'fixed_permutation_all_channels.csv');checked=[]
for gi,(ds,m) in enumerate(meta.groupby('dataset',sort=False)):
 y=pd.read_csv(root/f'{ds}_expression_contrasts.tsv.gz',sep='\t',index_col=0);yv=((y-y.mean())/y.std()).to_numpy();op=[]
 for name,s in sigs.items():
  x=s.loc[y.index];x=((x-x.mean())/x.std()).to_numpy();op.append(np.linalg.solve(x.T@x+10000*np.eye(x.shape[1]),x.T))
 op=np.vstack(op);coef=op@yv;s1=np.zeros_like(coef);s2=np.zeros_like(coef);rng=np.random.default_rng(20261013+gi)
 for _ in range(1000):
  v=op@yv[rng.permutation(len(y))];s1+=v;s2+=v*v
 mu=s1/1000;z=(coef-mu)/np.sqrt(s2/1000-mu*mu)
 for name,zs in zip(sigs,np.split(z,3)):
  ref=reference[(reference.dataset==ds)&(reference.signature==name)];cols=['z_'+k for k in sigs[name].columns];assert len(ref)==len(m)
  delta=float(np.max(np.abs(zs.T-ref[cols].to_numpy())));assert np.allclose(zs.T,ref[cols].to_numpy(),rtol=1e-10,atol=1e-10),(ds,name,delta)
  checked.append({'dataset':ds,'signature':name,'contrasts':len(m),'max_abs_z_replay_vs_archived':delta})
print(json.dumps({'status':'PASS','comparisons':len(checked),'checked':checked},indent=2));(out/'replay_verification.json').write_text(json.dumps({'status':'PASS','comparisons':len(checked),'checked':checked},indent=2))
