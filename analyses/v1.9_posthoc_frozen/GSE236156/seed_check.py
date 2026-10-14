# Seed robustness of the post-hoc direction counts (5 alternative permutation seeds); both datasets.
import subprocess,re,sys
for ds,path,script in [('GSE236156','/workspace/ici/v19/gse236156','rescore.py'),('GSE237879','/workspace/ici/v19/reaudit/GSE237879','rescore_GSE237879.py')]:
    src=open(f'{path}/{script}').read()
    for s in [1,2,3,4,5]:
        code=src.replace('seed=20261014',f'seed={s}').replace(".to_csv('per_donor",".to_csv('/tmp/x_per_donor").replace("pd.concat(allz,axis=1).T.to_csv('","pd.concat(allz,axis=1).T.to_csv('/tmp/x_").replace("open('input_sha256.json','w')","open('/tmp/x.json','w')")
        out=subprocess.run([sys.executable,'-c',code],cwd=path,capture_output=True,text=True).stdout
        print(ds,'seed',s,' | '.join(l for l in out.splitlines() if 'correct' in l))
