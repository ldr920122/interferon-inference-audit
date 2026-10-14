# M4: patient-level sensitivity analyses with real per-sample outputs (official signature)
import pandas as pd,numpy as np,itertools,json
from scipy import stats
out={}
L=pd.read_csv('../lupus_repro/per_sample_results.csv')
pm=L.groupby('patient').agg(titer=('titer','median'),ifn1=('IFN1_official','median'),n=('gsm','size'),n_undet=('titer',lambda s:(s<=100).sum()))
out['GSE72754_patient_table']=pm.reset_index().to_dict('records')
def exact_spearman(x,y):
    r0=stats.spearmanr(x,y).correlation; c=0;t=0
    for p in itertools.permutations(range(len(y))):
        r=stats.spearmanr(x,np.asarray(y)[list(p)]).correlation; t+=1; c+= abs(r)>=abs(r0)-1e-12
    return r0,c/t
r,p=exact_spearman(pm.titer.values,pm.ifn1.values); out['A_patient_median_exact_perm']=dict(rho=r,P_two_sided=p,n=6)
# detection limit: titer 100 = undetected. (B) treat undetected as tied lowest rank (already) ; (C) use per-patient fraction of detected samples instead of median titer
pm2=L.assign(det=(L.titer>100).astype(int)).groupby('patient').agg(fdet=('det','mean'),ifn1=('IFN1_official','median'))
r,p=exact_spearman(pm2.fdet.values,pm2.ifn1.values); out['B_fraction_detected_exact_perm']=dict(rho=r,P_two_sided=p,n=6)
# (D) patient means instead of medians
pm3=L.groupby('patient').agg(titer=('titer',lambda s: np.mean(np.log10(s))),ifn1=('IFN1_official','mean'))
r,p=exact_spearman(pm3.titer.values,pm3.ifn1.values); out['C_patient_mean_log10titer_exact_perm']=dict(rho=r,P_two_sided=p,n=6)
# (E) within-patient (centered) per-sample association: removes between-patient differences; permutation within patient
Lc=L.copy(); Lc['ltit']=np.log10(Lc.titer); Lc['lt_c']=Lc.ltit-Lc.groupby('patient').ltit.transform('mean'); Lc['i_c']=Lc.IFN1_official-Lc.groupby('patient').IFN1_official.transform('mean')
r0=stats.spearmanr(Lc.lt_c,Lc.i_c).correlation; rng=np.random.default_rng(1); c=0;B=10000
for _ in range(B):
    perm=Lc.groupby('patient').i_c.transform(lambda s: rng.permutation(s.values))
    c+= abs(stats.spearmanr(Lc.lt_c,perm).correlation)>=abs(r0)
out['D_within_patient_centered_perm']=dict(rho=r0,P_two_sided=(c+1)/(B+1),n_samples=len(Lc),n_patients=6,permutations=B)
out['ties']=dict(n_samples_undetected=int((L.titer<=100).sum()),n_patients_median_at_limit=int((pm.titer<=100).sum()))
# GSE11903 TNFA
G=pd.read_csv('../txcyto_clinical_per_patient/GSE11903_per_sample.csv'); G=G[G.condition=='lesional'] if 'condition' in G else G
gm=G.groupby(['patient','arm']).TNFA.median().reset_index()
R=gm[gm.arm=='responder'].TNFA.values; N=gm[gm.arm=='non-responder'].TNFA.values
# exact sign-flip test, responders change vs 0
obs=R.sum(); allp=[np.dot(s,np.abs(R)) for s in itertools.product([1,-1],repeat=len(R))]
out['GSE11903_TNFA_responders_signflip_exact']=dict(n=len(R),median=float(np.median(R)),P_two_sided=float(np.mean(np.abs(allp)>=abs(obs)-1e-12)))
# exact permutation responders vs non-responders (difference in medians)
allv=np.r_[R,N]; d0=np.median(R)-np.median(N); c=0;t=0
for idx in itertools.combinations(range(len(allv)),len(R)):
    m=np.zeros(len(allv),bool); m[list(idx)]=True; t+=1; c+= abs(np.median(allv[m])-np.median(allv[~m]))>=abs(d0)-1e-12
out['GSE11903_TNFA_R_vs_NR_exact_perm']=dict(n=f'{len(R)} vs {len(N)}',median_difference=float(d0),P_two_sided=c/t)
# bootstrap CI of median difference over patients
b=[np.median(rng.choice(R,len(R)))-np.median(rng.choice(N,len(N))) for _ in range(10000)]
out['GSE11903_TNFA_R_vs_NR_patient_bootstrap_CI95']=[float(np.percentile(b,2.5)),float(np.percentile(b,97.5))]
bb=[np.median(rng.choice(R,len(R))) for _ in range(10000)]
out['GSE11903_TNFA_responders_median_patient_bootstrap_CI95']=[float(np.percentile(bb,2.5)),float(np.percentile(bb,97.5))]
json.dump(out,open('patient_level_sensitivity.json','w'),indent=1,default=float)
for k,v in out.items():
    if k!='GSE72754_patient_table': print(k,v)
print(pm)
