from pathlib import Path
import numpy as np,pandas as pd,itertools,json
from scipy import stats
r=Path(__file__).parent.resolve();src=r.parent;out=r/'clinical_statistics_recomputed';out.mkdir(exist_ok=True)
l=pd.read_csv(src/'GSE72754_lupus/per_sample_results.csv');pm=l.groupby('patient')[['titer','IFN1_official']].median()
def exactcorr(x,y):
 v=stats.spearmanr(x,y).statistic;z=[stats.spearmanr(x,np.asarray(y)[list(p)]).statistic for p in itertools.permutations(range(len(y)))];return float(v),float((np.abs(z)>=abs(v)-1e-12).mean())
a=stats.spearmanr(l.titer,l.IFN1_official);b,p=exactcorr(pm.titer,pm.IFN1_official)
results={'GSE72754':{'n_samples':len(l),'n_patients':len(pm),'per_sample_rho':float(a.statistic),'per_sample_P':float(a.pvalue),'patient_median_rho':b,'patient_exact_P':p,'n_exact_permutations':720,'bootstrap_definition_audit':'Frozen 95% interval -0.72 to 0.33 is a patient-cluster bootstrap of the per-sample correlation. It is not a confidence interval for the six-patient median correlation. No frozen numeric result was changed.'}}
assert len(l)==36 and len(pm)==6 and abs(p-.2)<1e-12
clin=src/'txcyto_clinical_per_patient';old=pd.read_csv(clin/'results_per_sample_vs_per_patient.csv');new=[]
config=[('GSE100093',['IFNG'],'arm','AMG811','Placebo'),('GSE130588',['IL4','IL13'],'treatment','Dupilumab','Placebo'),('GSE59294',['IL4','IL13'],'arm','Dupilumab',None),('GSE11903',['TNFA','IL17A'],'arm','responder','non-responder'),('GSE80060',['IL1B'],'arm','Canakinumab','Placebo'),('GSE99802',['IL22'],'arm','Drug','Placebo')]
for ds,targets,ac,t,c in config:
 d=pd.read_csv(clin/f'{ds}_per_sample.csv')
 for target in targets:
  x=d[d[ac]==t];xm=x.groupby('patient')[target].median()
  rr={'series':ds,'target':target,'test':'A treated change vs 0','p_per_sample':float(stats.wilcoxon(x[target]).pvalue),'p_per_patient':float(stats.wilcoxon(xm).pvalue),'n_samples':str(len(x)),'n_patients':str(len(xm))};new.append(rr)
  if c is not None:
   y=d[d[ac]==c];ym=y.groupby('patient')[target].median()
   if len(y)>1:new.append({'series':ds,'target':target,'test':f'B {t} vs {c}','p_per_sample':float(stats.mannwhitneyu(x[target],y[target]).pvalue),'p_per_patient':float(stats.mannwhitneyu(xm,ym).pvalue),'n_samples':f'{len(x)} vs {len(y)}','n_patients':f'{len(xm)} vs {len(ym)}'})
new=pd.DataFrame(new);cmp=old.merge(new,on=['series','target','test'],suffixes=('_old','_new'),validate='one_to_one');assert len(cmp)==len(old)==16
for col in ['p_per_sample','p_per_patient']:
 err=float((cmp[col+'_old']-cmp[col+'_new']).abs().max());assert err<1e-12,(col,err)
results['further_six_clinical_series']={'comparison_rows':len(cmp),'max_P_abs_difference':float(max((cmp[c+'_old']-cmp[c+'_new']).abs().max() for c in ['p_per_sample','p_per_patient'])),'n_samples_and_patient_counts_match':bool(all((cmp[c+'_old'].astype(str)==cmp[c+'_new'].astype(str)).all() for c in ['n_samples','n_patients'])),'scope':'Statistics recomputed from frozen per-sample CytoSig outputs; not a new raw-expression inference run.'}
d=pd.read_csv(clin/'GSE11903_per_sample.csv');g=d.groupby(['patient','arm']).TNFA.median().reset_index();x=g[g.arm=='responder'].TNFA.to_numpy();y=g[g.arm=='non-responder'].TNFA.to_numpy()
null=np.asarray(list(itertools.product([-1,1],repeat=len(x))))@np.abs(x);ps=float((np.abs(null)>=abs(x.sum())-1e-12).mean());obs=np.median(x)-np.median(y);v=np.r_[x,y];cc=0;tot=0
for idx in itertools.combinations(range(len(v)),len(x)):
 m=np.zeros(len(v),bool);m[list(idx)]=True;cc+=abs(np.median(v[m])-np.median(v[~m]))>=abs(obs)-1e-12;tot+=1
results['GSE11903_TNFA']={'responders_n':len(x),'nonresponders_n':len(y),'responders_exact_signflip_P':ps,'between_groups_exact_P':cc/tot,'responders_median':float(np.median(x)),'median_difference':float(obs)}
expected=json.loads((src/'v1.7/patient_level_sensitivity.json').read_text());assert abs(ps-expected['GSE11903_TNFA_responders_signflip_exact']['P_two_sided'])<1e-12;assert abs(cc/tot-expected['GSE11903_TNFA_R_vs_NR_exact_perm']['P_two_sided'])<1e-12
new.to_csv(out/'six_series_statistics_recomputed.csv',index=False);cmp.to_csv(out/'six_series_statistics_comparison.csv',index=False);pm.to_csv(out/'GSE72754_patient_medians.csv');(out/'clinical_statistics_verification.json').write_text(json.dumps(results,indent=2));print(json.dumps(results,indent=2))
