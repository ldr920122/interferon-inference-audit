# v1.9: Fig. 2 only, from v1.9 eligibility file. Renders Fig. 2-5 for manuscript_GenomeBiology_routeB v1.1 from saved outputs (no new inference).
from pathlib import Path
import argparse
import pandas as pd, numpy as np, matplotlib, re, os
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'pdf.fonttype':42,'font.size':7,'axes.spines.top':False,'axes.spines.right':False})
p=argparse.ArgumentParser();p.add_argument('--out',default=str(Path(__file__).resolve().parent));a=p.parse_args();O=str(Path(a.out).resolve())+'/';Path(O).mkdir(parents=True,exist_ok=True)

W=170/25.4
def save(f,n): f.tight_layout(); f.savefig(O+n+'.pdf'); f.savefig(O+n+'.png',dpi=300); plt.close(f)
cb=['#0072B2','#E69F00','#009E73','#CC79A7','#56B4E9','#D55E00','#F0E442','#999999','#000000']
# ---- Fig 2
e=pd.read_csv(Path(__file__).resolve().parent.parent/'eligibility_audit/bg_candidate_eligibility_final.csv')
def cat(r):
    v,s=str(r.verdict),str(r.reason).lower()
    if v.startswith('ELIGIBLE'): return 'Strictly eligible (Tier 1, post hoc)'
    if v.startswith('ALREADY'): return 'Already used in this project'
    if v.startswith('CLOSED'): return 'Not evaluable (GSE202515)'
    if v.startswith('DOWNGRADED'): return 'Downgraded (E-MEXP-2032)'
    if re.search(r'cell line|line:|ipsc|es line|immortal|line;|myoblast',s): return 'Cell line or single stem-cell line'
    if re.search(r'superseries|same as',s): return 'SuperSeries or duplicate'
    if re.search(r'cytosig|aybey|txcyto|hirisa',s) and not re.search(r'epitheli|kerat|kidney|endotheli',s): return 'Method-source overlap'
    if re.search(r'epitheli|kerat|kidney|endotheli|organoid|brain|non-immune|not immune|skin',s): return 'Non-immune tissue'
    if re.search(r'no ifn|no b/g|single|1 b and 1 g|no exogenous|pdc',s): return 'Single arm or no exogenous β+γ pair'
    return 'Other design failure'
e['cat']=e.apply(cat,axis=1); c=e.cat.value_counts(); c.to_csv(O+'fig2b_counts.csv')
f,ax=plt.subplots(1,2,figsize=(W,70/25.4),gridspec_kw={'width_ratios':[1,1.3]})
labs=['GEO series-level','GEO sample-level\n(226 sample-only)','ArrayExpress','CELLxGENE\n(13 of 398 flagged)','Merged GEO matrix','Screen-flagged β+γ\ncandidates','Strictly eligible\n(post hoc re-audit)']
vals=[1854,1092,173,13,2080,66,2]
ax[0].barh(range(len(vals))[::-1],vals,color=[cb[0]]*4+[cb[2],cb[1],cb[5]])
for i,v in enumerate(vals): ax[0].text(v+20,len(vals)-1-i,f'{v:,}',va='center',fontsize=6)
ax[0].set_yticks(range(len(vals))[::-1]); ax[0].set_yticklabels(labs,fontsize=6); ax[0].set_xlabel('Series, studies or collections'); ax[0].set_title('a  Search branches',loc='left')
ax[1].barh(range(len(c))[::-1],c.values,color=cb[1]); ax[1].set_yticks(range(len(c))[::-1]); ax[1].set_yticklabels(c.index,fontsize=6)
for i,v in enumerate(c.values): ax[1].text(v+0.2,len(c)-1-i,str(v),va='center',fontsize=6)
ax[1].set_xlabel('Screen-flagged candidates (n = 66)'); ax[1].set_title('b  Eligibility (2 strictly eligible)',loc='left')
save(f,'Fig2_search_and_eligibility')
