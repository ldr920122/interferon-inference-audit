# Renders Fig. 2-5 for manuscript_GenomeBiology_routeB v1.1 from saved outputs (no new inference).
import pandas as pd, numpy as np, matplotlib, re, os
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'pdf.fonttype':42,'font.size':7,'axes.spines.top':False,'axes.spines.right':False})
R='/workspace/ici/deliverables/'; O=R+'manuscript/figures/'
B='/workspace/ici/pkg/项目资料/Daily/2026-10-11-跨来源归因边界并表与独立成文评估/analysis-output/'
W=170/25.4
def save(f,n): f.tight_layout(); f.savefig(O+n+'.pdf'); f.savefig(O+n+'.png',dpi=300); plt.close(f)
cb=['#0072B2','#E69F00','#009E73','#CC79A7','#56B4E9','#D55E00','#F0E442','#999999','#000000']
# ---- Fig 2
e=pd.read_csv(R+'routeB/bg_candidate_eligibility_final.csv')
def cat(r):
    v,s=str(r.verdict),str(r.reason).lower()
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
labs=['GEO series-level','GEO sample-level\n(226 sample-only)','ArrayExpress','CELLxGENE\n(13 of 398 flagged)','Merged GEO matrix','Screen-flagged β+γ\ncandidates','Strictly eligible']
vals=[1854,1092,173,13,2080,66,0]
ax[0].barh(range(len(vals))[::-1],vals,color=[cb[0]]*4+[cb[2],cb[1],cb[5]])
for i,v in enumerate(vals): ax[0].text(v+20,len(vals)-1-i,f'{v:,}',va='center',fontsize=6)
ax[0].set_yticks(range(len(vals))[::-1]); ax[0].set_yticklabels(labs,fontsize=6); ax[0].set_xlabel('Series, studies or collections'); ax[0].set_title('a  Search branches',loc='left')
ax[1].barh(range(len(c))[::-1],c.values,color=cb[1]); ax[1].set_yticks(range(len(c))[::-1]); ax[1].set_yticklabels(c.index,fontsize=6)
for i,v in enumerate(c.values): ax[1].text(v+0.2,len(c)-1-i,str(v),va='center',fontsize=6)
ax[1].set_xlabel('Screen-flagged candidates (n = 66)'); ax[1].set_title('b  Eligibility (0 strictly eligible)',loc='left')
save(f,'Fig2_search_and_eligibility')
# ---- Fig 3
m=pd.read_csv(R+'routeB/ifn_public_source_overlap_matrix_v3.csv')
tiers=['T1_IFN_stim_immune','T2_IFN_stim_nonimmune','T3_sample_level_only','T0_mention_only']
roles={'CytoSig final':m.cytosig_status.eq('final_training'),'CytoSig pool only':m.cytosig_status.eq('candidate_pool_only'),
 'TxCyto SI2':m.txcyto_SI2.astype(str).eq('True'),'TxCyto SI3':m.txcyto_SI3.astype(str).eq('True'),
 'Aybey 2025':m.aybey2025.notna()&m.aybey2025.astype(str).ne('nan'),'HIRISA':m.hirisa_reference.astype(str).eq('True')}
H=pd.DataFrame({k:[int((v&m.screen_tier.eq(t)).sum()) for t in tiers] for k,v in roles.items()},index=[t.split('_',1)[1].replace('_',' ') for t in tiers])
H['Series in tier']=[int(m.screen_tier.eq(t).sum()) for t in tiers]; H.to_csv(O+'fig3_counts.csv')
f,ax=plt.subplots(1,2,figsize=(W,62/25.4),gridspec_kw={'width_ratios':[2.2,1]})
V=H.drop(columns='Series in tier'); im=ax[0].imshow(np.log10(V.values+1),cmap='Blues',aspect='auto')
for i in range(V.shape[0]):
    for j in range(V.shape[1]): ax[0].text(j,i,V.values[i,j],ha='center',va='center',fontsize=6,color='k' if V.values[i,j]<60 else 'w')
ax[0].set_xticks(range(V.shape[1])); ax[0].set_xticklabels(V.columns,rotation=30,ha='right',fontsize=6)
ax[0].set_yticks(range(V.shape[0])); ax[0].set_yticklabels([f'{a} (n = {n:,})' for a,n in zip(H.index,H['Series in tier'])],fontsize=6)
ax[0].set_title('a  GEO series by text tier and method role',loc='left')
ax[1].bar(['Candidate\npool','Final\ntraining*'],[137,134],color=[cb[7],cb[0]]); ax[1].set_ylabel('CytoSig candidate IFN-β/γ sources',fontsize=6); ax[1].tick_params(axis='x',labelsize=6)
for i,v in enumerate([137,134]): ax[1].text(i,v+2,str(v),ha='center',fontsize=6)
ax[1].set_title('b  134 of 137 CytoSig\n    candidate sources',loc='left',fontsize=7); ax[1].set_ylim(0,155)
f.tight_layout(); save(f,'Fig3_overlap')
# ---- Fig 4
z=pd.read_csv(R+'routeB/validation/tables/per_contrast_z.csv'); z=z[z.signature=='official']
pos=z[z.target.isin(['IFN1','IFNG'])].copy()
neg=z[(z.target=='nonIFN')&z.stim.str.upper().str.replace('-','').isin(['IL4','IL13','GMCSF','GCSF','IL1B','IL8','TNF'])].copy()
f,ax=plt.subplots(1,2,figsize=(W,72/25.4),gridspec_kw={'width_ratios':[1,1.6]})
for t,col in [('IFN1',cb[0]),('IFNG',cb[5])]:
    s=pos[pos.target==t]; ax[0].scatter(s.z_IFN1,s.z_IFNG,s=10,color=col,label=f'{t} stimulation (n = {len(s)})',edgecolor='k',linewidth=0.2)
lim=[-5,50]; ax[0].plot(lim,lim,'--',color='grey',lw=0.5); ax[0].set_xlabel('IFN1 z'); ax[0].set_ylabel('IFNG z'); ax[0].legend(fontsize=5.5,frameon=False)
ax[0].set_title('a  IFN stimulations (official signature)',loc='left')
order=['IL4','IL13','GMCSF','GCSF','IL1B','IL8','TNF']; neg['st']=neg.stim.str.upper().str.replace('-','')
for k,(ch,col,off) in enumerate([('z_IFN1',cb[0],-0.15),('z_IFNG',cb[5],0.15)]):
    for i,st in enumerate(order):
        v=neg[neg.st==st][ch]; ax[1].scatter(np.full(len(v),i+off),v,s=12,color=col,edgecolor='k',linewidth=0.2,marker='o' if k==0 else 's',label=(ch[2:] if i==0 else None))
ax[1].axhline(2,ls='--',color='k',lw=0.6); ax[1].set_xticks(range(len(order))); ax[1].set_xticklabels(['IL-4','IL-13','GM-CSF','G-CSF','IL-1β','IL-8','TNF'])
ax[1].set_ylabel('z'); ax[1].legend(fontsize=5.5,frameon=False); ax[1].set_title('b  Non-IFN controls; dashed line z = 2',loc='left')
save(f,'Fig4_independent_validation')
# ---- Fig 5
d=pd.read_csv(B+'source_donor_responses_long.csv',dtype={'donor':str}); fp=d[d.method=='FixedPair']
f,ax=plt.subplots(1,2,figsize=(W,70/25.4),gridspec_kw={'width_ratios':[1.2,1]})
k=0
for src,mk in [('GSE181897','o'),('GSE168710','s')]:
    for cond,col in [('beta',cb[0]),('gamma',cb[5])]:
        s=fp[(fp.source==src)&(fp.condition==cond)]
        for _,r in s.iterrows(): ax[0].plot([k,k+0.6],[r.I,r.II],color=col,lw=0.5,alpha=0.7)
        ax[0].scatter(np.full(len(s),k),s.I,s=6,color=col,marker=mk); ax[0].scatter(np.full(len(s),k+0.6),s.II,s=6,color=col,marker=mk)
        ax[0].text(k+0.3,ax[0].get_ylim()[1] if False else -1.6,f'{src}\nIFN-{"β" if cond=="beta" else "γ"} (n = {len(s)})',ha='center',fontsize=5)
        k+=1.5
ax[0].axhline(0,color='grey',lw=0.4,ls=':'); ax[0].set_xticks([]); ax[0].set_ylabel('Module response (log2CPM, vs control)')
ax[0].set_title('a  Type I (left) and IFN-γ-labelled (right) modules',loc='left',fontsize=6.5)
sm=pd.read_csv(B+'fixed_module_boundary_summary.csv'); sm=sm[(sm.method=='FixedPair')&sm.time_h.isna()]
sm=sm[sm.condition.isin(['beta','gamma','TNF','mixture','alpha'])]
lab=sm.source+' '+sm.condition
y=np.arange(len(sm))
ax[1].barh(y-0.2,sm.both_above_zero/sm.n_observations,0.4,color=cb[2],label='Both modules > 0')
ax[1].barh(y+0.2,sm.II_gt_I/sm.n_observations,0.4,color=cb[5],label='Type II > type I')
ax[1].set_yticks(y); ax[1].set_yticklabels([f'{a} (n = {n})' for a,n in zip(lab,sm.n_observations)],fontsize=5.5); ax[1].set_xlabel('Fraction of donors')
ax[1].legend(fontsize=5.5,frameon=False,loc='upper center',bbox_to_anchor=(0.4,-0.17),ncol=2); ax[1].set_title('b  Fixed modules by source and condition',loc='left',fontsize=6.5)
sm.to_csv(O+'fig5b_counts.csv',index=False)
save(f,'Fig5_boundary_cases')
print(c.to_string())
