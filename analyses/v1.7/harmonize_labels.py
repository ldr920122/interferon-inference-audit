# v1.7 (M3/M4): one controlled vocabulary for training identity across the leakage index, AF2 and AF3,
# corrected GSE60236 counts in AF3, and recomputed CytoSig counts json.
import pandas as pd,json,re
RP='../resource_package'; AF='../manuscript/additional_files'
L=pd.read_csv(f'{RP}/query_tool/leakage_index.csv')
cols=[l.strip() for l in open(f'{RP}/cytosig_reproduction/final_columns_resolved_1307.txt') if l.strip()]
acc_of=lambda c:c.split('@')[-1].split('.')[0]
final=set(map(acc_of,cols))
def ident(r):
    if r.accession in final: return 'cytosig_final_training'
    if r.cytosig_status=='candidate_pool_only': return 'cytosig_candidate_pool_only'
    if isinstance(r.superseries_contains_final_training,str) and r.superseries_contains_final_training: return 'superseries_contains_cytosig_final_training'
    return 'not_in_cytosig'
L['cytosig_training_identity']=L.apply(ident,axis=1)
def other(r):
    o=[]
    if isinstance(r.txcyto_SI2_cytokines,str) and r.txcyto_SI2_cytokines: o.append('txcyto_SI2')
    if isinstance(r.txcyto_SI3_targets,str) and r.txcyto_SI3_targets: o.append('txcyto_SI3')
    if isinstance(r.aybey2025_role,str) and r.aybey2025_role: o.append('aybey2025_'+r.aybey2025_role.split('/')[0])
    if r.hirisa_reference in (True,'True'): o.append('hirisa_reference')
    return ';'.join(o)
L['other_method_roles']=L.apply(other,axis=1)
L.to_csv(f'{RP}/query_tool/leakage_index.csv',index=False); L.to_csv('../leakage_tool/leakage_index.csv',index=False)
recs=L.fillna('').to_dict('records')
for q in [f'{RP}/query_tool/data.js','../leakage_tool/data.js']: open(q,'w').write('window.LEAK='+json.dumps(recs,ensure_ascii=False)+';\n')
I=L.set_index('accession')
# AF2
d=pd.read_csv(f'{AF}/Additional_file_2_eligibility_audit.csv')
d['cytosig_training_identity']=d.accession.map(I.cytosig_training_identity).fillna('not_in_cytosig')
d['other_method_roles']=d.accession.map(I.other_method_roles).fillna('')
extra={'GSE202515':'lembas_training','GSE168710':'this_project_rule_forming','GSE181897':'this_project_rule_forming'}
for a,v in extra.items(): d.loc[d.accession==a,'other_method_roles']=(d.loc[d.accession==a,'other_method_roles']+';'+v).str.strip(';')
d=d.rename(columns={'method_overlap':'method_overlap_free_text_v1.4'})
both=(pd.to_numeric(d.gsm_beta,errors='coerce')>0)&(pd.to_numeric(d.gsm_gamma,errors='coerce')>0)
d['gsm_both_beta_and_gamma']=both
fix={'CytoSig candidate pool':'CytoSig final training','CytoSig pool; array':'CytoSig final training; array','CytoSig pool':'CytoSig final training',
     'SuperSeries of CytoSig-pool series':'SuperSeries containing CytoSig final-training SubSeries','SuperSeries of CytoSig series':'SuperSeries containing CytoSig final-training SubSeries',
     'SuperSeries of method-used series':'SuperSeries containing CytoSig final-training SubSeries (GSE131990)','CytoSig; no IFN-γ arm':'CytoSig final training; no IFN-γ arm'}
d['tier_note']=d.apply(lambda r: fix.get(r.tier_note,r.tier_note) if r.cytosig_training_identity!='not_in_cytosig' else r.tier_note,axis=1)
d.to_csv(f'{AF}/Additional_file_2_eligibility_audit.csv',index=False); d.to_csv(f'{RP}/eligibility_audit/eligibility_tiers_66.csv',index=False)
# AF3
M=pd.read_csv(f'{AF}/Additional_file_3_overlap_matrix.csv')
M['cytosig_training_identity']=M.accession.map(I.cytosig_training_identity).fillna('not_in_cytosig')
M['superseries_contains_final_training']=M.accession.map(I.superseries_contains_final_training).fillna('')
i=M.accession=='GSE60236'; M.loc[i,['gsm_ifnb','gsm_ifng','sample_level_beta_and_gamma']]=[411,0,False]
M.to_csv(f'{AF}/Additional_file_3_overlap_matrix.csv',index=False); M.to_csv(f'{RP}/overlap_matrix/ifn_public_source_overlap_matrix_v3.csv',index=False)
# counts json
pool=L[L.cytosig_status.notna()]
fl=[re.split('[&@]',c)[0] for c in cols]
ifn=lambda s:isinstance(s,str) and 'IFN' in s; bg=lambda s:isinstance(s,str) and bool(re.search('IFNB|IFNG',s))
j=dict(candidate_accessions_total=len(pool),candidate_GEO=int(pool.accession.str.startswith('GSE').sum()),candidate_ArrayExpress=int((~pool.accession.str.startswith('GSE')).sum()),
  final_training_accessions_total=len(final),final_GEO=sum(a.startswith('GSE') for a in final),final_ArrayExpress=sum(not a.startswith('GSE') for a in final),
  final_columns=len(cols),IFN_pool=int(pool.cytosig_pool_ligands.apply(ifn).sum()),IFN_final=int(pool.cytosig_final_ligands.apply(ifn).sum()),
  BG_pool=int(pool.cytosig_pool_ligands.apply(bg).sum()),BG_final=int(pool.cytosig_final_ligands.apply(bg).sum()),
  note='v1.7 recomputed from leakage_index.csv and final_columns_resolved_1307.txt; v1.6 file mislabelled the 842 candidates as GEO and carried a stale final count (609)')
json.dump(j,open(f'{RP}/search_logs/cytosig_final_vs_pool_counts.json','w'),indent=1)
print(j); print(L.cytosig_training_identity.value_counts().to_dict()); print(d.cytosig_training_identity.value_counts().to_dict(), 'both arms',int(both.sum()))
print(M.cytosig_training_identity.value_counts().to_dict()); print(M.sample_level_beta_and_gamma.sum())
