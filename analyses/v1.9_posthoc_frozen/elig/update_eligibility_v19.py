# v1.9 eligibility corrections from GSM-level re-audit (2026-10-14). Edits only listed rows; writes new files.
import pandas as pd
S='/workspace/ici/deliverables/resource_package/eligibility_audit/'; O='/workspace/ici/v19/elig/'
fix={
 'GSE236156':dict(gsm_beta=5,gsm_gamma=5,cell_type='MDM (5 donors; AM/MDM Mtb arms not used)',donors='5 (BB1, BB2, C352, TRIMA21, TRIMA22)',control='same-donor media',
   dose_time='IFNB 10 ng/mL, IFNG 10 ng/mL; 6 h',count_type='counts (GEO protein-coding; author featureCounts 60,664 rows)',
   verdict='ELIGIBLE (Tier 1; identified post hoc, v1.9)',
   reason='v1.9 correction: earlier exclusion relied on the GEO Series overall_design (IFNA1/A8/E/L1); GSM records and the paper show exogenous IFNB and IFNG arms with same-donor media controls in 5 MDM donors'),
 'GSE237879':dict(gsm_beta=4,gsm_gamma=4,cell_type='primary human plasmacytoid dendritic cells',donors='4 (d6, d8, d9, d10)',control='same-donor medium',
   dose_time='IFN-β 100 ng/mL, IFN-γ 100 ng/mL; 12 h',count_type='Salmon estimated counts (57,084 genes)',
   verdict='ELIGIBLE (Tier 1; identified post hoc, v1.9)',
   reason='v1.9 correction: earlier exclusion relied on the Series title/summary (IFN-III); GSM fields use labels "beta"/"gamma" (missed by automatic GSM count) and the author metadata gives 4 donors each with medium, IFN-β and IFN-γ'),
 'GSE188847':dict(reason='non-immune: post-mortem frontal cortex plus primary neurons exposed to IFNB/IFNG/TNF (v1.9 reason corrected from "no exogenous stimulation")'),
}
for fn in ['bg_candidate_eligibility_final.csv','eligibility_tiers_66.csv']:
    e=pd.read_csv(S+fn); 
    for acc,kv in fix.items():
        i=e.index[e.accession==acc][0]
        for k,v in kv.items():
            if k in e.columns: e[k]=e[k].astype(object); e.at[i,k]=v
        if 'tier' in e.columns and acc!='GSE188847':
            e.at[i,'tier']='Tier 1 (strict donor-paired β-versus-γ attribution)'; e.at[i,'tier1_strict_donor_paired_beta_vs_gamma']='yes (post hoc)'
            e.at[i,'tier_note']='identified in v1.9 GSM-level re-audit after the frozen analysis; scored post hoc (Figure 6)'
            e.at[i,'gsm_both_beta_and_gamma']=True
        if 'tier_note' in e.columns and acc=='GSE188847': e.at[i,'tier_note']=kv['reason']
    e.to_csv(O+fn,index=False)
    if 'tier' in e.columns: print(e.tier.str.split(r' \(').str[0].value_counts().to_string()); print('both arms',e.gsm_both_beta_and_gamma.sum())
