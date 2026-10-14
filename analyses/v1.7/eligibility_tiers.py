# v1.6 (GSE60236 moved from Tier 2 to Tier 3 after SuperSeries/GSM check; parent-child relations added). M5: three-tier eligibility for the 66 IFN-beta+gamma candidates (per-candidate reason retained from Additional file 2)
import pandas as pd
d=pd.read_csv('bg_candidate_eligibility_final.csv')
T2={'GSE168710':'primary MDM, 4 donors, counts, no matched-method use; already explored by this project (rule-forming), so exploratory only'}
T3={'GSE131990':'CytoSig candidate pool','GSE145647':'Aybey validation; whole blood; donor mapping unresolved','GSE145648':'SuperSeries of method-used series',
 'GSE16251':'CytoSig pool; array','GSE17762':'SuperSeries of CytoSig-pool series','GSE181897':'Aybey validation; used in this project',
 'GSE202515':'dose and culture records conflict; used by LEMBAS','GSE46903':'CytoSig pool','GSE47189':'SuperSeries of CytoSig series',
 'GSE85176':'targeted NanoString panel','GSE146338':'targeted NanoString panel; disease cohort','GSE297735':'targeted panel; n = 1',
 'E-MEXP-2032':'donor identity and pairing unresolved; two-colour array','E-MTAB-1326':'CytoSig pool','GSE306664':'HIRISA reference source',
 'GSE60236':'SuperSeries of GSE60234, GSE60235 and GSE60341 (all 2,054 GSMs are in these SubSeries); all 26 whole-transcriptome IFN-β GSMs are in GSE60234 (11) and GSE60235 (15), both CytoSig IFNB final-training sources; the remaining 385 IFN-β GSMs (GSE60341) are a 250-gene NanoString panel, and 5 of their 376 donors are also GSE60235 training donors (v1.6 sample-level check)','GSE60234':'CytoSig; no IFN-γ arm','GSE60235':'CytoSig; no IFN-γ arm','GSE60341':'targeted NanoString panel; IFN-β only',
 'GSE42924':'one IFN-β and one IFN-γ sample; no replicate design','GSE42926':'same samples as GSE42924',
 'GSE16251':'CytoSig pool; array','GSE224351':'in vivo PEG-IFN-β in patients; no IFN-γ arm; endogenous context'}
def tier(a):
    if a in T2: return 'Tier 2 (single-arm independent response detection)',T2[a]
    if a in T3: return 'Tier 3 (technical or development case only)',T3[a]
    return 'Not usable for IFN-β/γ validation',None
d[['tier','tier_note']]=d.accession.apply(lambda a: pd.Series(tier(a)))
d['tier_note']=d.tier_note.fillna(d.reason)
rel=pd.read_csv('superseries_relations_66.csv').set_index('accession')
d['geo_superseries_of']=d.accession.map(rel.superseries_of).fillna('')
d['geo_subseries_of']=d.accession.map(rel.subseries_of).fillna('')
d['related_series_in_cytosig_final_training']=d.accession.map(rel.related_in_final_training).fillna('')
d.loc[d.accession=='GSE60236','method_overlap']='cytosig (via SubSeries GSE60234, GSE60235)'
d.insert(0,'tier1_strict_donor_paired_beta_vs_gamma','no')
d.to_csv('eligibility_tiers_66.csv',index=False)
print(d.tier.value_counts()); print(d[d.tier.str.startswith('Tier')][['accession','tier','tier_note']].to_string())
