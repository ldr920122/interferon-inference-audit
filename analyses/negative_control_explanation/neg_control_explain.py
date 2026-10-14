"""Explain negative-control IFN false positives (manuscript v1.1, Results 'On independent data...').
Inputs: official CytoSig signature.centroid; per_contrast_z.csv from the independent validation (official signature rows).
Outputs: signature_correlations.csv, negative_control_ranks.csv."""
import sys, pandas as pd
SIG = sys.argv[1] if len(sys.argv) > 1 else 'signature.centroid'
Z = sys.argv[2] if len(sys.argv) > 2 else '../independent_validation/tables/per_contrast_z.csv'
s = pd.read_csv(SIG, sep='\t', index_col=0)
c = s.corr(method='pearson')
c.loc[['IFNG', 'IFN1'], ['IL4', 'IL13', 'GMCSF', 'GCSF', 'TNFA', 'IL1B', 'IFN1', 'IFNG']].round(3).to_csv('signature_correlations.csv')
z = pd.read_csv(Z); z = z[z.signature == 'official']
cols = [k for k in z.columns if k.startswith('z_')]
rows = []
for st in ['IL4', 'IL13', 'GMCSF', 'GCSF', 'IL1B', 'IL8', 'TNF']:
    sub = z[(z.target == 'nonIFN') & (z.stim.str.upper().str.replace('-', '') == st)]
    for _, r in sub.iterrows():
        v = r[cols].astype(float).sort_values(ascending=False); idx = [k[2:] for k in v.index]
        rows.append(dict(dataset=r.dataset, stim=st, donor=r.donor, n_cytokines_z_gt2=int((v > 2).sum()),
                         top1=idx[0], rank_IFNG=idx.index('IFNG') + 1, rank_IFN1=idx.index('IFN1') + 1,
                         z_IFNG=round(r.z_IFNG, 2), z_IFN1=round(r.z_IFN1, 2)))
pd.DataFrame(rows).to_csv('negative_control_ranks.csv', index=False)
print(pd.DataFrame(rows).to_string())
