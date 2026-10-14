"""Build leakage_index.csv / data.js from CytoSig diff.merge header, reproduced final columns,
overlap matrix v3, TxCyto SI2/SI3, Aybey 2025, HIRISA. No network."""
import gzip,json,re,collections,pandas as pd
R='/workspace/ici/deliverables/routeB'
hdr=gzip.open('/workspace/cytosig_repro/CytoSig_data_filter/data/diff.merge.gz','rt').readline().rstrip('\n').split('\t')
cols=[c.strip('"') for c in hdr if '@' in c]
final=set(l.strip('"') for l in open(f'{R}/cytosig_repro/final_columns_in_families_ge5.txt').read().split('\n'))
def parse(c):
    lig=re.split(r'[&@]',c)[0]; acc=c.split('@')[-1].split('.')[0]; return lig,acc
rec=collections.defaultdict(lambda:{'pool_cols':0,'final_cols':0,'pool_lig':set(),'final_lig':set()})
for c in cols:
    l,a=parse(c); r=rec[a]; r['pool_cols']+=1; r['pool_lig'].add(l)
    if c in final: r['final_cols']+=1; r['final_lig'].add(l)
m=pd.read_csv(f'{R}/ifn_public_source_overlap_matrix_v3.csv')
tx2=pd.read_csv(f'{R}/tx_SI2.csv'); tx3=pd.read_csv(f'{R}/tx_SI3.csv')
t2=tx2.groupby('GSE').cytokine.apply(lambda s:';'.join(sorted(set(s)))).to_dict()
t3=tx3.groupby('Dataset').Target.apply(lambda s:';'.join(sorted(set(map(str,s))))).to_dict()
ay={r.accession:r.aybey2025 for r in m.itertuples() if isinstance(r.aybey2025,str)}
hir={'GSE306664'}
el=pd.read_csv(f'{R}/bg_candidate_eligibility_final.csv').set_index('accession')
meta={r.accession:r for r in m.itertuples()}
accs=sorted(set(rec)|set(meta)|set(t2)|set(t3)|set(ay)|hir|set(el.index))
rows=[]
for a in accs:
    r=rec.get(a); x=meta.get(a)
    status='final_training' if r and r['final_cols'] else ('candidate_pool_only' if r else '')
    rows.append(dict(accession=a,title=(x.title if x is not None and isinstance(x.title,str) else (el.loc[a,'title'] if a in el.index else '')),
      year=(str(x.pdat)[:4] if x is not None and isinstance(x.pdat,str) else ''),
      ifn_screen_tier=(x.screen_tier if x is not None and isinstance(x.screen_tier,str) else ''),
      cytosig_status=status, cytosig_pool_columns=r['pool_cols'] if r else 0,
      cytosig_final_columns=r['final_cols'] if r else 0,
      cytosig_pool_ligands=';'.join(sorted(r['pool_lig'])) if r else '',
      cytosig_final_ligands=';'.join(sorted(r['final_lig'])) if r else '',
      txcyto_SI2_cytokines=t2.get(a,''), txcyto_SI3_targets=t3.get(a,''),
      aybey2025_role=ay.get(a,''), hirisa_reference=a in hir,
      bg_audit_verdict=(el.loc[a,'verdict'] if a in el.index else ''),bg_audit_reason=(el.loc[a,'reason'] if a in el.index else '')))
df=pd.DataFrame(rows); df.to_csv('leakage_index.csv',index=False)
open('data.js','w').write('window.LEAK='+df.to_json(orient='records')+';\nwindow.LEAK_META='+json.dumps(dict(
  built='2026-10-13',n=len(df),cytosig_pool_columns=len(cols),cytosig_final_columns=len(final-{''}),
  final_matched=sum(c in final for c in cols)))+';')
print(len(df),len(cols),sum(c in final for c in cols)); print(df.cytosig_status.value_counts())
