import json,re,csv,collections
rows=json.load(open('geo_raw.json'))
cy=open('cytosig_header.txt').read().split('\t')
cyg=collections.defaultdict(set)
for c in cy:
  m=re.search(r'@(GSE\d+)',c)
  if m: cyg[m.group(1)].add(c.split('@')[0])
tx2=set(re.findall(r'GSE\d+',open('tx_SI2.csv').read())); tx3=set(re.findall(r'GSE\d+',open('tx_SI3.csv').read()))
ay=set(re.findall(r'GSE\d+',open('aybey.xml',errors='ignore').read()))
ay_disc={'GSE161916','GSE120844','GSE133218','GSE148058','GSE74863'}  # per prior stage report (JTM2025 discovery sets)
hirisa={'GSE306664'}
pat={'b':r'(IFN[- ]?(beta|β|b\b)|interferon[- ]beta|IFNB1?\b)','g':r'(IFN[- ]?(gamma|γ|g\b|y\b)|interferon[- ]gamma|IFNG\b)','a':r'(IFN[- ]?(alpha|α|a\b)|interferon[- ]alpha|IFNA\d*\b)','l':r'(IFN[- ]?(lambda|λ|l\b)|interferon[- ]lambda|IFNL\d\b|IL-?29)'}
stim=re.compile(r'(stimulat|treated with|treatment with|exposed to|exposure to|incubated with|in response to|primed with|cultured with)',re.I)
imm=re.compile(r'(monocyt|macrophag|PBMC|peripheral blood|dendritic|myeloid|neutrophil|T cell|B cell|NK cell|whole blood)',re.I)
out=[]
for r in rows:
  t=r['title']+' '+r['summary']
  types=''.join(k for k,p in pat.items() if re.search(p,t,re.I))
  if not types or not stim.search(t): tier='T0_mention_only'
  elif imm.search(t): tier='T1_IFN_stim_immune'
  else: tier='T2_IFN_stim_nonimmune'
  a=r['accession']
  out.append(dict(accession=a,n_samples=r['n_samples'],pdat=r['pdat'],platform_type=r['gdstype'],ifn_types_in_text=types,screen_tier=tier,
   cytosig_pool_hit=';'.join(sorted(cyg[a])) if a in cyg else '',
   txcyto_SI2=a in tx2,txcyto_SI3=a in tx3,aybey2025=('discovery' if a in ay_disc else 'validation/other' if a in ay else ''),
   hirisa_reference=a in hirisa,title=r['title'][:200]))
out.sort(key=lambda x:(x['screen_tier'],x['accession']))
w=csv.DictWriter(open('ifn_public_source_overlap_matrix.csv','w'),fieldnames=list(out[0]));w.writeheader();w.writerows(out)
C=collections.Counter
print('total',len(out),C(o['screen_tier'] for o in out))
for t in ['T1_IFN_stim_immune','T2_IFN_stim_nonimmune']:
  s=[o for o in out if o['screen_tier']==t]
  ov=[o for o in s if o['cytosig_pool_hit'] or o['txcyto_SI2'] or o['txcyto_SI3'] or o['aybey2025'] or o['hirisa_reference']]
  print(t,'n',len(s),'any_overlap',len(ov),'cytosig',sum(bool(o['cytosig_pool_hit']) for o in s),'tx',sum(o['txcyto_SI2'] or o['txcyto_SI3'] for o in s),'aybey',sum(bool(o['aybey2025']) for o in s))
  bg=[o for o in s if 'b' in o['ifn_types_in_text'] and 'g' in o['ifn_types_in_text']]
  print('  beta+gamma in text',len(bg),'of which no overlap',sum(1 for o in bg if o not in ov))
  print('  by pdat<=2020',sum(1 for o in s if o['pdat']<'2021'),'cytosig hits among those',sum(1 for o in s if o['pdat']<'2021' and o['cytosig_pool_hit']))
# known sources check
for a in ['GSE131990','GSE181897','GSE168710','GSE190594','GSE226572','GSE202515','GSE297735','GSE236156','GSE145647','GSE306664']:
  m=[o for o in out if o['accession']==a]; print(a,'in search' if m else 'NOT retrieved', m[0]['screen_tier'] if m else '', bool(cyg.get(a)))
json.dump({'cytosig_GSE':len(cyg),'tx_SI2':len(tx2),'tx_SI3':len(tx3),'aybey_all':sorted(ay)},open('reference_pools.json','w'),indent=1)
