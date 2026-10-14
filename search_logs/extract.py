import json,re,csv,collections
cand=json.load(open('candidates.json'))
s1={r['accession']:r for r in csv.DictReader(open('../step1/ifn_public_source_overlap_matrix.csv'))}
B=re.compile(r'(IFN[- ]?(beta|β|b\b|b1|B1a|b1a)|interferon[- ]?beta|IFNB1?\b)',re.I)
G=re.compile(r'(IFN[- ]?(gamma|γ|g\b|y\b)|interferon[- ]?gamma|IFNG\b|IFNy)',re.I)
C=re.compile(r'(unstim|untreated|control|\bPBS\b|mock|medium|vehicle|buffer|naive|naïve|baseline|\bNT\b|none)',re.I)
D=re.compile(r'^(donor|individual|subject|patient|donor id|donor_id|participant|person|volunteer)\s*:\s*(.+)$',re.I)
rows=[]
for a in cand:
  st=open(f'soft/{a}_self.txt',errors='ignore').read(); gt=open(f'soft/{a}_gsm.txt',errors='ignore').read()
  def f(k,t=st): return [l.split('=',1)[1].strip() for l in t.splitlines() if l.startswith(k)]
  samples=gt.split('^SAMPLE')[1:]
  nb=ng=nc=0; donors=set(); cells=collections.Counter(); org=collections.Counter(); doses=set()
  for s in samples:
    ch=[l.split('=',1)[1].strip() for l in s.splitlines() if l.startswith(('!Sample_characteristics_ch1','!Sample_title','!Sample_source_name_ch1'))]
    txt=' | '.join(ch)
    org.update(f('!Sample_organism_ch1',s))
    if B.search(txt): nb+=1
    if G.search(txt): ng+=1
    if C.search(txt) and not (B.search(txt) or G.search(txt)): nc+=1
    for c in ch:
      m=D.match(c)
      if m: donors.add(m.group(2).strip())
      if re.match(r'^(cell type|cell_type|celltype|tissue)\s*:',c,re.I): cells[c.split(':',1)[1].strip()]+=1
    for l in s.splitlines():
      if l.startswith('!Sample_treatment_protocol_ch1'): doses.update(re.findall(r'\d+\s*(?:U|IU|ng|µg|ug|pg)\s*/\s*m[lL]',l))
  supp=f('!Series_supplementary_file')
  rawc=any(re.search(r'count|raw|RAW\.tar|\.h5|mtx',x,re.I) for x in supp)
  r=s1.get(a,{})
  rows.append(dict(accession=a,title=(f('!Series_title') or [''])[0][:150],n_gsm=len(samples),organism=';'.join(k for k,_ in org.most_common(2)),
    gsm_beta=nb,gsm_gamma=ng,gsm_control=nc,n_donor_labels=len(donors),cell_types=';'.join(k for k,_ in cells.most_common(4))[:150],
    dose_strings=';'.join(sorted(doses))[:120],supp_files=';'.join(x.rsplit('/',1)[-1] for x in supp)[:200],raw_like_supp=rawc,
    platform=(f('!Series_type') or [''])[0][:40],pdat=(f('!Series_submission_date') or [''])[0],
    cytosig=r.get('cytosig_pool_hit',''),txcyto=r.get('txcyto_SI2',''),aybey=r.get('aybey2025',''),hirisa=r.get('hirisa_reference','')))
cy=open('../step1/cytosig_header.txt').read()
tx=open('../step1/tx_SI2.csv').read()+open('../step1/tx_SI3.csv').read()
ay=open('../step1/aybey.xml',errors='ignore').read()
for r in rows:
  a=r['accession']; r['cytosig']= 'yes' if a+'.' in cy else ''; r['txcyto']='yes' if a in tx else ''; r['aybey']='yes' if re.search(a+r'\b',ay) else ''; r['hirisa']='yes' if a=='GSE306664' else ''
rows.sort(key=lambda x:x['accession'])
w=csv.DictWriter(open('bg_candidates_auto.csv','w'),fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
for r in rows: print(r['accession'],r['n_gsm'],r['organism'][:12],'b',r['gsm_beta'],'g',r['gsm_gamma'],'c',r['gsm_control'],'d',r['n_donor_labels'],'raw' if r['raw_like_supp'] else '-',('CY' if r['cytosig'] else '')+('TX' if r['txcyto'] else '')+('AY' if r['aybey'] else '')+('H' if r['hirisa'] else ''),'|',r['cell_types'][:50],'|',r['title'][:70])
