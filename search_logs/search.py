import urllib.request,urllib.parse,json,time,datetime
E='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
base='"Homo sapiens"[Organism] AND gse[Entry Type] AND ("expression profiling by high throughput sequencing"[DataSet Type] OR "expression profiling by array"[DataSet Type])'
Q={
'IFNB':'("interferon beta" OR "IFN-beta" OR "IFN beta" OR IFNb OR IFNB1 OR "IFN-b" OR "interferon-beta")',
'IFNG':'("interferon gamma" OR "IFN-gamma" OR "IFN gamma" OR IFNg OR IFNy OR IFNG OR "IFN-g" OR "interferon-gamma")',
'IFNA':'("interferon alpha" OR "IFN-alpha" OR "IFN alpha" OR IFNa OR IFNA2 OR "IFN-a" OR "interferon-alpha")',
'IFNL':'("interferon lambda" OR "IFN-lambda" OR IFNL1 OR IFNL2 OR IFNL3 OR "IL-29" OR "IFN-l")'}
stim='(stimulat* OR treated OR treatment OR "in vitro")'
log={'date':datetime.datetime.now().isoformat(),'endpoint':E,'queries':{}}
ids={}
for k,q in Q.items():
  term=f'{q} AND {stim} AND {base}'
  p=urllib.parse.urlencode({'db':'gds','term':term,'retmax':5000})
  r=json.loads(urllib.request.urlopen(E+'esearch.fcgi?'+p+'&retmode=json').read())['esearchresult']
  log['queries'][k]={'term':term,'count':int(r['count']),'retrieved':len(r['idlist'])}
  for i in r['idlist']: ids.setdefault(i,set()).add(k)
  time.sleep(0.5)
uids=list(ids); summ={}
for j in range(0,len(uids),300):
  d=urllib.parse.urlencode({'db':'gds','id':','.join(uids[j:j+300]),'retmode':'json'}).encode()
  r=json.loads(urllib.request.urlopen(E+'esummary.fcgi',d).read())['result']
  for u in r.get('uids',[]): summ[u]=r[u]
  time.sleep(0.5)
rows=[]
for u,s in summ.items():
  rows.append({'accession':s['accession'],'n_samples':s['n_samples'],'pdat':s['pdat'],'gdstype':s['gdstype'],'title':s['title'],'summary':s['summary'],'query_hits':'|'.join(sorted(ids[u]))})
json.dump(rows,open('geo_raw.json','w'));log['unique_GSE']=len(rows);json.dump(log,open('search_log.json','w'),indent=1)
print(json.dumps(log,indent=1)[:2000])
