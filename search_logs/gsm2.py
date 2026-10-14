import urllib.request,urllib.parse,json,time
E='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
res=json.load(open('gsm_ids.json'));out={}
allids=sorted(set(res['IFNB'])|set(res['IFNG']))
for j in range(0,len(allids),400):
  d=urllib.parse.urlencode({'db':'gds','id':','.join(allids[j:j+400]),'retmode':'json'}).encode()
  for t in range(3):
    try: r=json.loads(urllib.request.urlopen(E+'esummary.fcgi',d,timeout=60).read())['result'];break
    except Exception as e: time.sleep(3)
  for u in r.get('uids',[]): s=r[u]; out[u]={'gsm':s['accession'],'gse':s.get('gse',''),'title':s['title'],'tax':s.get('taxon','')}
  time.sleep(0.4)
for k in res: 
  for u in res[k]:
    if u in out: out[u].setdefault('hits',[]).append(k)
json.dump(out,open('gsm_summ.json','w')); print(len(out))
