import urllib.request,urllib.parse,json,time,datetime
E='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
Q={'IFNB':'("interferon beta" OR "IFN-beta" OR "IFN beta" OR IFNb OR IFNB1 OR "IFN-b" OR "IFNb1a" OR "IFN-B1a")',
'IFNG':'("interferon gamma" OR "IFN-gamma" OR "IFN gamma" OR IFNg OR IFNy OR IFNG OR "IFN-g" OR "IFN-y")'}
base='"Homo sapiens"[Organism] AND gsm[Entry Type]'
log={'date':datetime.datetime.now().isoformat(),'queries':{}};res={}
for k,q in Q.items():
  term=f'{q} AND {base}'
  r=json.loads(urllib.request.urlopen(E+'esearch.fcgi?'+urllib.parse.urlencode({'db':'gds','term':term,'retmax':100000,'retmode':'json'})).read())['esearchresult']
  log['queries'][k]={'term':term,'count':int(r['count']),'retrieved':len(r['idlist'])}; res[k]=r['idlist']; print(k,r['count'])
json.dump(res,open('gsm_ids.json','w'));json.dump(log,open('gsm_search_log.json','w'),indent=1)
