"""Rebuild data.js from the current leakage_index.csv, without network or legacy build paths."""
from pathlib import Path
import pandas as pd,json
r=Path(__file__).parent
old=(r/'data.js').read_text();prefix=old.split('=',1)[0]
d=pd.read_csv(r/'leakage_index.csv',keep_default_na=False)
(r/'data.js').write_text(prefix+'='+json.dumps(d.to_dict('records'),ensure_ascii=False,separators=(',',':'))+';\n')
