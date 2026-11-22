import os
from pathlib import Path
import pandas as pd,pickle,hashlib
O=Path(os.environ.get('DIAGNOSTIC_OUTPUT_DIR','private_outputs/diagnostics/tcga_screening'))
m=pd.read_csv(O/'feature_modality_map.csv').set_index('feature').modality
rows=[]
for p in Path(os.environ.get('TCGA_PYTHON_FITS_DIR','private_outputs/tcga/fits_python')).rglob('*.pkl'):
 if not any(s in p.name for s in ['_CC_SF_','_MI_SF_']):continue
 f=pickle.load(open(p,'rb'));cols=f['cols'];lv=m.reindex(cols);assert lv.notna().all();rows.append(dict(engine='corrected Python saved objects',method=f['model'],style=f['style'],imputation=f.get('imp') or 0,features=len(cols),mrna=int((lv=='mrna').sum()),methylation=int((lv=='methylation').sum()),rppa=int((lv=='rppa').sum()),file=p.name,feature_set_signature='|'.join(sorted(cols))))
r=pd.read_csv(O/'saved_full_fit_feature_counts_internal.csv');a=pd.concat([r,pd.DataFrame(rows)],ignore_index=True)
a['feature_set_sha256']=a.feature_set_signature.map(lambda x:hashlib.sha256(x.encode()).hexdigest());a.drop(columns='feature_set_signature').to_csv(O/'saved_full_fit_feature_counts.csv',index=False);print(a.groupby(['engine','style','mrna','methylation','rppa']).size().to_string());print(a.groupby(['engine','style']).feature_set_sha256.nunique())
