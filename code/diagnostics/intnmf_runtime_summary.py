"""Original intNMF runtime aggregation, kept separate from primary masks; no fitting."""
from pathlib import Path
import argparse
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[2];SRC=ROOT
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'build/intnmf_runtime');a=p.parse_args();OUT=a.output.resolve()
if not OUT.is_relative_to(ROOT/'build'):p.error('Choose output under build/')
OUT.mkdir(parents=True,exist_ok=True)
def read(p):return pd.read_csv(p,keep_default_na=False,na_values=[''])
def stats(v):
 v=pd.to_numeric(v,errors='coerce');v=v[np.isfinite(v)&(v>0)]
 return dict(n_timed=len(v),median_seconds=v.median(),q25_seconds=v.quantile(.25),q75_seconds=v.quantile(.75),p025_seconds=v.quantile(.025),p975_seconds=v.quantile(.975),mean_seconds=v.mean())
# Separate intNMF sensitivity timing; do not pool repaired TCGA masks with main experiment.
sp=SRC/'results/semi_synthetic/intnmf_sensitivity_metrics_by_repetition.csv';sens=read(sp)
sr=[]
for key,g in sens.groupby(['dataset','experiment','model','dropped_modality','target_rate','mechanism']):
 sr.append(dict(zip(['dataset','experiment','model','dropped_modality','target_rate','mechanism'],key),n_recorded=len(g),n_success=int((g.status=='success').sum()),n_failed=int((g.status!='success').sum()),**stats(g.loc[g.status=='success','fit_elapsed_seconds'])))
pd.DataFrame(sr).to_csv(OUT/'intnmf_sensitivity_runtime_summary.csv',index=False)
