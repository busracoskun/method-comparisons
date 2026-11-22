"""Recalculate feature-recovery summaries, pooled text results, and S4/S7 checks."""
from pathlib import Path
import argparse

import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[2];PUB=ROOT
parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=ROOT/'build/feature_recovery');args=parser.parse_args();OUT=args.output.resolve()
if not OUT.is_relative_to(ROOT/'build'):parser.error('Choose output under build/')
OUT.mkdir(parents=True,exist_ok=True)
def read(p):return pd.read_csv(p,keep_default_na=False,na_values=['','NaN','nan','NULL'])
p=PUB/'results/semi_synthetic/recovery_metrics_by_repetition.csv';d=read(p)
d['target_rate']=d.target_rate.round(1)
assert not d.duplicated(['dataset','model','dropped_modality','mechanism','target_rate','repetition']).any()
# The manuscript narrative averages finite repetition values. Figure S4
# instead gives each mechanism equal weight inside a plotted cell.
text_metrics = ["rv", "recon_cor", "wgene_rankcor", "wgene_jaccard"]
successful = d.loc[d.status.eq("success")].copy()
successful[text_metrics] = successful[text_metrics].where(np.isfinite(successful[text_metrics]))
paired = successful.dropna(subset=["rv", "wgene_jaccard"])
base = ["dataset", "target_rate"]
groups = {
    "rate": (base, successful, text_metrics),
    "modality": (base + ["dropped_modality"], successful, text_metrics),
    "mechanism": (base + ["mechanism"], successful, text_metrics),
    "method_modality": (base + ["model", "dropped_modality"], successful, text_metrics),
    "paired_method_modality": (base + ["model", "dropped_modality"], paired, ["rv", "wgene_jaccard"]),
}
text_tables = []
for group, (keys, records, columns) in groups.items():
    for metric in columns:
        table = records.groupby(keys)[metric].agg(n="count", mean="mean").reset_index()
        text_tables.append(table.assign(group=group, metric=metric))
text_columns = ["group", *base, "model", "dropped_modality", "mechanism", "metric", "n", "mean"]
pd.concat(text_tables, ignore_index=True)[text_columns].to_csv(
    OUT / "pooled_recovery_summary.csv", index=False, float_format="%.17g"
)

metrics=['wgene_rankcor','wgene_jaccard','rv','rv_W']
rows=[]
for keys,g in d.groupby(['dataset','model','dropped_modality','target_rate','mechanism']):
 r=dict(zip(['dataset','model','modality','rate','mechanism'],keys));r.update(n_recorded=len(g),n_success=int(g.status.eq('success').sum()))
 for metric in metrics:
  v=g.loc[g.status.eq('success'),metric].dropna();v=v[np.isfinite(v)]
  r.update({metric:v.mean(),f'n_{metric}':len(v),f'q025_{metric}':v.quantile(.025),f'q975_{metric}':v.quantile(.975)})
 rows.append(r)
by=pd.DataFrame(rows);by.to_csv(OUT/'recovery_by_mechanism.csv',index=False)
pooled=[]
for keys,g in by.groupby(['dataset','model','modality','rate']):
 r=dict(zip(['dataset','model','modality','rate'],keys));r.update(n_recorded=int(g.n_recorded.sum()),n_success=int(g.n_success.sum()))
 for metric in metrics:
  r[metric]=g[metric].mean() if len(g)==3 and g[metric].notna().all() else np.nan
  r[f'n_{metric}']=g[f'n_{metric}'].sum()
 pooled.append(r)
pool=pd.DataFrame(pooled);pool.to_csv(OUT/'feature_recovery_summary.csv',index=False)
wide=pool.pivot(index=['dataset','model','modality'],columns='rate',values=metrics)
wide.columns=[f'{m}_{int(rate*100)}' for m,rate in wide.columns]
for m in metrics:wide[f'{m}_change_70_minus_20']=wide[f'{m}_70']-wide[f'{m}_20']
wide.reset_index().to_csv(OUT/'feature_recovery_20_vs_70.csv',index=False)
# Independent reconstruction of every figure table value and its denominators.
s4=read(PUB/'results/semi_synthetic/feature_recovery_plotted_values.csv');checks=[]
for r in s4.itertuples():
 g=by[(by.dataset==r.dataset)&(by.model==r.model)&(by.modality==r.modality)&np.isclose(by.rate,r.rate)]
 actual=g[r.metric].mean() if len(g)==3 and g[r.metric].notna().all() else np.nan
 for field,expected,actual in [('value',r.value,actual),('mechanisms',r.mechanisms,int((g[f'n_{r.metric}']>0).sum())),('successful_metric_records',r.successful_metric_records,int(g[f'n_{r.metric}'].sum()))]+[(f'mean_{mech}',getattr(r,f'mean_{mech}'),g.loc[g.mechanism==mech,r.metric].iloc[0]) for mech in ['MCAR','MAR','MNAR']]+[(f'n_{mech}',getattr(r,f'n_{mech}'),g.loc[g.mechanism==mech,f'n_{r.metric}'].iloc[0]) for mech in ['MCAR','MAR','MNAR']]:
  checks.append(dict(figure='S4',dataset=r.dataset,model=r.model,modality=r.modality,rate=r.rate,metric=r.metric,field=field,expected=expected,recalculated=actual,difference=actual-expected,passed=bool(np.isclose(expected,actual,atol=1e-12,rtol=0,equal_nan=True))))
s7=read(PUB/'results/semi_synthetic/patient_feature_recovery_plotted_values.csv')
for r in s7.itertuples():
 g=d[(d.dataset==r.dataset)&(d.model==r.model)&(d.dropped_modality==r.modality)&np.isclose(d.target_rate,r.rate)&(d.mechanism==r.mechanism)&(d.status=='success')].dropna(subset=['rv','wgene_jaccard'])
 for field,actual in [('rv',g.rv.mean()),('jaccard',g.wgene_jaccard.mean()),('n',len(g))]:
  expected=getattr(r,field);checks.append(dict(figure='S7',dataset=r.dataset,model=r.model,modality=r.modality,rate=r.rate,mechanism=r.mechanism,metric=field,field=field,expected=expected,recalculated=actual,difference=actual-expected,passed=bool(np.isclose(expected,actual,atol=1e-12,rtol=0))))
check=pd.DataFrame(checks);check.to_csv(OUT/'feature_recovery_numerical_verification.csv',index=False);assert check.passed.all()
