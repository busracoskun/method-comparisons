#!/usr/bin/env python3
"""Select the same repetition at both rates, nearest median RV by rank distance.
Selection depends only on numerical results, before any embedding is inspected.
"""
from pathlib import Path
import pandas as pd
root=Path(__file__).resolve().parents[2]
rows=[]
for ds,mods in [('ADNI',['gene','meta']),('TCGA',['methylation','rppa'])]:
    d=pd.read_csv(root/f'results_semi_synthetic_{ds.lower()}/all_results_{ds.lower()}.csv')
    for mod in mods:
        g=d[(d.model=='tica')&(d.mechanism=='MCAR')&(d.dropped_modality==mod)&(d.status=='success')]
        p=g.pivot(index='repetition',columns='target_rate',values='rv').dropna()
        score=(p.rank(pct=True)-.5).abs().sum(axis=1)
        rep=score.sort_values(kind='stable').index[0]
        z=g[g.repetition==rep].copy();z['selection_distance']=score.loc[rep];z['plot_dataset']=ds;rows.append(z)
pd.concat(rows).to_csv(root/'figures_revised/tables/semisynthetic_UMAP_selection.csv',index=False)
