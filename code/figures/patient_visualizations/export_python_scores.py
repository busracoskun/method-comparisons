#!/usr/bin/env python3
"""Export the latest TCGA Python full-run scores, preserving patient IDs.
Only repository-owned model artifacts are unpickled.
"""
from pathlib import Path
import pickle
import numpy as np
import pandas as pd
root=Path(__file__).resolve().parents[2]
for strategy in ['NA','MI_SF']:
    suffix='_imp1' if strategy=='MI_SF' else ''
    source=root/f'FullRuns_TCGA_python/mofa/mofa_{strategy}_0{suffix}.pkl'
    if not source.exists(): continue
    with source.open('rb') as f: fit=pickle.load(f)
    z=fit['Z']
    assert isinstance(z,pd.DataFrame) and z.index.is_unique and np.isfinite(z.to_numpy()).all()
    assert z.index.astype(str).tolist()==list(map(str,fit['rows']))
    z.rename_axis('patient_id').to_csv(root/f'figures_revised/cache/TCGA_mofa_{strategy}_latest_Z.csv')
