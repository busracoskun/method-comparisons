"""TCGA adapter applies the historical TCGA SIS and imputation loaders."""
import os
from pathlib import Path
import python_methods_tcga as tcga
pm=tcga.base
mode=os.environ['OBSERVED_MODE']
work=Path(os.environ['BENCH_WORK_DIR'])
settings=pm.Settings(data_dir=tcga.TCGA_INPUT_DIR,
    out_dir=work/('full_python' if mode=='full' else 'bootstrap_python'),
    model_dir=work/'fits_python',save_models=os.environ['SAVE_BOOT_MODELS']=='1',
    models=tuple(os.environ['PY_METHODS'].split(',')),
    bootstrap_B=int(os.environ['B_BOOT']),bootstrap_start=int(os.environ['BOOT_START']),
    n_jobs_bootstrap=int(os.environ['BOOT_CORES']))
if mode=='full':pm.run_full_benchmark(settings)
else:pm.run_benchmark(settings)
