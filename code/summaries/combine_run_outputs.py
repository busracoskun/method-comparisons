"""Combine completed outputs, preserving failures, literal NA style, and provenance."""
from pathlib import Path
import argparse,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from project_paths import local_path,check_output_tree
import pandas as pd

def combine(files,keys):
    frames=[]
    for p in files:
        d=pd.read_csv(p,keep_default_na=False,na_values=['','NaN','nan'])
        d['source_file']=p.name
        frames.append(d)
    if not frames:raise ValueError('No completed results found')
    data=pd.concat(frames,ignore_index=True,sort=False)
    if data.duplicated(keys).any():raise ValueError('Overlapping result keys; resolve duplicates explicitly before combining')
    return data.sort_values(keys)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('work',type=Path,help='Dataset output directory, e.g. private_outputs/adni')
    p.add_argument('analysis',choices=['observed','semisynthetic'])
    p.add_argument('--replace-intnmf-from',type=Path,help='Explicitly use matching rows from a separate clean intNMF repair output')
    p.add_argument('--publication',action='store_true',help='Also write the manuscript semi-synthetic method subset')
    a=p.parse_args();work=check_output_tree(a.work)
    if a.publication and a.analysis!='semisynthetic':p.error('--publication applies to semi-synthetic method exclusions only')
    root=Path(__file__).resolve().parents[2]
    if not work.is_relative_to(root):p.error('Outputs must be inside the clean project')
    def collect(folder):
        if a.analysis=='observed':
            files=list((folder/'bootstrap_r').glob('*.csv'))
            for name in ['results_boot_python','bootstrap_python']:
                files+=list((folder/name).glob('boot_ALL_styles*.csv'))
            return files,['boot','style','model']
        files=[]
        for name in ['semisynthetic_r','semisynthetic_python']:
            files+=list((folder/name/'tasks').glob('task_*_results.csv'))
        return files,['task_id','dataset','dropped_modality','model']
    files,keys=collect(work)
    data=combine(sorted(files),keys)
    if a.replace_intnmf_from:
        repair_work=check_output_tree(a.replace_intnmf_from)
        if repair_work.name!=work.name:raise ValueError('Main and repair directories must name the same dataset')
        repair_files,_=collect(repair_work)
        repair=combine(sorted(repair_files),keys)
        if set(repair['model'])!={'intnmf'}:
            raise ValueError('The replacement directory must contain intNMF only')
        replacement_keys=pd.MultiIndex.from_frame(repair[keys])
        original_keys=pd.MultiIndex.from_frame(data[keys])
        # Only keys actually present in the repair are replaced; no missing results invented.
        data=pd.concat([data.loc[~original_keys.isin(replacement_keys)],repair],ignore_index=True,sort=False).sort_values(keys)
    out=local_path(work/'combined');out.mkdir(parents=True,exist_ok=True)
    data.to_csv(out/f'{a.analysis}_all.csv',index=False)
    data[data['status'].ne('success')].to_csv(out/f'{a.analysis}_failures.csv',index=False)
    counts=data.groupby(['model']+(['style'] if a.analysis=='observed' else ['dropped_modality','mechanism','target_rate']),dropna=False)['status'].value_counts().rename('records').reset_index()
    counts.to_csv(out/f'{a.analysis}_counts.csv',index=False)
    if a.publication:
        reported=data.loc[data['model'].ne('intnmf')].copy()
        unsupported=reported['dataset'].str.upper().str.startswith('TCGA') & reported['model'].isin(['cia','o2pls'])
        reported=reported.loc[~unsupported]
        reported.to_csv(out/'semisynthetic_reported_methods.csv',index=False)
    print(f'Combined {len(data)} rows into {out}; missing replicates are not filled')
if __name__=='__main__':main()
