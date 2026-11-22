import os
"""Read-only audit. Public outputs contain aggregate counts, never subject identifiers."""
import csv, pathlib, collections, hashlib, json, statistics
ROOT=pathlib.Path(os.environ.get('DIAGNOSTIC_OUTPUT_DIR','private_outputs/diagnostics/gipca_bidifac'))
for part in ['private','public']:(ROOT/part).mkdir(parents=True,exist_ok=True)
P=pathlib.Path(os.environ['BENCH_ARCHIVE_DIR'])
def read(p):
 with open(p,newline='') as f:return list(csv.DictReader(f))
def write(p,rows):
 if not rows:return
 with open(p,'w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
summary=[];types=[];manifest=[];assoc=[];checks=[];metrics=[]
for ds,method,combined in [('adni','gipca',P/'results_boot_adni/all_results.csv'),('tcga','bidifac',P/'results_boot_tcga/results_all_tcga.csv')]:
 rows=[r for r in read(combined) if r['model']==method and r['style']=='NA']
 ids=[int(r['boot']) for r in rows]; assert len(ids)==len(set(ids))==1000 and set(ids)==set(range(1,1001))
 original=[];files=[]
 for f in sorted((P/('results_boot' if ds=='adni' else 'results_boot_tcga')).glob('boot_NA_*.csv' if ds=='adni' else 'boot_TCGA_NA_*.csv')):
  rr=[r for r in read(f) if r['model']==method];original+=rr; files.append(f)
 # TCGA style-specific naming is inspected rather than guessed when absent.
 if not original:
  for f in sorted((P/'results_boot_tcga').glob('boot_TCGA_ALL_styles_????-????.csv')):
   original += [r for r in read(f) if r['model']==method and r['style']=='NA'];files.append(f)
 orig={int(r['boot']):r for r in original};assert len(original)==len(orig)==1000
 fields=['status','error_message','n_attempted','n_successful','n_failed']
 norm=lambda x:'' if x in ('NA','') else x
 mismatch=sum(any(norm(r[k])!=norm(orig[int(r['boot'])][k]) for k in fields) for r in rows)
 assert mismatch==0
 checks.append(dict(dataset=ds,method=method,records=1000,unique_bootstrap_ids=1000,missing_ids=0,duplicate_ids=0,original_chunk_files=len(files),chunk_combined_status_mismatches=mismatch))
 for f in [combined,*files,P/('final.R' if ds=='adni' else 'final_tcga.R'),P/'vendor/bidifac/BIDIFAC.R',P/'variables'/('analysis_variables.rds' if ds=='adni' else 'tcga_bootstrap_analysis_variables.rds')]:
  manifest.append(dict(path=str(f),sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
 c=collections.Counter((r['status'],norm(r['error_message'])) for r in rows)
 for (status,err),n in c.items():types.append(dict(dataset=ds,method=method,status=status,count=n,error_message=err))
 summary.append(dict(dataset=ds,method=method,style='NA',recorded=1000,success=sum(r['status']=='success' for r in rows),failed=sum(r['status']=='failed' for r in rows),missing=0,convergence_verified='not recorded'))
 write(ROOT/'private'/f'{ds}_status_records.csv',rows)
 for status in ['success','failed']:
  rr=[r for r in rows if r['status']==status];times=[float(r['elapsed_seconds']) for r in rr]
  metrics.append(dict(dataset=ds,status=status,n=len(rr),median_elapsed_seconds=statistics.median(times),success_rows_missing_any_metric=sum(any(r.get(k,'') in ('','NA','NaN') for k in ['rand_diag','ari_diag','chi_diag','db_diag','pve_diag','R2','rand_miss','ari_miss','chi_miss','db_miss','pve_miss']) for r in rr)))
 inp=ROOT/'private'/f'{ds}_bootstrap_input_checks.csv'
 if inp.exists():
  ir={int(r['boot']):r for r in read(inp)}
  for status in ['success','failed']:
   rr=[ir[int(r['boot'])] for r in rows if r['status']==status]
   for k in rr[0]:
    if k in ('dataset','boot'):continue
    a=[float(r[k]) for r in rr];assoc.append(dict(dataset=ds,status=status,property=k,n=len(a),minimum=min(a),median=statistics.median(a),maximum=max(a),nonzero_count=sum(x!=0 for x in a)))
write(ROOT/'public/failure_signatures.csv',types);write(ROOT/'public/reliability_counts.csv',summary);write(ROOT/'public/record_validation.csv',checks);write(ROOT/'public/timing_and_metric_availability.csv',metrics);write(ROOT/'public/input_associations.csv',assoc)
(ROOT/'private/source_sha256.json').write_text(json.dumps(manifest,indent=2))
print(json.dumps(checks,indent=2))

# Reproducible sparse-support cross-tabulation, aggregate only.
support_path=ROOT/'private/adni_unique_support.csv'
if support_path.exists():
    support={r['boot']:r for r in read(support_path)}
    records=read(ROOT/'private/adni_status_records.csv')
    tab=collections.Counter((r['status'],int(support[r['boot']]['min_unique_observations'])) for r in records)
    write(ROOT/'public/sparse_support_by_status.csv',[dict(status=status,minimum_distinct_observed_subjects_per_feature=n,bootstrap_count=count) for (status,n),count in sorted(tab.items())])
