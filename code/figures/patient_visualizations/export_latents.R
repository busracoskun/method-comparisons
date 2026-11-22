#!/usr/bin/env Rscript
# Read saved fits only; export explicitly aligned scores and deterministic UMAPs.
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) normalizePath(args[1]) else normalizePath('.')
out <- file.path(root, 'figures_revised')
.libPaths(c(file.path(Sys.getenv('HOME'),'R','4.5'), .libPaths()))
stopifnot(requireNamespace('uwot', quietly=TRUE))
manifest <- list(); failures <- list()
for(ds in c('ADNI','TCGA')) {
  vars <- readRDS(file.path(root,'variables',if(ds=='ADNI') 'analysis_variables.rds' else 'tcga_bootstrap_analysis_variables.rds'))
  # Clinical eligibility is applied before UMAP, using the saved fitted scores.
  ids_all <- rownames(vars$X_stand)
  lv <- vars$layer_vec[colnames(vars$X_stand)]
  available <- sapply(unique(unname(lv)),function(m) rowSums(!is.na(vars$X_stand[,lv==m,drop=FALSE]))>0)
  colnames(available) <- unique(unname(lv))
  patterns <- apply(available,1,function(a) if(all(a)) 'CC' else paste0('missing_',paste(names(a)[!a],collapse='+')))
  metadata <- data.frame(patient_id=ids_all,
    diagnosis=as.character(vars$diag_labels[ids_all]),missingness_pattern=patterns)
  valid_labels <- if(ds=='ADNI') c('CN','MCI','AD') else c('GBM','LGG')
  eligible <- !is.na(metadata$diagnosis) & toupper(metadata$diagnosis) %in% valid_labels
  excluded <- metadata[!eligible,,drop=FALSE]
  metadata <- metadata[eligible,,drop=FALSE]
  stopifnot(nrow(metadata)==if(ds=='ADNI') 889L else 841L)
  write.csv(excluded,file.path(out,'tables',paste0(ds,'_observed_UMAP_excluded_patients.csv')),row.names=FALSE)
  write.csv(metadata,file.path(out,'cache',paste0(ds,'_clinical.csv')),row.names=FALSE)
  folder <- if(ds=='ADNI') 'FullRuns' else 'FullRuns_TCGA'
  for(spec in list(c('mofa','NA'),c('mofa','MI_SF'),c('gpcca','NA'),c('tica','MI_SF'))) {
    method <- spec[1]; strategy <- spec[2]; suffix <- if(strategy=='MI_SF') '_imp1' else ''
    name <- paste(ds,method,strategy,sep='_')
    # TCGA Python reruns used corrected SIS and supersede older Python fits.
    candidate <- file.path(root,folder,method,paste0(method,'_',strategy,'_0',suffix,'.rds'))
    if(ds=='TCGA' && method=='mofa') {
      latest <- file.path(root,'FullRuns_TCGA_python',method,paste0(method,'_',strategy,'_0',suffix,'.rds'))
      if(file.exists(latest)) candidate <- latest
    }
    tryCatch({
      latest_scores <- file.path(out,'cache',paste0('TCGA_mofa_',strategy,'_latest_Z.csv'))
      if(ds=='TCGA' && method=='mofa' && !file.exists(latest_scores)) stop('Run python figures_revised/code/export_python_scores.py first.')
      if(ds=='TCGA' && method=='mofa' && file.exists(latest_scores)) {
        Z <- as.matrix(read.csv(latest_scores,row.names=1,check.names=FALSE))
        candidate <- file.path(root,'FullRuns_TCGA_python','mofa',paste0('mofa_',strategy,'_0',suffix,'.pkl'))
        fit <- list(rows=rownames(Z))
      } else {
        fit <- readRDS(candidate); Z <- as.matrix(fit$Z)
      }
      ids <- rownames(Z)
      if(is.null(ids)) ids <- fit$rows
      stopifnot(length(ids)==nrow(Z), !anyDuplicated(ids), all(ids %in% ids_all), all(is.finite(Z)), ncol(Z)>1)
      rownames(Z) <- ids
      source_n <- nrow(Z)
      keep <- ids %in% metadata$patient_id
      ids <- ids[keep]; Z <- Z[keep,,drop=FALSE]
      stopifnot(length(ids)==nrow(metadata),setequal(ids,metadata$patient_id))
      # Preserve fitted factor geometry: center, but do not rescale each factor.
      set.seed(20261018)
      xy <- uwot::umap(Z,n_neighbors=15,min_dist=.15,metric='euclidean',
        n_components=2,n_epochs=500,scale=FALSE,n_threads=1,n_sgd_threads=1,verbose=FALSE)
      write.csv(data.frame(patient_id=ids,Z,check.names=FALSE),file.path(out,'cache',paste0(name,'_Z.csv')),row.names=FALSE)
      write.csv(data.frame(patient_id=ids,UMAP1=xy[,1],UMAP2=xy[,2]),file.path(out,'cache',paste0(name,'_UMAP.csv')),row.names=FALSE)
      manifest[[length(manifest)+1]] <- data.frame(dataset=ds,model=method,style=strategy,imputation=if(strategy=='MI_SF') 1 else NA,
        n=nrow(Z),source_n=source_n,excluded_n=source_n-nrow(Z),factors=ncol(Z),source=sub(paste0(root,'/'),'',candidate,fixed=TRUE),seed=20261018)
    },error=function(e) {failures[[length(failures)+1]] <<- data.frame(example=name,reason=conditionMessage(e))})
  }
}
write.csv(do.call(rbind,manifest),file.path(out,'tables','observed_UMAP_manifest.csv'),row.names=FALSE)
if(length(failures)) {
  write.csv(do.call(rbind,failures),file.path(out,'tables','observed_UMAP_unavailable.csv'),row.names=FALSE)
  stop('Some observed-score exports failed; see tables/observed_UMAP_unavailable.csv.')
}
writeLines(capture.output(sessionInfo()),file.path(out,'logs','R_sessionInfo.txt'))
