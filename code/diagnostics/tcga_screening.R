args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(args[1]);cores <- if(length(args)>1)as.integer(args[2]) else 4L
out <- Sys.getenv('DIAGNOSTIC_OUTPUT_DIR','private_outputs/diagnostics/gpcca');dir.create(file.path(out,'screening_replay'),recursive=TRUE,showWarnings=FALSE)
source('code/methods/tcga.R',local=TRUE)
v <- readRDS(file.path(root,'tcga/variables/tcga_bootstrap_analysis_variables.rds'))
# Identify exact duplicate features in each complete-case modality without fitting.
duplicates <- list()
for(mod in unique(v$layer_vec[colnames(v$X_cc)])) {
 Z <- v$X_cc[,v$layer_vec[colnames(v$X_cc)]==mod,drop=FALSE]
 for(j in which(duplicated(t(Z)))) {
  matches <- which(vapply(seq_len(j-1L),function(k)identical(as.numeric(Z[,k]),as.numeric(Z[,j])),logical(1)))
  duplicates[[length(duplicates)+1]] <- data.frame(modality=mod,feature_keep=colnames(Z)[matches[1]],feature_duplicate=colnames(Z)[j],equal_complete_cases=TRUE,equal_full_observed=identical(as.numeric(v$X_stand[,colnames(Z)[matches[1]]]),as.numeric(v$X_stand[,colnames(Z)[j]])))
 }
}
dup <- do.call(rbind,duplicates);write.csv(dup,file.path(out,'duplicate_features.csv'),row.names=FALSE)
replay_one <- function(b) {
 path <- file.path(out,'screening_replay',sprintf('boot_%04d.rds',b))
 if(file.exists(path))return(readRDS(path))
 tryCatch({
  pilot <- file.path(out,paste0('pilot_boot_',b,'.rds'))
  if(file.exists(pilot))r <- readRDS(pilot) else {
   start <- proc.time()[['elapsed']]
   set.seed(100L+b);ids <- sample(rownames(v$X_stand),nrow(v$X_stand),replace=TRUE)
   X <- v$X_stand[ids,,drop=FALSE];rownames(X) <- paste0(ids,'_b',b,'_',seq_along(ids))
   lab <- v$diag_labels[ids];names(lab) <- rownames(X)
   is_cc <- ids %in% v$common_ids
   X <- X[is_cc,,drop=FALSE];X <- X[,colSums(is.na(X))==0,drop=FALSE]
   # Exact preserved sel_SIS call. Capture console only; do not change arguments.
   logs <- capture.output(selected <- sel_SIS(X,lab[rownames(X)]))
   r <- list(boot=b,sampling_seed=100L+b,sampled_ids=ids,cc_ids=ids[is_cc],selected=selected,n_cc=nrow(X),n_cc_unique=length(unique(ids[is_cc])),counts=table(factor(v$layer_vec[selected],levels=c('mrna','methylation','rppa'))),screening_seconds=proc.time()[['elapsed']]-start)
  }
  r$status <- 'reconstructed';r$error <- NA_character_
  saveRDS(r,path);r
 },error=function(e) {
  r <- list(boot=b,status='replay_failed',error=conditionMessage(e),selected=character(),counts=rep(NA_integer_,3),n_cc=NA_integer_,n_cc_unique=NA_integer_,sampling_seed=100L+b,screening_seconds=NA_real_)
  saveRDS(r,path);r
 })
}
# Batch checkpointing and bounded independent parallel screening. No integration
# model, imputation or original benchmark result is run/written in this replay.
for(start in seq.int(1L,1000L,by=40L)) {
 ids <- seq.int(start,min(1000L,start+39L))
 values <- parallel::mclapply(ids,replay_one,mc.cores=cores,mc.preschedule=FALSE,mc.set.seed=FALSE)
 cat('SCREENING REPLAY',min(ids),'-',max(ids),'complete\n');flush.console()
}
all <- lapply(1:1000,replay_one)
rows <- lapply(all,function(r) data.frame(boot=r$boot,provenance=r$status,sampling_seed=r$sampling_seed,n_complete_rows=r$n_cc,n_unique_complete_patients=r$n_cc_unique,mrna=as.integer(r$counts[1]),methylation=as.integer(r$counts[2]),rppa=as.integer(r$counts[3]),total_features=if(r$status=='reconstructed')length(r$selected) else NA_integer_,gpcca_number_nonempty_blocks=if(r$status=='reconstructed')length(unique(v$layer_vec[r$selected])) else NA_integer_,duplicate_methylation_pair_retained=if(r$status=='reconstructed')all(unlist(dup[dup$modality=='methylation',c('feature_keep','feature_duplicate')]) %in% r$selected) else NA,screening_seconds=r$screening_seconds,error=r$error))
d <- do.call(rbind,rows);write.csv(d,file.path(out,'bootstrap_feature_counts_reconstructed.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'replay_sessionInfo.txt'))
cat('COMPLETE',sum(d$provenance=='reconstructed'),'reconstructions; RPPA-zero',sum(d$rppa==0,na.rm=TRUE),'duplicate pair retained',sum(d$duplicate_methylation_pair_retained,na.rm=TRUE),'\n')
print(summary(d[c('mrna','methylation','rppa','n_complete_rows','n_unique_complete_patients')]))
