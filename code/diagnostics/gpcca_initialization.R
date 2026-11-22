args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(args[1]);case <- args[2];boot <- if(length(args)>2L)as.integer(args[3]) else 0L;id <- if(boot)paste0(case,'_boot',boot) else case;out <- file.path(Sys.getenv('DIAGNOSTIC_OUTPUT_DIR','private_outputs/diagnostics/gpcca'),'fits');dir.create(out,recursive=TRUE,showWarnings=FALSE)
v <- readRDS(file.path(root,'tcga/variables/tcga_bootstrap_analysis_variables.rds'))
if(boot) {
 r <- readRDS(file.path(Sys.getenv('DIAGNOSTIC_OUTPUT_DIR','private_outputs/diagnostics/gpcca'),'screening_replay',sprintf('boot_%04d.rds',boot)))
 stopifnot(r$status=='reconstructed')
 Xall <- v$X_stand[r$sampled_ids,,drop=FALSE]
 rownames(Xall) <- paste0(r$sampled_ids,'_b',boot,'_',seq_along(r$sampled_ids))
 Xcc <- Xall[r$sampled_ids %in% v$common_ids,,drop=FALSE]
 Xcc <- Xcc[,colSums(is.na(Xcc))==0,drop=FALSE]
 v$X_cc <- Xcc;v$X_cc_sf <- Xcc[,r$selected,drop=FALSE]
}
lv <- v$layer_vec;X <- v$X_cc_sf;changed_model <- FALSE
if(case %in% c('cc_no_sis','cc_no_sis_deduplicated'))X <- v$X_cc
if(case %in% c('ccsf_add_rppa','ccsf_add_rppa_deduplicated')) {
 add <- setdiff(colnames(v$X_cc)[lv[colnames(v$X_cc)]=='rppa'],colnames(X))
 X <- cbind(X,v$X_cc[rownames(X),add,drop=FALSE])
}
dropped <- character()
if(grepl('deduplicated$',case)) {
 keep <- unlist(lapply(unique(lv[colnames(X)]),function(mod) {
  columns <- colnames(X)[lv[colnames(X)]==mod]
  columns[!duplicated(t(X[,columns,drop=FALSE]))]
 }),use.names=FALSE)
 dropped <- setdiff(colnames(X),keep);X <- X[,colnames(X) %in% keep,drop=FALSE]
}
X_list <- lapply(unique(lv[colnames(X)]),function(mod)t(X[,lv[colnames(X)]==mod,drop=FALSE]))
names(X_list) <- unique(lv[colnames(X)])
model <- GPCCA::GPCCAmodel
if(case=='ccsf_initial_ridge_diagnostic') {
 # Isolated function copies: apply the routine's existing M-step ridge formula
 # at initialization as well. Do not modify the package namespace or benchmark.
 changed_model <- TRUE
 env <- new.env(parent=asNamespace('GPCCA'))
 body0 <- body(get('EM_completeDAT',asNamespace('GPCCA')))
 parts <- as.list(body0)
 index <- which(vapply(parts,function(e)is.call(e)&&identical(e[[1]],as.name('<-'))&&identical(e[[2]],as.name('Psi0')),logical(1)))
 stopifnot(length(index)==1)
 extra <- quote(Psi0 <- Psi0 + (1/lambda - 1) * diag(diag(Psi0)))
 em <- get('EM_completeDAT',asNamespace('GPCCA'));body(em) <- as.call(append(parts,list(extra),after=index));environment(em) <- env
 env$EM_completeDAT <- em;environment(model) <- env
 writeLines(deparse(em),file.path(out,'isolated_initial_ridge_function.R'))
}
messages <- warnings <- calls <- character();stage <- list();fit <- NULL;status <- 'success';error <- NA_character_
set.seed(20261110L);start <- proc.time()[['elapsed']]
console <- capture.output(tryCatch(withCallingHandlers({
 fit <- model(X.list=X_list,d=3L,missVal=FALSE,plot_RMSE=FALSE,verbose=1)
 if(!all(is.finite(fit$Z)) || !all(is.finite(fit$W)))stop('Non-finite fitted Z/W')
},message=function(e){messages <<- c(messages,conditionMessage(e));invokeRestart('muffleMessage')},warning=function(e){warnings <<- c(warnings,conditionMessage(e));invokeRestart('muffleWarning')},error=function(e){
 calls <<- vapply(sys.calls(),function(z)paste(deparse(z),collapse=' '),character(1))
 for(frame in sys.frames()) if(exists('Psi',frame,inherits=FALSE)&&exists('iter',frame,inherits=FALSE)) {
  psi <- get('Psi',frame)
  stage <<- list(iteration=get('iter',frame),Psi_dimensions=dim(psi),Psi_zero_diagonal=sum(Matrix::diag(psi)==0),covariance_initialization=if(get('iter',frame)==0)'empirical covariance, before any M-step ridge' else 'after M-step')
 }
}),error=function(e){status <<- 'failed';error <<- conditionMessage(e)}))
elapsed <- proc.time()[['elapsed']]-start
if(!is.null(fit))saveRDS(fit,file.path(out,paste0(id,'_fit.rds')))
res <- list(case=id,base_case=case,boot=boot,status=status,error=error,seed=20261110L,settings=list(d=3L,missVal=FALSE,lambda=.5,maxiter=100L,tol=1e-4,EarlyStop=TRUE,W_init='PCA',model_changed_for_diagnostic=changed_model),patients=nrow(X),features=ncol(X),block_features=sapply(X_list,nrow),dropped_features=dropped,elapsed_seconds=elapsed,warnings=warnings,messages=messages,error_stage=stage,error_calls=calls,score_dimensions=if(!is.null(fit))dim(fit$Z) else NULL,finite_scores=if(!is.null(fit))all(is.finite(fit$Z)) else NA)
saveRDS(res,file.path(out,paste0(id,'_result.rds')))
jsonlite::write_json(res,file.path(out,paste0(id,'_result.json')),pretty=TRUE,auto_unbox=TRUE,na='null')
writeLines(c(messages,console,warnings,error,calls),file.path(out,paste0(id,'_trace.log')))
cat(id,status,'patients',nrow(X),'features',ncol(X),'seconds',elapsed,'error',error,'\n')
