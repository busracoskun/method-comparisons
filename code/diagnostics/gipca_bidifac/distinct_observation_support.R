# Path-only adaptation of the recorded diagnostic; fitting settings are unchanged.
# BENCH_ARCHIVE_DIR: original benchmark directory, including variables/ and scalar result CSVs.
# Diagnostic outputs stay under ignored private_outputs by default.
diag_out <- Sys.getenv('DIAGNOSTIC_OUTPUT_DIR','private_outputs/diagnostics/gipca_bidifac')
for(part in c('private','public')) dir.create(file.path(diag_out,part),recursive=TRUE,showWarnings=FALSE)
root<-Sys.getenv('DIAGNOSTIC_OUTPUT_DIR','private_outputs/diagnostics/gipca_bidifac');p<-Sys.getenv('BENCH_ARCHIVE_DIR')
res<-list();features<-list()
for(ds in c('adni','tcga')){
v<-readRDS(file.path(p,'variables',if(ds=='adni')'analysis_variables.rds' else 'tcga_bootstrap_analysis_variables.rds')); X<-v$X_stand;lv<-v$layer_vec[colnames(X)]
for(bl in unique(lv)){
 A<-X[,lv==bl,drop=FALSE];nobs<-colSums(!is.na(A));features[[paste(ds,bl)]]<-data.frame(dataset=ds,modality=bl,features=ncol(A),min_observations=min(nobs),median_observations=median(nobs),max_observations=max(nobs),features_under_10_observations=sum(nobs<10),features_with_partial_missingness=sum(nobs>0 & nobs<nrow(A)))
}
rr<-lapply(1:1000,function(b){set.seed(100+b);ix<-sample(seq_len(nrow(X)),nrow(X),TRUE);ux<-unique(ix);nobs<-colSums(!is.na(X[ux,,drop=FALSE]));data.frame(dataset=ds,boot=b,min_unique_observations=min(nobs),features_below_4_unique=sum(nobs<4),features_exactly_4_unique=sum(nobs==4),features_0_unique=sum(nobs==0))})
res[[ds]]<-do.call(rbind,rr)
write.csv(res[[ds]],file.path(root,'private',paste0(ds,'_unique_support.csv')),row.names=FALSE)
}
write.csv(do.call(rbind,features),file.path(root,'public/feature_support.csv'),row.names=FALSE)
