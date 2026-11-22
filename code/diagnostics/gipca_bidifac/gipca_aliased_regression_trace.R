# Path-only adaptation of the recorded diagnostic; fitting settings are unchanged.
# BENCH_ARCHIVE_DIR: original benchmark directory, including variables/ and scalar result CSVs.
# Diagnostic outputs stay under ignored private_outputs by default.
diag_out <- Sys.getenv('DIAGNOSTIC_OUTPUT_DIR','private_outputs/diagnostics/gipca_bidifac')
for(part in c('private','public')) dir.create(file.path(diag_out,part),recursive=TRUE,showWarnings=FALSE)
root<-Sys.getenv('DIAGNOSTIC_OUTPUT_DIR','private_outputs/diagnostics/gipca_bidifac')
for(f in list.files(Sys.getenv('GIPCA_SOURCE_DIR','external_sources/gipca/Functions'),full.names=TRUE,pattern='\\.R$'))source(f)
v<-readRDS(file.path(Sys.getenv('BENCH_ARCHIVE_DIR'),'variables/analysis_variables.rds'))
# Trace returns only; coefficients/feature names and patient-level values are never exported.
trace(stats::lm,exit=quote({z<-returnValue();if(anyNA(z$coefficients)){cat('ALIASED LM rank',z$rank,'coefficients',length(z$coefficients),'observations',nobs(z),'\n');for(e in sys.frames())if(exists('V_pre.joint',e,inherits=FALSE)&&exists('l',e,inherits=FALSE))cat('iteration',get('l',e),'\n')}}),print=FALSE)
for(b in c(3L,865L)){
set.seed(100+b);ix<-sample(rownames(v$X_stand),nrow(v$X_stand),TRUE);X<-v$X_stand[ix,,drop=FALSE];bl<-split(colnames(X),v$layer_vec[colnames(X)]);dat<-lapply(bl,function(cc)X[,cc,drop=FALSE]);cat('BOOT',b,'min feature obs',min(colSums(!is.na(X))),'all missing rows',sum(rowSums(!is.na(X))==0),'\n')
tryCatch(EPCAJIVEMissbio(dat,3,rep(0L,length(dat)),vapply(dat,ncol,integer(1)),rep('gaussian',length(dat)),.1,500,NULL,0),error=function(e)cat('ERROR',conditionMessage(e),'\n'))
}
untrace(stats::lm)
