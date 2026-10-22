suppressPackageStartupMessages(library(jsonlite))
Sys.setenv(OMP_NUM_THREADS=2,OPENBLAS_NUM_THREADS=2,MKL_NUM_THREADS=2)
p<-'/root/five_new_studies_20261020/projects/N02_keloid_genetic_cell_state';args<-commandArgs(trailingOnly=TRUE);fit<-('--fit-root-approved'%in%args);contract<-file.path(p,'reviews/coloc_contract_v1');inp<-file.path(p,'results/coloc_harmonized_inputs_v1');o<-file.path(p,if(fit)'results/coloc_complete_family_ABF_v1' else 'results/coloc_input_contract_checks_v1');dir.create(o,recursive=TRUE,showWarnings=FALSE)
sha<-function(f){z<-system2('sha256sum',shQuote(f),stdout=TRUE);sub(' .*','',z[1])}
script<-file.path(p,'scripts/coloc_complete_family_driver_v1.R');manifest<-file.path(inp,'complete_gene_context_region_family_harmonized_manifest.tsv');stopifnot(file.exists(file.path(inp,'COMPLETE.json')));info<-fromJSON(file.path(inp,'COMPLETE.json'));stopifnot(info$status=='COMPLETE_TECHNICAL_HARMONIZATION_NO_ASSOCIATION_FIT',!info$association_fitted)
backend<-file.path(contract,'official_source');expected<-c('R_check.R'='26ff6491f8d722c0fd916b4acbb7854607848db1781f025f5625f5a8f08b644f','R_claudia.R'='428f68c70d9351668c2eb61b0c870a31877dae731ff34b83d7c3096edd57b2cc');for(n in names(expected))stopifnot(sha(file.path(backend,n))==expected[[n]])
if(file.exists(file.path(o,'COMPLETE.json')))stop('alreadycomplete,no duplicate')
if(fit){approval<-fromJSON(file.path(contract,'ROOT_APPROVAL_TO_FIT.json'));stopifnot(isTRUE(approval$root_reviewed),approval$driver_SHA256==sha(script),approval$manifest_SHA256==sha(manifest))}
guard<-function(){
 disk<-system2('df',c('-PB1',shQuote('/root'),shQuote('/root/autodl-tmp')),stdout=TRUE)
 rows<-strsplit(trimws(disk[-1]),'[[:space:]]+');free<-vapply(rows,function(x)as.numeric(x[4]),numeric(1));stopifnot(length(free)==2,all(is.finite(free)),all(free>=5*1024^3))
 rss<-grep('^VmRSS:',readLines('/proc/self/status'),value=TRUE);rss<-as.numeric(strsplit(trimws(rss),'[[:space:]]+')[[1]][2])*1024;stopifnot(is.finite(rss),rss<=8*1024^3)
 invisible(TRUE)
}
guard()
core<-new.env(parent=globalenv());sys.source(file.path(backend,'R_check.R'),envir=core);sys.source(file.path(backend,'R_claudia.R'),envir=core)
stopifnot(is.function(core$check_dataset),is.function(core$coloc.abf));fam<-read.delim(manifest,check.names=FALSE,stringsAsFactors=FALSE);stopifnot(nrow(fam)==2391)
prior<-data.frame(prior_id=c('sensitivity_low','primary','sensitivity_high'),p1=1e-4,p2=1e-4,p12=c(1e-6,1e-5,1e-4));write.table(prior,file.path(o,'frozen_priors.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
state<-list(status=if(fit)'RUNNING_ROOT_APPROVED_ABF' else 'RUNNING_INPUT_CONTRACT_ONLY_NO_ABF',PID=Sys.getpid(),driver_SHA256=sha(script),manifest_SHA256=sha(manifest),backend_commit='27776b1bb14e19d1542798ac7477fa378fc2fc0d',backend_DESCRIPTION_version='6.0.3',official_source_unmodified=TRUE,model_fitted=fit,shared_family_rows=nrow(fam),completed_units=0L);write_json(state,file.path(o,'run_status_v1.json'),pretty=TRUE,auto_unbox=TRUE)
checks<-list();outputs<-list();warnlog<-list();ncheck<-0L;nout<-0L
if(fit){snps<-gzfile(file.path(o,'all_primary_per_SNP_ABF_results.tsv.gz'),'wt');writeLines('context\tregion_id\tgene_id\tvariant_id_GRCh38\tlABF_GWAS\tlABF_QTL\tSNP_PP_H4_conditional_on_H4',snps)}
paths<-unique(fam$path[fam$path!='NA'&!is.na(fam$path)]);loaded<-0L
for(f in paths){
 guard()
 subfam<-fam[!is.na(fam$path)&fam$path==f,,drop=FALSE];stopifnot(length(unique(subfam$input_SHA256))==1,sha(f)==subfam$input_SHA256[1]);a<-read.delim(gzfile(f),check.names=FALSE,stringsAsFactors=FALSE);groups<-split(a,a$gene_id)
 for(j in seq_len(nrow(subfam))){
  z<-subfam[j,];d<-groups[[z$gene_id]];n<-if(is.null(d))0L else nrow(d);stopifnot(n==z$common_allele_unique_SNPs);rec<-list(context=z$context,region_id=z$region_id,gene_id=z$gene_id,n_shared=n,nominal_original_rows=z$nominal_original_rows,input_SHA256=z$input_SHA256,model_status='INPUT_NA_NO_SHARED_SNP',error=NA_character_,warnings=NA_character_)
  warnings<-character();err<-NA_character_;ok<-FALSE
  if(n>0){
   stopifnot(!anyDuplicated(d$variant_id_GRCh38),all(is.finite(d$GWAS_beta_raw)),all(is.finite(d$GWAS_SE)),all(d$GWAS_SE>0),all(is.finite(d$QTL_beta_ALT)),all(is.finite(d$QTL_SE)),all(d$QTL_SE>0),all(is.finite(d$QTL_MAF)),all(d$QTL_MAF>0&d$QTL_MAF<1),length(unique(d$QTL_N_context))==1)
   D1<-list(beta=d$GWAS_beta_raw,varbeta=d$GWAS_SE^2,snp=d$variant_id_GRCh38,position=d$position_GRCh38_1based,type='cc',N=1282582,s=4086/1282582)
   D2<-list(beta=d$QTL_beta_ALT,varbeta=d$QTL_SE^2,snp=d$variant_id_GRCh38,position=d$position_GRCh38_1based,type='quant',N=d$QTL_N_context[1],MAF=d$QTL_MAF)
   ok<-tryCatch(withCallingHandlers({core$check_dataset(D1);core$check_dataset(D2);TRUE},warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e){err<<-conditionMessage(e);FALSE})
   rec$model_status<-if(ok)'CONTRACT_PASS_NO_ASSOCIATION' else 'CONTRACT_FAILED_INPUT_NA';rec$error<-err;rec$warnings<-paste(unique(warnings),collapse=' | ');rec$QTL_N<-d$QTL_N_context[1];rec$sdY_mode<-'official_coloc_estimate_from_MAF_N';rec$n_direction_unresolved<-sum(!is.finite(d$GWAS_beta_ALT_if_direction_resolved));rec$position_min<-min(d$position_GRCh38_1based);rec$position_max<-max(d$position_GRCh38_1based)
   if(fit&&ok){
    for(k in seq_len(nrow(prior))){
     pr<-prior[k,];fw<-character();fe<-NA_character_;res<-tryCatch(withCallingHandlers({invisible(capture.output(rr<-core$coloc.abf(D1,D2,p1=pr$p1,p2=pr$p2,p12=pr$p12)));rr},warning=function(w){fw<<-c(fw,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(e){fe<<-conditionMessage(e);NULL});out<-list(context=z$context,region_id=z$region_id,gene_id=z$gene_id,prior_id=pr$prior_id,n_shared=n,p1_requested=pr$p1,p2_requested=pr$p2,p12_requested=pr$p12,fit_status=if(is.null(res))'MODEL_ERROR_NA_RETAINED' else 'COMPLETE_SINGLE_CAUSAL_ABF',error=fe,warnings=paste(unique(fw),collapse=' | '),input_SHA256=z$input_SHA256)
     for(h in 0:4)out[[paste0('PP_H',h)]]<-if(is.null(res))NA_real_ else unname(res$summary[paste0('PP.H',h,'.abf')])
     out$p1_effective<-if(is.null(res))NA_real_ else res$priors['p1'];out$p2_effective<-if(is.null(res))NA_real_ else res$priors['p2'];out$p12_effective<-if(is.null(res))NA_real_ else res$priors['p12'];nout<-nout+1L;outputs[[nout]]<-out
     if(!is.null(res)&&pr$prior_id=='primary'){
      r<-res$results;stopifnot(nrow(r)==n);table<-data.frame(context=z$context,region_id=z$region_id,gene_id=z$gene_id,variant_id_GRCh38=r$snp,lABF_GWAS=r$lABF.df1,lABF_QTL=r$lABF.df2,SNP_PP_H4_conditional_on_H4=r$SNP.PP.H4);write.table(table,snps,sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
     }
    }
   }
  }
  ncheck<-ncheck+1L;checks[[ncheck]]<-rec
 }
 loaded<-loaded+1L;state$completed_units<-ncheck;state$loaded_region_files<-loaded;write_json(state,file.path(o,'run_status_v1.json'),pretty=TRUE,auto_unbox=TRUE)
}
for(i in which(fam$path=='NA'|is.na(fam$path))){z<-fam[i,];ncheck<-ncheck+1L;checks[[ncheck]]<-list(context=z$context,region_id=z$region_id,gene_id=z$gene_id,n_shared=0L,nominal_original_rows=0L,input_SHA256='NA',model_status=z$input_status,error=NA_character_,warnings=NA_character_)}
bind<-function(rows){namesall<-unique(unlist(lapply(rows,names)));do.call(rbind,lapply(rows,function(x){x[setdiff(namesall,names(x))]<-NA;for(nm in names(x))if(is.character(x[[nm]]))x[[nm]]<-gsub('[\r\n\t]+',' ',x[[nm]]);as.data.frame(x[namesall],stringsAsFactors=FALSE)}))};tab<-bind(checks);write.table(tab,file.path(o,'complete_all_units_input_contract_checks.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
if(fit){
 close(snps);base<-fam[rep(seq_len(nrow(fam)),each=nrow(prior)),c('context','region_id','gene_id'),drop=FALSE];base$prior_id<-rep(prior$prior_id,nrow(fam));joined<-if(length(outputs))merge(base,bind(outputs),by=c('context','region_id','gene_id','prior_id'),all.x=TRUE,sort=FALSE) else base
 key<-function(z)paste(z$context,z$region_id,z$gene_id,sep='|');ck<-match(key(joined),key(tab));if(!('fit_status'%in%names(joined)))joined$fit_status<-NA_character_;missing<-is.na(joined$fit_status);joined$fit_status[missing]<-paste0('NOT_FITTED_NA_',tab$model_status[ck[missing]])
 for(h in 0:4)if(!(paste0('PP_H',h)%in%names(joined)))joined[[paste0('PP_H',h)]]<-NA_real_
 write.table(joined,file.path(o,'all_units_all_three_priors_PP0_to_PP4.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='NA');stopifnot(nrow(joined)==nrow(fam)*nrow(prior))
}

state$status<-if(fit)'COMPLETE_ROOT_APPROVED_SINGLE_CAUSAL_ABF_ALL_FAMILY_RETAINED' else 'COMPLETE_INPUT_CONTRACT_ONLY_NO_ABF';state$model_fitted<-fit;state$family_rows_checked<-nrow(tab);state$input_status_counts<-as.list(table(tab$model_status));state$posterior_rows<-nout;state$CC_sample_size_role<-'metadata_only_when_beta_varbeta_supplied_no_P_to_beta_reconstruction';state$QTL_sdY_mode<-'official_estimate_from_MAF_N_not_assumed_1';state$LD_model<-'none_ABF_singlecausal_no_SuSiE';state$end_UTC<-format(Sys.time(),tz='UTC',usetz=TRUE);write_json(state,file.path(o,'run_status_v1.json'),pretty=TRUE,auto_unbox=TRUE);write_json(state,file.path(o,'COMPLETE.json'),pretty=TRUE,auto_unbox=TRUE);writeLines(capture.output(sessionInfo()),file.path(o,'sessionInfo.txt'))
