# One descriptive estimator reuse, not new limma/CAMERA/FRY or patient P.
suppressPackageStartupMessages({library(edgeR);library(jsonlite)})
args<-commandArgs(TRUE);stopifnot(length(args)==1);H<-args[1];OUT<-file.path(H,'results')
orig<-file.path(dirname(dirname(H)),'literature_intervention_v1','run_paired_intervention_analysis_v1.R')
# Evaluate only pinned original input loaders/metadata and CI helper, never run_one.
ee<-new.env(parent=.GlobalEnv);allowed<-c('D','meta','hg','probe','load_counts','load_HC','ci')
for(ex in parse(file=orig)){if(is.call(ex)&&as.character(ex[[1]])=='<-'&&is.symbol(ex[[2]])&&as.character(ex[[2]])%in%allowed)eval(ex,envir=ee)}
stopifnot(all(vapply(allowed,exists,logical(1),envir=ee,inherits=FALSE)))
members<-read.delim(file.path(OUT,'fixed16_program_members.tsv'),check.names=FALSE,colClasses='character')
catalog<-read.delim(file.path(OUT,'all15_native_state_matching_catalog.tsv'),check.names=FALSE,colClasses='character')
panels<-c('ORIGINAL_SHARED_FB01_25',catalog$program)
old<-read.delim(gzfile(file.path(ee$D,'results_v1','allPS_complete.tsv.gz')),check.names=FALSE)
old<-old[old$program=='fixed_shared_nonGS_FB01_25',]
allscore<-list();alldelta<-list();allsum<-list();allcoverage<-list();repro<-list()
writegz<-function(x,n){cc<-gzfile(file.path(OUT,n),'wt');write.table(x,cc,sep='\t',row.names=FALSE,quote=FALSE,na='NA');close(cc)}
for(ds in c('GSE232079','GSE232435','GSE7890')){
 raw<-if(ds=='GSE7890')ee$load_HC()else ee$load_counts(ds);x<-raw$x;m<-raw$m
 L<-if(ds=='GSE7890')x else cpm(calcNormFactors(DGEList(x),method='TMM'),log=TRUE,prior.count=2)
 ranks<-apply(L,2,rank,ties.method='average')/nrow(L);rownames(ranks)<-rownames(L)
 scores<-list()
 for(pr in panels){
  genes<-unique(members$ensembl_gene_id[members$program==pr]);ix<-which(rownames(ranks)%in%genes)
  ok<-if(pr=='ORIGINAL_SHARED_FB01_25')length(ix)>0 else length(genes)==25&&length(ix)==25
  val<-if(ok)colMeans(ranks[ix,,drop=FALSE])else rep(NA_real_,ncol(ranks))
  scores[[pr]]<-val
  allscore[[paste(ds,pr)]]<-data.frame(dataset=ds,program=pr,analysis_sample=m$analysis_sample,source_strain=m$source_strain,material=m$material,treatment=m$treatment,score=val,assayed_members=length(ix),fixed_members=length(genes),whole_assay_HGNC_background=nrow(L),status=if(ok)'DESCRIPTIVE_COMPLETE'else'INSUFFICIENT_FIXED_MEMBERS_NA')
  allcoverage[[paste(ds,pr)]]<-data.frame(dataset=ds,program=pr,fixed_members=length(genes),assayed_members=length(ix),missing_members=paste(setdiff(genes,rownames(ranks)),collapse='|'),status=if(ok)'DESCRIPTIVE_COMPLETE'else'INSUFFICIENT_FIXED_MEMBERS_NA')
 }
 rr<-allscore[[paste(ds,'ORIGINAL_SHARED_FB01_25')]];oo<-old[old$dataset==ds,];oo<-oo[match(rr$analysis_sample,oo$analysis_sample),];err<-max(abs(rr$score-oo$score));stopifnot(all(rr$analysis_sample==oo$analysis_sample),is.finite(err),err<1e-12)
 repro[[ds]]<-data.frame(dataset=ds,max_abs_original25_reproduction_error=err,samples=nrow(rr),no_original_result_overwrite=TRUE)
 specs<-if(ds=='GSE232079')list(list(name='Corin_minus_DMSO_primary4',from='DMSO',to='corin',strains=paste0('PK',1:4),material='primary_keloid'))else if(ds=='GSE232435')list(list(name='TGFbeta_minus_Control',from='Control',to='TGFbeta',strains=c('N5','N6','N8'),material='primary_foreskin'),list(name='HS276_TGFbeta_minus_TGFbeta',from='TGFbeta',to='HS276_TGFbeta',strains=c('N5','N6','N8'),material='primary_foreskin'),list(name='HS276_TGFbeta_minus_Control',from='Control',to='HS276_TGFbeta',strains=c('N5','N6','N8'),material='primary_foreskin'))else list(list(name='HC_minus_NoHC_keloid5',from='NoHC',to='HC',strains=c('33','50','145','261','124'),material='keloid'),list(name='HC_minus_NoHC_normalScar4',from='NoHC',to='HC',strains=c('116','131','170','21'),material='normal_scar'))
 for(pr in panels)for(s in specs){
  vals<-scores[[pr]];dd<-sapply(s$strains,function(st){a<-which(m$source_strain==st&m$treatment==s$to);b<-which(m$source_strain==st&m$treatment==s$from);stopifnot(length(a)==1,length(b)==1);vals[a]-vals[b]})
  ref<-scores[['ORIGINAL_SHARED_FB01_25']];rd<-sapply(s$strains,function(st)ref[m$source_strain==st&m$treatment==s$to]-ref[m$source_strain==st&m$treatment==s$from])
  alldelta[[paste(ds,pr,s$name)]]<-data.frame(dataset=ds,program=pr,contrast=s$name,source_strain=s$strains,material=s$material,delta=dd,delta_minus_original25=dd-rd,not_patient_P_or_clinical_selectivity=TRUE)
  a<-ee$ci(dd);b<-ee$ci(dd-rd)
  allsum[[paste(ds,pr,s$name)]]<-data.frame(dataset=ds,program=pr,contrast=s$name,n=a['n'],mean_delta=a['mean'],CI95_lower=a['lower'],CI95_upper=a['upper'],positive=sum(dd>0,na.rm=TRUE),negative=sum(dd<0,na.rm=TRUE),NA_sources=sum(!is.finite(dd)),mean_delta_minus_original25=b['mean'],relative_CI95_lower=b['lower'],relative_CI95_upper=b['upper'],description='Same original within-assay rank estimator; smalln descriptive paired-source intervals, no inferential P')
 }
 if(ds=='GSE7890')for(pr in panels){
  a<-alldelta[[paste(ds,pr,'HC_minus_NoHC_keloid5')]]$delta;b<-alldelta[[paste(ds,pr,'HC_minus_NoHC_normalScar4')]]$delta
  ra<-alldelta[[paste(ds,pr,'HC_minus_NoHC_keloid5')]]$delta_minus_original25;rb<-alldelta[[paste(ds,pr,'HC_minus_NoHC_normalScar4')]]$delta_minus_original25
  wi<-function(a,b){if(any(!is.finite(c(a,b))))return(c(mean=NA,lower=NA,upper=NA));se<-sqrt(var(a)/length(a)+var(b)/length(b));df<-(var(a)/length(a)+var(b)/length(b))^2/((var(a)/length(a))^2/(length(a)-1)+(var(b)/length(b))^2/(length(b)-1));d<-mean(a)-mean(b);if(se==0)return(c(mean=d,lower=d,upper=d));c(mean=d,lower=d-qt(.975,df)*se,upper=d+qt(.975,df)*se)}
  z<-wi(a,b);w<-wi(ra,rb)
  allsum[[paste(ds,pr,'interaction')]]<-data.frame(dataset=ds,program=pr,contrast='HC_response_keloid_minus_normalScar',n=9,mean_delta=z['mean'],CI95_lower=z['lower'],CI95_upper=z['upper'],positive=NA,negative=NA,NA_sources=sum(!is.finite(c(a,b))),mean_delta_minus_original25=w['mean'],relative_CI95_lower=w['lower'],relative_CI95_upper=w['upper'],description='Welch descriptive interval for5KL-minus4ordinary-scar source deltas; no P or clinical response claim')
 }
}
S<-do.call(rbind,allscore);D<-do.call(rbind,alldelta);M<-do.call(rbind,allsum)
stopifnot(nrow(M)==16*7,nrow(D)==16*22)
writegz(S,'all16_program_all_assay_sample_rank_scores.tsv.gz');writegz(D,'all16_program_all22_paired_source_deltas.tsv.gz');writegz(M,'all16_program_all7_contrast_descriptive_intervals.tsv.gz');writegz(do.call(rbind,allcoverage),'all16_program_assay_coverage.tsv.gz')
write.table(do.call(rbind,repro),file.path(OUT,'original25_estimator_exact_reproduction.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write_json(list(status='COMPLETE_DESCRIPTIVE_ALL15_NATIVE_STATE_COMPARATORS_NO_NEW_FITS_OR_P',programs=16,contrasts=7,paired_source_rows=nrow(D),summary_rows=nrow(M),sessionInfo=capture.output(sessionInfo())),file.path(OUT,'R_DESCRIPTIVE_COMPLETE.json'),auto_unbox=TRUE,pretty=TRUE)
