# Fixed, donor-level paired public-intervention analysis. No environment changes.
suppressPackageStartupMessages({library(edgeR);library(limma);library(jsonlite)})
D <- '/root/five_new_studies_20261020/projects/N02_keloid_genetic_cell_state/mechanisms_20261021/literature_intervention_v1'
args <- commandArgs(TRUE);stopifnot(length(args)==1);OUT<-args[1];dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
options(stringsAsFactors=FALSE);set.seed(20261021)
meta<-read.delim(file.path(D,'all56_GSM_sample_metadata.tsv'),check.names=FALSE)
hg<-read.delim(file.path(D,'all_unique_approved_HGNC_gene_axis.tsv'),check.names=FALSE,colClasses='character')
program<-read.delim(file.path(D,'fixed4_program_members_and_weights.tsv'),check.names=FALSE)
cand<-read.delim(file.path(D,'fixed9_candidates.tsv'),check.names=FALSE,colClasses='character')
probe<-read.delim(gzfile(file.path(D,'GPL570_fixed_all_probe_to_HGNC_Ensembl.tsv.gz')),check.names=FALSE,colClasses='character')
writegz<-function(x,n){con<-gzfile(file.path(OUT,n),'wt');write.table(x,con,sep='\t',row.names=FALSE,quote=FALSE,na='NA');close(con)}
allDE<-list();allCAM<-list();allFRY<-list();allRaw<-list();allPS<-list();allDelta<-list();allSummary<-list();coverage<-list();QC<-list();allCand<-list();fitDiag<-list()
# Each program remains fixed. Missing genes are reported; there is no backfilling.
sets<-split(program$ensembl_gene_id,program$program);sets<-lapply(sets,function(x)unique(x[!is.na(x)&nzchar(x)]))
ci<-function(x){x<-x[is.finite(x)];n<-length(x);if(n<2)return(c(n=n,mean=if(n)mean(x)else NA,lower=NA,upper=NA));se<-sd(x)/sqrt(n);c(n=n,mean=mean(x),lower=mean(x)-qt(.975,n-1)*se,upper=mean(x)+qt(.975,n-1)*se)}
load_counts<-function(ds){
 fn<-if(ds=='GSE232079')'GSE232079_counts.txt.gz'else'GSE232435_counts_ready.txt.gz'
 x<-as.matrix(read.delim(gzfile(file.path(D,'sources',fn)),row.names=1,check.names=FALSE));storage.mode(x)<-'double'
 stopifnot(!anyDuplicated(rownames(x)),all(is.finite(x)),all(x>=0),all(x==round(x)))
 m<-meta[meta$dataset==ds,];stopifnot(setequal(colnames(x),m$count_column));m<-m[match(colnames(x),m$count_column),]
 # Same donor/condition culture repeats are collapsed before primary inference.
 keys<-paste(m$source_strain,m$treatment,sep='::');ord<-unique(keys)
 a<-vapply(ord,function(k)rowSums(x[,keys==k,drop=FALSE]),numeric(nrow(x)));rownames(a)<-rownames(x);colnames(a)<-ord
 am<-m[match(ord,keys),];am$analysis_sample<-ord;am$cultures_collapsed<-as.integer(table(factor(keys,levels=ord)))
 a<-a[rownames(a)%in%hg$ensembl_gene_id,,drop=FALSE]
 list(x=a,m=am,raw_library=colSums(x),raw_gene_axis=nrow(x))
}
load_HC<-function(){
 z<-readLines(gzfile(file.path(D,'sources','GSE7890_series_matrix.txt.gz')));b<-which(z=='!series_matrix_table_begin');e<-which(z=='!series_matrix_table_end');stopifnot(length(b)==1,length(e)==1)
 x<-as.matrix(read.delim(text=paste(z[(b+1):(e-1)],collapse='\n'),row.names=1,check.names=FALSE));stopifnot(all(is.finite(x)),all(x>0))
 # Actual positive normalized export scale is frozen by full-matrix preflight.
 stopifnot(max(x)>1000);x<-log2(x)
 a<-probe[probe$status=='UNIQUE_APPROVED_MAPPING'&probe$probe_id%in%rownames(x),];a<-a[!duplicated(a$probe_id),]
 v<-x[match(a$probe_id,rownames(x)),,drop=FALSE];s<-rowsum(v,a$ensembl_gene_id,reorder=FALSE);n<-table(factor(a$ensembl_gene_id,levels=rownames(s)));s<-s/as.numeric(n)
 m<-meta[meta$dataset=='GSE7890',];m<-m[match(colnames(s),m$count_column),];m$analysis_sample<-m$count_column;m$cultures_collapsed<-3L
 list(x=s,m=m,raw_library=rep(NA_real_,ncol(s)),raw_gene_axis=nrow(x))
}
run_one<-function(ds){
 raw<-if(ds=='GSE7890')load_HC()else load_counts(ds);x<-raw$x;m<-raw$m
 if(ds=='GSE7890')L<-x else{y0<-calcNormFactors(DGEList(x),method='TMM');L<-cpm(y0,log=TRUE,prior.count=2)}
 # Rank within the same assay's complete measured HGNC expression background.
 # This is RNA-program expression, not an individual genetic PRS/scDRS score.
 ranks<-apply(L,2,rank,ties.method='average')/nrow(L);rownames(ranks)<-rownames(L)
 ps<-list()
 for(pr in names(sets)){
  pm<-program[program$program==pr&program$ensembl_gene_id%in%rownames(ranks),];pm<-pm[!duplicated(pm$ensembl_gene_id),];stopifnot(nrow(pm)>0)
  w<-as.numeric(pm$weight);sc<-colSums(ranks[pm$ensembl_gene_id,,drop=FALSE]*w)/sum(w)
  r<-data.frame(dataset=ds,program=pr,analysis_sample=m$analysis_sample,source_strain=m$source_strain,material=m$material,treatment=m$treatment,cultures_collapsed=m$cultures_collapsed,score=as.numeric(sc),assayed_program_genes=nrow(pm),HGNC_background_genes=nrow(L),not_individual_PRS_or_drug_efficacy=TRUE)
  allPS[[paste(ds,pr)]]<<-r;ps[[pr]]<-r
 }
 for(g in cand$gene_id){for(j in seq_len(ncol(x))){ix<-match(g,rownames(L));allRaw[[length(allRaw)+1]]<<-data.frame(dataset=ds,analysis_sample=m$analysis_sample[j],source_strain=m$source_strain[j],material=m$material[j],treatment=m$treatment[j],gene_id=g,approved_symbol=cand$approved_symbol[match(g,cand$gene_id)],normalized_log_expression=if(!is.na(ix))L[ix,j]else NA_real_,assayed=!is.na(ix),no_new_independent_patient=TRUE)}}
 # QC only; preserve all samples and report every material. No PCA-based removal.
 qc<-data.frame(dataset=ds,analysis_sample=m$analysis_sample,source_strain=m$source_strain,material=m$material,treatment=m$treatment,HGNC_input_library=if(ds=='GSE7890')NA_real_ else colSums(x),nonzero_HGNC_genes=colSums(x>0),cultures_collapsed=m$cultures_collapsed)
 pca<-prcomp(t(L),center=TRUE,scale.=FALSE);qc$PC1<-pca$x[,1];qc$PC2<-pca$x[,2];QC[[ds]]<<-qc
 pdf(file.path(OUT,paste0(ds,'_sample_QC.pdf')),width=9,height=5);par(mfrow=c(1,2));boxplot(L,outline=FALSE,las=2,cex.axis=.45,main=paste(ds,'all HGNC log-expression'));plot(qc$PC1,qc$PC2,pch=19,col=as.integer(factor(m$treatment)),xlab='PC1',ylab='PC2',main='No sample exclusions');text(qc$PC1,qc$PC2,labels=paste(m$source_strain,m$treatment),cex=.45,pos=3);dev.off()
 if(ds=='GSE232079'){
  use<-m$source_strain%in%paste0('PK',1:4);mm<-m[use,];mm$source_strain<-factor(mm$source_strain);mm$treatment<-factor(mm$treatment,levels=c('DMSO','corin'));design<-model.matrix(~source_strain+treatment,mm);cm<-matrix(0,ncol(design),1,dimnames=list(colnames(design),'Corin_minus_DMSO_primary4'));cm['treatmentcorin',1]<-1
  specs<-list(list(name=colnames(cm)[1],material='primary_keloid',from='DMSO',to='corin',strains=paste0('PK',1:4)))
 }else if(ds=='GSE232435'){
  use<-rep(TRUE,nrow(m));mm<-m;mm$source_strain<-factor(mm$source_strain);mm$treatment<-factor(mm$treatment,levels=c('Control','TGFbeta','HS276_TGFbeta'));design<-model.matrix(~source_strain+treatment,mm);cm<-matrix(0,ncol(design),3,dimnames=list(colnames(design),c('TGFbeta_minus_Control','HS276_TGFbeta_minus_TGFbeta','HS276_TGFbeta_minus_Control')));cm['treatmentTGFbeta',1]<-1;cm['treatmentHS276_TGFbeta',2]<-1;cm['treatmentTGFbeta',2]<- -1;cm['treatmentHS276_TGFbeta',3]<-1
  specs<-list(list(name=colnames(cm)[1],material='primary_foreskin',from='Control',to='TGFbeta',strains=c('N5','N6','N8')),list(name=colnames(cm)[2],material='primary_foreskin',from='TGFbeta',to='HS276_TGFbeta',strains=c('N5','N6','N8')),list(name=colnames(cm)[3],material='primary_foreskin',from='Control',to='HS276_TGFbeta',strains=c('N5','N6','N8')))
 }else{
  use<-m$pairing_status=='PAIRED_WITHIN_STRAIN';mm<-m[use,];mm$source_strain<-factor(mm$source_strain);mm$HC<-as.integer(mm$treatment=='HC');mm$KHC<-as.integer(mm$material=='keloid')*mm$HC;design<-model.matrix(~source_strain+HC+KHC,mm);cm<-matrix(0,ncol(design),3,dimnames=list(colnames(design),c('HC_minus_NoHC_keloid5','HC_minus_NoHC_normalScar4','HC_response_keloid_minus_normalScar')));cm['HC',1]<-1;cm['KHC',1]<-1;cm['HC',2]<-1;cm['KHC',3]<-1
  specs<-list(list(name=colnames(cm)[1],material='keloid',from='NoHC',to='HC',strains=c('33','50','145','261','124')),list(name=colnames(cm)[2],material='normal_scar',from='NoHC',to='HC',strains=c('116','131','170','21')))
 }
 stopifnot(qr(design)$rank==ncol(design),nrow(design)>ncol(design));write.table(data.frame(analysis_sample=m$analysis_sample[use],design),file.path(OUT,paste0(ds,'_paired_design.tsv')),sep='\t',row.names=FALSE,quote=FALSE)
 write.table(cm,file.path(OUT,paste0(ds,'_all_contrast_coefficients.tsv')),sep='\t',quote=FALSE)
 if(ds=='GSE7890'){
  exp<-x[,use,drop=FALSE];keep<-apply(exp,1,function(z)all(is.finite(z)));vv<-exp[keep,,drop=FALSE]
 }else{
  y<-DGEList(x[,use,drop=FALSE]);keep<-filterByExpr(y,design=design,min.count=10,min.total.count=15);y<-calcNormFactors(y[keep,,keep.lib.sizes=FALSE],method='TMM');vv<-voom(y,design,plot=FALSE)
 }
 fit<-eBayes(contrasts.fit(lmFit(vv,design),cm),robust=TRUE)
 ids<-if(ds=='GSE7890')rownames(vv)else rownames(vv$E);ann<-hg[match(rownames(x),hg$ensembl_gene_id),];stopifnot(!anyNA(ann$ensembl_gene_id))
 fitDiag[[ds]]<<-list(samples=nrow(design),paired_design_rank=ncol(design),residual_df=nrow(design)-ncol(design),input_HGNC_genes=nrow(x),DE_expression_admitted_genes=length(ids),dropped_low_expression_genes=nrow(x)-length(ids),method=if(ds=='GSE7890')'limma on log2 provider normalized expression; no new RMA'else'TMM-voom paired limma moderated t',no_patient_or_sample_exclusions_by_effect=TRUE)
 for(k in seq_len(ncol(cm))){
  nm<-colnames(cm)[k];tt<-topTable(fit,coef=k,number=Inf,sort.by='none');tt$gene_id<-rownames(tt);base<-data.frame(dataset=ds,contrast=nm,gene_id=ann$ensembl_gene_id,HGNC_symbol=ann$symbol,entrez_id=ann$entrez_id,locus_type=ann$locus_type,test_status=ifelse(rownames(x)%in%ids,'TESTED','LOW_EXPRESSION_NOT_TESTED_NA'));tab<-merge(base,tt,by='gene_id',all.x=TRUE,sort=FALSE);writegz(tab,paste0(ds,'_',nm,'_complete_HGNC_DE.tsv.gz'));allDE[[paste(ds,nm)]]<<-data.frame(dataset=ds,contrast=nm,gene_id=tab$gene_id,HGNC_symbol=tab$HGNC_symbol,t=tab$t,logFC=tab$logFC,P.Value=tab$P.Value,adj.P.Val=tab$adj.P.Val,test_status=tab$test_status)
  for(g in cand$gene_id){r<-tab[match(g,tab$gene_id),];allCand[[length(allCand)+1]]<<-data.frame(dataset=ds,contrast=nm,gene_id=g,approved_symbol=cand$approved_symbol[match(g,cand$gene_id)],logFC=r$logFC,t=r$t,P=r$P.Value,whole_HGNC_FDR=r$adj.P.Val,status=if(is.na(r$gene_id))'NOT_ASSAYED_NA'else r$test_status)}
  ind<-lapply(sets,function(g)which(ids%in%g));cam<-camera(vv,ind,design=design,contrast=cm[,k],inter.gene.cor=.01,allow.neg.cor=FALSE,sort=FALSE);cam$program<-rownames(cam);cam$dataset<-ds;cam$contrast<-nm;cam$gene_set_P_not_patient_P<-TRUE;allCAM[[paste(ds,nm)]]<<-cam
  fr<-fry(vv,ind,design=design,contrast=cm[,k],sort="none");fr$program<-rownames(fr);fr$dataset<-ds;fr$contrast<-nm;fr$gene_set_P_not_patient_P<-TRUE;allFRY[[paste(ds,nm)]]<<-fr
  coverage[[paste(ds,nm)]]<<-data.frame(dataset=ds,contrast=nm,program=names(ind),fixed_assayed_members=sapply(sets,function(g)sum(g%in%rownames(x))),DE_tested_members=sapply(ind,length),weighted_MAGMA_only_in_rank_score_not_CAMERA_or_FRY=TRUE)
 }
 for(pr in names(ps))for(s in specs){
  sr<-ps[[pr]];dlt<-sapply(s$strains,function(st){a<-sr$score[sr$source_strain==st&sr$treatment==s$to];b<-sr$score[sr$source_strain==st&sr$treatment==s$from];stopifnot(length(a)==1,length(b)==1);a-b})
  allDelta[[paste(ds,pr,s$name)]]<<-data.frame(dataset=ds,contrast=s$name,program=pr,source_strain=s$strains,material=s$material,delta=as.numeric(dlt),patient_or_source_not_culture_repeat=TRUE,descriptive_no_patient_P=TRUE)
  z<-ci(dlt);allSummary[[paste(ds,pr,s$name)]]<<-data.frame(dataset=ds,contrast=s$name,program=pr,n=z['n'],mean=z['mean'],CI95_lower=z['lower'],CI95_upper=z['upper'],positive=sum(dlt>0),negative=sum(dlt<0),zero=sum(dlt==0),CI_method='t-interval over independent source deltas; smalln descriptive',no_test_or_direction_based_stop=TRUE)
 }
 if(ds=='GSE7890')for(pr in names(ps)){
  a<-allDelta[[paste(ds,pr,specs[[1]]$name)]]$delta;b<-allDelta[[paste(ds,pr,specs[[2]]$name)]]$delta;se<-sqrt(var(a)/length(a)+var(b)/length(b));df<-(var(a)/length(a)+var(b)/length(b))^2/((var(a)/length(a))^2/(length(a)-1)+(var(b)/length(b))^2/(length(b)-1));mean<-mean(a)-mean(b)
  allSummary[[paste(ds,pr,'interaction')]]<<-data.frame(dataset=ds,contrast='HC_response_keloid_minus_normalScar',program=pr,n=length(a)+length(b),mean=mean,CI95_lower=mean-qt(.975,df)*se,CI95_upper=mean+qt(.975,df)*se,positive=NA,negative=NA,zero=NA,CI_method='Welch interval on5KLvs4normalstrain paired deltas; descriptive',no_test_or_direction_based_stop=TRUE)
 }
}
for(ds in c('GSE232079','GSE232435','GSE7890'))run_one(ds)
cam<-do.call(rbind,allCAM);fr<-do.call(rbind,allFRY);stopifnot(nrow(cam)==28,nrow(fr)==28);cam$BH28_joint<-p.adjust(cam$PValue,method='BH');fr$BH28_directional<-p.adjust(fr$PValue,method='BH');fr$BH28_mixed<-p.adjust(fr$PValue.Mixed,method='BH')
for(obj in c('cam','fr')){z<-get(obj);primary<-z$program%in%c('fixed_shared_nonGS_FB01_25','original_MAGMA1000_Z_RNA');p<-z$PValue;z$BH14_primary2<-NA_real_;z$BH14_primary2[primary]<-p.adjust(p[primary],method='BH');z$BH14_ECM_growth_comparator<-NA_real_;z$BH14_ECM_growth_comparator[!primary]<-p.adjust(p[!primary],method='BH');assign(obj,z)}
cc<-do.call(rbind,allCand);stopifnot(nrow(cc)==63);cc$BH63_candidate_family<-NA_real_;ok<-is.finite(cc$P);cc$BH63_candidate_family[ok]<-p.adjust(replace(cc$P,!ok,1),method='BH')[ok]
writegz(do.call(rbind,allDE),'all7_contrasts_full_HGNC_t_logFC_for_fixed_prior_footprints.tsv.gz')
for(k in c('allRaw','allPS','allDelta','allSummary','coverage','QC'))writegz(do.call(rbind,get(k)),paste0(k,'_complete.tsv.gz'))
writegz(cam,'fixed4programs_all7contrasts_CAMERA_BH28.tsv.gz');writegz(fr,'fixed4programs_all7contrasts_FRY_BH28.tsv.gz');writegz(cc,'all9candidates_all7contrasts_complete63.tsv.gz')
write_json(list(status='COMPLETE_PAIRED_INTERVENTION_STATISTICS_ROOT_REVIEW_PENDING',datasets=fitDiag,n_contrasts=7,n_program_tests=28,n_candidate_rows=63,method_scope='Public human fibroblast transcriptional response, not clinical efficacy, target engagement or causal mediation. Corin not NEDD4-selective; HS276 targetsTAK1 notTAB2; HC is not TAC clinical response.',normal_scar_not_healthy=TRUE,program_rank_not_individual_PRS=TRUE,all_non_significant_opposite_and_NA_preserved=TRUE,sessionInfo=capture.output(sessionInfo())),file.path(OUT,'R_STATISTICS_COMPLETE.json'),pretty=TRUE,auto_unbox=TRUE,na='null')
