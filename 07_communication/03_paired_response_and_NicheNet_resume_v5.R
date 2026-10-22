#!/usr/bin/env Rscript
# Engineering-only eight-worker single-ligand parallel resume; unchanged official single=TRUE metrics and full families.
if(Sys.info()[['sysname']]!='Linux'||Sys.getenv('FIVE_STUDIES_COMPUTE')!='autodl')stop('AutoDL only')
suppressPackageStartupMessages({library(edgeR);library(limma);library(jsonlite);library(digest)})
root<-'/root/five_new_studies_20261020/projects/N02_keloid_genetic_cell_state';old<-'/root/autodl-tmp/Five_Studies_Execution_20261020';here<-file.path(root,'communication_mechanism_20261021_v1')
sha<-function(p)digest(file=p,algo='sha256');f<-fromJSON(file.path(here,'FREEZE.json'),simplifyVector=FALSE);a<-fromJSON(file.path(here,'ROOT_RUN_APPROVAL.json'),simplifyVector=FALSE)
stopifnot(a$status=='APPROVED','activity'%in%unlist(a$stages),a$freeze_sha256==sha(file.path(here,'FREEZE.json')))
for(n in names(a$source_sha256))stopifnot(identical(a$source_sha256[[n]],sha(file.path(here,'scripts',n))))
prep<-file.path(here,'profiles_v1');stopifnot(file.exists(file.path(prep,'COMPLETE.json')));pc<-fromJSON(file.path(prep,'COMPLETE.json'),simplifyVector=FALSE);stopifnot(pc$status=='COMPLETE_ALL_PATIENT_PROFILES_NO_S1_EFFECTS',pc$freeze_sha256==sha(file.path(here,'FREEZE.json')),pc$source_sha256==sha(file.path(here,'scripts/02_aggregate_profiles.py')));for(study_entry in pc$studies){for(rr in study_entry$files){stopifnot(file.exists(rr$path),unname(file.info(rr$path)$size)==rr$bytes,sha(rr$path)==rr$sha256)}};
resume_source<-file.path(here,'gpu_switch_20261021_v1','03_paired_response_and_NicheNet_resume_v5.R')
original_source<-file.path(here,'scripts','03_paired_response_and_NicheNet_v3.R')
original_source_sha<-'4297aca11616a9b6b7a473c95ec0f817060993bb0bb5ee4fedc67149b7ed4494'
grant_path<-file.path(here,'gpu_switch_20261021_v1','ROOT_GPU_PARALLEL_RESUME_APPROVAL_v1.json')
stopifnot(file.exists(grant_path));g<-fromJSON(grant_path,simplifyVector=FALSE)
stopifnot(g$status=='APPROVED_FIXED8_PARALLEL_RESUME',g$source_sha256==sha(resume_source),g$original_source_sha256==original_source_sha,sha(original_source)==original_source_sha,g$freeze_sha256==sha(file.path(here,'FREEZE.json')),g$profile_complete_sha256==sha(file.path(prep,'COMPLETE.json')))
stopifnot(is.character(g$equivalence_fixture_path),file.exists(g$equivalence_fixture_path),sha(g$equivalence_fixture_path)==g$equivalence_fixture_sha256)
eq<-fromJSON(g$equivalence_fixture_path,simplifyVector=FALSE)
stopifnot(eq$status=='PASS_REAL_FIXED8_SINGLE_LIGAND_EQUIVALENCE',eq$parallel_source_sha256==sha(resume_source),eq$fixture_source_sha256==g$equivalence_fixture_source_sha256,eq$workspace_sha256==g$workspace_sha256,eq$freeze_sha256==g$freeze_sha256,eq$ligand_target_sha256==g$ligand_target_sha256,isTRUE(eq$all_five_fields_exact_equal),eq$mc_cores==8L,identical(eq$single,TRUE),identical(eq$mc_set_seed,FALSE))
stopifnot(length(g$workspace_path)==1L,is.character(g$workspace_path),file.exists(g$workspace_path),normalizePath(g$workspace_path)==g$workspace_path,unname(file.info(g$workspace_path)$size)==g$workspace_bytes,sha(g$workspace_path)==g$workspace_sha256)
out<-file.path(here,'response_activity_GPU_v1');stopifnot(g$output_path==out)
if(file.exists(file.path(out,'COMPLETE.json')))stop('Already complete; no duplicate')
if(!dir.exists(out))dir.create(out)
atom_json<-function(x,p){tmp<-paste0(p,'.incoming');stopifnot(!file.exists(tmp));write_json(x,tmp,pretty=TRUE,auto_unbox=TRUE,na='null');stopifnot(file.rename(tmp,p))}
binding<-list(schema='N02_GPU_FIXED8_PARALLEL_RESUME_v1',source_sha256=sha(resume_source),original_source_sha256=original_source_sha,freeze_sha256=sha(file.path(here,'FREEZE.json')),profile_complete_sha256=sha(file.path(prep,'COMPLETE.json')),workspace_path=g$workspace_path,workspace_bytes=g$workspace_bytes,workspace_sha256=g$workspace_sha256)
binding_path<-file.path(out,'RESUME_EXECUTION_BINDING.json')
if(file.exists(binding_path)){oldbind<-fromJSON(binding_path,simplifyVector=FALSE);stopifnot(identical(oldbind,binding),isTRUE(g$allow_resume_existing_output))}else{stopifnot(length(list.files(out,all.files=TRUE,no..=TRUE))==0L);atom_json(binding,binding_path)}
cache_dir<-file.path(out,'score_cache');if(!dir.exists(cache_dir))dir.create(cache_dir)
ws<-new.env(parent=emptyenv());loaded_names<-load(g$workspace_path,envir=ws)
required_workspace<-c('f','a','pc','LT','lrraw','allact','allpatient','allcov','allsets','alllr')
stopifnot(all(required_workspace%in%loaded_names),identical(ws$f,f),identical(ws$a,a),identical(ws$pc,pc))
legacy<-list(allact=ws$allact,allpatient=ws$allpatient,allcov=ws$allcov,allsets=ws$allsets,alllr=ws$alllr)
for(n in names(legacy)){stopifnot(is.list(legacy[[n]]));if(length(legacy[[n]]))stopifnot(!is.null(names(legacy[[n]])),!anyNA(names(legacy[[n]])),!anyDuplicated(names(legacy[[n]])))}
legacy_used<-list(allact=character(),allpatient=character());cache_log<-list();score_call_count<-0L

options(error=function(){write_json(list(status='FAILED_PRESERVED',message=geterrmessage()),file.path(out,'FAILED.json'),pretty=TRUE,auto_unbox=TRUE);quit(save='no',status=1)})
write<-function(x,p)write.table(x,p,sep='\t',quote=FALSE,row.names=FALSE,na='NA')
wgz<-function(x,p){tmp<-paste0(p,'.incoming');stopifnot(!file.exists(tmp));con<-gzfile(tmp,'wt');tryCatch(write(x,con),finally=close(con));stopifnot(file.rename(tmp,p))}
rds<-file.path(old,'metadata/mechanism_primary_sources/nichenet_priors');ltpath<-file.path(rds,'ligand_target_matrix_nsga2r_final.rds');lrpath<-file.path(rds,'lr_network_human_21122021.rds')
for(x in f$inputs)if(x$path%in%c(ltpath,lrpath))stopifnot(identical(x$sha256,sha(x$path)))
stopifnot(unname(tools::md5sum(ltpath))=='b09606b04b2d4490418d9028c0e58b9f')
manifest<-fromJSON(file.path(old,'metadata/nichenetr_optional_source/manifest.json'))
for(i in seq_len(nrow(manifest$files))){z<-manifest$files[i,];stopifnot(sha(file.path(old,'metadata/nichenetr_optional_source',z$path))==z$sha256)}
stopifnot(g$ligand_target_sha256==sha(ltpath),g$LR_prior_sha256==sha(lrpath));nn<-source(file.path(old,'shared/load_nichenet_activity.R'))$value;LT<-readRDS(ltpath);lrraw<-readRDS(lrpath);stopifnot(identical(ws$LT,LT),identical(ws$lrraw,lrraw));rm(list=loaded_names,envir=ws);rm(ws);gc();stopifnot(is.matrix(LT),!anyDuplicated(rownames(LT)),!anyDuplicated(colnames(LT)),all(c('from','to')%in%names(lrraw)))
write(lrraw,file.path(out,'official_complete_LR_prior_unfiltered.tsv'));lr<-unique(lrraw[,c('from','to')]);allact<-list();allcov<-list();alllr<-list();allsets<-list();alltargets<-list();allpatient<-list();k<-0L;priorAllTargets<-vapply(seq_len(ncol(LT)),function(j)sum(LT[,j]>0),numeric(1));priorBGCache<-new.env(parent=emptyenv())
expected_keys<-list(allact=character(),allpatient=character())
for(ss in names(f$cohorts)){
 pp<-unlist(f$cohorts[[ss]]$patients);scopes<-c('main',paste0('LODO_',pp))
 for(sc in scopes)for(vv in names(f$views_excluded_symbols))for(sg in c('Up','Down')){
  kk<-paste(ss,sc,vv,sg,sep='__');expected_keys$allact<-c(expected_keys$allact,kk)
  if(sc=='main')expected_keys$allpatient<-c(expected_keys$allpatient,paste(kk,pp))
 }
}
for(kind in c('allact','allpatient'))stopifnot(all(names(legacy[[kind]])%in%expected_keys[[kind]]))
predict_official_fixed8<-function(nn,LT,target,bg,use){
 stopifnot(Sys.info()[['sysname']]=='Linux',length(use)>0L,!anyDuplicated(use),all(use%in%colnames(LT)))
 for(n in c('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS'))stopifnot(Sys.getenv(n)=='1')
 cpu<-strsplit(trimws(readLines('/sys/fs/cgroup/cpu.max',warn=FALSE)),' +')[[1]];stopifnot(length(cpu)==2L,cpu[1]!='max',as.numeric(cpu[1])/as.numeric(cpu[2])>=8)
 mlim<-readLines('/sys/fs/cgroup/memory.max',warn=FALSE);mcur<-as.numeric(readLines('/sys/fs/cgroup/memory.current',warn=FALSE));stopifnot(mlim!='max',as.numeric(mlim)-mcur>=24*2^30)

 # Contiguous deterministic blocks in original use order. At n<8, empty blocks do not call NicheNet.
 cuts<-floor(seq(0,length(use),length.out=9L))
 blocks<-lapply(seq_len(8L),function(i)if(cuts[i+1L]>cuts[i])use[seq.int(cuts[i]+1L,cuts[i+1L])]else character())
 stopifnot(identical(unlist(blocks,use.names=FALSE),use))
 parts<-parallel::mclapply(blocks,function(block){
  Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',NUMBA_NUM_THREADS='1')
  if(!length(block))return(NULL)
  # The frozen original function remains unchanged; single-ligand evaluation, not multi-ligand learning.
  ac<-nn$predict_ligand_activities(geneset=target,background_expressed_genes=bg,ligand_target_matrix=LT,potential_ligands=block,single=TRUE)
  stopifnot(is.data.frame(ac),identical(names(ac),c('test_ligand','auroc','aupr','aupr_corrected','pearson')),nrow(ac)==length(block),!anyDuplicated(ac$test_ligand),setequal(ac$test_ligand,block))
  ac<-as.data.frame(ac,stringsAsFactors=FALSE);ac<-ac[match(block,ac$test_ligand),,drop=FALSE];rownames(ac)<-NULL;ac
 },mc.cores=8L,mc.set.seed=FALSE,mc.preschedule=TRUE)
 if(any(vapply(parts,function(x)inherits(x,'try-error'),logical(1))))stop('Official single-ligand worker failed; no key checkpoint published')
 parts<-parts[!vapply(parts,is.null,logical(1))];stopifnot(length(parts)>0L)
 ac<-do.call(rbind,parts);stopifnot(nrow(ac)==length(use),!anyDuplicated(ac$test_ligand),setequal(ac$test_ligand,use));ac<-ac[match(use,ac$test_ligand),,drop=FALSE];rownames(ac)<-NULL;stopifnot(identical(ac$test_ligand,use));ac
}
score<-function(bg,target,eligible,study,view,sign,scope,cache_kind,cache_key){
 stopifnot(cache_kind%in%c('allact','allpatient'),length(cache_key)==1L,cache_key%in%expected_keys[[cache_kind]],!anyDuplicated(bg),!anyDuplicated(target),all(target%in%bg))
 universe<-colnames(LT);use<-intersect(universe,eligible);can_score<-length(target)>=10L&&length(bg)>=50L&&length(target)<length(bg)&&length(use)>0L
 bgkey<-digest(bg,algo='sha256');if(!exists(bgkey,envir=priorBGCache,inherits=FALSE)){bi<-match(bg,rownames(LT));assign(bgkey,if(length(bg))vapply(seq_len(ncol(LT)),function(j)sum(LT[bi,j]>0),numeric(1))else rep(0,ncol(LT)),envir=priorBGCache)}
 expected_bg_prior<-get(bgkey,envir=priorBGCache,inherits=FALSE)
 context<-list(cache_kind=cache_kind,cache_key=cache_key,study=study,view=view,response_sign=sign,scope=scope,bg=bg,target=target,eligible=eligible,universe=universe,freeze_sha256=binding$freeze_sha256,profile_complete_sha256=binding$profile_complete_sha256,original_source_sha256=original_source_sha,ligand_target_sha256=g$ligand_target_sha256,LR_prior_sha256=g$LR_prior_sha256)
 validate_table<-function(ans){
  want_names<-c('test_ligand','auroc','aupr','aupr_corrected','pearson','status','study','view','response_sign','scope','background_n','target_n','prior_all_targets_n','prior_expressed_targets_n')
  stopifnot(is.data.frame(ans),identical(names(ans),want_names),nrow(ans)==length(universe),identical(ans$test_ligand,universe),!anyDuplicated(ans$test_ligand))
  for(n in c('auroc','aupr','aupr_corrected','pearson'))stopifnot(is.numeric(ans[[n]]))
  stopifnot(!anyNA(ans$status),all(ans$study==study),all(ans$view==view),all(ans$response_sign==sign),all(ans$scope==scope),all(ans$background_n==length(bg)),all(ans$target_n==length(target)),identical(as.numeric(ans$prior_all_targets_n),as.numeric(priorAllTargets)),identical(as.numeric(ans$prior_expressed_targets_n),as.numeric(expected_bg_prior)))
  unsupported<-!universe%in%eligible
  stopifnot(all(ans$status[unsupported]=='NO_SAME_PATIENT_LR_SUPPORT'))
  for(n in c('auroc','aupr','aupr_corrected','pearson'))stopifnot(all(is.na(ans[[n]][unsupported])))
  if(can_score){ii<-match(use,universe);stopifnot(identical(ans$status[ii],ifelse(is.finite(ans$aupr_corrected[ii]),'SCORED','OFFICIAL_METRIC_NONFINITE_PRESERVED')))}else{
   stopifnot(all(ans$status[!unsupported]=='RESPONSE_TOO_SMALL_OR_CONSTANT'))
   for(n in c('auroc','aupr','aupr_corrected','pearson'))stopifnot(all(is.na(ans[[n]])))
  }
  invisible(TRUE)
 }
 cp<-file.path(cache_dir,paste0(cache_kind,'_',digest(cache_key,algo='sha256',serialize=FALSE),'.rds'))
 origin<-'NEW_OFFICIAL_SCORE'
 if(file.exists(cp)){
  cc<-readRDS(cp);stopifnot(identical(cc$schema,'N02_SMALL_SCORE_CHECKPOINT_v1'),identical(cc$context,context),identical(cc$table_sha256,digest(cc$table,algo='sha256')));ans<-cc$table;validate_table(ans);origin<-'REUSED_SMALL_CHECKPOINT'
 }else if(cache_key%in%names(legacy[[cache_kind]])){
  ans<-legacy[[cache_kind]][[cache_key]];validate_table(ans);origin<-'RECOVERED_ORIGINAL_V3_WORKSPACE'
  if(cache_kind=='allact'&&length(target)){
   stopifnot(cache_key%in%names(legacy$allsets));mm<-legacy$allsets[[cache_key]]
   stopifnot(identical(mm$gene_symbol,target),identical(as.integer(mm$position),seq_along(target)),all(mm$study==study),all(mm$scope==scope),all(mm$view==view),all(mm$sign==sign))
  }
  legacy_used[[cache_kind]]<<-c(legacy_used[[cache_kind]],cache_key)
 }else{
  ans<-data.frame(test_ligand=universe,auroc=NA_real_,aupr=NA_real_,aupr_corrected=NA_real_,pearson=NA_real_,status=ifelse(universe%in%eligible,'RESPONSE_TOO_SMALL_OR_CONSTANT','NO_SAME_PATIENT_LR_SUPPORT'),stringsAsFactors=FALSE)
  if(can_score){
   ac<-predict_official_fixed8(nn,LT,target,bg,use)
   stopifnot(!anyDuplicated(ac$test_ligand),setequal(ac$test_ligand,use));ix<-match(ac$test_ligand,ans$test_ligand)
   for(n in intersect(names(ac),names(ans)))ans[ix,n]<-ac[[n]]
   ans$status[ix]<-ifelse(is.finite(ans$aupr_corrected[ix]),'SCORED','OFFICIAL_METRIC_NONFINITE_PRESERVED')
  }
  ans$study<-study;ans$view<-view;ans$response_sign<-sign;ans$scope<-scope;ans$background_n<-length(bg);ans$target_n<-length(target);ans$prior_all_targets_n<-priorAllTargets;ans$prior_expressed_targets_n<-expected_bg_prior
  validate_table(ans)
 }
 if(!file.exists(cp)){
  tmp<-paste0(cp,'.incoming');stopifnot(!file.exists(tmp))
  saveRDS(list(schema='N02_SMALL_SCORE_CHECKPOINT_v1',context=context,table=ans,table_sha256=digest(ans,algo='sha256'),origin=origin,legacy_context_boundary=if(origin=='RECOVERED_ORIGINAL_V3_WORKSPACE')'Original v3 stored counts and main response members, not per-call full bg/eligible/patient-target hashes; full contexts deterministically reconstructed from unchanged inputs and exact original source.'else NA_character_),tmp,compress='gzip')
  chk<-readRDS(tmp);stopifnot(identical(chk$context,context),identical(chk$table,ans));stopifnot(file.rename(tmp,cp))
 }
 score_call_count<<-score_call_count+1L
 cache_log[[score_call_count]]<<-data.frame(kind=cache_kind,key=cache_key,origin=origin,checkpoint_path=cp,checkpoint_bytes=unname(file.info(cp)$size),checkpoint_sha256=sha(cp),target_n=length(target),background_n=length(bg),eligible_n=length(use),stringsAsFactors=FALSE)
 atom_json(list(status='RUNNING_SMALL_CHECKPOINTS',completed_score_calls=score_call_count,latest_kind=cache_kind,latest_key=cache_key,latest_origin=origin,source_sha256=sha(resume_source),freeze_sha256=binding$freeze_sha256,UTC=format(Sys.time(),tz='UTC',usetz=TRUE)),file.path(out,'PROGRESS.json'))
 cat(format(Sys.time(),tz='UTC',usetz=TRUE),cache_kind,cache_key,origin,'complete\n');flush.console()
 ans
}
compare_original_effect<-function(whole,path){
 if(!file.exists(path))return(invisible(NULL))
 original<-read.delim(gzfile(path),check.names=FALSE,stringsAsFactors=FALSE)
 stopifnot(identical(names(original),names(whole)),setequal(original$gene_symbol,whole$gene_symbol),!anyDuplicated(original$gene_symbol));original<-original[match(whole$gene_symbol,original$gene_symbol),,drop=FALSE]
 for(n in names(whole)){
  x<-whole[[n]];y<-original[[n]]
  if(is.numeric(x)){stopifnot(identical(is.na(x),is.na(y)),identical(is.finite(x),is.finite(y)));bad<-!is.finite(x)&!is.na(x);stopifnot(identical(x[bad],y[bad]));ok<-is.finite(x);stopifnot(all(abs(x[ok]-y[ok])<=1e-10*pmax(1,abs(x[ok]))))}else stopifnot(identical(as.character(x),as.character(y)))
 }
 invisible(TRUE)
}

for(st in names(f$cohorts)){
 d<-file.path(prep,st);dest<-file.path(out,st);if(!dir.exists(dest))dir.create(dest);M<-read.delim(file.path(d,'profiles.tsv'),check.names=FALSE,stringsAsFactors=FALSE);known_raw<-tolower(as.character(M$known_patient));stopifnot(!anyNA(known_raw),all(known_raw%in%c('true','false','1','0')));M$known_patient<-known_raw%in%c('true','1');C<-as.matrix(read.delim(file.path(d,'paired_FB_counts.tsv.gz'),row.names=1,check.names=FALSE));B<-read.delim(file.path(d,'paired_FB_profiles.tsv'),check.names=FALSE,stringsAsFactors=FALSE);D<-as.matrix(read.delim(file.path(d,'all_profile_detection.tsv.gz'),row.names=1,check.names=FALSE));E<-as.matrix(read.delim(file.path(d,'all_profile_meanlog_CP10K.tsv.gz'),row.names=1,check.names=FALSE));stopifnot(identical(colnames(C),B$profile_id),identical(colnames(D),M$profile_id),identical(rownames(D),rownames(E)))
 known<-M$known_patient;Q<-M[M$group=='FB01'&known&M$n_cells>=f$profile_min_cells,,drop=FALSE];bg0<-intersect(rownames(D)[rowSums(D[,Q$profile_id,drop=FALSE]>=f$detection_fraction)>0],rownames(LT))
 # All source-specific LR patient rows; absent genes/support are explicit, no cross-patient mixing.
 ledger<-list();j<-0L
 for(pat in unlist(f$cohorts[[st]]$patients)){
  rr<-M[M$patient==pat&M$group=='FB01',,drop=FALSE];send<-M[M$patient==pat,,drop=FALSE];sender_groups<-sort(unique(M$group[M$known_patient]))
  for(sender_group in sender_groups){ss<-send[send$group==sender_group,,drop=FALSE];if(!nrow(ss))ss<-data.frame(group=sender_group,n_cells=0L,profile_id=NA_character_);rrOK<-nrow(rr)==1L&&rr$n_cells>=f$profile_min_cells;ssOK<-ss$n_cells>=f$profile_min_cells
   lv<-lr$from%in%rownames(D);rv<-lr$to%in%rownames(D);ld<-rep(NA_real_,nrow(lr));rd<-ld;lm<-ld;rm<-ld
   if(!is.na(ss$profile_id)){ld[lv]<-D[lr$from[lv],ss$profile_id];lm[lv]<-E[lr$from[lv],ss$profile_id]}
   if(nrow(rr)==1L){rd[rv]<-D[lr$to[rv],rr$profile_id];rm[rv]<-E[lr$to[rv],rr$profile_id]}
   ok<-ssOK&rrOK&lv&rv&!is.na(ld)&!is.na(rd)&ld>=f$detection_fraction&rd>=f$detection_fraction;j<-j+1L
   ledger[[j]]<-data.frame(study=st,patient=pat,sender=ss$group,receiver='FB01',ligand=lr$from,receptor=lr$to,sender_cells=ss$n_cells,receiver_cells=if(nrow(rr))rr$n_cells else 0L,ligand_assayed=lv,receptor_assayed=rv,ligand_detection=ld,receptor_detection=rd,ligand_meanlogCP10K=lm,receptor_meanlogCP10K=rm,same_patient_supported=ok,complex_completeness='NOT_ESTABLISHED_BY_MONOMERIC_LR_PRIOR',stringsAsFactors=FALSE)
  }
 }
 L<-do.call(rbind,ledger);wgz(L,file.path(dest,'all_LR_sender_patient_support.tsv.gz'));alllr[[st]]<-aggregate(same_patient_supported~study+sender+ligand+receptor,data=L,FUN=sum)
 complete<-unlist(f$cohorts[[st]]$patients);complete<-complete[sapply(complete,function(z){b<-B[B$patient==z,,drop=FALSE];nrow(b)==2L&&setequal(b$group,c('FB01','otherFB'))&&all(b$n_cells>=f$profile_min_cells)})]
 if(length(complete)<2L){allcov[[st]]<-data.frame(study=st,scope='main',view=NA_character_,sign=NA_character_,background_n=NA_integer_,target_n=NA_integer_,eligible_ligands_n=NA_integer_,n_pairs=length(complete),status='INSUFFICIENT_QUALIFIED_PATIENT_PAIRS');next}
 response_models<-c(list(main=complete),setNames(lapply(complete,function(z)setdiff(complete,z)),paste0('LODO_',complete)))
 for(scope in names(response_models)){
  pats<-response_models[[scope]];if(length(pats)<2L)next;b<-B[B$patient%in%pats,,drop=FALSE];counts<-C[,b$profile_id,drop=FALSE];b$patient<-factor(b$patient);b$group<-factor(b$group,levels=c('otherFB','FB01'));design<-model.matrix(~patient+group,b);coef<-which(colnames(design)=='groupFB01');stopifnot(qr(design)$rank==ncol(design))
  y<-DGEList(counts=counts,lib.size=b$full_library_UMI);keep<-filterByExpr(y,design=design,min.count=10,min.total.count=15);y<-calcNormFactors(y[keep,,keep.lib.sizes=TRUE],method='TMM');v<-voom(y,design,plot=FALSE);fit<-eBayes(lmFit(v,design),robust=TRUE);T<-topTable(fit,coef=coef,number=Inf,sort.by='none');T$gene_symbol<-rownames(T);whole<-merge(data.frame(gene_symbol=rownames(C)),T,by='gene_symbol',all.x=TRUE,sort=FALSE);whole$tested<-whole$gene_symbol%in%rownames(T);whole$study<-st;whole$scope<-scope;whole$patient_pairs<-length(pats);compare_original_effect(whole,file.path(here,'response_activity_v1',st,paste0(scope,'__complete_gene_effects_BH.tsv.gz')));wgz(whole,file.path(dest,paste0(scope,'__complete_gene_effects_BH.tsv.gz')))
  if(scope=='main'){
   lc<-cpm(y,log=TRUE,prior.count=0.25);delta<-sapply(pats,function(z)lc[,b$patient==z&b$group=='FB01']-lc[,b$patient==z&b$group=='otherFB']);colnames(delta)<-pats;write(data.frame(gene_symbol=rownames(delta),delta,check.names=FALSE),file.path(dest,'patient_FB01_minus_otherFB_TMM_logCPM.tsv'))
  }
  for(view in names(f$views_excluded_symbols)){
   scopeQ<-Q[Q$patient%in%pats,,drop=FALSE];scopeBG<-intersect(rownames(D)[rowSums(D[,scopeQ$profile_id,drop=FALSE]>=f$detection_fraction)>0],rownames(LT));bg<-setdiff(scopeBG,unlist(f$views_excluded_symbols[[view]]));eligible<-unique(L$ligand[L$patient%in%pats&L$same_patient_supported]);t<-T[T$gene_symbol%in%bg,,drop=FALSE]
   for(sign in c('Up','Down')){
    z<-if(sign=='Up')t[t$t>0&t$logFC>0,,drop=FALSE]else t[t$t<0&t$logFC<0,,drop=FALSE];z<-z[order(if(sign=='Up')-z$t else z$t,z$gene_symbol),,drop=FALSE];target<-head(z$gene_symbol,f$response_top_n);key<-paste(st,scope,view,sign,sep='__');allact[[key]]<-score(bg,target,eligible,st,view,sign,scope,'allact',key);allcov[[key]]<-data.frame(study=st,scope=scope,view=view,sign=sign,background_n=length(bg),target_n=length(target),eligible_ligands_n=length(intersect(eligible,colnames(LT))),n_pairs=length(pats),status='RESPONSE_AND_COVERAGE_RECORDED');if(length(target))allsets[[key]]<-data.frame(study=st,scope=scope,view=view,sign=sign,gene_symbol=target,position=seq_along(target))
    if(scope=='main'){
     # Patient-specific compatibility is descriptive; not four independent clinical validations.
     for(pat in pats){dd<-delta[,pat];names(dd)<-rownames(delta);dd<-dd[names(dd)%in%bg];dd<-dd[if(sign=='Up')dd>0 else dd<0];dd<-sort(dd,decreasing=sign=='Up');tg<-head(names(dd),f$response_top_n);elig<-unique(L$ligand[L$patient==pat&L$same_patient_supported]);allpatient[[paste(key,pat)]]<-score(bg,tg,elig,st,view,sign,paste0('PATIENT_',pat),'allpatient',paste(key,pat))}
    }
   }
  }
 }
 cat(st,'response/activity complete\\n');gc()
}
for(kind in c('allact','allpatient')){stopifnot(setequal(names(get(kind)),expected_keys[[kind]]));stopifnot(all(names(legacy[[kind]])%in%names(get(kind))))}
write(do.call(rbind,cache_log),file.path(out,'score_checkpoint_provenance.tsv'))
atom_json(list(status='COMPLETE_CACHED_KEYS_ACCOUNTED',workspace_sha256=g$workspace_sha256,legacy_allact_available=length(legacy$allact),legacy_allpatient_available=length(legacy$allpatient),score_calls=score_call_count,expected_allact=length(expected_keys$allact),expected_allpatient=length(expected_keys$allpatient),legacy_input_boundary='Original v3 did not store per-call bg/eligible/patient-target vectors; exact original source and all immutable input identities plus deterministic context reconstruction used. New small checkpoints store complete call contexts.'),file.path(out,'WORKSPACE_RECOVERY_ACCOUNTING.json'))
bind<-function(x)if(length(x))do.call(rbind,x)else data.frame(status='NO_ELIGIBLE_RESULT')
wgz(bind(allact),file.path(out,'all_ligands_all_views_main_LODO.tsv.gz'));write(bind(allcov),file.path(out,'all_response_and_coverage.tsv'));wgz(bind(allsets),file.path(out,'all_response_members.tsv.gz'));wgz(bind(allpatient),file.path(out,'all_ligands_patient_descriptive.tsv.gz'));wgz(bind(alllr),file.path(out,'all_LR_patient_support_counts.tsv.gz'));writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
files<-list.files(out,recursive=TRUE,full.names=TRUE);files<-files[!dir.exists(files)];rows<-lapply(files,function(p)list(path=p,bytes=unname(file.info(p)$size),sha256=sha(p)));write_json(list(status='COMPLETE_S1_S2_TECHNICAL_ROOT_SCIENTIFIC_REVIEW_PENDING',freeze_sha256=sha(file.path(here,'FREEZE.json')),profile_complete_sha256=sha(file.path(prep,'COMPLETE.json')),activity_source_sha256=sha(resume_source),original_activity_source_sha256=original_source_sha,resume_execution_binding_sha256=sha(binding_path),workspace_sha256=g$workspace_sha256,equivalence_fixture_sha256=g$equivalence_fixture_sha256,parallel_execution=list(single=TRUE,mc_cores=8L,mc_set_seed=FALSE,BLAS_threads_per_worker=1L),outputs=rows,clinical_boundary='Candidate communication prior compatibility, not causal or clinical validation'),file.path(out,'COMPLETE.json.tmp'),pretty=TRUE,auto_unbox=TRUE,na='null');stopifnot(file.rename(file.path(out,'COMPLETE.json.tmp'),file.path(out,'COMPLETE.json')))
