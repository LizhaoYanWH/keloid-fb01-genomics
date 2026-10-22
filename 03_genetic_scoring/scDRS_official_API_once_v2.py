"""One whole-cohort official scDRS API invocation. Root SHA approval required; never batch-score cells."""
from pathlib import Path
import os,sys,json,time,hashlib,shutil,argparse,fcntl,multiprocessing as mp,signal
R=Path('/root/five_new_studies_20261020/projects/N02_keloid_genetic_cell_state');P=R/'reviews/scDRS_production_preparation_v1';F=P/'RNA_freeze_v1';G=P/'official_MAGMA_GS_v1';LOCK=R.parents[1]/'metadata/assess01_CPU_worker.lock';EXPECTED_PY=R/'tools/scdrs104_env/bin/python'
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb')as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()
def atom(p,j):
 q=p.with_suffix(p.suffix+'.tmp');q.write_text(json.dumps(j,indent=2,allow_nan=False)+'\n');q.replace(p)
def need(b,s):
 if not b:raise RuntimeError(s)
def resources():
 m={x.split(':')[0]:int(x.split()[1])*1024 for x in Path('/proc/meminfo').read_text().splitlines()};mx=int(Path('/sys/fs/cgroup/memory.max').read_text());used=int(Path('/sys/fs/cgroup/memory.current').read_text());return {'effective_free_bytes':min(m['MemAvailable'],mx-used),'cgroup_current_bytes':used,'cgroup_limit_bytes':mx,'system_free_bytes':shutil.disk_usage('/root').free,'data_free_bytes':shutil.disk_usage('/root/autodl-tmp').free}
def disk_guard(x):need(min(x['system_free_bytes'],x['data_free_bytes'])>=5*2**30,'Both disks must retain5GiB')
def process_identity(pid):
 p=Path('/proc')/str(pid);s=(p/'stat').read_text().split();return {'pid':pid,'start_ticks':int(s[21]),'boot_id':Path('/proc/sys/kernel/random/boot_id').read_text().strip()}
def score_worker(st,out,approval):
 try:
  # Imports and precision are inherited from the pinned isolated environment; no source changes.
  import numpy as np,pandas as pd,h5py,scipy.sparse as sp,scipy.io as sio,gzip,resource,gc
  import scanpy as sc,anndata as ad,scdrs
  from threadpoolctl import threadpool_limits
  phases={}
  def tick():
   u=resource.getrusage(resource.RUSAGE_SELF);return(time.monotonic(),u.ru_utime+u.ru_stime)
  def mark(name,t,accumulate=False):
   u=resource.getrusage(resource.RUSAGE_SELF);w=time.monotonic()-t[0];cpu=u.ru_utime+u.ru_stime-t[1]
   if accumulate and name in phases:
    q=phases[name];q['wall_seconds']+=w;q['CPU_seconds_all_threads']+=cpu;q['calls']+=1;q['cumulative_peak_RSS_bytes_at_phase_end']=int(u.ru_maxrss*1024);q['mean_CPU_cores']=q['CPU_seconds_all_threads']/q['wall_seconds']
   else:phases[name]={'wall_seconds':w,'CPU_seconds_all_threads':cpu,'mean_CPU_cores':cpu/w if w else 0,'cumulative_peak_RSS_bytes_at_phase_end':int(u.ru_maxrss*1024),'calls':1}
  t=tick();d=F/st;freeze=json.loads((F/'COMPLETE.json').read_text());c=freeze['cohorts'][st];inp=c['input'];counts=Path(inp['counts_path']);need(sha(counts)==inp['counts_sha256'],'Counts source SHA mismatch')
  obs=pd.read_csv(d/'all_final_cell_metadata.tsv.gz',sep='\t',dtype=str,keep_default_na=False).set_index('cell_id');genes=pd.read_csv(d/'admitted_gene_axis.tsv',sep='\t',dtype={'canonical_entrez':str});cov=pd.read_csv(d/'technical_covariates.tsv.gz',sep='\t',index_col=0);need(list(obs.index)==list(cov.index),'Covariate order mismatch')
  if st=='GSE241132_D0':
   with gzip.open(counts,'rb')as f,threadpool_limits(limits=approval['threads']):X=sio.mmread(f).T.tocsr()
  else:
   with h5py.File(counts)as h:
    z=h['layers/counts'];X=sp.csr_matrix((z['data'][:],z['indices'][:],z['indptr'][:]),shape=tuple(z.attrs['shape']))
  need(X.shape==(c['n_cells'],c['input_genes']),'Counts schema changed');X.sum_duplicates();X.eliminate_zeros();need(np.issubdtype(X.dtype,np.integer)and(X.data>=0).all(),'Raw count unit changed')
  U=np.asarray(X.sum(1)).ravel();N=np.diff(X.indptr);need(np.array_equal(U,obs.full_input_UMI.astype('int64').to_numpy())and np.array_equal(N,obs.full_input_nGenes.astype('int64').to_numpy()),'Full-input technical cov source changed');need(np.allclose(cov.iloc[:,0],np.log1p(U),rtol=0,atol=1e-14)and np.allclose(cov.iloc[:,1],np.log1p(N),rtol=0,atol=1e-14),'Covariates changed')
  X=X[:,genes.source_index.to_numpy()].tocsr();a=ad.AnnData(X=X.astype(np.float64),obs=obs,var=pd.DataFrame(index=genes.canonical_entrez));del X;gc.collect();need(a.obs_names.is_unique and a.var_names.is_unique,'Nonunique admitted axes');need((np.asarray(a.X.sum(1)).ravel()>0).all(),'Zero library after fixed mapping; preserve failure, no cell deletion')
  mark('load_count_and_identity_check',t);t=tick()
  # Exact integer counts cast float64 avoids legacy Scanpy int64 divide casting failure; then unchanged official normalization recipe.
  before=list(a.obs_names);sc.pp.normalize_per_cell(a,counts_per_cell_after=1e4);sc.pp.log1p(a);need(before==list(a.obs_names),'Official normalization dropped frozen cells');mark('normalize_log_once',t);t=tick()
  scdrs.preprocess(a,cov=cov,adj_prop=None,n_mean_bin=20,n_var_bin=20,n_chunk=128);mark('official_preprocess',t);t=tick();a.uns['SCDRS_PARAM']['GENE_STATS'].to_csv(out/'actual_gene_mean_variance_bins.tsv.gz',sep='\t');a.uns['SCDRS_PARAM']['COV_MAT'].to_csv(out/'actual_covariate_design.tsv.gz',sep='\t')
  gs=pd.read_csv(G/'Keloid_EUR.top1000_Z.gs',sep='\t');member=[x.rsplit(':',1)for x in gs.GENESET.iloc[0].split(',')];gl=[x[0]for x in member];gw=[float(x[1])for x in member];need(len(gl)==1000,'Frozen GS must contain1000 original members')
  # Transparent capture of the one actual official selector call; arguments/returns remain unaltered.
  selector=scdrs.method._select_ctrl_geneset;selection_calls=[];raw_fn=scdrs.method._compute_raw_score;background_fn=scdrs.method._correct_background;pooled_fn=scdrs.method._get_p_from_empi_null;raw_calls=[0]
  def raw_timed(*args,**kwargs):
   t0=tick();v=raw_fn(*args,**kwargs);mark('raw_trait_kernel'if raw_calls[0]==0 else'raw1000_control_kernels',t0,accumulate=True);raw_calls[0]+=1;return v
  def background_timed(*args,**kwargs):
   t0=tick();v=background_fn(*args,**kwargs);mark('whole_cohort_global_background_normalization',t0);return v
  def pooled_timed(*args,**kwargs):
   t0=tick();v=pooled_fn(*args,**kwargs);mark('whole_cohort_pooled_empirical_P',t0);return v
  def capture(*args,**kwargs):
   t0=tick();lists,weights=selector(*args,**kwargs);mark('official_control_gene_selection',t0)
   with gzip.open(out/'actual_control_gene_sets_and_input_weights.tsv.gz','wt')as w:
    w.write('control_id\tgene_position\tcanonical_entrez\tinput_weight\n')
    for k in sorted(lists):
     for j,(g,z)in enumerate(zip(lists[k],weights[k])):w.write(f'{k}\t{j}\t{g}\t{z:.17g}\n')
   selection_calls.append({'n_ctrl':len(lists),'disease_members_measurable':len(args[1]),'seed':args[6],'mean_var_matching':args[3]});return lists,weights
  scdrs.method._select_ctrl_geneset=capture;scdrs.method._compute_raw_score=raw_timed;scdrs.method._correct_background=background_timed;scdrs.method._get_p_from_empi_null=pooled_timed
  try:y=scdrs.score_cell(a,gl,gene_weight=gw,ctrl_match_key='mean_var',n_ctrl=1000,weight_opt='vs',return_ctrl_raw_score=False,return_ctrl_norm_score=True,random_seed=20261020)
  finally:scdrs.method._select_ctrl_geneset=selector;scdrs.method._compute_raw_score=raw_fn;scdrs.method._correct_background=background_fn;scdrs.method._get_p_from_empi_null=pooled_fn
  need(len(selection_calls)==1,'Exactly one official whole-cohort control selection required');mark('official_score1000_controls_and_background',t);t=tick();need(y.index.equals(obs.index)and y.shape==(len(obs),1006),'Whole-cohort output schema mismatch');need(str(y.values.dtype)=='float32','Pinned official score dtype changed')
  cache=np.lib.format.open_memmap(out/'full_score_with1000_norm_controls.float32.npy',mode='w+',dtype='float32',shape=y.shape)
  for start in range(0,len(y),4096):cache[start:start+4096]=y.iloc[start:start+4096].to_numpy(dtype='float32',copy=False)
  cache.flush();del cache;y.iloc[:,:6].to_csv(out/'six_cell_scores.tsv.gz',sep='\t');pd.Series(y.index,name='cell_id').to_csv(out/'cell_axis.tsv',sep='\t',index=False);atom(out/'cache_columns.json',list(y.columns))
  # Reuse the same cache. Official association-only q95 MC formula, without unnecessary Geary/KNN.
  control=[f'ctrl_norm_score_{i}'for i in range(1000)];group_rows=[]
  for family,col,levels in [('lineage','lineage',c['coarse_lineage_family']),('FB_state','FB_state',c['FB_state_family'])]:
   for level in levels:
    ix=np.flatnonzero(obs[col].eq(level).to_numpy());r={'cohort':st,'family':family,'group':level,'n_cells':len(ix),'n_ctrl':1000,'FDR_pending_full_three_cohort_family':True}
    if len(ix):
     q=float(np.quantile(y.iloc[ix].norm_score,0.95));z=np.quantile(y.iloc[ix][control],.95,axis=0);sd=float(z.std());r.update(score_q95=q,assoc_mcp=float(((z>=q).sum()+1)/1001),assoc_mcz=float((q-z.mean())/sd)if sd>0 else None)
    else:r.update(score_q95=None,assoc_mcp=None,assoc_mcz=None)
    group_rows.append(r)
  pd.DataFrame(group_rows).to_csv(out/'all_frozen_lineage_FBstate_association_MC_unadjusted.tsv',sep='\t',index=False)
  meta=obs.join(y[['norm_score','mc_pval','pval']],how='left');summ=[]
  for view,col,keep in [('original_GEO','region_original',np.ones(len(meta),dtype=bool)),('P2_label_swap','region_P2_swap',np.ones(len(meta),dtype=bool)),('P2_omit','region_original',meta.include_P2_omit.isin(['True','true','1']).to_numpy())]:
   z=meta.loc[keep].copy();z['region_view']=z[col]
   for keys,sub in z.groupby(['independent_patient_id','capture_id','region_view','lineage','FB_state'],dropna=False,observed=True,sort=True):summ.append({'cohort':st,'view':view,**dict(zip(['independent_patient_id','capture_id','region','lineage','FB_state'],keys)),'n_cells':len(sub),'mean_norm_score':float(sub.norm_score.mean()),'median_norm_score':float(sub.norm_score.median()),'q95_norm_score':float(sub.norm_score.quantile(.95)),'descriptive_only_not_independent_cell_replication':True})
  pd.DataFrame(summ).to_csv(out/'patient_capture_lineage_FBstate_P2_views_descriptive.tsv',sep='\t',index=False);mark('cache_and_all_fixed_summary',t)
  diagnostics={'actual_raw_trait_plus_control_kernel_calls':raw_calls[0],'actual_cpu_quota_file':Path('/sys/fs/cgroup/cpu.max').read_text().strip(),'affinity_cpus':sorted(os.sched_getaffinity(0)),'n_scores':int(y.size),'nonfinite_total':int((~np.isfinite(y.values)).sum()),'nonfinite_by_column':{str(k):int(v)for k,v in (~np.isfinite(y)).sum().items()},'GS_requested':1000,'GS_measurable':selection_calls[0]['disease_members_measurable'],'selection_calls':selection_calls,'cell_count':len(y),'no_cell_batch_scoring':True,'no_FDR_until_complete_frozen_family':True}
  rec={'status':'INTERNAL_WHOLE_COHORT_SCORE_RETURNED_PARENT_VALIDATION_PENDING','cohort':st,'phases':phases,'peak_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,'diagnostics':diagnostics,'normalization':'official normalize_per_cell counts_after1e4 thenlog1p; after fixed min50/one-to-one gene mapping; all final cells kept','internal_dtype':str(a.X.dtype),'returned_dtype':str(y.values.dtype),'outputs':{f.name:{'bytes':f.stat().st_size,'sha256':sha(f)}for f in out.iterdir()if f.is_file()and f.name not in['RUNNING.json','ROOT_APPROVAL_USED.json']}}
  atom(out/'INTERNAL_COMPLETE.json',rec)
 except BaseException as e:atom(out/'INTERNAL_FAILURE.json',{'error':repr(e),'cohort':st});raise

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--approval',required=True);ap.add_argument('--approval-sha256',required=True);a=ap.parse_args();af=Path(a.approval);need(sha(af)==a.approval_sha256,'Root approval SHA mismatch');j=json.loads(af.read_text());need(j['action']=='SCORE_ONE_WHOLE_COHORT_OFFICIAL_API','Approval action mismatch');need(sha(Path(__file__))==j['executor_sha256'],'Executor SHA mismatch');need(sys.executable==str(EXPECTED_PY),'Use exact isolated interpreter');st=j['cohort'];need(st in['GSE241132_D0','GSE220300','GSE335482'],'Unexpected cohort');need(j['threads']in[4,8],'Root-approved4or8threads required')
 for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]=str(j['threads'])
 allowed=sorted(os.sched_getaffinity(0));quota,period=Path('/sys/fs/cgroup/cpu.max').read_text().split();need(len(allowed)>=j['threads']and(quota=='max'or int(quota)/int(period)>=j['threads']),'CPU quota cannot satisfy approved threads');os.sched_setaffinity(0,allowed[-j['threads']:])
 required=[(P/'environment_install_receipt_v1.json','environment_install_receipt_sha256'),(P/'synthetic_compatibility_receipt_v1.json','synthetic_receipt_sha256'),(F/'COMPLETE.json','RNA_freeze_sha256'),(G/'COMPLETE.json','GS_freeze_sha256')]
 for p,k in required:need(sha(p)==j[k],k+' mismatch')
 rf=json.loads((F/'COMPLETE.json').read_text());need(rf['status']=='COMPLETE_INPUT_FREEZE_NO_SCDRS_FIT','RNA incomplete');gf=json.loads((G/'COMPLETE.json').read_text());need(gf['status']=='COMPLETE_OFFICIAL_TOP1000_GS_NO_RNA_SELECTION_NO_SCORE','GS incomplete');need(gf['RNA_freeze_sha256']==j['RNA_freeze_sha256'],'GS/RNA cross-bind mismatch');need(gf['native_complete_sha256']==j['native_complete_sha256'],'Native cross-bind mismatch')
 for name,x in rf['cohorts'][st]['files'].items():p=F/st/name;need(p.stat().st_size==x['bytes']and sha(p)==x['sha256'],'RNA artifact changed '+name)
 for name,x in gf['files'].items():p=G/name;need(p.stat().st_size==x['bytes']and sha(p)==x['sha256'],'GS artifact changed '+name)
 ei=json.loads((P/'environment_install_receipt_v1.json').read_text());site=R/'tools/scdrs104_env/lib/python3.12/site-packages'
 for x in ei['verified_source_files']:
  if x['path'].startswith('scdrs/')and len(Path(x['path']).parts)==2:need(sha(site/x['path'])==x['sha256'],'Installed pinned source differs '+x['path'])
 out=R/f'results/scDRS_whole_cohort_v1/{st}';out.parent.mkdir(exist_ok=True,parents=True);need(not out.exists(),'Never overwrite prior partial/complete cohort');claim=out.parent/(st+'.once_launch_claim.json')
 with claim.open('x')as f:json.dump({**process_identity(os.getpid()),'executor_sha256':j['executor_sha256'],'approval_sha256':a.approval_sha256},f)
 with LOCK.open('a+')as lock:
  fcntl.flock(lock,fcntl.LOCK_EX);x=resources();disk_guard(x);need(x['system_free_bytes']>=5*2**30+j['output_reserve_bytes'],'Output reserve unavailable');need(x['effective_free_bytes']>=j['job_RAM_cap_bytes']+j['effective_memory_floor_bytes']+j['peer_reserved_bytes'],'Memory ownership admission failed; no automatic retry');out.mkdir();atom(out/'ROOT_APPROVAL_USED.json',j);atom(out/'RUNNING.json',{**process_identity(os.getpid()),'start_resources':x,'UTC':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())});start=time.monotonic();proc=mp.get_context('spawn').Process(target=score_worker,args=(st,out,j));proc.start();peak=0;failure=None
  while proc.is_alive():
   proc.join(5);x=resources()
   try:
    status=Path(f'/proc/{proc.pid}/status').read_text();rss=next((int(z.split()[1])*1024 for z in status.splitlines()if z.startswith('VmRSS:')),0);peak=max(peak,rss);disk_guard(x);need(rss<=j['job_RAM_cap_bytes'],'Owned RSS cap exceeded');need(x['effective_free_bytes']>=j['effective_memory_floor_bytes'],'Effective cgroup/host memory floor violated');need(time.monotonic()-start<=j['max_wall_seconds'],'Approved wall cap exceeded')
   except FileNotFoundError:pass
   except BaseException as e:failure=repr(e);proc.terminate();proc.join(30);break
  rc=proc.exitcode
  if failure or rc!=0:atom(out/'FAILED_PRESERVED_NO_AUTOMATIC_RETRY.json',{'error':failure,'child_returncode':rc,'peak_owned_RSS_bytes':peak,'last_resources':resources()});raise RuntimeError('Scientific process failed/guarded; retained partials')
  p=out/'INTERNAL_COMPLETE.json';need(p.exists(),'Missing whole-cohort internal result');result=json.loads(p.read_text())
  for name,z in result['outputs'].items():q=out/name;need(q.stat().st_size==z['bytes']and sha(q)==z['sha256'],'Final artifact validation failed '+name)
  need(result['diagnostics']['cell_count']==rf['cohorts'][st]['n_cells'],'Frozen cell coverage lost');need(result['diagnostics']['selection_calls'][0]['n_ctrl']==1000,'Control count changed');need(sha(Path(__file__))==j['executor_sha256'],'Source changed during score');atom(out/'COMPLETE.json',{'status':'COMPLETE_ONE_WHOLE_COHORT_OFFICIAL_API_NO_BATCH_SPLITTING','approval_sha256':a.approval_sha256,'executor_sha256':j['executor_sha256'],'internal_complete_sha256':sha(p),'frozen_input_receipts':{k:j[k]for _,k in required},'native_complete_sha256':j['native_complete_sha256'],'elapsed_seconds':time.monotonic()-start,'max_owned_RSS_bytes':max(peak,result['peak_RSS_bytes']),'returncode':rc,'diagnostics':result['diagnostics'],'outputs':result['outputs'],'root_scientific_and_full_family_FDR_review_pending':True,'final_resources':resources()});print(json.dumps({'COMPLETE':str(out/'COMPLETE.json'),'SHA256':sha(out/'COMPLETE.json')}))
if __name__=='__main__':main()
