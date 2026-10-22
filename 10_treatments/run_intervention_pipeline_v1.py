"""One input freeze, then one root-approved bounded paired-public-RNA run."""
from pathlib import Path
import argparse,os,sys,json,hashlib,subprocess,shutil,time,gzip,traceback,resource,io
D=Path('/root/five_new_studies_20261020/projects/N02_keloid_genetic_cell_state/mechanisms_20261021/literature_intervention_v1')
P=D.parents[1];OLD=Path('/root/autodl-tmp/Five_Studies_Execution_20261020');RS=Path('/root/five_studies_envs/r/bin/Rscript');R=D/'run_paired_intervention_analysis_v1.R';CAP=128*2**20

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()
def rec(p):return {'path':str(p),'bytes':p.stat().st_size,'sha256':sha(p)}
def atom(p,j):
 q=Path(str(p)+'.tmp');q.write_text(json.dumps(j,indent=2,ensure_ascii=False,allow_nan=False)+'\n');q.replace(p)
def verify(r):
 p=Path(r['path']);need(p.stat().st_size==r['bytes'] and sha(p)==r['sha256'],'Input changed '+str(p))
def usage():
 mem={x.split(':')[0]:int(x.split()[1])*1024 for x in Path('/proc/meminfo').read_text().splitlines()};free=min(mem['MemAvailable'],int(Path('/sys/fs/cgroup/memory.max').read_text())-int(Path('/sys/fs/cgroup/memory.current').read_text()))
 return {'effective_free_bytes':free,'system_free_bytes':shutil.disk_usage('/root').free,'data_free_bytes':shutil.disk_usage('/root/autodl-tmp').free,'owned_directory_bytes':sum(p.stat().st_size for p in D.rglob('*') if p.is_file())}
def guard(job=0,peer=12):
 r=usage();need(r['effective_free_bytes']>=(job+8+peer)*2**30,'Job+8floor+peer RAM guard')
 need(r['data_free_bytes']>=5*2**30 and r['system_free_bytes']>=6*2**30+max(0,CAP-r['owned_directory_bytes']),'Two5GiB plus original02future1GiB and remaining128MiB output budget')
 need(r['owned_directory_bytes']<CAP,'128MiB owned source+result cap');return r

def freeze(args):
 import pandas as pd,numpy as np
 need(args.expected_source_sha256==sha(Path(__file__)) and not(D/'FREEZE.json').exists(),'Source SHA/fresh freeze')
 os.sched_setaffinity(0,[min(os.sched_getaffinity(0))]);guard(job=2,peer=16)
 admission=D/'ADMISSION_COMPLETE.json';a=json.loads(admission.read_text());need(a['no_differential_expression_or_drug_response_computed'],'Do not freeze after viewing drug response')
 for r in a['inputs']+a['outputs']:verify(r)
 bg=P/'extensions_20261021/specificity_locus_v1/FREEZE.json';need(sha(bg)=='d7d8d74a641078f35baa3b0cfbdee3a23ffff8daf9f9c11bd82d6c1e7c1d83ab','Original frozen CellCycle changed')
 b=json.loads(bg.read_text());growth=b['fixed_modules']['CellCycle']['canonical_entrez'];ecm=b['fixed_modules']['ECM']['canonical_entrez'];need(len(growth)==689 and len(ecm)==321,'Original689CellCycle/321ECM required')
 h=pd.read_csv(D/'all_unique_approved_HGNC_gene_axis.tsv',sep='\t',dtype=str,keep_default_na=False);ct=h.loc[h.entrez_id.ne(''),'entrez_id'].value_counts();mp={r.entrez_id:r.ensembl_gene_id for r in h.itertuples() if r.entrez_id and ct[r.entrez_id]==1}
 pr=pd.read_csv(D/'fixed2_program_members_and_weights.tsv',sep='\t',keep_default_na=False);more=[]
 for name,members in [('fixed_Reactome_CellCycle_comparator',growth),('fixed_Reactome_ECM_comparator',ecm)]:
  more.extend({'program':name,'source_id':str(x),'ensembl_gene_id':mp.get(str(x),''),'weight':1.,'mapping_status':'UNIQUE_HGNC_ENTREZ_TO_ENSEMBL' if str(x) in mp else 'UNRESOLVED_NA_NO_BACKFILL'} for x in members)
 pd.concat([pr,pd.DataFrame(more)],ignore_index=True).to_csv(D/'fixed4_program_members_and_weights.tsv',sep='\t',index=False)
 # Read actual complete HC export scale; no effects, contrasts or P-values.
 path=D/'sources/GSE7890_series_matrix.txt.gz';inside=False;body=[]
 for l in gzip.open(path,'rt'):
  if l.startswith('!series_matrix_table_begin'):inside=True;continue
  if l.startswith('!series_matrix_table_end'):inside=False;continue
  if inside:body.append(l)
 x=pd.read_csv(io.StringIO(''.join(body)),sep='\t',index_col=0);quant=[0,.01,.25,.5,.75,.9,.99,1];vv=x.to_numpy();q=np.quantile(vv,quant);need((vv>0).all() and q[-2]>100 and q[-1]>1000,'HC exported scale is not demonstrated positive nonlog; stop before transformation')
 scale={'source':rec(path),'shape':list(x.shape),'whole_matrix_quantiles':{str(k):float(v) for k,v in zip(quant,q)},'each_column_quantiles':{c:{str(k):float(v) for k,v in zip(quant,np.quantile(x[c],quant))} for c in x.columns},'provider_method':'GeneSpring7 RMA; actual stored values are positive linear normalized signal export, not log2 numeric scale','chosen_transform':'log2 each positive provider intensity, then average all uniquely HGNC-mapped probes; no new RMA/quantile normalization','no_response_statistics':True};atom(D/'HC_COMPLETE_MATRIX_SCALE_PREFLIGHT_v1.json',scale)
 gene_axes={}
 for ds,fn in [('GSE232079','GSE232079_counts.txt.gz'),('GSE232435','GSE232435_counts_ready.txt.gz')]:
  with gzip.open(D/'sources'/fn,'rt') as handle:
   next(handle);ids=[line.split('\t',1)[0] for line in handle]
  gene_axes[ds]={'rows':len(ids),'version_suffixed_rows':sum(bool(__import__('re').search(r'\.[0-9]+$',z)) for z in ids),'version_stripping_performed':False,'mapping_rule':'Same exact stable Ensembl IDs as input-admission, direct HGNC membership in R'}
 atom(D/'SOURCE_GENE_AXIS_SEMANTICS_v1.json',gene_axes)
 rp=OLD/'shared_inputs/programs_01_03/sources/NCBI2Reactome_All_Levels.txt';prior=OLD/'shared/programs/curated_regulatory_LR_v1/CollecTRI_human_interactions.tsv.gz';need(sha(prior)=='c3703341dde2254eee4e9c4a5ba65179489eef905e90afd4cc9ed52ccaba2ae4','Shared CollecTRI source mismatch')
 cmd=[str(RS),'-e','invisible(parse(file='+json.dumps(str(R))+'));for(p in c("edgeR","limma","statmod","jsonlite")){stopifnot(requireNamespace(p,quietly=TRUE));cat(p,as.character(packageVersion(p)),find.package(p));cat(intToUtf8(10))}']
 rr=subprocess.run(cmd,capture_output=True,text=True);(D/'R_PARSE_AND_EXISTING_PACKAGES_PREFLIGHT.txt').write_text(rr.stdout+rr.stderr);need(rr.returncode==0,'Existing R parse/environment preflight failed')
 inp=[Path(r['path']) for r in a['inputs']]+[Path(r['path']) for r in a['outputs']]+[admission,bg,rp,prior,D/'fixed4_program_members_and_weights.tsv',D/'HC_COMPLETE_MATRIX_SCALE_PREFLIGHT_v1.json',D/'SOURCE_GENE_AXIS_SEMANTICS_v1.json',D/'R_PARSE_AND_EXISTING_PACKAGES_PREFLIGHT.txt',RS,R]
 for line in rr.stdout.splitlines():
  vals=line.split()
  if len(vals)>=3:
   pkg=Path(vals[2]);inp += [pkg/'DESCRIPTION']+[q for q in (pkg/'R').glob('*.rdb')]
 j={'status':'INPUT_AND_SOURCE_FROZEN_NO_DRUG_RESPONSE_RESULTS','source_sha256':sha(Path(__file__)),'R_source_sha256':sha(R),'input_admission_sha256':sha(admission),'inputs':[rec(p) for p in dict.fromkeys(inp)],'sample_units':'Corin4donor paired after2culture repeats pooled withinpatient/condition; other3material sources descriptive. TAK1N5,N6,N8 biological source-lines paired3conditions. HC5KL and4ordinaryscar pairedstrains; untreated58 preserved descriptive only.','programs':['fixed_shared_nonGS_FB01_25','original_MAGMA1000_Z_RNA','fixed_Reactome_CellCycle_comparator','fixed_Reactome_ECM_comparator'],'program_semantics':'MAGMAZweights are association strength not eQTL effect signs or individualPRS. RNA programme decrease cannot be called geneticrisk decrease.','nuisance_sources':'Original already frozen689canonicalEntrez CellCycle R-HSA-1640170 and321ECM R-HSA-1474244, officialReactome; no newHallmark or effect-selected sets','contrasts7':['Corin_minus_DMSO_primary4','TGFbeta_minus_Control','HS276_TGFbeta_minus_TGFbeta','HS276_TGFbeta_minus_Control','HC_minus_NoHC_keloid5','HC_minus_NoHC_normalScar4','HC_response_keloid_minus_normalScar'],'RNA_counts_model':'uniqueapprovedHGNCallgene axis; preserve all rawrows; sum culture repeats within samepatient/condition; edgeR filterByExpr(design,min.count10,min.total15), TMM,voom,paired limma moderatedt robust; eachcontrast geneBHwholeHGNC; all filtered/unassayed rowsNA retained','array_model':'OriginalGEO RMA positive normalized export log2 once; officialGPL570 uniqueEntrez->HGNC; average all unambiguous probes per gene; no effect-based probe choice; pairedstrain+HC+KHC interaction','programme_score':'Within assay rankpercentile over complete measuredHGNClogexpression background, mean for25/growth or originalZweightedmean for1000; missing genes no backfill. Sourcepaired delta and95tCI descriptive; HCinteraction WelchCI; no programme patientP.','competitive_selfcontained':'CAMERA rho.01 andFRY on all expressed HGNC gene universe, unweighted fixed membership even for1000 (Zweights used only rankRNAprogramme). Original2x7 BH14 primary, ECM/growth2x7 BH14, joint4x7 BH28 and FRYmixed separately; gene-setP is notpatientP.','candidate_family':'All9x7=63 rows; percontrastwholegeneFDR and jointBH63 with missingP at1 for multiplicity butNA outputs.','regulatory_handoff':'Save complete7contrast t/logFC/HGNCtables for otherowner same frozenCollecTRI fullnetULM/MLM; no locally reinvented TF scores','technical_outputs':'All samplesQC/library/PCA PDFs and complete tables; no PCA/effect-driven exclusions; resource/format failures preserved notbiologicalnegative','clinical_boundary':'Human fibroblast transcriptional response only, not clinical efficacy/targetengagement/causal mediation; CorinCoREST multitarget notNEDD4 drug; HS276TAK1 notTAB2; hydrocortisone notTACclinicalresponse; ordinaryscar nothealthy','runtime_limits':{'RAM_GiB':8,'CPU':4,'wall_seconds_max':1800,'owned_total_cap_bytes':CAP,'data_writes':0,'disk_floors_GiB_each':5,'old02_future_system_GiB_reserved':1},'resources':usage()};atom(D/'FREEZE.json',j);print(json.dumps({'status':j['status'],'source_sha256':j['source_sha256'],'R_source_sha256':j['R_source_sha256'],'FREEZE_sha256':sha(D/'FREEZE.json'),'input_count':len(j['inputs']),'HC_all_matrix_quantiles':scale['whole_matrix_quantiles'],'resources':usage()}))

def run(args):
 ap=Path(args.approval);need(sha(ap)==args.approval_sha256,'Root approvalSHA');a=json.loads(ap.read_text());need(a['action']=='RUN_FIXED_PAIRED_PUBLIC_INTERVENTION_RNA','Root action')
 need(a['source_sha256']==sha(Path(__file__)) and a['R_source_sha256']==sha(R) and a['FREEZE_sha256']==sha(D/'FREEZE.json'),'Root source/freeze binding')
 f=json.loads((D/'FREEZE.json').read_text());need(f['source_sha256']==sha(Path(__file__)) and f['R_source_sha256']==sha(R),'Frozen source changed')
 for r in f['inputs']:verify(r)
 need(a['RAM_GiB']==8 and a['CPU']==4 and a['wall_seconds']<=1800,'Finite resource recipe')
 cp=a['CPU_ids'];need(len(cp)==4 and set(cp)<=os.sched_getaffinity(0),'Four grantedCPUaffinity');os.sched_setaffinity(0,cp)
 for n in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:os.environ[n]='4'
 guard(job=8,peer=a.get('peer_RAM_GiB',12));out=D/'results_v1';need(not out.exists() and not(D/'ONCE_CLAIM_v1.json').exists(),'No duplicate science or overwrite')
 with(D/'ONCE_CLAIM_v1.json').open('x')as h:json.dump({'PID':os.getpid(),'boot_id':Path('/proc/sys/kernel/random/boot_id').read_text().strip(),'source_sha256':sha(Path(__file__)),'approval_sha256':sha(ap)},h)
 out.mkdir();tmp=D/'owned_R_tmp';tmp.mkdir(exist_ok=False);env=os.environ.copy();env['TMPDIR']=str(tmp);start=time.monotonic();problem=None
 try:
  with(out/'R_run.log').open('x')as log:
   child=subprocess.Popen([str(RS),'--vanilla',str(R),str(out)],stdout=log,stderr=subprocess.STDOUT,env=env,cwd=D)
   st=Path(f'/proc/{child.pid}/stat').read_text().rsplit(')',1)[1].split();atom(D/'ACTUAL_R_CHILD_v1.json',{'PID':child.pid,'start_ticks':int(st[19]),'parent_PID':os.getpid(),'CPU_ids':cp,'command':[str(RS),'--vanilla',str(R),str(out)]})
   while child.poll() is None:
    time.sleep(1);rr=usage();need(rr['effective_free_bytes']>=8*2**30,'Runtime8GiB RAMfloor');need(rr['data_free_bytes']>=5*2**30 and rr['system_free_bytes']>=6*2**30,'Both5GiB plus original02future1GiB');need(rr['owned_directory_bytes']<CAP,'Runtime128MiB cap');need(time.monotonic()-start<=a['wall_seconds'],'Runtimewallcap')
    try:
     status=Path(f'/proc/{child.pid}/status').read_text();rss=next((int(z.split()[1])*1024 for z in status.splitlines() if z.startswith('VmRSS:')),0);need(rss<=8*2**30,'Owned8GiB RAMcap')
    except FileNotFoundError:pass
   need(child.returncode==0,'R scientific run failed; preserve outputs, no automatic retry')
  need((out/'R_STATISTICS_COMPLETE.json').exists(),'R actual complete marker missing')
  for q in out.iterdir():
   if q.is_file() and q.suffix=='.gz':
    with gzip.open(q,'rb')as z:
     for _ in iter(lambda:z.read(2**20),b''):pass
  atom(out/'COMPLETE.json',{'status':'COMPLETE_FIXED7_PAIRED_PUBLIC_INTERVENTION_RNA_ROOT_REVIEW_PENDING','returncode':0,'source_sha256':sha(Path(__file__)),'R_source_sha256':sha(R),'FREEZE_sha256':sha(D/'FREEZE.json'),'approval_sha256':sha(ap),'elapsed_seconds':time.monotonic()-start,'all_outputs':{q.name:rec(q) for q in out.iterdir() if q.is_file()},'resources':usage(),'no_other_worker_signals_or_raw_mutation':True})
 except BaseException as e:
  problem=repr(e)
  if 'child' in locals() and child.poll() is None:child.terminate();child.wait(timeout=30)
  atom(D/'FAILED_PRESERVED_v1.json',{'error':problem,'traceback':traceback.format_exc(),'source_sha256':sha(Path(__file__))});raise
 finally:atom(D/'ACTUAL_RETURN_v1.json',{'returncode':child.returncode if 'child'in locals() else None,'error':problem,'elapsed_seconds':time.monotonic()-start,'source_sha256':sha(Path(__file__)),'approval_sha256':sha(ap)})

def main():
 p=argparse.ArgumentParser();p.add_argument('mode',choices=['freeze','run']);p.add_argument('--expected-source-sha256');p.add_argument('--approval');p.add_argument('--approval-sha256');a=p.parse_args();need(Path('/proc').exists() and D.exists(),'ExistingAutoDLonly')
 if a.mode=='freeze':freeze(a)
 else:run(a)
if __name__=='__main__':main()
