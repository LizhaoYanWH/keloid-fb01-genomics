"""Conditional S3/S4 evidence audit from COMPLETE profiles and frozen priors.
No new clustering, normalization, differential model, TF fit, pathway search or KO.
Prepare-only source: requires a root-approved actual candidate plan before AutoDL run.
"""
from pathlib import Path
import argparse,collections,csv,gzip,hashlib,json,math,os,resource,shutil,time
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS'):os.environ[k]='1'

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb')as f:
  for z in iter(lambda:f.read(2**20),b''):h.update(z)
 return h.hexdigest()
def verify(r):
 p=Path(r['path']);need(p.is_file() and p.stat().st_size==r['bytes'] and sha(p)==r['sha256'],'Frozen input changed: '+str(p))
def atom(p,j):
 q=p.with_name(p.name+'.tmp');q.write_text(json.dumps(j,indent=2,allow_nan=False,ensure_ascii=False)+'\n');q.replace(p)
def tsread(p):
 with (gzip.open(p,'rt')if str(p).endswith('.gz')else Path(p).open())as f:return list(csv.DictReader(f,delimiter='\t'))
def writetsv(path,rows,fields):
 with gzip.open(path,'wt')as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='raise',lineterminator='\n');w.writeheader();w.writerows({k:('NA'if v is None else v)for k,v in r.items()}for r in rows)
def truth(x):return str(x).lower()in ('true','1')
def resources():
 m={x.split(':')[0]:int(x.split()[1])*1024 for x in Path('/proc/meminfo').read_text().splitlines()};v=m['MemAvailable'];lim=Path('/sys/fs/cgroup/memory.max').read_text().strip()
 if lim!='max':v=min(v,int(lim)-int(Path('/sys/fs/cgroup/memory.current').read_text()))
 rss=int(next(x.split()[1]for x in Path('/proc/self/status').read_text().splitlines()if x.startswith('VmRSS:')))*1024
 return {'effective_free_bytes':v,'RSS_bytes':rss,'system_free_bytes':shutil.disk_usage('/root').free,'data_free_bytes':shutil.disk_usage('/root/autodl-tmp').free}
def guard(out,start=False):
 r=resources();need(r['effective_free_bytes']>=256*2**20+(int(1.25*2**30)if start else 0),'Own1.25GiB plus runtime256MiB floor');need(r['RSS_bytes']<=1.25*2**30,'Owned RSS1.25GiB bound')
 need(r['system_free_bytes']>=5*2**30+64*2**20 and r['data_free_bytes']>=5*2**30,'Both5GiB floors plus output reserve')
 if out.exists():need(sum(p.stat().st_size for p in out.rglob('*')if p.is_file())<64*2**20,'Owned64MiB output bound')
 return r

def main():
 pa=argparse.ArgumentParser();pa.add_argument('--plan',required=True);pa.add_argument('--approval',required=True);pa.add_argument('--approval-sha256',required=True);a=pa.parse_args()
 need(os.getenv('FIVE_STUDIES_COMPUTE')=='autodl' and Path('/proc/self/status').exists(),'AutoDL only')
 need(sha(a.approval)==a.approval_sha256,'Root grant changed');g=json.loads(Path(a.approval).read_text());need(g['action']=='AUDIT_FROZEN_CANDIDATE_COMPLEX_AND_TF_PRIOR_EDGES','Wrong grant')
 need(g['source_sha256']==sha(__file__) and g['plan_sha256']==sha(a.plan),'Source/plan not root-bound');need(g['CPU']==1 and g['own_RSS_GiB']==1.25 and g['wall_seconds']<=900,'Bounded1CPU audit grant')
 plan=json.loads(Path(a.plan).read_text());need(not plan.get('schema_only',False),'No execution of template');need(plan['scientific_fits_permitted'] is False,'No new fitting');need(0<len(plan['chains'])<=3,'At most3 chains')
 out=Path(plan['output']).resolve();base=Path('/root/five_new_studies_20261020/projects/N02_keloid_genetic_cell_state/communication_mechanism_20261021_v1/reuse_bridge').resolve();need(out.is_relative_to(base)and not out.exists(),'Fresh owned output only');guard(out,start=True)
 os.sched_setaffinity(0,[min(os.sched_getaffinity(0))]);resource.setrlimit(resource.RLIMIT_AS,(int(1.5*2**30),int(1.5*2**30)))
 for r in plan['inputs']:verify(r)
 by={r['key']:r for r in plan['inputs']};need(len(by)==len(plan['inputs']),'Duplicate input keys');need(all(k in by for k in ['LIANA_consensus','canonical_CollecTRI','profile_COMPLETE','S2_COMPLETE','S1_S2_FREEZE']),'Exact priors and actual producer records required')
 frozen=json.loads(Path(by['S1_S2_FREEZE']['path']).read_text());need(plan['profile_min_cells']==frozen['profile_min_cells'] and plan['detection_fraction']==frozen['detection_fraction'],'Original expression-support gates changed')
 for st in plan['studies']:need(sorted(plan['known_patients'][st])==sorted(frozen['cohorts'][st]['patients']),'Original known-patient membership changed')
 li=tsread(by['LIANA_consensus']['path']);need(set(li[0])=={'ligand','receptor'},'Observed LIANA v1.10.0 consensus ligand/receptor schema required')
 net=tsread(by['canonical_CollecTRI']['path']);need(all(c in net[0]for c in ['source','target','weight','resources','references','sign_decision','source_original','target_original']),'Frozen canonical CollecTRI ledger required')
 # Existing ledger preserves each source edge; do not recompute TF footprint or resolve signs by desired outcome.
 need(plan['complex_subunit_separator']=='_' and plan['complex_annotation_scope']=='LISTED_RESOURCE_ENTITY_NOT_ALL_BIOLOGICALLY_REQUIRED_SUBUNITS','Explicit LIANA entity semantics required')
 import numpy as np
 profiles={};profile_meta={};axes={};prod=json.loads(Path(by['profile_COMPLETE']['path']).read_text())
 need(prod['status']=='COMPLETE_ALL_PATIENT_PROFILES_NO_S1_EFFECTS','Profiles not true complete')
 for st in plan['studies']:
  found=[x for x in prod['studies']if x['study']==st];need(len(found)==1,'Study producer missing')
  for r in found[0]['files']:
   if Path(r['path']).name in ['all_profile_expression.npz','profiles.tsv']:verify(r)
  d=Path(by['profile_COMPLETE']['path']).parent/st;meta=tsread(d/'profiles.tsv');z=np.load(d/'all_profile_expression.npz',allow_pickle=False)
  ids=z['profile_ids'].tolist();genes=z['gene_symbols'].tolist();need(ids==[r['profile_id']for r in meta],'Profile axis/order changed');need(len(set(genes))==len(genes),'Duplicate gene symbols')
  need(z['raw_counts'].shape==z['detected_cells'].shape==z['mean_log_CP10K'].shape==(len(meta),len(genes)),'Profile matrix shape changed')
  profiles[st]={k:z[k]for k in ['raw_counts','detected_cells','mean_log_CP10K']};profile_meta[st]={r['profile_id']:dict(r,index=i)for i,r in enumerate(meta)};axes[st]={s:i for i,s in enumerate(genes)}
 def expr(st,patient,group,symbol):
  m=profile_meta[st].get(patient+'::'+group);j=axes[st].get(symbol)
  if m is None:return {'status':'PROFILE_ABSENT_NA','n_cells':0,'assayed':j is not None,'detection':None,'meanlog':None,'raw_count':None,'pass':False}
  if j is None:return {'status':'GENE_NOT_ASSAYED_NA','n_cells':int(m['n_cells']),'assayed':False,'detection':None,'meanlog':None,'raw_count':None,'pass':False}
  i=m['index'];n=int(m['n_cells']);a=profiles[st];det=float(a['detected_cells'][i,j])/n;val=float(a['mean_log_CP10K'][i,j]);need(math.isfinite(det)and math.isfinite(val),'Nonfinite original profile')
  return {'status':'MEASURED','n_cells':n,'assayed':True,'detection':det,'meanlog':val,'raw_count':int(a['raw_counts'][i,j]),'pass':n>=plan['profile_min_cells'] and det>=plan['detection_fraction']}
 # All exact original edge rows are retained. Tokens are resource-defined subunits, not aliases or causal proof.
 parts=lambda x:x.split('_')
 sub=[];entities=[];sender_delta=[];tfedges=[];missing=[];start=time.monotonic();out.mkdir();atom(out/'RUNNING.json',{'PID':os.getpid(),'source_sha256':sha(__file__),'plan_sha256':sha(a.plan)})
 try:
  need(len({c['chain_id']for c in plan['chains']})==len(plan['chains']),'Duplicate chains')
  for chain in plan['chains']:
   guard(out);need(time.monotonic()-start<g['wall_seconds'],'Audit wall bound');cid=chain['chain_id'];need(chain['source_S2_COMPLETE_sha256']==by['S2_COMPLETE']['sha256'] and chain['actual_selection_reason'],'Real selected-chain producer binding required');need(set(chain['studies'])<=set(plan['studies']),'Unadmitted study')
   need(chain['direction']in ['FB01_input','FB01_output'],'Unknown estimand');L,R=chain['ligand'],chain['receptor'];send,recv=chain['sender'],chain['receiver']
   need(recv=='FB01'if chain['direction']=='FB01_input'else send=='FB01','FB01 direction differs')
   matches=[dict(x,row_index=i+1)for i,x in enumerate(li)if(L==x['ligand']or L in parts(x['ligand']))and(R==x['receptor']or R in parts(x['receptor']))]
   if not matches:missing.append({'chain_id':cid,'study':'ALL','kind':'LIANA_ENTITY','status':'NO_MATCHING_LIANA_ENTITY_COMPLEX_UNRESOLVED_NOT_NEGATIVE'})
   # Original monomer support remains available even when LIANA annotation is missing.
   for st in chain['studies']:
    for patient in plan['known_patients'][st]:
     original=[('ligand',L,send),('receptor',R,recv)]
     for role,symbol,group in original:
      e=expr(st,patient,group,symbol);sub.append(dict(chain_id=cid,study=st,patient=patient,resource_row='ORIGINAL_S2_MONOMER',entity=symbol,role=role,group=group,symbol=symbol,**e))
     for hit in matches:
      vals={}
      for role,col,group in [('ligand','ligand',send),('receptor','receptor',recv)]:
       vv=[]
       for symbol in parts(hit[col]):
        e=expr(st,patient,group,symbol);vv.append(e);sub.append(dict(chain_id=cid,study=st,patient=patient,resource_row=str(hit['row_index']),entity=hit[col],role=role,group=group,symbol=symbol,**e))
       vals[role]=vv
      allv=vals['ligand']+vals['receptor'];entities.append(dict(chain_id=cid,study=st,patient=patient,LIANA_row=hit['row_index'],ligand_entity=hit['ligand'],receptor_entity=hit['receptor'],all_listed_subunits_assayed=all(e['assayed']for e in allv),all_listed_subunits_detection_supported=all(e['pass']for e in allv),n_listed_ligand_subunits=len(vals['ligand']),n_listed_receptor_subunits=len(vals['receptor']),evidence='RESOURCE_ENTITY_PLUS_SAME_PATIENT_RNA_NOT_COMPLEX_ASSEMBLY_OR_PHYSICAL_SIGNAL',not_universal_required_components=True))
     if chain['direction']=='FB01_output':
      for symbol in sorted({L}|{s for row in matches for s in parts(row['ligand'])}):
       x,y=expr(st,patient,'FB01',symbol),expr(st,patient,'otherFB',symbol);ok=x['status']==y['status']=='MEASURED'
       sender_delta.append(dict(chain_id=cid,study=st,patient=patient,symbol=symbol,FB01_status=x['status'],otherFB_status=y['status'],FB01_meanlog=x['meanlog'],otherFB_meanlog=y['meanlog'],FB01_minus_otherFB_meanlog=x['meanlog']-y['meanlog']if ok else None,FB01_detection=x['detection'],otherFB_detection=y['detection'],receiver_response_status='NOT_ESTIMABLE_WITHOUT_INDEPENDENT_RECEIVER_CONTRAST',no_new_patient_P=True))
   if chain['direction']=='FB01_output':
    missing.append({'chain_id':cid,'study':'ALL','kind':'INDEPENDENT_RECEIVER_RESPONSE','status':'NO_INDEPENDENT_RECEIVER_CONTRAST_IN_CURRENT_POOLED_PROFILES_NO_S1_FB01_RESPONSE_REUSED_AS_OTHER_RECEIVER_RESPONSE'});continue
   need('response_members' in by,'Exact S1 response membership required for conditional TF-to-response edge audit')
   members=tsread(by['response_members']['path'])
   for st in chain['studies']:
    targets={x['gene_symbol']for x in members if x['study']==st and x['scope']=='main'and x['view']==chain['response_view']and x['sign']==chain['response_sign']}
    # All TF->response edges retained separately by study. No borrowed target or receptor->TF claim.
    for edge in net:
     if edge['target']in targets:tfedges.append(dict(chain_id=cid,study=st,response_view=chain['response_view'],response_sign=chain['response_sign'],TF=edge['source'],response_target=edge['target'],prior_weight=edge['weight'],prior_resources=edge['resources'],prior_references=edge['references'],prior_sign_decision=edge['sign_decision'],source_original=edge['source_original'],target_original=edge['target_original'],edge_type='COLLECTRI_PRIOR_TF_TO_TARGET',receptor_to_TF='NOT_ESTABLISHED_BY_THIS_EDGE_AUDIT',observed_response='MEMBER_OF_EXPLORATORY_S1_RESPONSE_NOT_INDEPENDENT_CONFIRMATION'))
    if not targets:missing.append({'chain_id':cid,'study':st,'kind':'RESPONSE','status':'NO_EXACT_SELECTED_RESPONSE_MEMBERS_NA_NO_REPLACEMENT'})
  fields=['chain_id','study','patient','resource_row','entity','role','group','symbol','status','n_cells','assayed','detection','meanlog','raw_count','pass'];writetsv(out/'all_candidate_subunit_patient_RNA.tsv.gz',sub,fields)
  fields=['chain_id','study','patient','LIANA_row','ligand_entity','receptor_entity','all_listed_subunits_assayed','all_listed_subunits_detection_supported','n_listed_ligand_subunits','n_listed_receptor_subunits','evidence','not_universal_required_components'];writetsv(out/'all_matching_LIANA_entities_support.tsv.gz',entities,fields)
  fields=['chain_id','study','patient','symbol','FB01_status','otherFB_status','FB01_meanlog','otherFB_meanlog','FB01_minus_otherFB_meanlog','FB01_detection','otherFB_detection','receiver_response_status','no_new_patient_P'];writetsv(out/'conditional_outgoing_sender_descriptions.tsv.gz',sender_delta,fields)
  fields=['chain_id','study','response_view','response_sign','TF','response_target','prior_weight','prior_resources','prior_references','prior_sign_decision','source_original','target_original','edge_type','receptor_to_TF','observed_response'];writetsv(out/'all_TF_response_prior_edges_NOT_receptor_paths.tsv.gz',tfedges,fields)
  atom(out/'all_missing_bridges_and_NA.json',missing);atom(out/'plan_snapshot.json',plan);guard(out)
  atom(out/'COMPLETE_AUDIT_NOT_CAUSAL.json',{'status':'COMPLETE_BOUNDED_COMPLEX_AND_PRIOR_EDGE_AUDIT_ROOT_REVIEW_PENDING','source_sha256':sha(__file__),'plan_sha256':sha(a.plan),'approval_sha256':sha(a.approval),'elapsed_seconds':time.monotonic()-start,'resources':resources(),'files':{p.name:{'bytes':p.stat().st_size,'sha256':sha(p)}for p in out.iterdir()if p.is_file()},'not_performed':['receptor_to_TF_path_search','new_TF_fit','new_RNA_DE','new_communication_model','new_KO','complex_protein_assembly_validation','clinical_validation']})
 except BaseException as e:atom(out/'FAILED_PRESERVED.json',{'error':repr(e),'no_automatic_retry':True});raise
if __name__=='__main__':main()
