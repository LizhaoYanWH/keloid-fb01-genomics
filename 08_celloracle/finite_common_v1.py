"""Finite CellOracle IO and resource contracts. Existing environments only. Python 3.10 compatible."""
import os
for _k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[_k]='1'
import datetime,hashlib,json,platform,resource,shutil
from pathlib import Path
SEED=20261021
OFFICIAL_COMMIT='83c8f6aedf91c3406ec49145e76e9debbc3ae4d0'
PRIOR_BYTES=5520100
PRIOR_GIT_BLOB_SHA1='3ec976bc176e57ba5ab6b164e7b56b15891ea3ea'
STUDIES=('GSE220300','GSE335482')
def need(x,msg):
    if not x: raise RuntimeError(msg)
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(2**20),b''):h.update(b)
    return h.hexdigest()
def utc():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def rec(p):
    p=Path(p);return {'path':str(p),'bytes':p.stat().st_size,'sha256':sha(p)}
def checkref(r):
    p=Path(r['path']);need(p.is_file() and p.stat().st_size==r['bytes'] and sha(p)==r['sha256'],'Input identity mismatch: '+str(p));return p
def atom(p,x):
    p=Path(p);t=p.with_name(p.name+'.incoming');need(not t.exists(),'Prior interrupted write retained: '+str(t));t.write_text(json.dumps(x,indent=2,ensure_ascii=False,allow_nan=False)+'\n');t.replace(p)
def load(p):return json.loads(Path(p).read_text())
def grant(path,status,script,extra=None):
    g=load(path);need(g.get('status')==status,'Root stage grant required');need(g['host']==platform.node(),'Host ownership mismatch');need(g['source_sha256']==sha(script),'Root source binding changed');need(g['common_sha256']==sha(__file__),'Common source binding changed')
    for k,v in (extra or {}).items():need(g.get(k)==v,'Root binding mismatch: '+k)
    return g
def guard(own_GiB=4,remaining_GiB=8,output_reserve_MiB=256,start=True):
    need(platform.system()=='Linux' and os.getenv('FIVE_STUDIES_COMPUTE')=='autodl','Scientific/real input work AutoDL only')
    lim=Path('/sys/fs/cgroup/memory.max').read_text().strip();cur=int(Path('/sys/fs/cgroup/memory.current').read_text());host={x.split(':')[0]:int(x.split()[1])*1024 for x in Path('/proc/meminfo').read_text().splitlines()};free=host['MemAvailable']
    if lim!='max':free=min(free,int(lim)-cur)
    rss=int(next(x.split()[1] for x in Path('/proc/self/status').read_text().splitlines() if x.startswith('VmRSS:')))*1024
    need(free >= (own_GiB+remaining_GiB if start else remaining_GiB)*2**30,'Actual memory reserve unmet')
    if not start:need(rss<=own_GiB*2**30,'Actual own RSS exceeds approved bound')
    need(shutil.disk_usage('/root').free>=5*2**30+output_reserve_MiB*2**20,'System 5GiB plus own output reserve unmet');need(shutil.disk_usage('/root/autodl-tmp').free>=5*2**30,'Data 5GiB floor unmet')
    return {'actual_available_bytes':free,'self_RSS_bytes':rss,'cgroup_current_bytes':cur,'memory_limit':lim,'system_free_bytes':shutil.disk_usage('/root').free,'data_free_bytes':shutil.disk_usage('/root/autodl-tmp').free}
def harray(n):
    import numpy as np,h5py
    if isinstance(n,h5py.Dataset):
        a=n[:];return np.array([x.decode() if isinstance(x,bytes) else str(x) for x in a]) if a.dtype.kind in 'OSU' else a
    cats=harray(n['categories']);return np.array([cats[i] if i>=0 else '' for i in n['codes'][:]])
def hindex(h,key):return harray(h[key][h[key].attrs.get('_index','_index')])
def tf_catalog(spec):
    import pandas as pd
    p=checkref(spec)
    if spec['format']=='one_symbol_per_line':
        vals=p.read_text().splitlines();need(all(x.strip()==x and '\t' not in x and ' ' not in x for x in vals if x),'Catalog is not explicit one-symbol-per-line format')
    elif spec['format']=='tsv':
        d=pd.read_csv(p,sep='\t',dtype=str,keep_default_na=False);need(spec['symbol_column'] in d,'Explicit TF catalog column missing');vals=d[spec['symbol_column']].tolist()
    else:raise RuntimeError('Unsupported TF catalog format; root must inspect real schema')
    return set(x for x in vals if x)
def prior_frame(spec):
    import pandas as pd,numpy as np
    p=checkref(spec);need(p.stat().st_size==PRIOR_BYTES,'Unexpected pinned prior bytes');payload=p.read_bytes();need(hashlib.sha1(b'blob '+str(len(payload)).encode()+b'\0'+payload).hexdigest()==PRIOR_GIT_BLOB_SHA1,'Pinned publisher Git blob identity mismatch');d=pd.read_parquet(p);need({'peak_id','gene_short_name'}<=set(d.columns),'Actual promoter prior schema is incompatible');need(not d.columns.duplicated().any(),'Duplicated prior columns')
    for n in d.columns.difference(['peak_id','gene_short_name']):need(pd.api.types.is_numeric_dtype(d[n].dtype) or pd.api.types.is_bool_dtype(d[n].dtype),'Non-numeric/non-boolean promoter TF column: '+n)
    return d

def check_plan(plan_path):
    p=load(plan_path);need(p['status']=='ROOT_FROZEN_CELLORACLE_CANDIDATE_PLAN','Candidate scientific selection must be frozen by root');need(len(p['chains'])>=1,'At least one actual candidate chain required');need(len({x['chain_id'] for x in p['chains']})==len(p['chains']),'Duplicate chain IDs');tfs=[x['TF'] for x in p['chains']];need(1<=len(set(tfs))<=3,'Exactly 1-3 distinct finite TFs; overlapping chains do not create duplicate KO fits');need(all(x.get('selection_reason') and x.get('specific_question') for x in p['chains']),'Specific non-forced question/evidence required')
    f=load(checkref(p['communication_freeze']));need(set(f['cohorts'])==set(STUDIES),'Unexpected frozen studies')
    comp=load(checkref(p['S2_complete']));need(comp['status']=='COMPLETE_S1_S2_TECHNICAL_ROOT_SCIENTIFIC_REVIEW_PENDING' and comp['freeze_sha256']==p['communication_freeze']['sha256'],'Actual complete S2 binding required')
    members=checkref(p['S2_response_members']);need(any(r['path']==str(members) and r['sha256']==p['S2_response_members']['sha256'] and r['bytes']==p['S2_response_members']['bytes'] for r in comp['outputs']),'Response membership not bound to S2 producer')
    for r in p['supporting_evidence']:checkref(r)
    rt=load(checkref(p['runtime_adoption']));need(rt['status']=='ROOT_ACCEPTED_CELLORACLE_RUNTIME_AND_MINI_FIXTURE','Actual import/API/prior/small-fixture adoption required before building production input')
    need(rt['celloracle_version']=='0.20.0','Different CellOracle version needs new review')
    pc=load(checkref(p['profile_complete']));need(pc['status']=='COMPLETE_ALL_PATIENT_PROFILES_NO_S1_EFFECTS' and pc['freeze_sha256']==p['communication_freeze']['sha256'],'Frozen profile producer mismatch')
    reg=load(checkref(p['regulatory_freeze']));need(p['regulatory_freeze']['sha256']==next(x['sha256'] for x in f['inputs'] if x['path']==p['regulatory_freeze']['path']),'Regulatory program freeze is not original communication input')
    return p,f,reg

def files_receipt(directory):return [rec(p) for p in sorted(Path(directory).rglob('*')) if p.is_file() and p.name not in ('COMPLETE.json','COMPLETE.json.incoming')]
