"""Common IO and guards. New communication branch; existing results read-only."""
from pathlib import Path
import hashlib,json,os,platform,resource,shutil,time
for _k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[_k]='1'
OLD=Path('/root/autodl-tmp/Five_Studies_Execution_20261020')
P=Path('/root/five_new_studies_20261020/projects/N02_keloid_genetic_cell_state')
HERE=P/'communication_mechanism_20261021_v1'
RNA=P/'reviews/scDRS_production_preparation_v1/RNA_freeze_v1'
HGNC=OLD/'shared_inputs/programs_01_03/sources/hgnc_complete_set.txt'
PRIOR=OLD/'metadata/mechanism_primary_sources/nichenet_priors'
STUDIES=('GSE220300','GSE335482')
def need(x,m):
 if not x:raise RuntimeError(m)
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()
def receipt(p):
 p=Path(p);return dict(path=str(p),bytes=p.stat().st_size,sha256=sha(p))
def atom(p,x):
 p=Path(p);q=p.with_suffix(p.suffix+'.tmp');q.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n');q.replace(p)
def resources():
 m={l.split(':')[0]:int(l.split()[1])*1024 for l in Path('/proc/meminfo').read_text().splitlines()};f=m['MemAvailable'];lim=Path('/sys/fs/cgroup/memory.max').read_text().strip()
 if lim!='max':f=min(f,int(lim)-int(Path('/sys/fs/cgroup/memory.current').read_text()))
 return dict(effective_free_bytes=f,system_free_bytes=shutil.disk_usage('/root').free,data_free_bytes=shutil.disk_usage('/root/autodl-tmp').free)
def guard(start=True):
 need(platform.system()=='Linux' and os.getenv('FIVE_STUDIES_COMPUTE')=='autodl','Scientific run AutoDL only')
 r=resources();budget=int(1.25*2**30);floor=256*2**20
 if start:need(r['effective_free_bytes']>=budget+floor,'Own1.25GiB plus256MiB remaining required at admission')
 else:
  rss=int(next(l.split()[1] for l in Path('/proc/self/status').read_text().splitlines() if l.startswith('VmRSS:')))*1024
  need(r['effective_free_bytes']>=floor and rss<=budget,'Running guard: remaining256MiB plus ownRSS<=1.25GiB; no double reservation')
 need(r['system_free_bytes']>=5*2**30+256*2**20,'System needs5GiB plus own256MiB outputs');need(r['data_free_bytes']>=5*2**30,'Read-only data volume needs5GiB')
 return r
def authorize(stage):
 a=json.loads((HERE/'ROOT_RUN_APPROVAL.json').read_text());need(a.get('status')=='APPROVED','Root approval required')
 need(stage in a['stages'],'Stage not approved')
 for name,h in a['source_sha256'].items():need(sha(HERE/'scripts'/name)==h,'Approved source changed: '+name)
 return a
def harray(n):
 import h5py,numpy as np
 if isinstance(n,h5py.Dataset):
  a=n[:];return np.array([x.decode() if isinstance(x,bytes) else str(x) for x in a]) if a.dtype.kind in 'OSU' else a
 cats=harray(n['categories']);return np.array([cats[i] if i>=0 else '' for i in n['codes'][:]])
def hindex(h,key):return harray(h[key][h[key].attrs.get('_index','_index')])
def utc():return time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())
