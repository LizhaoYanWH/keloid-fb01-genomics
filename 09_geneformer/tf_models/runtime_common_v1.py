"""Shared finite perturbation runtime; existing AutoDL only."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import json,hashlib,time,shutil,resource,sys
from pathlib import Path
P=Path('/root/five_new_studies_20261020/projects/N02_keloid_genetic_cell_state')
OLD=Path('/root/autodl-tmp/Five_Studies_Execution_20261020')
L=Path('/root/autodl-tmp/Wound_Target_Finetuning_Execution_20261019')
D=P/'mechanisms_20261021/TF_followup_v1';BASE=P/'mechanisms_20261021/perturbation_multimodel_v1';E=P/'extensions_20261021/virtual_KO_v1'
MODEL=Path('/root/autodl-tmp/NRP2_virtual_knockout_20261016/models/Geneformer-V2-104M')
GIB=1024**3;PROTEINS=['RELA','NFKB1','NFKB2','REL','RELB','JUN','JUNB','JUND','FOS','FOSB','FOSL1','FOSL2','SMAD2','SMAD3','SMAD4'];CONTROLS=['DYM','SLC25A39','UTP11','QARS1'];PANEL=PROTEINS+CONTROLS

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for z in iter(lambda:f.read(8388608),b''):h.update(z)
 return h.hexdigest()
def atom(p,j):
 p=Path(p);assert not p.exists(),str(p)
 t=Path(str(p)+'.tmp');t.write_text(json.dumps(j,indent=2,allow_nan=False)+'\n');t.replace(p)
def effective_free():
 vals=[int(next(s for s in Path('/proc/meminfo').read_text().splitlines() if s.startswith('MemAvailable:')).split()[1])*1024]
 for b,a,c in [('/sys/fs/cgroup','memory.max','memory.current'),('/sys/fs/cgroup/memory','memory.limit_in_bytes','memory.usage_in_bytes')]:
  try:
   cap=(Path(b)/a).read_text().strip()
   if cap!='max':vals.append(int(cap)-int((Path(b)/c).read_text()))
  except FileNotFoundError:pass
 return min(vals)
def own_rss():
 return int(next(s for s in Path('/proc/self/status').read_text().splitlines() if s.startswith('VmRSS:')).split()[1])*1024

def guard(grant,stage,start=False):
 budget=6 if stage=='heat' else 16
 peer=int(grant['peer_unmaterialized_GiB'][stage])
 assert effective_free()>=(8+peer+(budget if start else 0))*GIB,'RAM_HOLD_FROZEN_BUDGET'
 assert own_rss()<=budget*GIB,'OWN_RSS_LIMIT'
 assert shutil.disk_usage('/root/autodl-tmp').free>=5*GIB,'DATA_FIVE_GIB_FLOOR'
 assert shutil.disk_usage('/root').free>=5*GIB,'SYSTEM_FIVE_GIB_FLOOR'
 used=sum(p.stat().st_size for p in D.rglob('*') if p.is_file())
 assert used<128*1024**2,'OWN_128MIB_OUTPUT_CAP'

def approved(stage):
 assert sys.platform=='linux' and os.getenv('FIVE_STUDIES_COMPUTE')=='autodl','EXISTING_AUTODL_ONLY'
 f=json.loads((D/'source_freeze_v1.json').read_text());g=json.loads((D/'ROOT_APPROVAL_TO_RUN.json').read_text())
 assert g['approved'] is True and stage in g['stages'] and sha(D/'source_freeze_v1.json')==g['freeze_sha256']
 for n,h in g['source_sha256'].items():assert sha(D/n)==h,n
 guard(g,stage,True)
 for x in f['files']:assert Path(x['path']).stat().st_size==x['bytes'] and sha(x['path'])==x['sha256'],x['path']
 os.sched_setaffinity(0,set(sorted(os.sched_getaffinity(0))[:2]))
 return f,g

def claim(stage):
 atom(D/(stage+'_ONCE_CLAIM.json'),{'pid':os.getpid(),'boot_id':Path('/proc/sys/kernel/random/boot_id').read_text().strip(),'UTC':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())})
def file_receipt(folder):
 return {str(p.relative_to(folder)):{'bytes':p.stat().st_size,'sha256':sha(p)} for p in folder.rglob('*') if p.is_file() and p.name!='COMPLETE.json'}
