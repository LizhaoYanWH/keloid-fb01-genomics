#!/usr/bin/env python3
"""One locked engineering audit of existing mapping; never rejoins GWAS or fits."""
import collections,csv,fcntl,gc,gzip,hashlib,json,math,os,platform,resource,shutil,stat,subprocess,sys,time
from pathlib import Path
sys.dont_write_bytecode=True
ROOT=Path('/root/five_new_studies_20261020')
P=ROOT/'projects/N02_keloid_genetic_cell_state'
S=P/'engineering/MAGMA_inputs_v1'
M=S/'whole_EUR_mapping_v1'
R=P/'reviews/MAGMA_native_association_v1'
GIB=1024**3
NROWS=9759875
ADMITTED=8862187

def need(ok,why):
 if not ok:raise RuntimeError(why)

def guard():
 free={str(p):shutil.disk_usage(p).free for p in (Path('/'),Path('/root/autodl-tmp'))}
 need(min(free.values())>=5*GIB,'Both disks require5GiB')
 h=int(next(x.split()[1]for x in Path('/proc/meminfo').read_text().splitlines()if x.startswith('MemAvailable:')))*1024
 c=Path('/sys/fs/cgroup');lim=(c/'memory.max').read_text().strip();current=int((c/'memory.current').read_text())
 available=h if lim=='max' else min(h,int(lim)-current)
 rss=int(next(x.split()[1]for x in Path('/proc/self/status').read_text().splitlines()if x.startswith('VmRSS:')))*1024
 need(available>=4*GIB and rss<8*GIB,'Effective memory4GiB / own RSS8GiB guard')
 return dict(disk_free_bytes=free,effective_available_bytes=available,RSS_bytes=rss,cgroup_current_bytes=current)

def identity(p):
 s=p.stat();a=dict(bytes=s.st_size,mtime_ns=s.st_mtime_ns,device=s.st_dev,inode=s.st_ino)
 md=hashlib.md5();hs=hashlib.sha256()
 with p.open('rb')as f:
  for block in iter(lambda:f.read(4*1024**2),b''):md.update(block);hs.update(block);guard()
 need(p.stat().st_size==a['bytes'] and p.stat().st_mtime_ns==a['mtime_ns'],'Input changed while hashing')
 return dict(path=str(p),**a,MD5=md.hexdigest(),SHA256=hs.hexdigest())

def js(p,x):
 tmp=p.with_name(p.name+'.tmp')
 with tmp.open('x')as f:json.dump(x,f,ensure_ascii=False,indent=2);f.write('\n');f.flush();os.fsync(f.fileno())
 need(not p.exists(),'Existing metadata retained')
 os.link(tmp,p);tmp.unlink()

def table(p,cols,data):
 with p.open('x',newline='')as f:
  w=csv.DictWriter(f,fieldnames=cols,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)

def main():
 need(platform.system()=='Linux' and os.getenv('FIVE_STUDIES_COMPUTE')=='autodl','Existing AutoDL only')
 need(Path(__file__).resolve()==R/'audit_and_freeze_inputs_v1.py','Exact isolated source required')
 os.sched_setaffinity(0,set(sorted(os.sched_getaffinity(0))[-2:]));os.nice(10);os.umask(0o077)
 lock=(ROOT/'metadata/assess01_CPU_worker.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 out=R/'input_audit_v1';need(not out.exists(),'One audit only; preserve prior attempt')
 out.mkdir();started=time.monotonic()
 status=dict(status='RUNNING_ENGINEERING_AUDIT_NO_ASSOCIATION',PID=os.getpid(),start_ticks=Path('/proc/self/stat').read_text().split()[21],
  UTC_start=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),source_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),admission=guard())
 js(out/'RUNNING.json',status)
 try:
  official=json.loads((S/'COMPLETE.json').read_text());mapping=json.loads((M/'COMPLETE.json').read_text())
  need(mapping['mapping_status_counts']=={'ADMITTED_UNDIRECTED_P':ADMITTED,'NO_REFERENCE_COORDINATE_ALLELE_MATCH':896968,'AMBIGUOUS_SOURCE_DUPLICATE_ESTIMATES':720},'Original mapping family differs')
  need(official['status']=='COMPLETE_OFFICIAL_INPUTS_AND_TOOLS_UNPACKED_NO_FIT','Official input stage incomplete')
  paths=[P/'inputs/gwas/GCST90652488.tsv.gz',S/'COMPLETE.json',M/'COMPLETE.json',P/'metadata/EUR_GWAS_download_recovery_v2.json',
   P/'metadata/coloc_production_launcher_v1/COMPLETE.json',P/'results/coloc_complete_family_ABF_v1/ROOT_ACCEPTANCE.json',
   P/'scripts/prepare_complete_EUR_MAGMA_mapping_v1.py',S/'全EUR输入与映射冻结协议_v1.md']
  for item in official['items']:
   paths.append(Path(item['file']))
   dest=S/(item['kind']+'_verified_unpacked')
   paths += [dest/x['member']for x in item['members']]
  paths += [M/n for n in mapping['output_files']]
  files={str(p):identity(p)for p in paths}
  gw=files[str(paths[0])];need(gw['MD5']=='5eddd7bdd0d98cf7b49615c1df7ea0bd'and gw['SHA256']==mapping['source_full_GWAS_SHA256'],'GWAS official MD5/SHA mismatch')
  for item in official['items']:
   a=files[item['file']];need(a['MD5']==item['MD5'] and a['SHA256']==item['actual_SHA256'],'Official archive mismatch')
   for x in item['members']:
    a=files[str(S/(item['kind']+'_verified_unpacked')/x['member'])]
    need(a['bytes']==x['bytes']and a['SHA256']==x['SHA256'],'Unpacked member mismatch')
  for name,x in mapping['output_files'].items():
   a=files[str(M/name)];need(a['bytes']==x['bytes']and a['SHA256']==x['SHA256'],'Committed mapping output mismatch')
  # Independent full BIM duplicate/chr/allele-format metadata pass, not a GWAS join.
  seen=set();dup=set();chr_ref=collections.Counter();bim_alleles=collections.Counter();nr=0
  for line in (S/'EUR_LD_verified_unpacked/g1000_eur.bim').open():
   z=line.split();need(len(z)==6,'BIM schema');nr+=1;chr_ref[z[0]]+=1
   token=int(z[1][2:])if z[1].startswith('rs')and z[1][2:].isdigit()else z[1]
   if token in seen:dup.add(z[1])
   seen.add(token);bim_alleles['SNV'if len(z[4])==len(z[5])==1 else 'indel_or_other_sequence']+=1
   if nr%250000==0:guard()
  need(nr==mapping['reference_rows']and not dup,'Reference duplicate IDs changed')
  del seen;gc.collect()
  fam=[x.split()for x in (S/'EUR_LD_verified_unpacked/g1000_eur.fam').read_text().splitlines()if x.strip()]
  need(all(len(x)==6 for x in fam),'FAM schema')
  need(len({(x[0],x[1])for x in fam})==len(fam),'Duplicate FAM samples')
  aliases={};alias_duplicates=0;groups=0
  for line in (S/'EUR_LD_verified_unpacked/g1000_eur.synonyms').open():
   if line.startswith('#')or not line.strip():continue
   z=[int(x)for x in line.split()];groups+=1
   for x in z:
    if x in aliases and aliases[x]!=z[0]:alias_duplicates+=1
    aliases[x]=z[0]
   if groups%250000==0:guard()
  def number(x):return int(x[2:])if x.startswith('rs')and x[2:].isdigit()else None
  categories=collections.Counter();cross=collections.Counter();duphist=collections.Counter();rsclass=collections.Counter();seen_source=set();dupsource=set();seen_ref=set()
  bits=bytearray(NROWS+1);nr=0;id_hash=hashlib.sha256();loc_hash=hashlib.sha256()
  examples=collections.defaultdict(list)
  with gzip.open(M/'all_original_rows_mapping_audit.tsv.gz','rt',newline='')as f:
   for row in csv.DictReader(f,delimiter='\t'):
    nr+=1;i=int(row['original_source_row']);need(1<=i<=NROWS and not bits[i],'Duplicate/missing original row identity');bits[i]=1
    src=row['source_rsid'];ref=row['LD_reference_ID'];status0=row['status'];categories[status0]+=1
    source_token=number(src);source_key=source_token if source_token is not None else src
    if source_key in seen_source:dupsource.add(src)
    seen_source.add(source_key)
    duphist[row['exact_duplicate_source_rows']]+=1
    if src in ('','NA','.'):relation='SOURCE_RSID_MISSING'
    elif ref=='NA':relation='NO_REFERENCE_ID_ASSIGNED'
    elif src==ref:relation='IDENTICAL_SOURCE_AND_LD_ID'
    elif source_token is not None and number(ref)is not None and source_token in aliases and number(ref)in aliases and aliases[source_token]==aliases[number(ref)]:relation='OFFICIAL_DBSNP151_SYNONYM'
    else:relation='DISTINCT_ID_NOT_CONFIRMED_BY_OFFICIAL_SYNONYMS'
    rsclass[(status0,relation)]+=1
    cross[(row['chromosome_GRCh37'],status0,row['mapped_strand_relation'],row['palindromic_SNV'],row['indel_sequence'],relation)]+=1
    if len(examples[relation])<8:examples[relation].append(row)
    if status0=='ADMITTED_UNDIRECTED_P':
     rt=number(ref);rt=rt if rt is not None else ref
     need(rt not in seen_ref,'Duplicate admitted LD ID');seen_ref.add(rt)
     id_hash.update((ref+'\n').encode());loc_hash.update((ref+'\t'+row['chromosome_GRCh37']+'\t'+str(int(row['position_GRCh37']))+'\n').encode())
    if nr%100000==0:guard()
  need(nr==NROWS and sum(bits)==NROWS and dict(categories)==mapping['mapping_status_counts'],'Incomplete original-row audit')
  del bits,seen_source,seen_ref,aliases;gc.collect()
  # Compare pre-existing P/location order and all eligible counts, never recompute mapping.
  ids=hashlib.sha256();locs=hashlib.sha256();np=0
  with (M/'complete_mapped_EUR_SNP_P.tsv').open()as f,(M/'complete_mapped_EUR_SNP_location_GRCh37.tsv').open()as h:
   pr=csv.DictReader(f,delimiter='\t');lr=csv.DictReader(h,delimiter='\t')
   import itertools
   for pv,lo in itertools.zip_longest(pr,lr):
    need(pv is not None and lo is not None and pv['SNP']==lo['SNP'],'P/location axis differs');np+=1
    need(math.isfinite(float(pv['P']))and 0<float(pv['P'])<=1,'Invalid P in original admitted input')
    ids.update((pv['SNP']+'\n').encode());locs.update((lo['SNP']+'\t'+lo['CHR']+'\t'+str(int(lo['BP']))+'\n').encode())
    if np%250000==0:guard()
  need(np==ADMITTED and ids.hexdigest()==id_hash.hexdigest()and locs.hexdigest()==loc_hash.hexdigest(),'Admitted IDs/location differ from audit')
  cov=list(csv.DictReader((M/'complete_official_gene_window10kb_coverage_before_association.tsv').open(),delimiter='\t'))
  official_genes=[x.split()for x in (S/'gene37_verified_unpacked/NCBI37.3.gene.loc').read_text().splitlines()if x.strip()]
  need(len(cov)==len(official_genes)==19427 and len({r['official_gene_id']for r in cov})==19427,'Official gene family differs')
  genes=[];gene_status=collections.Counter();assignment_tokens=0
  for row,z in zip(cov,official_genes):
   need([row[k]for k in ('official_gene_id','chromosome','official_start','official_stop')]==z[:4],'Gene coordinate/order identity differs')
   n=int(row['n_unique_admitted_SNP']);assignment_tokens+=n
   reason='UNSUPPORTED_Y_NO_NATIVE_ANALYSIS'if row['chromosome']in('Y','24')else 'NO_ADMITTED_SNP_FIXED_WINDOW'if n==0 else 'ELIGIBLE_FOR_NATIVE_MODEL_INTERNAL_QC_PENDING'
   gene_status[reason]+=1;genes.append({**row,'gene_symbol':z[5],'native_preanalysis_status':reason})
  table(out/'all19427_gene_admission_ledger.tsv',list(genes[0]),genes)
  cr=[dict(zip(('chromosome','mapping_status','strand_relation','palindromic_SNV','indel_sequence','rsID_relation'),k),n_rows=v)for k,v in sorted(cross.items())]
  table(out/'all_chr_mapping_and_rsid_diagnostics.tsv',list(cr[0]),cr)
  js(out/'bounded_examples_by_rsid_relation.json',dict(examples))
  # Exact official binary is made executable in its isolated tool directory; bytes unchanged.
  binary=S/'MAGMALinux1.10_verified_unpacked/magma';mode=stat.S_IMODE(binary.stat().st_mode)
  binary.chmod(mode|stat.S_IXUSR)
  fx=out/'native_annotation_synthetic_fixture';fx.mkdir()
  (fx/'snps.loc').write_text('SNP CHR BP\ns1 1 95000\ns2 1 200000\n')
  (fx/'genes.loc').write_text('SYNTH_HIT 1 100000 101000 + SYNTH_HIT\nSYNTH_ZERO 1 1000000 1001000 + SYNTH_ZERO\nSYNTH_X X 100000 101000 + SYNTH_X\n')
  version=subprocess.run([str(binary),'--version'],cwd=fx,text=True,capture_output=True)
  cmd=[str(binary),'--annotate','window=10,10','--snp-loc',str(fx/'snps.loc'),'--gene-loc',str(fx/'genes.loc'),'--out',str(fx/'tiny')]
  rr=subprocess.run(cmd,cwd=fx,text=True,capture_output=True)
  (fx/'version.stdout_stderr.txt').write_text(version.stdout+version.stderr)
  (fx/'annotate.stdout_stderr.txt').write_text(rr.stdout+rr.stderr)
  need(rr.returncode==0 and (fx/'tiny.genes.annot').exists(),'Native synthetic annotate failed; no statistical execution')
  js(fx/'tool_probe_receipt.json',dict(version_returncode=version.returncode,version_stdout=version.stdout,version_stderr=version.stderr,
   original_mode=oct(mode),new_mode=oct(stat.S_IMODE(binary.stat().st_mode)),binary_SHA256_unchanged=identity(binary)['SHA256']==files[str(binary)]['SHA256'],
   annotation_argv=cmd,annotation_returncode=rr.returncode,synthetic_only=True,trait_annotation_or_association=False))
  for name,a in files.items():
   f=Path(name);need(f.stat().st_size==a['bytes']and f.stat().st_mtime_ns==a['mtime_ns'],'Input changed during audit')
  result={**status,'status':'COMPLETE_ENGINEERING_AUDIT_FROZEN_INPUTS_NO_ASSOCIATION','files':list(files.values()),
   'source_rows':nr,'status_counts':dict(categories),'reference_rows':mapping['reference_rows'],'reference_duplicate_ID_count':len(dup),
   'reference_chr_counts':dict(chr_ref),'reference_allele_class_counts':dict(bim_alleles),'FAM_individuals':len(fam),'FAM_sex_counts':dict(collections.Counter(x[4]for x in fam)),
   'official_synonym_groups':groups,'official_synonym_ID_conflicts':alias_duplicates,'source_rsid_duplicate_ID_count':len(dupsource),
   'source_exact_duplicate_size_histogram':dict(duphist),'rsID_relations':[dict(status=k[0],relation=k[1],n=v)for k,v in sorted(rsclass.items())],
   'admitted_P_rows':np,'admitted_ID_order_SHA256':ids.hexdigest(),'admitted_location_order_SHA256':locs.hexdigest(),
   'official_gene_family':len(genes),'gene_admission_counts':dict(gene_status),'total_gene_SNP_annotation_assignments':assignment_tokens,
   'constant_total_N':1282582,'window_kb_each_side':10,'native_model':'snp-wise=mean','seed':20261020,'reference_synonyms_native':'official automatic default synonym-dup=skip',
   'MAF_note':'Source audit has no per-SNP MAF; no MAF filter or invented per-SNP MAF audit. Native reference QC remains default and logs retained.',
   'no_new_GWAS_join':True,'input_family_unchanged':True,'MAGMA_association_run':False,'scDRS_run':False,'native_annotation_fixture_only':True,
   'max_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,'elapsed_seconds':time.monotonic()-started,'final_resources':guard()}
  js(out/'frozen_inputs_and_audit_v1.json',result)
  print(json.dumps({k:result[k]for k in ('status','source_rows','status_counts','rsID_relations','FAM_individuals','gene_admission_counts','total_gene_SNP_annotation_assignments','max_RSS_bytes','elapsed_seconds')},indent=2),flush=True)
 except BaseException as e:
  js(out/'FAILED_PRESERVED.json',dict(error=type(e).__name__,reason=str(e),PID=os.getpid()));raise

if __name__=='__main__':main()
