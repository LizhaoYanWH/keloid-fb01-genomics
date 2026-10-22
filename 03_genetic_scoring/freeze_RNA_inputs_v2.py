"""Input-only CPU admission; full sparse counts QA/covariates, never scDRS scoring."""
from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='2'
os.sched_setaffinity(0,sorted(os.sched_getaffinity(0))[:2])
import json,hashlib,shutil,time,gzip,gc,resource
import numpy as np,pandas as pd,h5py,scipy.sparse as sp,scipy.io as sio
from threadpoolctl import threadpool_limits
ROOT=Path('/root/five_new_studies_20261020/projects/N02_keloid_genetic_cell_state');OLD=Path('/root/autodl-tmp/Five_Studies_Execution_20261020');O=ROOT/'reviews/scDRS_production_preparation_v1/RNA_freeze_v1'
def sha(f):
 h=hashlib.sha256()
 with open(f,'rb')as q:
  for b in iter(lambda:q.read(2**20),b''):h.update(b)
 return h.hexdigest()
def atom(p,j):
 q=p.with_suffix(p.suffix+'.tmp');q.write_text(json.dumps(j,indent=2)+'\n');q.replace(p)
def disk():
 assert all(shutil.disk_usage(m).free>=5*2**30 for m in ['/root','/root/autodl-tmp'])
def guard():
 disk();limit=int(Path('/sys/fs/cgroup/memory.max').read_text());cur=int(Path('/sys/fs/cgroup/memory.current').read_text());mem={r.split(':')[0]:int(r.split()[1])*1024 for r in Path('/proc/meminfo').read_text().splitlines()};free=min(limit-cur,mem['MemAvailable']);assert free>=(32+8+8)*2**30,('native32+job8+floor8 resource guard',free)
def dec(v):return v.decode()if isinstance(v,bytes)else str(v)
def ser(h,k):
 z=h[k]
 if isinstance(z,h5py.Group):
  cats=[dec(x)for x in z['categories'][:]];return[cats[int(i)]if i>=0 else''for i in z['codes'][:]]
 return[dec(x)for x in z[:]]
def idx(h,k):return ser(h,k+'/'+dec(h[k].attrs.get('_index','_index')))
def save(df,p):df.to_csv(p,sep='\t',index=False,compression='gzip'if str(p).endswith('.gz')else None)
assert not O.exists();guard();O.mkdir();claim={'PID':os.getpid(),'UTC':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'source_sha256':sha(Path(__file__))};atom(O/'RUNNING.json',claim)
hgnc=OLD/'shared_inputs/programs_01_03/sources/hgnc_complete_set.txt';hg=pd.read_csv(hgnc,sep='\t',dtype=str).fillna('');hg=hg.loc[hg.status.eq('Approved')&hg.entrez_id.ne('')&hg.ensembl_gene_id.ne('')].copy();hg=hg.loc[~hg.entrez_id.duplicated(keep=False)&~hg.ensembl_gene_id.duplicated(keep=False)];emap=hg.set_index('ensembl_gene_id');imap=hg.set_index('entrez_id')
rec={'status':'INPUT_ONLY_RUNNING','source':claim,'canonical_axis':'HGNC-approved one-to-one Ensembl/Entrez; disease GS original MAGMA Entrez; never RNA backfill','HGNC_path':str(hgnc),'HGNC_sha256':sha(hgnc),'technical_gene_min_cells':50,'all_final_cells_preserved':True,'covariates':['log1p_full_input_UMI','log1p_full_input_nGenes'],'normalization_order':'cov/UMI on full original counts; feature one-to-one plus detection>=50; normalize admitted sparse counts to1e4 per cell THEN log1p, equivalent to official CLI gene-filter then normalize; full-input library totals are cov only','adj_prop':None,'cohorts':{}}
try:
 for st in ['GSE241132_D0','GSE220300','GSE335482']:
  guard();t=time.monotonic();d=O/st;d.mkdir();inp={}
  if st=='GSE241132_D0':
   src=OLD/'projects/03_susceptibility/inputs/Bisque_reference/D0';cp=src/'actual_cells_counts.mtx.gz';gp=src/'genes.tsv';mp=src/'actual_cell_manifest.tsv'
   m=pd.read_csv(mp,sep='\t',dtype=str,keep_default_na=False);g=pd.read_csv(gp,sep='\t',dtype=str,keep_default_na=False);assert m.cell_id.is_unique and g.entrez_id.is_unique
   with gzip.open(cp,'rb')as h,threadpool_limits(limits=2):X=sio.mmread(h).T.tocsr()
   assert X.shape==(len(m),len(g));g['source_gene_id']=g.entrez_id;g['source_symbol']=g.entrez_id.map(imap.symbol).fillna('');g['canonical_entrez']=g.entrez_id.where(g.entrez_id.isin(imap.index),'')
   m['lineage']=m.cellType;m['FB_state']='NA';m['region_original']='D0';m['study_id']=st;m['independent_patient_id']=m.patient_id;m['region_P2_swap']=m.region_original;m['include_P2_omit']=True
   inp={'counts_path':str(cp),'genes_path':str(gp),'cell_manifest_path':str(mp),'counts_sha256':sha(cp),'genes_sha256':sha(gp),'cell_manifest_sha256':sha(mp),'input_counts_unit':'MTX integer gene_by_cell transposed; no expression scaling before QA'}
  else:
   cp=OLD/f'projects/01_regional/results/atlas/{st}/atlas_provisional.h5ad';ap=OLD/f'projects/01_regional/results/lineage_annotation_v2/{st}_all_QC_final_mutually_exclusive_lineages.tsv.gz';vp=OLD/f'projects/01_regional/versions/lineage_v2/results/pseudobulk/{st}/all_QC_cell_annotations.tsv.gz';fp=OLD/f'projects/01_regional/versions/lineage_v2/results/pseudobulk/{st}/FB_states.h5ad'
   m=pd.read_csv(ap,sep='\t',dtype=str,keep_default_na=False);v=pd.read_csv(vp,sep='\t',dtype=str,keep_default_na=False).set_index('cell_id');assert m.cell_id.is_unique;assert set(m.cell_id)==set(v.index);assert (m.set_index('cell_id').major_type_final==v.reindex(m.cell_id).major_type_frozen).all()
   with h5py.File(cp)as h:
    ci=idx(h,'obs');assert ci==list(m.cell_id),'Frozen final row order must match original counts, no implicit reorder';z=h['layers/counts'];X=sp.csr_matrix((z['data'][:],z['indices'][:],z['indptr'][:]),shape=tuple(z.attrs['shape']));g=pd.DataFrame({'source_gene_id':ser(h,'var/source_gene_id'),'source_symbol':ser(h,'var/gene_symbol')})
   with h5py.File(fp)as h:states=dict(zip(idx(h,'obs'),ser(h,'obs/FB_state')))
   m['lineage']=m.major_type_final;m['FB_state']=m.cell_id.map(states).fillna('NA');assert set(m.loc[m.lineage.eq('fibroblast'),'cell_id'])==set(states);assert len(states)==(6774 if st=='GSE335482'else 26498)
   m['independent_patient_id']=m.patient_id.where(~m.patient_id.str.contains(r'unknown|^NA$|^nan$|^$',case=False,regex=True),'');m['region_P2_swap']=m.region_original;m['include_P2_omit']=True
   if st=='GSE335482':
    p2=m.patient_id.eq('GSE335482::P2');assert p2.any();m.loc[p2,'region_P2_swap']=m.loc[p2,'region_original'].map({'hypercellular':'infiltrating','infiltrating':'hypercellular'});assert m.region_P2_swap.notna().all();m.loc[p2,'include_P2_omit']=False
   en=g.source_gene_id.str.replace(r'\.[0-9]+$','',regex=True);g['canonical_entrez']=en.map(emap.entrez_id).fillna('')
   inp={'counts_path':str(cp),'counts_layer':'layers/counts','counts_sha256':sha(cp),'authoritative_lineage_path':str(ap),'authoritative_lineage_sha256':sha(ap),'version_lineage_path':str(vp),'version_lineage_sha256':sha(vp),'FB_state_path':str(fp),'FB_state_sha256':sha(fp),'input_counts_unit':'original raw UMI layer; full unit scan in this freeze'}
  assert np.issubdtype(X.dtype,np.integer)and X.data.min()>=0;X.sum_duplicates();X.eliminate_zeros();assert X.shape==(len(m),len(g));assert X.data.min()>0
  U=np.asarray(X.sum(1)).ravel();N=np.diff(X.indptr);D=np.bincount(X.indices,minlength=X.shape[1]);assert(U>0).all()and(N>0).all()
  if 'total_counts'in m:assert np.array_equal(U,m.total_counts.astype(float).to_numpy())and np.array_equal(N,m.n_genes_by_counts.astype(float).to_numpy())
  g['source_index']=np.arange(len(g));g['detected_cells']=D;g['gene_id_ambiguous_within_input']=g.canonical_entrez.ne('')&g.canonical_entrez.duplicated(keep=False);g['passes_min50']=D>=50;g['admitted']=g.canonical_entrez.ne('')&~g.gene_id_ambiguous_within_input&g.passes_min50;g['reason']=np.select([g.canonical_entrez.eq(''),g.gene_id_ambiguous_within_input,~g.passes_min50],['HGNC_missing_or_non_one_to_one','duplicate_canonical_Entrez','detected_in_less_than50_cells'],default='ADMITTED')
  m['full_input_UMI']=U;m['full_input_nGenes']=N;m['log1p_full_input_UMI']=np.log1p(U);m['log1p_full_input_nGenes']=np.log1p(N);cov=m[['cell_id','log1p_full_input_UMI','log1p_full_input_nGenes']];assert np.linalg.matrix_rank(np.column_stack([np.ones(len(cov)),cov.iloc[:,1:].to_numpy()]))==3
  for df,name in [(g,'complete_gene_admission.tsv.gz'),(g.loc[g.admitted,['source_index','source_gene_id','source_symbol','canonical_entrez','detected_cells']],'admitted_gene_axis.tsv'),(m,'all_final_cell_metadata.tsv.gz'),(cov,'technical_covariates.tsv.gz')]:save(df,d/name)
  rec['cohorts'][st]={'input':inp,'n_cells':len(m),'input_genes':len(g),'admitted_genes':int(g.admitted.sum()),'lineages':m.lineage.value_counts().to_dict(),'FB_states':m.FB_state.value_counts().to_dict(),'independent_patients':sorted(m.loc[m.independent_patient_id.ne(''),'independent_patient_id'].unique()),'patientNA_context_cells':int(m.independent_patient_id.eq('').sum()),'coarse_lineage_family':sorted(m.lineage.unique()),'FB_state_family':sorted(m.loc[m.FB_state.ne('NA'),'FB_state'].unique()),'cell_IDs_unique':True,'nonnegative_integer_full_scan':True,'full_input_QC_matches_prior_if_present':True,'full_input_UMI_range':[int(U.min()),int(U.max())],'covariate_rank_with_intercept':3,'files':{f.name:{'bytes':f.stat().st_size,'sha256':sha(f)}for f in d.iterdir()},'elapsed_seconds':time.monotonic()-t}
  del X,g,m,cov,U,N,D;gc.collect();atom(O/'status.json',rec)
 sets={st:set(pd.read_csv(O/st/'admitted_gene_axis.tsv',sep='\t',dtype=str).canonical_entrez)for st in rec['cohorts']};rec['common_measurable_Entrez_count_report_only']=len(set.intersection(*sets.values()));rec['status']='COMPLETE_INPUT_FREEZE_NO_SCDRS_FIT';rec['peak_RSS_bytes']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024;atom(O/'COMPLETE.json',rec);print(json.dumps({'status':rec['status'],'cohorts':{k:{x:v[x]for x in ['n_cells','input_genes','admitted_genes','patientNA_context_cells']}for k,v in rec['cohorts'].items()},'common_report_only':rec['common_measurable_Entrez_count_report_only'],'peak_RSS_bytes':rec['peak_RSS_bytes'],'receipt_sha256':sha(O/'COMPLETE.json')}))
except BaseException as e:rec['status']='FAILED_PRESERVED_NO_RETRY';rec['error']=repr(e);atom(O/'FAILED.json',rec);raise
