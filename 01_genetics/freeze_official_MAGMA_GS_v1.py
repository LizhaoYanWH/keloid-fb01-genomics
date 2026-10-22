"""Real complete MAGMA -> unchanged official munge_gs top1000 Z weights; never RNA selection/score."""
from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:os.environ[k]='2'
os.sched_setaffinity(0,sorted(os.sched_getaffinity(0))[:2])
import json,hashlib,runpy,time,shutil
import pandas as pd,numpy as np
R=Path('/root/five_new_studies_20261020/projects/N02_keloid_genetic_cell_state');N=R/'results/MAGMA_complete_EUR_native_v1';O=R/'reviews/scDRS_production_preparation_v1';D=O/'official_MAGMA_GS_v1'
def sha(f):return hashlib.sha256(Path(f).read_bytes()).hexdigest()
assert all(shutil.disk_usage(x).free>=5*2**30 for x in ['/root','/root/autodl-tmp']);assert not D.exists();D.mkdir()
j=json.loads((N/'COMPLETE.json').read_text());assert j['status']=='COMPLETE_NATIVE_MAGMA_FULL_FAMILY_TECHNICAL_ONLY'
checked=[]
for x in j['outputs']:
 p=Path(x['path']);assert p.stat().st_size==x['bytes']and sha(p)==x['SHA256'];checked.append({'path':str(p),'sha256':x['SHA256']})
full=pd.read_csv(N/'all19427_gene_results_and_NA_ledger.tsv',sep='\t',dtype=str,keep_default_na=False);assert len(full)==19427
native=pd.read_csv(N/'complete_EUR.genes.out',sep=r'\s+',comment='#',dtype={'GENE':str});assert native.GENE.is_unique;assert len(native)==18334
pv=native[['GENE','P']].rename(columns={'P':'Keloid_EUR'}).copy();pv['Keloid_EUR']=pd.to_numeric(pv.Keloid_EUR);assert np.isfinite(pv.Keloid_EUR).all()and pv.Keloid_EUR.between(0,1).all();pv.to_csv(D/'complete_MAGMA_gene_P.tsv',sep='\t',index=False)
api=runpy.run_path(str(R/'tools/scdrs104_env/bin/scdrs'),run_name='official_scdrs_import_not_main');api['munge_gs'](out_file=str(D/'Keloid_EUR.top1000_Z.gs'),pval_file=str(D/'complete_MAGMA_gene_P.tsv'),weight='zscore',fdr=None,fwer=None,n_min=1000,n_max=1000)
gs=pd.read_csv(D/'Keloid_EUR.top1000_Z.gs',sep='\t');assert len(gs)==1;members=[x.rsplit(':',1)for x in gs.GENESET.iloc[0].split(',')];assert len(members)==1000 and len({x[0]for x in members})==1000;g=pd.DataFrame(members,columns=['canonical_entrez','official_GS_weight_string']);g['rank_in_official_GS']=np.arange(1,1001);g=g.merge(native.rename(columns={'GENE':'canonical_entrez'}),on='canonical_entrez',validate='1:1',how='left');g.to_csv(D/'fixed_top1000_genetic_members_and_native_values.tsv',sep='\t',index=False)
rna=O/'RNA_freeze_v1/COMPLETE.json';f=json.loads(rna.read_text());assert f['status']=='COMPLETE_INPUT_FREEZE_NO_SCDRS_FIT';coverage=[]
for st,c in f['cohorts'].items():
 ax=pd.read_csv(O/f'RNA_freeze_v1/{st}/admitted_gene_axis.tsv',sep='\t',dtype=str);present=set(ax.canonical_entrez);z=g[['canonical_entrez','rank_in_official_GS','official_GS_weight_string']].copy();z['cohort']=st;z['present_after_fixed_mapping_and_min50']=z.canonical_entrez.isin(present);coverage.append(z)
pd.concat(coverage).to_csv(D/'all1000_members_x3_RNA_coverage_no_backfill.tsv',sep='\t',index=False)
r={'status':'COMPLETE_OFFICIAL_TOP1000_GS_NO_RNA_SELECTION_NO_SCORE','source_sha256':sha(Path(__file__)),'native_complete_sha256':sha(N/'COMPLETE.json'),'native_full_gene_family':19427,'native_output_genes':len(native),'native_outputs_all_SHA_checked':checked,'official_installed_bin_sha256':sha(R/'tools/scdrs104_env/bin/scdrs'),'RNA_freeze_sha256':sha(rna),'official_parameters':{'n_min':1000,'n_max':1000,'weight':'zscore','fdr':None,'fwer':None},'actual_native_genes_used':len(pv),'top1000_frozen_before_RNA_coverage':True,'no_ABF_or_RNA_filtering_or_gene_backfill':True,'coverage':{st:int(z.present_after_fixed_mapping_and_min50.sum())for st,z in pd.concat(coverage).groupby('cohort')},'files':{p.name:{'bytes':p.stat().st_size,'sha256':sha(p)}for p in D.iterdir()}}
q=D/'COMPLETE.json.tmp';q.write_text(json.dumps(r,indent=2)+'\n');q.replace(D/'COMPLETE.json');print(json.dumps({'status':r['status'],'coverage':r['coverage'],'receipt_sha256':sha(D/'COMPLETE.json')}))
