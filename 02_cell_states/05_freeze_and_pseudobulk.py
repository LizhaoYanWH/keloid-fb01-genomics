#!/usr/bin/env python3
"""Freeze reviewed broad labels, preserve counts, construct FB states and pseudobulk."""
import argparse,hashlib,inspect,json,os,platform
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

def main():
 p=argparse.ArgumentParser();p.add_argument('--project-root',type=Path,required=True);p.add_argument('--study',required=True);p.add_argument('--reviewed-map',type=Path,required=True);p.add_argument('--all-cell-annotations',type=Path);a=p.parse_args()
 if platform.system()!='Linux' or os.getenv('FIVE_STUDIES_COMPUTE')!='autodl':raise SystemExit('AutoDL only')
 sc.settings.n_jobs=4;np.random.seed(20261020)
 atlas=a.project_root/'results/atlas'/a.study;out=a.project_root/'results/pseudobulk'/a.study
 if (out/'pseudobulk_receipt.json').exists():raise RuntimeError('Existing pseudobulk; create a new version')
 mapping=json.loads(a.reviewed_map.read_text());assert mapping['program_effects_seen_at_annotation_freeze'] is False
 obj=ad.read_h5ad(atlas/'atlas_provisional.h5ad');cluster=json.loads((atlas/'atlas_annotation_receipt.json').read_text())['cluster_key']
 labels=obj.obs[cluster].astype(str).map(mapping['cluster_to_major_type'])
 if a.all_cell_annotations:
  final=pd.read_csv(a.all_cell_annotations,sep='\t',index_col=0,low_memory=False);labels=final.major_type_final.reindex(obj.obs_names)
 if labels.isna().any():raise ValueError('Incomplete reviewed cluster map')
 obj.obs['major_type_frozen']=labels.astype('category');counts=obj.layers['counts'].tocsr()
 out.mkdir(parents=True,exist_ok=True);obj.obs.to_csv(out/'all_QC_cell_annotations.tsv.gz',sep='\t',compression='gzip',index_label='cell_id')
 composition=obj.obs.groupby(['capture_id','patient_id','region_original','major_type_frozen'],observed=True).size().rename('n_cells').reset_index()
 composition['capture_QC_cells']=composition.capture_id.map(obj.obs.groupby('capture_id',observed=True).size()).astype(np.int64)
 composition['captured_cell_fraction']=composition.n_cells/composition.capture_QC_cells
 composition.to_csv(out/'all_lineage_composition.tsv',sep='\t',index=False)
 fb=obj[labels.eq('fibroblast').to_numpy()].copy()
 if fb.n_obs==0:raise ValueError('No marker-supported fibroblast compartment')
 # Fixed state construction precedes program scores; neither inflammatory genes nor state support are selected by effects.
 sc.pp.highly_variable_genes(fb,n_top_genes=min(3000,fb.n_vars),flavor='seurat',batch_key='capture_id')
 sc.pp.pca(fb,n_comps=40,mask_var='highly_variable',random_state=20261020)
 import harmonypy as hm
 kwargs={'max_iter_harmony':20,'random_state':20261020,'verbose':False}
 if 'ncores' in inspect.signature(hm.run_harmony).parameters:kwargs['ncores']=4
 patient=fb.obs.patient_id.astype(str)
 fb.obs['patient_batch']=patient.where(~patient.isin(['NA','']),fb.obs.capture_id.astype(str))
 h=hm.run_harmony(fb.obsm['X_pca'],fb.obs,['patient_batch'],**kwargs);z=h.Z_corr;fb.obsm['X_pca_harmony']=z.T if z.shape[1]==fb.n_obs else z
 sc.pp.neighbors(fb,use_rep='X_pca_harmony',n_pcs=30,n_neighbors=20,random_state=20261020)
 sc.tl.leiden(fb,resolution=.4,key_added='FB_state',random_state=20261020,flavor='igraph',n_iterations=2)
 fb.obs['FB_state']='FB'+fb.obs.FB_state.astype(str).str.zfill(2)
 sc.tl.umap(fb,random_state=20261020)
 # Raw counts remain the only expression input for pseudobulk.
 genes=fb.var[['source_gene_id','gene_symbol']].copy();genes.to_csv(out/'genes.tsv',sep='\t',index=False)
 for mode,keys in [('all_FB',['capture_id']),('FB_state',['capture_id','FB_state'])]:
  rows=[];metadata=[]
  for key,indices in fb.obs.groupby(keys,observed=True).indices.items():
   first=fb.obs.iloc[indices[0]];raw=fb.layers['counts'][indices,:]
   count=np.asarray(raw.sum(axis=0)).ravel();rows.append(count)
   record={key:str(first[key]) for key in ['study_id','gsm','patient_id','capture_id','region_original','activity','role']}
   record['FB_state']=str(first.FB_state) if mode=='FB_state' else 'all_FB'
   record['pseudobulk_id']=str(first.capture_id)+'::'+record['FB_state'];record['n_FB_cells']=len(indices);record['total_gene_UMI']=int(count.sum())
   record['median_cell_UMI']=float(np.median(np.asarray(raw.sum(axis=1)).ravel()));record['mean_mt_pct']=float(fb.obs.iloc[indices].pct_counts_mt.mean())
   metadata.append(record)
  matrix=pd.DataFrame(np.vstack(rows).T,index=genes.source_gene_id,columns=[row['pseudobulk_id'] for row in metadata])
  matrix.to_csv(out/(mode+'_counts.tsv.gz'),sep='\t',compression='gzip',index_label='source_gene_id')
  pd.DataFrame(metadata).to_csv(out/(mode+'_samples.tsv'),sep='\t',index=False)
 fb.write_h5ad(out/'FB_states.h5ad',compression='gzip')
 receipt={'study':a.study,'all_QC_cells':obj.n_obs,'FB_cells':fb.n_obs,'FB_states':sorted(fb.obs.FB_state.unique().tolist()),'reviewed_map_path':str(a.reviewed_map),
          'reviewed_map_sha256':hashlib.sha256(a.reviewed_map.read_bytes()).hexdigest(),'annotation_program_effects_seen':bool(a.all_cell_annotations),
          'lineage_revision':str(a.all_cell_annotations) if a.all_cell_annotations else 'none; initial marker-reviewed broad labels','lineage_revision_timing':'old regional results seen, revised effects not seen when independent-reference/multi-marker rule frozen' if a.all_cell_annotations else 'initial labels before regional effects',
          'state_construction':'FB-only 3000 capture-aware HVGs/PCA40, patient-only Harmony, neighbors20, Leiden resolution0.4; unchanged raw counts',
          'state_coverage_line_for_secondary_comparison':30,'counts_status':'unintegrated raw UMI pseudobulk; no minimum program nonzero gate'}
 (out/'pseudobulk_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)

if __name__=='__main__':main()
