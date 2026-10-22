#!/usr/bin/env python
"""Root-gated complete identity-only 32patient pools/token/query freeze; no model."""
from runtime_common_v1 import *
import gzip,collections
import h5py,numpy as np,pandas as pd
from scipy import sparse

def arr(x):
 if isinstance(x,h5py.Dataset):
  a=x[:];return np.asarray([s.decode() if isinstance(s,bytes) else str(s) for s in a]) if a.dtype.kind in 'OSU' else a
 c=arr(x['categories']);return np.array([c[int(v)] if v>=0 else '' for v in x['codes'][:]])
def main():
 f,g=approved('prepare');claim('prepare');out=D/'prepared';out.mkdir();start=time.time()
 sys.path.insert(0,str(L/'scripts'));from gf_encoding_reference import load_assets,build_axis,encode_sparse_column,canonical_id
 a,am=load_assets(L/'inputs/gf_reference')
 rf=json.loads((P/'reviews/scDRS_production_preparation_v1/RNA_freeze_v1/COMPLETE.json').read_text())
 prep=json.loads((E/'results/PREPARATION_COMPLETE.json').read_text())
 candidates=pd.read_csv(P/'extensions_20261021/gene_state_bridge_v1/results/all_primary_PP4_ge08_unique_genes.tsv',sep='\t')
 records=[]
 for r in candidates.to_dict('records'):
  _,canon=canonical_id(r['gene_id'],a['names'],a['mapping'],a['tokens'])
  records.append({**r,'model_canonical_ensembl':canon,'model_token':a['tokens'].get(canon),'role':'complete_genetic_candidate','forward_family':r['approved_symbol'] in PROTEINS})
 for symbol in CONTROLS:
  ens,canon=canonical_id(symbol,a['names'],a['mapping'],a['tokens'])
  records.append({'approved_symbol':symbol,'gene_id':ens,'model_canonical_ensembl':canon,'model_token':a['tokens'].get(canon),'role':'existing_expression_background','forward_family':True})
 C=pd.DataFrame(records);assert set(C.loc[C.forward_family,'approved_symbol'])==set(PANEL)
 assert C.loc[C.forward_family,'model_token'].notna().all()
 C['measurement_status']=np.where(C.model_token.notna(),'OFFICIAL_V2_VOCABULARY_MEASURABLE','NO_OFFICIAL_V2_VOCABULARY_NA')
 C.to_csv(out/'complete9_candidate_and4_control_token_audit.tsv',sep='\t',index=False,na_rep='NA')
 memberships=[]
 for cx in prep['contexts']:
  m=pd.read_csv(Path(cx['path'])/'cells.tsv',sep='\t',dtype=str)
  assert len(m)==488 and m.independent_patient_id.value_counts().to_list()==[122]*4
  for patient,z in m.groupby('independent_patient_id',sort=True):
   z=z.copy();z['blind_hash']=[hashlib.sha256(('N02_GF32_v1|'+cx['study']+'|'+cx['context']+'|'+str(cx['seed'])+'|'+patient+'|'+c).encode()).hexdigest() for c in z.cell_id]
   chosen=z.sort_values(['blind_hash','cell_id']).head(32)
   for r in chosen.to_dict('records'):memberships.append({'study':cx['study'],'context':cx['context'],'seed':int(cx['seed']),'patient':patient,'cell_id':r['cell_id'],'capture_id':r['capture_id'],'original_FB_state':r['FB_state'],'blind_selection_hash':r['blind_hash']})
 M=pd.DataFrame(memberships);assert len(M)==1536 and not M.duplicated(['study','context','seed','cell_id']).any()
 M.to_csv(out/'all12_patient_balanced32_pool_memberships.tsv',sep='\t',index=False)
 cells=[];axes=[]
 for st in ['GSE220300','GSE335482']:
  guard(g,'prepare');keep=set(M[M.study.eq(st)].cell_id);inp=rf['cohorts'][st]['input']
  meta=pd.read_csv(P/'reviews/scDRS_production_preparation_v1/RNA_freeze_v1'/st/'all_final_cell_metadata.tsv.gz',sep='\t',dtype=str,keep_default_na=False).set_index('cell_id')
  with h5py.File(inp['counts_path'],'r') as h:
   obs=arr(h['obs'][h['obs'].attrs['_index']]);genes=arr(h['var'][h['var'].attrs['_index']]);axis,srcmap,audit=build_axis(genes,a)
   tok=np.array([a['tokens'][x] for x in axis]);med=np.array([a['medians'][x] for x in axis]);ix={c:i for i,c in enumerate(obs)};assert keep<=set(ix)
   axes.extend({'study':st,**r}for r in audit)
   for c in sorted(keep):
    i=ix[c];x=h['layers/counts'];lo,hi=map(int,x['indptr'][i:i+2]);jj=x['indices'][lo:hi];v=x['data'][lo:hi];total=float(v.sum())
    assert int(total)==int(float(meta.loc[c,'full_input_UMI'])) and total>0
    seq,_,_=encode_sparse_column(jj,v,total,srcmap,tok,med,a['tokens']['<cls>'],a['tokens']['<eos>'])
    cells.append({'study':st,'cell_id':c,'patient':meta.loc[c,'independent_patient_id'],'capture_id':meta.loc[c,'capture_id'],'FB_state':meta.loc[c,'FB_state'],'total_UMI':int(total),'input_ids':seq})
  # Every available pre-existing335 token identity must agree; no old forward is used.
  if st=='GSE335482':
   previous=OLD/'projects/02_selectivity/versions/lineage_v2/prepared/GSE335482/cells.jsonl.gz'
   with gzip.open(previous,'rt')as z:old={r['cell_id']:r['input_ids']for r in map(json.loads,z)}
   assert all(r['input_ids']==old[r['cell_id']]for r in cells if r['study']==st),'EXACT_PRIOR335_TOKEN_IDENTITY'
 pd.DataFrame(axes).to_csv(out/'all_original_count_gene_axis_official_token_mapping.tsv.gz',sep='\t',index=False,compression='gzip')
 with gzip.open(out/'unique_cells_official_tokens.jsonl.gz','wt')as z:
  for r in cells:z.write(json.dumps(r)+'\n')
 q=[];visibility=[];bycell={(r['study'],r['cell_id']):r for r in cells}
 for r in C.to_dict('records'):
  gene=r['approved_symbol'];token=int(r['model_token']) if pd.notna(r['model_token']) else None
  for m in M.to_dict('records'):
   c=bycell[(m['study'],m['cell_id'])];visible=token is not None and token in c['input_ids'][1:-1]
   visibility.append({**m,'gene':gene,'role':r['role'],'vocab_measurable':token is not None,'visible':visible,'source_full_library_UMI':c['total_UMI'],'original_token_length':len(c['input_ids']),'target_rank':c['input_ids'].index(token)if visible else None,'coverage_status':'VISIBLE_FORWARD_ADMITTED'if visible and r['forward_family'] else 'NOT_VISIBLE_OR_NO_MODEL_VOCAB_NA'})
   if visible and r['forward_family']:
    q.append({'study':m['study'],'cell_id':m['cell_id'],'gene':gene,'token':token,'rank':c['input_ids'].index(token),'length':len(c['input_ids'])})
 Q=pd.DataFrame(q).drop_duplicates(['study','cell_id','gene']).sort_values(['length','study','cell_id','gene']).reset_index(drop=True)
 Q['query_id']=np.arange(len(Q));Q.to_csv(out/'unique_single_gene_queries_before_any_forward.tsv.gz',sep='\t',index=False,compression='gzip')
 V=pd.DataFrame(visibility);V.to_csv(out/'all9candidate4control_visibility_and_NA_before_forward.tsv.gz',sep='\t',index=False,compression='gzip',na_rep='NA')
 V.groupby(['study','context','seed','patient','gene'],dropna=False).agg(n_identity_selected=('cell_id','size'),visible_cells=('visible','sum')).reset_index().to_csv(out/'all_gene_patient_context_coverage_before_forward.tsv',sep='\t',index=False)
 assert len(Q)*2<=33792
 bound=(len(cells)+2*len(Q))*768*4
 assert bound<128*1024**2,'RAW_EMBEDDING_BUDGET'
 atom(out/'COMPLETE.json',{'status':'COMPLETE_IDENTITY_BLIND_INPUT_QUERY_FREEZE_NO_FORWARD','n_memberships':len(M),'n_unique_cells':len(cells),'n_unique_gene_cell_queries_per_operator':len(Q),'n_delete_queries_all2operators':len(Q)*2,'embedding_uncompressed_bytes_upper':bound,'complete_candidate9':list(candidates.approved_symbol),'protein_forward_family7':PROTEINS,'control4':CONTROLS,'no_target_specific_matched_controls_for_other5proteins':True,'position_preserved_not_independent_model':True,'files':file_receipt(out),'elapsed_seconds':time.time()-start})
 print((out/'COMPLETE.json').read_text())
if __name__=='__main__':main()
