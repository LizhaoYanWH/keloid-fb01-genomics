#!/usr/bin/env python
"""Discovery-followup complete fixed15TF coverage, same1450 token cells; no model."""
from runtime_common_v1 import *
import gzip,numpy as np,pandas as pd

def main():
 f,g=approved('prepare');claim('prepare');out=D/'prepared';out.mkdir()
 sys.path.insert(0,str(L/'scripts'));from gf_encoding_reference import load_assets,canonical_id
 a,_=load_assets(L/'inputs/gf_reference')
 with gzip.open(BASE/'prepared/unique_cells_official_tokens.jsonl.gz','rt')as z:cells=list(map(json.loads,z))
 M=pd.read_csv(BASE/'prepared/all12_patient_balanced32_pool_memberships.tsv',sep='\t');assert len(M)==1536 and len(cells)==1450
 index={(c['study'],c['cell_id']):c for c in cells};genes=[];visibility=[];queries=[]
 for target in PROTEINS:
  raw,can=canonical_id(target,a['names'],a['mapping'],a['tokens']);token=a['tokens'].get(can)
  genes.append({'gene':target,'official_original_ensembl':raw,'canonical_ensembl':can,'token':token,'status':'VOCAB_MEASURABLE'if token is not None else'NO_OFFICIAL_MAPPING_OR_VOCAB_NA','followup_status':'fixed15_family_discovery_followup_not_independent_validation'})
  for m in M.to_dict('records'):
   c=index[(m['study'],m['cell_id'])];seq=c['input_ids'];yes=token is not None and token in seq[1:-1]
   visibility.append({**m,'gene':target,'visible':yes,'vocab_measurable':token is not None,'input_length':len(seq),'original_token_length':len(seq),'source_full_library_UMI':c['total_UMI'],'role':'fixed15_discovery_followup_TF','target_rank':seq.index(token)if yes else None,'status':'VISIBLE_ADMITTED'if yes else'NOT_VISIBLE_OR_NO_VOCAB_NA'})
   if yes:queries.append({'study':m['study'],'cell_id':m['cell_id'],'gene':target,'token':int(token),'rank':seq.index(token),'length':len(seq)})
 pd.DataFrame(genes).to_csv(out/'all15_TF_official_vocabulary.tsv',sep='\t',index=False,na_rep='NA')
 V=pd.DataFrame(visibility);assert len(V)==23040;V.to_csv(out/'all15_TF_fixed_membership_visibility_and_NA.tsv.gz',sep='\t',index=False,compression='gzip',na_rep='NA')
 Q=pd.DataFrame(queries).drop_duplicates(['study','cell_id','gene']).sort_values(['length','study','cell_id','gene']).reset_index(drop=True);Q['query_id']=np.arange(len(Q));Q.to_csv(out/'unique_TF_query_family_before_forward.tsv.gz',sep='\t',index=False,compression='gzip')
 V.groupby(['study','context','seed','patient','gene']).agg(n_selected=('cell_id','size'),n_visible=('visible','sum')).reset_index().to_csv(out/'all15_TF_patient_context_coverage.tsv',sep='\t',index=False)
 atom(out/'COMPLETE.json',{'status':'COMPLETE_FIXED15_TF_IDENTITY_COVERAGE_NO_MODEL','TF15':PROTEINS,'n_unique_cells_reused':1450,'n_memberships':1536,'n_query_per_operator':len(Q),'n_query_all2operators':len(Q)*2,'after_embedding_upper_bytes':len(Q)*2*768*4,'existing4controls_are_diagnostic_not_matched_to15TF':True,'no_new_encoding_or_sampling':True,'files':file_receipt(out)})
 print((out/'COMPLETE.json').read_text())
if __name__=='__main__':main()
