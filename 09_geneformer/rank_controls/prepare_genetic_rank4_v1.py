#!/usr/bin/env python
"""Complete original9-family rank-matched input controls, blind to reference effect."""
from runtime_common_v1 import *
import gzip,numpy as np,pandas as pd
TF15=['RELA','NFKB1','NFKB2','REL','RELB','JUN','JUNB','JUND','FOS','FOSB','FOSL1','FOSL2','SMAD2','SMAD3','SMAD4']

def main():
 f,g=approved('prepare');claim('prepare');out=D/'prepared';out.mkdir()
 sys.path.insert(0,str(L/'scripts'));from gf_encoding_reference import load_assets,canonical_id
 a,_=load_assets(L/'inputs/gf_reference')
 with gzip.open(BASE/'prepared/unique_cells_official_tokens.jsonl.gz','rt')as z:cells=list(map(json.loads,z))
 cells={(c['study'],c['cell_id']):c for c in cells}
 C=pd.read_csv(BASE/'prepared/complete9_candidate_and4_control_token_audit.tsv',sep='\t');nine=list(C[C.role.eq('complete_genetic_candidate')].approved_symbol);assert len(nine)==9
 exclude=nine+TF15+CONTROLS;forbidden=set();audit=[]
 for gene in exclude:
  _,can=canonical_id(gene,a['names'],a['mapping'],a['tokens']);token=a['tokens'].get(can)
  if token is not None:forbidden.add(token)
  audit.append({'gene':gene,'canonical_ensembl':can,'token':token,'excluded_from_background':True})
 pd.DataFrame(audit).to_csv(out/'all9_plus15TF_plus4oldcontrol_exclusion_audit.tsv',sep='\t',index=False,na_rep='NA')
 T=pd.read_csv(BASE/'prepared/unique_single_gene_queries_before_any_forward.tsv.gz',sep='\t');T=T[T.gene.isin(PROTEINS)]
 bind=[];control={}
 for r in T.itertuples(index=False):
  c=cells[(r.study,r.cell_id)];seq=c['input_ids'];rank=int(r.rank);assert seq[rank]==r.token
  candidates=[]
  for j in range(max(1,rank-32),min(len(seq)-1,rank+33)):
   token=seq[j]
   if token in forbidden:continue
   h=hashlib.sha256(('N02rank4_v1|'+r.study+'|'+r.cell_id+'|'+str(int(r.token))+'|'+str(token)).encode()).hexdigest();candidates.append((abs(j-rank),h,j,token))
  selected=sorted(candidates)[:4];assert len(selected)==4,'NO_RESULT_BASED_REPLACEMENT_IF_NO4_RANK_CONTROLS'
  for ordinal,(dist,h,j,token)in enumerate(selected):
   key=(r.study,r.cell_id,token);control[key]={'study':r.study,'cell_id':r.cell_id,'gene':'rankcontrol_token_'+str(token),'token':token,'rank':j,'length':len(seq)}
   bind.append({'study':r.study,'cell_id':r.cell_id,'target_gene':r.gene,'target_token':int(r.token),'target_original_query_id':int(r.query_id),'target_rank':rank,'control_ordinal':ordinal,'control_token':token,'control_rank':j,'absolute_rank_distance':dist,'tie_identity_hash':h})
 Q=pd.DataFrame(control.values()).sort_values(['length','study','cell_id','token']).reset_index(drop=True);Q['query_id']=np.arange(len(Q));lookup={(r.study,r.cell_id,int(r.token)):int(r.query_id)for r in Q.itertuples(index=False)}
 B=pd.DataFrame(bind);B['reference_query_id']=[lookup[(r.study,r.cell_id,int(r.control_token))]for r in B.itertuples(index=False)]
 assert B.groupby(['study','cell_id','target_gene']).size().eq(4).all()
 Q.to_csv(out/'unique_rankcontrol_queries_before_forward.tsv.gz',sep='\t',index=False,compression='gzip');B.to_csv(out/'all_target_to4_rankcontrol_bindings_before_forward.tsv.gz',sep='\t',index=False,compression='gzip')
 V=pd.read_csv(BASE/'prepared/all9candidate4control_visibility_and_NA_before_forward.tsv.gz',sep='\t');V=V[V.gene.isin(nine)];assert len(V)==13824;V.to_csv(out/'all9_original_candidate_memberships_visibility_NA.tsv.gz',sep='\t',index=False,compression='gzip',na_rep='NA')
 assert 2*len(Q)<=30000 and 2*len(Q)*768*4<96*1024**2
 atom(out/'COMPLETE.json',{'status':'COMPLETE_FIXED_ORIGINAL9_FAMILY_RANK4_INPUT_CONTROL_FREEZE_NO_FORWARD','family9':nine,'protein7':PROTEINS,'n_target_unique_queries':len(T),'n_target_reference_bindings':len(B),'n_reference_query_per_operator':len(Q),'n_new_forward_two_operators':2*len(Q),'full_after_FP32_upper_bytes':2*len(Q)*768*4,'max_rank_distance':int(B.absolute_rank_distance.max()),'all_targetqueries_have4_visible_refs':True,'excludes9_plus15TF_plus4oldcontrols':True,'target_responses_previously_seen_not_independent_validation':True,'not_expression_matching_or_empirical_FDR':True,'no_target_baseline_profile_forward_rerun':True,'files':file_receipt(out)})
 print((out/'COMPLETE.json').read_text())
if __name__=='__main__':main()
