#!/usr/bin/env python
"""Original9-family fixed4 same-cell rank input controls; cached targets/baseline."""
from runtime_common_v1 import *
import gzip,argparse,collections
import numpy as np,pandas as pd

def unit(x):
 x=np.asarray(x,np.float64);n=np.linalg.norm(x,axis=-1,keepdims=True);assert np.all(n>0);return x/n

def read_prepared():
 p=D/'prepared';c=json.loads((p/'COMPLETE.json').read_text())
 for n,z in c['files'].items():assert sha(p/n)==z['sha256'] and (p/n).stat().st_size==z['bytes']
 with gzip.open(BASE/'prepared/unique_cells_official_tokens.jsonl.gz','rt')as z:cells=list(map(json.loads,z))
 q=pd.read_csv(p/'unique_rankcontrol_queries_before_forward.tsv.gz',sep='\t');return cells,q,c

def forward():
 f,g=approved('gf');claim('gf');cells,Q,prep=read_prepared();out=D/'geneformer_forward';out.mkdir();started=time.time()
 sys.path.insert(0,str(OLD/'shared'));from gpu_admission import acquire_gpu
 lease=acquire_gpu(32,'N02_original9_rank4_input_controls',trace_path=out/'gpu_admission.jsonl',stop_path=out/'STOP_AFTER_CHUNK')
 import torch
 torch.set_num_threads(2);torch.backends.cuda.matmul.allow_tf32=False;torch.backends.cudnn.allow_tf32=False
 total=torch.cuda.get_device_properties(0).total_memory;torch.cuda.set_per_process_memory_fraction(min(1.,30*GIB/total),0)
 sys.path.insert(0,str(L/'scripts'));from model_utils import load_bert,embed_sequences
 sys.path.insert(0,str(OLD/'projects/02_selectivity/scripts'));from forward_GF_equivalent import embed_cls_equivalent
 model=load_bert(MODEL);model.eval();assert next(model.parameters()).dtype==torch.float32
 assert model.config.hidden_size==768 and all(z.attention.self.position_embedding_type=='absolute'for z in model.bert.encoder.layer[:11])
 def embed(seqs,positions=None):
  if positions is None:return embed_cls_equivalent(model,seqs,batch=8)
  bert=model.bert;ans=[]
  with torch.inference_mode():
   for begin in range(0,len(seqs),8):
    ss=seqs[begin:begin+8];pos=positions[begin:begin+8];length=max(map(len,ss));ids=torch.zeros((len(ss),length),dtype=torch.long,device='cuda');pp=torch.zeros_like(ids);mask=torch.zeros_like(ids)
    for j,s in enumerate(ss):ids[j,:len(s)]=torch.tensor(s,device='cuda');pp[j,:len(s)]=torch.tensor(pos[j],device='cuda');mask[j,:len(s)]=1
    am=bert.get_extended_attention_mask(mask,ids.shape);h=bert.embeddings(input_ids=ids,position_ids=pp)
    for layer in bert.encoder.layer[:10]:h=layer(h,attention_mask=am)[0]
    layer=bert.encoder.layer[10];a=layer.attention.self;query=a.transpose_for_scores(a.query(h[:,:1,:]));key=a.transpose_for_scores(a.key(h));value=a.transpose_for_scores(a.value(h))
    z=torch.nn.functional.scaled_dot_product_attention(query,key,value,attn_mask=am,dropout_p=0.,is_causal=False).transpose(1,2).contiguous().view(h.shape[0],1,a.all_head_size)
    z=layer.attention.output(z,h[:,:1,:]);ans.append(layer.output(layer.intermediate(z),z)[:,0,:].float().cpu().numpy())
  return np.concatenate(ans)
 def check(stage):
  guard(g,'gf');assert not(out/'STOP_AFTER_CHUNK').exists(),'ROOT_STOP_AFTER_CHUNK'
  assert time.time()-started<1800,'THIRTY_MINUTE_RANKCONTROL_ENGINEERING_LIMIT_KEEP_ATOMIC_OUTPUTS'
  assert torch.cuda.max_memory_reserved()<=30*GIB,'OWN_ALLOCATOR_BUDGET'
  with(out/'resource_by_chunk.jsonl').open('a')as z:z.write(json.dumps({'UTC':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'stage':stage,'pid':os.getpid(),'allocated':torch.cuda.memory_allocated(),'reserved':torch.cuda.memory_reserved(),'peak_reserved':torch.cuda.max_memory_reserved(),'RSS':own_rss(),'elapsed':time.time()-started})+'\n')
 index={(c['study'],c['cell_id']):i for i,c in enumerate(cells)}
 seqs=[c['input_ids']for c in cells]
 # Reuse the completed same-model32x4096 profile, exact original1450 tokens,
 # baseline and four background responses. No original candidate/control rerun.
 prior_profile=json.loads((BASE/'geneformer_forward/PROFILE_COMPLETE.json').read_text())
 assert prior_profile['status']=='PASS32_LONG4096_BATCH8_AND_ORIGINAL_DELTA_EQUIVALENCE'
 assert prior_profile['FP32']is True and prior_profile['TF32']is False and prior_profile['batch']==8
 assert prior_profile['prepared_sha256']==sha(BASE/'prepared/COMPLETE.json')
 assert prior_profile['delta_relative_L2_max']<1e-3
 with np.load(BASE/'geneformer_forward/all_unique_cell_baseline.npz')as z:
  B=z['embeddings'];assert list(z['cell_ids'].astype(str))==[c['cell_id']for c in cells]and list(z['studies'].astype(str))==[c['study']for c in cells]
 assert np.isfinite(B).all()
 torch.cuda.reset_peak_memory_stats();check('reuse_verified_profile_baseline_no_new_profile_or_baseline')
 for operator in ['delete_natural','delete_position_preserved']:
  folder=out/operator;folder.mkdir()
  for begin in range(0,len(Q),128):
   part=Q.iloc[begin:begin+128];ss=[];pp=[]
   for r in part.itertuples():
    s=seqs[index[(r.study,r.cell_id)]];rank=int(r.rank);assert s[rank]==r.token and s[1:-1].count(int(r.token))==1
    ss.append(s[:rank]+s[rank+1:]);pp.append(list(range(rank))+list(range(rank+1,len(s))))
   Eafter=embed(ss,pp if operator=='delete_position_preserved'else None);assert Eafter.shape==(len(part),768)and np.isfinite(Eafter).all()
   path=folder/f'chunk_{begin:07d}.npz';tmp=path.with_suffix('.tmp.npz');np.savez_compressed(tmp,query_id=part.query_id.to_numpy(),embedding_after=Eafter);tmp.replace(path)
   check(operator+'_'+str(begin));print(json.dumps({'operator':operator,'begin':begin,'total_queries':len(Q),'elapsed_seconds':time.time()-started}),flush=True)
 atom(out/'COMPLETE.json',{'status':'COMPLETE_ORIGINAL9_FAMILY_NEW_RANK4_CONTROLS_TWO_OPERATORS_CACHED_TARGETS_BASELINE','n_unique_cells':len(cells),'query_count_per_operator':len(Q),'prepared_sha256':sha(D/'prepared/COMPLETE.json'),'elapsed_seconds':time.time()-started,'peak_reserved_bytes':int(torch.cuda.max_memory_reserved()),'files':file_receipt(out),'operators_are_same_model_sensitivity':True,'not_RNA_KO_or_therapy':True})

def readouts():
 f,g=approved('readouts');claim('readouts');cells,Q,prep=read_prepared();src=D/'geneformer_forward';fc=json.loads((src/'COMPLETE.json').read_text())
 for n,x in fc['files'].items():assert (src/n).stat().st_size==x['bytes']and sha(src/n)==x['sha256']
 out=D/'rank_readouts';out.mkdir();M=pd.read_csv(BASE/'prepared/all12_patient_balanced32_pool_memberships.tsv',sep='\t');V=pd.read_csv(D/'prepared/all9_original_candidate_memberships_visibility_NA.tsv.gz',sep='\t');binding=pd.read_csv(D/'prepared/all_target_to4_rankcontrol_bindings_before_forward.tsv.gz',sep='\t')
 with np.load(BASE/'geneformer_forward/all_unique_cell_baseline.npz')as z:B=unit(z['embeddings'])
 idx={(c['study'],c['cell_id']):i for i,c in enumerate(cells)};anchors={}
 for(st,seed),part in M.groupby(['study','seed']):
  pts=sorted(part.patient.unique())
  for held in pts:
   v={}
   for ct in ['FB01','otherFB']:
    means=[]
    for patient in pts:
     if patient==held:continue
     ids=part[part.patient.eq(patient)&part.context.eq(ct)].cell_id;assert len(ids)==32
     means.append(B[[idx[(st,c)]for c in ids]].mean(0))
    v[ct]=unit(np.mean(means,axis=0))
   anchors[(st,int(seed),held)]=v
 raw=pd.read_csv(BASE/'geneformer_readouts/all9candidate4control_cell_responses_and_NA.tsv.gz',sep='\t');keys=['study','context','seed','patient','cell_id','gene','operator'];raw=raw.set_index(keys,verify_integrity=True)
 groups={k:np.asarray(z.sort_values('control_ordinal').reference_query_id,dtype=int)for k,z in binding.groupby(['study','cell_id','target_gene'])};assert all(len(x)==4 for x in groups.values())
 rows=[]
 for op in ['delete_natural','delete_position_preserved']:
  A=np.full((len(Q),768),np.nan,np.float32)
  for p in sorted((src/op).glob('chunk_*.npz')):
   with np.load(p)as z:A[z['query_id']]=z['embedding_after']
  assert np.isfinite(A).all();A=unit(A)
  for v in V.to_dict('records'):
   base={k:v[k]for k in ['study','context','seed','patient','cell_id','capture_id','gene','vocab_measurable','visible']};base['operator']=op
   original=raw.loc[(v['study'],v['context'],v['seed'],v['patient'],v['cell_id'],v['gene'],op)]
   record={**base,'original_raw_target_S':original.S_toward_otherFB_away_FB01,'original_matched2_S_unchanged':original.matched2_centered_S,'n_rank_controls':0,'rank4_mean_reference_S':np.nan,'rank4_median_reference_S':np.nan,'rank4_centered_mean_S':np.nan,'rank4_centered_median_S':np.nan,'status':'NO_TARGET_VISIBLE_OR_MODEL_VOCAB_NA'}
   k=(v['study'],v['cell_id'],v['gene'])
   if k in groups:
    ids=groups[k];before=B[idx[(v['study'],v['cell_id'])]];an=anchors[(v['study'],int(v['seed']),v['patient'])];axis=an['otherFB']-an['FB01'];effects=(A[ids]-before)@axis;assert np.isfinite(effects).all()and pd.notna(record['original_raw_target_S'])
    record.update({'n_rank_controls':4,'rank4_mean_reference_S':float(np.mean(effects)),'rank4_median_reference_S':float(np.median(effects)),'rank4_min_reference_S':float(min(effects)),'rank4_max_reference_S':float(max(effects)),'rank4_centered_mean_S':float(record['original_raw_target_S']-np.mean(effects)),'rank4_centered_median_S':float(record['original_raw_target_S']-np.median(effects)),'status':'COMPLETE_TARGET_WITH4_VISIBLE_SAMECELL_RANK_INPUT_CONTROLS'})
   rows.append(record)
 R=pd.DataFrame(rows);assert len(R)==27648
 R.to_csv(out/'all9_candidate_cell_raw_and_rank4_centered_NA.tsv.gz',sep='\t',index=False,compression='gzip',na_rep='NA')
 patient=[]
 for k,z in R.groupby(['study','context','seed','patient','gene','operator']):
  n=int(z.original_raw_target_S.notna().sum());nr=int(z.rank4_centered_mean_S.notna().sum());assert n==nr
  patient.append({**dict(zip(['study','context','seed','patient','gene','operator'],k)),'n_fixed_selected':len(z),'n_target_visible':n,'n_all4_rank_supported':nr,'patient_raw_mean_S':z.original_raw_target_S.mean(),'patient_rank4_centered_mean_S':z.rank4_centered_mean_S.mean(),'patient_rank4_median_centered_mean_S':z.rank4_centered_median_S.mean(),'support_at_least3_visible_targetcells':n>=3})
 Pp=pd.DataFrame(patient);Pp.to_csv(out/'all9_patient_seed_raw_rank4_means_and_unchanged_support.tsv',sep='\t',index=False,na_rep='NA');sums=[]
 for k,z in Pp.groupby(['study','context','seed','gene','operator']):
  ok=len(z)==4 and z.support_at_least3_visible_targetcells.all()
  sums.append({**dict(zip(['study','context','seed','gene','operator'],k)),'n_supported_patients':int(z.support_at_least3_visible_targetcells.sum()),'all4_patient_raw_S':z.patient_raw_mean_S.mean()if ok else np.nan,'all4_patient_rank4_centered_S':z.patient_rank4_centered_mean_S.mean()if ok else np.nan,'all4_patient_rank4_median_centered_S':z.patient_rank4_median_centered_mean_S.mean()if ok else np.nan,'positive_rank4_patient_means':int(z.patient_rank4_centered_mean_S.gt(0).sum()),'negative_rank4_patient_means':int(z.patient_rank4_centered_mean_S.lt(0).sum())})
 S=pd.DataFrame(sums);S.to_csv(out/'complete_all9_patient_context_seed_operator_rank4_summary.tsv',sep='\t',index=False,na_rep='NA')
 atom(out/'COMPLETE.json',{'status':'COMPLETE_ORIGINAL9_FAMILY_POSTHOC_RANK4_INPUT_CONTROL_SENSITIVITY','n_cell_slots':len(R),'n_patient_seed_slots':len(Pp),'n_group_summaries':len(S),'target_visibility_and_patient_eligibility_unchanged':True,'all_visible_targetqueries_have4_reference_responses':True,'prior_matched2_NA_unchanged':True,'not_expression_matching_biological_null_or_empirical_FDR':True,'targets_previously_seen_posthoc_bounded_extension_not_independent_validation':True,'no_new_target_WT_baseline_or_profile':True,'files':file_receipt(out)})
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--stage',choices=['forward','readouts'],required=True);a=p.parse_args();forward()if a.stage=='forward'else readouts()
