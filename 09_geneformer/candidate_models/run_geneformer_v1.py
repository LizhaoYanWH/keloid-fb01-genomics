#!/usr/bin/env python
"""Frozen all7candidate+4control Geneformer perturbation, FP32, two operators."""
from runtime_common_v1 import *
import gzip,argparse,collections
import numpy as np,pandas as pd

def unit(x):
 x=np.asarray(x,np.float64);n=np.linalg.norm(x,axis=-1,keepdims=True);assert np.all(n>0);return x/n

def read_prepared():
 p=D/'prepared';c=json.loads((p/'COMPLETE.json').read_text())
 for n,z in c['files'].items():assert sha(p/n)==z['sha256'] and (p/n).stat().st_size==z['bytes']
 with gzip.open(p/'unique_cells_official_tokens.jsonl.gz','rt')as z:cells=list(map(json.loads,z))
 q=pd.read_csv(p/'unique_single_gene_queries_before_any_forward.tsv.gz',sep='\t');return cells,q,c

def forward():
 f,g=approved('gf');claim('gf');cells,Q,prep=read_prepared();out=D/'geneformer_forward';out.mkdir();started=time.time()
 sys.path.insert(0,str(OLD/'shared'));from gpu_admission import acquire_gpu
 lease=acquire_gpu(32,'N02_finite_multimodel_geneformer',trace_path=out/'gpu_admission.jsonl',stop_path=out/'STOP_AFTER_CHUNK')
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
  assert time.time()-started<7200,'TWO_HOUR_ENGINEERING_LIMIT_KEEP_ATOMIC_OUTPUTS'
  assert torch.cuda.max_memory_reserved()<=30*GIB,'OWN_ALLOCATOR_BUDGET'
  with(out/'resource_by_chunk.jsonl').open('a')as z:z.write(json.dumps({'UTC':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'stage':stage,'pid':os.getpid(),'allocated':torch.cuda.memory_allocated(),'reserved':torch.cuda.memory_reserved(),'peak_reserved':torch.cuda.max_memory_reserved(),'RSS':own_rss(),'elapsed':time.time()-started})+'\n')
 index={(c['study'],c['cell_id']):i for i,c in enumerate(cells)}
 seqs=[c['input_ids']for c in cells]
 long=sorted(range(len(cells)),key=lambda i:(-len(seqs[i]),cells[i]['cell_id']))[:32]
 assert len(long)==32 and all(len(seqs[i])==4096 for i in long),'NO_32_FROZEN_LONG4096_PROFILE_NO_ALTERED_CELLS'
 torch.cuda.reset_peak_memory_stats();t=time.time();profile=embed([seqs[i]for i in long]);seconds=time.time()-t
 n4=[seqs[i]for i in long[:4]];fast=profile[:4];hook=embed_sequences(model,n4,batch=2);pos=embed(n4,[list(range(len(s)))for s in n4])
 assert np.allclose(fast,hook,rtol=2e-5,atol=2e-5)and np.allclose(fast,pos,rtol=2e-5,atol=2e-5)
 test=Q.sort_values(['length','study','cell_id','gene'],ascending=[False,True,True,True]).head(4);before=[seqs[index[(r.study,r.cell_id)]]for r in test.itertuples()];after=[s[:int(r.rank)]+s[int(r.rank)+1:]for s,r in zip(before,test.itertuples())]
 fast_b=embed(before);fast_a=embed(after);hook_b=embed_sequences(model,before,batch=2);hook_a=embed_sequences(model,after,batch=2)
 errors=np.linalg.norm((fast_a-fast_b)-(hook_a-hook_b),axis=1)/np.maximum(np.linalg.norm(hook_a-hook_b,axis=1),1e-12)
 assert np.max(errors)<1e-3,'PRESPECIFIED_DELTA_RELATIVE_L2_ENGINEERING_QA'
 np.savez_compressed(out/'profile4_legacy_vs_equivalent_delta.npz',fast_before=fast_b,fast_after=fast_a,legacy_before=hook_b,legacy_after=hook_a,relative_L2=errors)
 atom(out/'PROFILE_COMPLETE.json',{'status':'PASS32_LONG4096_BATCH8_AND_ORIGINAL_DELTA_EQUIVALENCE','prepared_sha256':sha(D/'prepared/COMPLETE.json'),'n_profile':32,'input_length':4096,'batch':8,'forward_seconds':seconds,'queries_per_second':32/seconds,'peak_reserved_bytes':int(torch.cuda.max_memory_reserved()),'delta_relative_L2_max':float(max(errors)),'FP32':True,'TF32':False,'auto_continue_allowed_by_root':g.get('profile_then_continue')is True})
 assert g.get('profile_then_continue')is True,'PROFILE_ONLY_NO_FORWARD_APPROVAL'
 check('after_profile')
 B=np.empty((len(cells),768),np.float32);B[:]=np.nan
 # Reuse completed original335 baseline exactly: same source hashes and cell axes.
 old=OLD/'projects/02_selectivity/versions/lineage_v2/results/31_regional_GF/base';man=json.loads((old/'manifest.json').read_text());oc=json.loads((old/'COMPLETE.json').read_text());spec=man['spec']
 assert oc['signature']==man['signature']and oc['n_cells']==6774 and spec['base_weight_SHA256']==sha(MODEL/'model.safetensors')and spec['prepared_SHA256']==sha(OLD/'projects/02_selectivity/versions/lineage_v2/prepared/GSE335482/cells.jsonl.gz')
 assert spec['representation']=='penultimateCLS' and spec['precision']=='FP32 TF32off' and spec['forward_helper_SHA256']==sha(OLD/'projects/02_selectivity/scripts/forward_GF_equivalent.py')and spec['model_helper_SHA256']==sha(L/'scripts/model_utils.py')
 with gzip.open(OLD/'projects/02_selectivity/versions/lineage_v2/prepared/GSE335482/cells.jsonl.gz','rt')as z:oldcells=list(map(json.loads,z))
 for p in sorted(old.glob('chunk_*.npz')):
  with np.load(p)as z:
   assert str(z['signature'])==man['signature']
   for ci,v in zip(z['cell_index'],z['embeddings']):
    key=('GSE335482',oldcells[int(ci)]['cell_id'])
    if key in index:B[index[key]]=v
 assert np.isfinite(B[[i for i,c in enumerate(cells)if c['study']=='GSE335482']]).all(),'ALL335_NEEDED_BASELINES_MUST_REUSE_VERIFIED_COMPLETE_BLOCKS'
 for begin in range(0,len(cells),128):
  ii=[i for i in range(begin,min(len(cells),begin+128))if not np.isfinite(B[i]).all()]
  if ii:B[ii]=embed([seqs[i]for i in ii])
  check('baseline_'+str(begin))
 assert np.isfinite(B).all();np.savez_compressed(out/'all_unique_cell_baseline.npz',embeddings=B,cell_ids=np.array([c['cell_id']for c in cells]),studies=np.array([c['study']for c in cells]))
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
 atom(out/'COMPLETE.json',{'status':'COMPLETE_FROZEN_GENEFORMER_ALL7_CANDIDATE4_CONTROL_TWO_OPERATOR_FAMILY','n_unique_cells':len(cells),'query_count_per_operator':len(Q),'prepared_sha256':sha(D/'prepared/COMPLETE.json'),'elapsed_seconds':time.time()-started,'peak_reserved_bytes':int(torch.cuda.max_memory_reserved()),'files':file_receipt(out),'operators_are_same_model_sensitivity':True,'not_RNA_KO_or_therapy':True})

def readouts():
 f,g=approved('readouts');claim('readouts');cells,Q,prep=read_prepared();src=D/'geneformer_forward';z=json.loads((src/'COMPLETE.json').read_text())
 for n,x in z['files'].items():assert (src/n).stat().st_size==x['bytes']and sha(src/n)==x['sha256']
 out=D/'geneformer_readouts';out.mkdir();M=pd.read_csv(D/'prepared/all12_patient_balanced32_pool_memberships.tsv',sep='\t');V=pd.read_csv(D/'prepared/all9candidate4control_visibility_and_NA_before_forward.tsv.gz',sep='\t')
 index={(c['study'],c['cell_id']):i for i,c in enumerate(cells)}
 with np.load(src/'all_unique_cell_baseline.npz')as a:B=unit(a['embeddings'])
 anchors={};coverage=[]
 for(st,seed),part in M.groupby(['study','seed']):
  patients=sorted(part.patient.unique());assert len(patients)==4
  for held in patients:
   sides={}
   for ct in ['FB01','otherFB']:
    means=[]
    for pt in patients:
     if pt==held:continue
     zz=part[part.patient.eq(pt)&part.context.eq(ct)];assert len(zz)==32
     means.append(B[[index[(st,c)]for c in zz.cell_id]].mean(0));coverage.append({'study':st,'seed':seed,'held_out_patient':held,'anchor':ct,'included_patient':pt,'n_cells':len(zz)})
    sides[ct]=unit(np.mean(means,axis=0))
   anchors[(st,int(seed),held)]=sides
 pd.DataFrame(coverage).to_csv(out/'all_LODO_anchor_actual_patient_coverage.tsv',sep='\t',index=False)
 records=[]
 for operator in ['delete_natural','delete_position_preserved']:
  A=np.empty((len(Q),768),np.float32);A[:]=np.nan
  for p in sorted((src/operator).glob('chunk_*.npz')):
   with np.load(p)as a:A[a['query_id']]=a['embedding_after']
  assert np.isfinite(A).all();A=unit(A);lookup={(r.study,r.cell_id,r.gene):int(r.query_id)for r in Q.itertuples()}
  for v in V.to_dict('records'):
   base={k:v[k]for k in ['study','context','seed','patient','cell_id','capture_id','gene','role','original_token_length','source_full_library_UMI','vocab_measurable','visible']};base['operator']=operator;key=(v['study'],v['cell_id'],v['gene'])
   if key not in lookup:records.append({**base,'status':'NO_OFFICIAL_VOCAB_OR_NOT_VISIBLE_NA'});continue
   i=index[(v['study'],v['cell_id'])];j=lookup[key];delta=A[j]-B[i];an=anchors[(v['study'],int(v['seed']),v['patient'])]
   records.append({**base,'status':'COMPLETE_MODEL_CELL_RESPONSE','global_embedding_cosine_distance':float(1-A[j]@B[i]),'delta_cos_FB01':float(delta@an['FB01']),'delta_cos_otherFB':float(delta@an['otherFB']),'S_toward_otherFB_away_FB01':float(delta@(an['otherFB']-an['FB01']))})
 R=pd.DataFrame(records);assert len(R)==39936
 key=['study','context','seed','patient','cell_id','operator'];matched={'NEDD4':['DYM','SLC25A39'],'TAB2':['UTP11','QARS1']};metric='S_toward_otherFB_away_FB01'
 R['matched2_centered_S']=np.nan;R['diagnostic4_centered_S']=np.nan;R['n_visible_matched_controls']=0;R['n_visible_diagnostic_controls']=0
 for _,ii in R.groupby(key,dropna=False).groups.items():
  part=R.loc[ii];vals=part.set_index('gene')[metric].to_dict();avail=[vals[g]for g in CONTROLS if pd.notna(vals.get(g))]
  for k in ii:
   gene=R.at[k,'gene'];val=R.at[k,metric]
   own=[vals[g]for g in matched.get(gene,[])if pd.notna(vals.get(g))];R.at[k,'n_visible_matched_controls']=len(own);R.at[k,'n_visible_diagnostic_controls']=len(avail)
   if gene in PROTEINS and pd.notna(val):
    if len(own)==2:R.at[k,'matched2_centered_S']=val-np.mean(own)
    if len(avail)>=2:R.at[k,'diagnostic4_centered_S']=val-np.median(avail)
 R.to_csv(out/'all9candidate4control_cell_responses_and_NA.tsv.gz',sep='\t',index=False,compression='gzip',na_rep='NA')
 keys=['study','context','seed','patient','gene','operator'];rows=[]
 for k,part in R.groupby(keys,dropna=False):
  rows.append({**dict(zip(keys,k)),'n_fixed_selected':len(part),'n_visible':int(part[metric].notna().sum()),'raw_mean_S':part[metric].mean(),'matched2_centered_mean_S':part.matched2_centered_S.mean(),'n_matched2_support':int(part.matched2_centered_S.notna().sum()),'diagnostic4_centered_mean_S':part.diagnostic4_centered_S.mean(),'n_diagnostic4_support':int(part.diagnostic4_centered_S.notna().sum()),'inference_support_ge3_visible':int(part[metric].notna().sum())>=3})
 A=pd.DataFrame(rows);A.to_csv(out/'all9candidate4control_patient_seed_means_and_support.tsv',sep='\t',index=False,na_rep='NA')
 sums=[]
 for k,part in A.groupby(['study','context','seed','gene','operator'],dropna=False):
  ok=len(part)==4 and part.inference_support_ge3_visible.all()
  sums.append({**dict(zip(['study','context','seed','gene','operator'],k)),'n_patient_slots':len(part),'n_supported_patients_ge3_visible':int(part.inference_support_ge3_visible.sum()),'all4_supported_raw_mean_S':part.raw_mean_S.mean()if ok else np.nan,'positive_patient_raw_mean':int(part.raw_mean_S.gt(0).sum()),'negative_patient_raw_mean':int(part.raw_mean_S.lt(0).sum()),'matched2_all4_patient_mean_S':part.matched2_centered_mean_S.mean()if len(part)==4 and part.n_matched2_support.ge(3).all()else np.nan,'diagnostic4_all4_patient_mean_S':part.diagnostic4_centered_mean_S.mean()if len(part)==4 and part.n_diagnostic4_support.ge(3).all()else np.nan})
 pd.DataFrame(sums).to_csv(out/'complete_all4_patient_context_seed_gene_summary.tsv',sep='\t',index=False,na_rep='NA')
 atom(out/'COMPLETE.json',{'status':'COMPLETE_ALL_FROZEN_CANDIDATE_PATIENT_CONTEXT_OPERATOR_READOUTS_NA_RETAINED','n_cell_slots':len(R),'n_patient_seed_slots':len(A),'n_context_seed_gene_operator_slots':len(sums),'no_new_patient_P_values_or_clinical_claim':True,'only_NEDD4_TAB2_have_target_specific_expression_matched2':True,'otherFB_is_not_healthy_or_therapy_outcome':True,'files':file_receipt(out)})
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--stage',choices=['forward','readouts'],required=True);args=a.parse_args();forward()if args.stage=='forward'else readouts()
