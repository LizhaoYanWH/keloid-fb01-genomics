#!/usr/bin/env python3
"""Symmetric cross-study state markers, lineage diagnostics and common-panel scores."""
import argparse,itertools,json,os,platform,re
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata,spearmanr,t,ttest_1samp

PANELS={'FB':['DCN','LUM','COL1A1','COL1A2','PDGFRA'],'epithelial':['KRT1','KRT10','KRT5','KRT14','EPCAM','DMKN'],'endothelial':['PECAM1','VWF','CLDN5','CDH5','PLVAP'],'mural':['RGS5','ACTA2','MYH11','NOTCH3','CSPG4'],'immune':['PTPRC','LYZ','CD3E','MS4A1','TYROBP']}
def effect(x):
 x=np.asarray(x);n=len(x);beta=x.mean();se=x.std(ddof=1)/np.sqrt(n);null=np.array(list(itertools.product([-1,1],repeat=n)))@x/n
 return {'n_patients':n,'beta':float(beta),'ci_low':float(beta-t.ppf(.975,n-1)*se),'ci_high':float(beta+t.ppf(.975,n-1)*se),'p_value':float(ttest_1samp(x,0).pvalue),'signflip_reference_p':float(np.mean(abs(null)>=abs(beta)-1e-12))}
def main():
 p=argparse.ArgumentParser();p.add_argument('--project-root',type=Path,required=True);p.add_argument('--program-root',type=Path,required=True);a=p.parse_args()
 if platform.system()!='Linux' or os.getenv('FIVE_STUDIES_COMPUTE')!='autodl':raise SystemExit('AutoDL only')
 out=a.project_root/'results/state_correspondence_v1';out.mkdir(exist_ok=True);studies=['GSE335482','GSE220300'];markers={};objs={};purity=[];centroids={};programs={}
 for study in studies:
  m=pd.read_csv(a.project_root/'results/expansion_readouts'/study/'FB_state_top100_annotation_markers.tsv',sep='\t')
  for state,g in m.groupby('group',sort=True):
   keep=[str(s) for s in g.gene_symbol if pd.notna(s) and not re.match(r'^(MT-|RPL|RPS|HB[ABDEGMQZ])',str(s))];markers[(study,state)]=set(keep[:50]);programs[study+'__'+state+'__top50']=set(keep[:50])
  fb=ad.read_h5ad(a.project_root/'results/pseudobulk'/study/'FB_states.h5ad');objs[study]=fb;X=fb.layers['counts'].tocsr();symbols=fb.var.gene_symbol.astype(str)
  det={name:np.asarray((X[:,symbols.isin(panel).to_numpy()]>0).sum(axis=1)).ravel() for name,panel in PANELS.items()}
  for (capture,state),idx in fb.obs.groupby(['capture_id','FB_state'],observed=True).indices.items():
   raw=np.asarray(X[idx,:].sum(axis=0)).ravel();cpm=raw/max(raw.sum(),1)*1e6;rec={'study':study,'capture_id':capture,'patient_id':str(fb.obs.iloc[idx[0]].patient_id),'region_original':str(fb.obs.iloc[idx[0]].region_original),'FB_state':state,'n_cells':len(idx),'mean_mt_pct':float(fb.obs.iloc[idx].pct_counts_mt.mean()),'median_UMI':float(fb.obs.iloc[idx].total_counts.median())}
   for name,panel in PANELS.items():
    mask=symbols.isin(panel).to_numpy();rec[name+'_atleast2_detected_fraction']=float(np.mean(det[name][idx]>=2));rec[name+'_panel_CPM']=float(cpm[mask].sum());rec[name+'_and_FB_atleast2_fraction']=float(np.mean((det[name][idx]>=2)&(det['FB'][idx]>=2)))
   purity.append(rec)
  for state,idx in fb.obs.groupby('FB_state',observed=True).indices.items():
   percap=[]
   for capture,ci in fb.obs.iloc[idx].groupby('capture_id',observed=True).indices.items():
    trueidx=np.array(idx)[ci]
    if len(trueidx)<10:continue
    raw=np.asarray(X[trueidx,:].sum(axis=0)).ravel();percap.append(np.log1p(raw/max(raw.sum(),1)*1e6))
   if percap:centroids[(study,state)]=pd.Series(np.mean(percap,axis=0),index=fb.var.source_gene_id)
 pd.DataFrame(purity).to_csv(out/'all_state_capture_multilineage_diagnostics.tsv',sep='\t',index=False)
 overlap=[]
 for key1,genes1 in markers.items():
  if key1[0]!=studies[0]:continue
  for key2,genes2 in markers.items():
   if key2[0]!=studies[1]:continue
   common=centroids[key1].index.intersection(centroids[key2].index);rho=float(spearmanr(centroids[key1].loc[common],centroids[key2].loc[common]).statistic)
   overlap.append({'study1':key1[0],'state1':key1[1],'study2':key2[0],'state2':key2[1],'top50_marker_Jaccard':len(genes1&genes2)/len(genes1|genes2),'shared_top50_markers':'|'.join(sorted(genes1&genes2)),'common_gene_centroid_spearman':rho,'centroid_weighting':'mean logCPM across available capture/state with>=10cells; annotation only, not patient regional inferential support'})
 pd.DataFrame(overlap).to_csv(out/'all_cross_study_state_correspondence.tsv',sep='\t',index=False)
 (out/'all_symmetric_state_marker_programs.gmt').write_text('\n'.join(name+'\tannotation-derived top50 non-MT/RP/HB markers\t'+'\t'.join(sorted(genes)) for name,genes in programs.items())+'\n')
 # Every state-derived marker program is scored in both studies; no focus-state-only screen.
 effects=[];scores=[];deltas=[]
 members=pd.read_csv(a.program_root/'gene_members.tsv',sep='\t',dtype=str);panel=pd.read_csv(a.project_root/'results/spatial/GSE341250/admission_v1/GSM9902814_full_feature_axis.tsv',sep='\t');spatial_genes=set(panel.loc[panel.feature_type.eq('Gene Expression'),'gene_symbol'])
 for study,fb in objs.items():
  source=a.project_root/'results/pseudobulk'/study;cts=pd.read_csv(source/'all_FB_counts.tsv.gz',sep='\t',index_col=0);meta=pd.read_csv(source/'all_FB_samples.tsv',sep='\t');genes=pd.read_csv(source/'genes.tsv',sep='\t');ranks=np.vstack([rankdata(cts[c],method='average') for c in cts]).T;ranks=(ranks-1)/(len(cts)-1)
  region0,region1=('hypercellular','infiltrating') if study=='GSE335482' else ('center','periphery')
  scored=dict(programs)
  for pid,g in members.groupby('program_id'):scored['XENIUM_covered__'+pid]=set(g.symbol)&spatial_genes
  for name,program in scored.items():
   mask=genes.gene_symbol.isin(program).to_numpy();score=2*ranks[mask,:].mean(axis=0)-1;table=meta.copy();table['study']=study;table['program_id']=name;table['rank_score']=score;table['assayed_members']=int(mask.sum());scores.append(table)
   g=table.loc[table.region_original.isin([region0,region1])&table.patient_id.ne('NA')];w=g.pivot(index='patient_id',columns='region_original',values='rank_score');d=(w[region1]-w[region0]).dropna()
   if len(d)!=4:raise ValueError('Whole FB pairing changed')
   effects.append({'study':study,'program_id':name,'program_source':'all symmetric state annotation programs' if not name.startswith('XENIUM') else 'official program intersected with360gene Xenium assay; score uses complete scRNA gene universe, matched-subset sensitivity only',**effect(d)})
   deltas.extend({'study':study,'program_id':name,'patient_id':patient,'delta':delta} for patient,delta in d.items())
 pd.concat(scores).to_csv(out/'all_state_marker_and_spatial_subset_scores.tsv',sep='\t',index=False);pd.DataFrame(effects).to_csv(out/'all_state_marker_and_spatial_subset_paired_effects.tsv',sep='\t',index=False);pd.DataFrame(deltas).to_csv(out/'all_state_marker_and_spatial_subset_patient_deltas.tsv',sep='\t',index=False)
 (out/'interpretation.txt').write_text('State IDs are not equated across studies. All marker sets compared symmetrically. Annotation-derived continuous scores combine state abundance and expression and are not state-intrinsic effects or independent preregistered validation. Co-expression and multiple markers, QC and capture distribution must be considered together: no single epithelial/endothelial/immune marker proves contamination. Existing frozen main results are unchanged.\n')
 print('state correspondence and diagnostics complete',flush=True)
if __name__=='__main__':main()
