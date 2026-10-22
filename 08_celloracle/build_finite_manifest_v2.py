"""Freeze sampled RNA, HVGs and program coverage before any new CellOracle fit/effect."""
from finite_common_v1 import *
import argparse,traceback,random

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--candidate-plan',required=True);ap.add_argument('--root-grant',required=True);ap.add_argument('--out',required=True);a=ap.parse_args();guard(own_GiB=4)
    grant(a.root_grant,'APPROVED_CELLORACLE_INPUT_MANIFEST_BUILD',__file__,{'candidate_plan_sha256':sha(a.candidate_plan)})
    p,f,reg=check_plan(a.candidate_plan);out=Path(a.out);need(not out.exists(),'Input output already exists; do not repeat build');out.mkdir(parents=True)
    atom(out/'RUNNING.json',{'PID':os.getpid(),'UTC':utc(),'candidate_plan':rec(a.candidate_plan),'source':rec(__file__)})
    try:
        import h5py,numpy as np,pandas as pd,scanpy as sc,anndata as ad
        from scipy import sparse
        need(sc.__version__=='1.9.8','Reviewed Scanpy version changed');random.seed(SEED);np.random.seed(SEED)
        prior=prior_frame(p['promoter_prior']);TFnames=tf_catalog(p['human_tf_catalog']);tfs=list(dict.fromkeys(c['TF'] for c in p['chains']))
        for tf in tfs:need(tf in TFnames and tf in prior and bool((prior[tf]>0).any()),'Candidate is not cataloged TF with promoter prior outedge: '+tf)
        response=pd.read_csv(checkref(p['S2_response_members']),sep='\t',dtype=str,keep_default_na=False);need({'study','scope','view','sign','gene_symbol','position'}<=set(response),'Actual S2 response schema mismatch')
        hgnc_ref=next(r for r in f['inputs'] if Path(r['path']).name=='hgnc_complete_set.txt');hg=pd.read_csv(checkref(hgnc_ref),sep='\t',dtype=str,keep_default_na=False);need({'status','entrez_id','symbol'}<=set(hg),'Actual HGNC schema mismatch');hg=hg[hg.status.eq('Approved')]
        ecmmodule=reg['fixed_modules']['ECM'];need(ecmmodule['id']=='R-HSA-1474244' and len(set(ecmmodule['canonical_entrez']))==321,'Original ECM321 changed');module_ecm=set(hg.loc[hg.entrez_id.isin(set(ecmmodule['canonical_entrez'])),'symbol']);e2s=f['ensembl_to_symbol'];fixed25=[e2s[x] for x in f['original_fixed25_Ensembl'] if x in e2s];need(len(set(fixed25))==25,'Original fixed25 mapping incomplete')
        # Original regulatory freeze includes 992 unambiguously HGNC-mapped Ensembl entries of the original 1000 Entrez set; no RNA backfill.
        original_gs=[e2s[x] for x in reg['GS_Ensembl'] if x in e2s];baseprogram={'fixed25':fixed25,'fixed25_no_ECM':[x for x in fixed25 if x not in module_ecm],'original1000_model_covered_only':list(dict.fromkeys(original_gs))}
        builds=[]
        for st,c in f['cohorts'].items():
            guard(own_GiB=4,start=False);d=out/st;d.mkdir();inp=[r for r in f['inputs'] if r['path'] in (c['raw'],c['metadata'])];need(len(inp)==2,'Original raw/metadata source receipts absent')
            for r in inp:checkref(r)
            m=pd.read_csv(c['metadata'],sep='\t',dtype=str,keep_default_na=False);need(m.cell_id.is_unique,'Duplicate final cell IDs');req={'cell_id','lineage','FB_state','independent_patient_id','full_input_UMI','capture_id','region_original','region_P2_swap','include_P2_omit'};need(req<=set(m),'Authoritative metadata schema mismatch')
            fb=m[m.lineage.eq('fibroblast')].copy();need(len(fb)==c['n_FB'],'Original final FB count changed');known=fb[fb.independent_patient_id.isin(c['patients'])].copy();need(set(known.independent_patient_id)==set(c['patients']),'Four known patients not covered')
            schema=p['input_schema'][st];need(schema['known_patients']==c['patients'],'Explicit patient order differs from original freeze');states=sorted(fb.FB_state.unique());need(states==sorted(schema['original_FB_states']) and all(x.startswith('FB') and x!='NA' for x in states),'Root actual state schema mismatch')
            selected=[];counts=[]
            for pat in c['patients']:
                for state in states:
                    ids=known.loc[known.independent_patient_id.eq(pat)&known.FB_state.eq(state),'cell_id'].tolist();ranked=sorted(ids,key=lambda x:hashlib.sha256(f'{SEED}|{st}|{pat}|{state}|{x}'.encode()).digest());keep=ranked[:50];selected+=keep;counts.append({'study':st,'independent_patient_id':pat,'FB_state':state,'before':len(ids),'after':len(keep),'cap':50})
            selpos=np.flatnonzero(m.cell_id.isin(selected).to_numpy());obs=m.iloc[selpos].copy().set_index('cell_id');need(len(obs)==len(selected) and len(obs)>=4,'Selected axis mismatch')
            with h5py.File(c['raw']) as h:
                need(hindex(h,'obs').tolist()==m.cell_id.tolist(),'Raw/final cell axis changed');source_genes=harray(h['var']['source_gene_id']);source_genes=np.asarray([x.split('.')[0] for x in source_genes]);need(len(set(source_genes))==len(source_genes),'Duplicate source gene IDs');ix=np.array([i for i,x in enumerate(source_genes) if x in e2s]);symbols=np.array([e2s[source_genes[i]] for i in ix]);need(len(set(symbols))==len(symbols),'Ambiguous HGNC symbol mapping');z=h[c['raw_node']];need(str(z.attrs['encoding-type'])=='csr_matrix','Expected frozen raw integer CSR')
                rows=[];libraries=[]
                for row in selpos:
                    lo,hi=map(int,z['indptr'][row:row+2]);v=z['data'][lo:hi];need(np.isfinite(v).all() and (v>=0).all() and np.equal(v,np.rint(v)).all(),'Raw is not nonnegative integer');X=sparse.csr_matrix((v,z['indices'][lo:hi],np.array([0,hi-lo])),shape=(1,len(source_genes)));X.sum_duplicates();X.eliminate_zeros();total=int(X.sum());need(total==int(m.iloc[row].full_input_UMI) and total>0,'Full original UMI conservation mismatch');libraries.append(total);rows.append(X[:,ix])
                raw=sparse.vstack(rows,format='csr');linear=(sparse.diags(1e4/np.array(libraries,dtype=np.float64))@raw).tocsr();embedding_key=schema['embedding_key'];need(embedding_key in h['obsm'],'Root-named existing embedding missing');em=h['obsm'][embedding_key];need(isinstance(em,h5py.Dataset) and em.shape==(len(m),2),'Actual embedding must be explicit original 2D dataset');embedding=em[selpos,:];need(np.isfinite(embedding).all(),'Nonfinite original embedding')
            var=pd.DataFrame({'source_gene_id':source_genes[ix]},index=pd.Index(symbols,name='gene_symbol'));hv=ad.AnnData(X=linear.copy(),obs=obs.copy(),var=var.copy());sc.pp.log1p(hv);sc.pp.highly_variable_genes(hv,flavor='seurat',n_top_genes=2000,inplace=True);hvg=hv.var.highly_variable.to_numpy()
            # Amendment before any KO/fit: retain every member of all four frozen response families.
            fixed_responses=response[response.study.eq(st)&response.scope.eq('main')&response.view.isin(['non_definition','non_definition_no_ECM'])&response.sign.isin(['Up','Down'])]
            need(len(fixed_responses.groupby(['view','sign']))==4,'All four original response families required')
            for _,rr in fixed_responses.groupby(['view','sign']):need(len(rr)==250 and not rr.gene_symbol.duplicated().any(),'Original full250 response family changed')
            must=set(fixed25+tfs+fixed_responses.gene_symbol.tolist());modelmask=hvg|np.array([x in must for x in symbols]);modelgenes=symbols[modelmask];need(len(modelgenes)>=2,'No model gene axis')
            ainput=ad.AnnData(X=linear[:,modelmask].toarray().astype(np.float64),obs=obs.copy(),var=var.loc[modelgenes].copy());ainput.obs['FB_state']=ainput.obs.FB_state.astype('category');ainput.obsm[embedding_key]=embedding;ainput.uns['normalization_contract']='full original UMI -> linear CP10K; CellOracle import natural_log applies log1p once; no new library normalization';ainput.write_h5ad(d/'sampled_linear_CP10K.h5ad',compression='gzip')
            priorcols=[x for x in modelgenes if x in prior.columns];ps=prior.loc[prior.gene_short_name.isin(modelgenes),['peak_id','gene_short_name']+priorcols].copy();ps.to_parquet(d/'model_promoter_prior.parquet',index=False)
            programs=dict(baseprogram)
            for view in ['non_definition','non_definition_no_ECM']:
                for sign in ['Up','Down']:
                    rr=response[response.study.eq(st)&response.scope.eq('main')&response.view.eq(view)&response.sign.eq(sign)].copy();need(not rr.gene_symbol.duplicated().any(),'Duplicated response members');rr['position']=rr.position.astype(int);programs['S1_'+view+'_'+sign]=rr.sort_values('position').gene_symbol.tolist()
            atom(d/'programs.json',{'programs':programs,'original1000_identifier_count':len(reg['original_GS_Entrez']),'original1000_unambiguous_Ensembl_count':len(reg['GS_Ensembl']),'original1000_mapped_symbols_count':len(baseprogram['original1000_model_covered_only']),'ECM_symbols':sorted(module_ecm),'boundary':'Means of model-covered fixed gene members, not scDRS or complete genetic risk recomputation'})
            modelpos={g:i for i,g in enumerate(modelgenes)};measured=set(symbols[np.asarray(raw.sum(axis=0)).ravel()>0]);prior_targets=set(ps.gene_short_name);coverage=[]
            for name,members in programs.items():
                for gene in members:coverage.append({'program':name,'gene_symbol':gene,'raw_HGNC_assayed':gene in set(symbols),'positive_in_sampled_raw':gene in measured,'in_model':gene in modelpos,'in_prior_targets':gene in prior_targets,'in_prior_regulators':gene in priorcols and bool((ps[gene]>0).any())})
            pd.DataFrame(coverage).to_csv(d/'pre_effect_program_coverage.tsv',sep='\t',index=False);pd.DataFrame(counts).to_csv(d/'sampling_group_counts.tsv',sep='\t',index=False);obs.to_csv(d/'selected_cell_metadata.tsv.gz',sep='\t',index_label='cell_id');pd.DataFrame({'gene_symbol':symbols,'source_gene_id':source_genes[ix],'HVG_scanpy_seurat':hvg,'fixed25':np.isin(symbols,fixed25),'selected_TF':np.isin(symbols,tfs),'in_model':modelmask}).to_csv(d/'whole_gene_selection.tsv.gz',sep='\t',index=False)
            tfcov=[]
            for tf in tfs:tfcov.append({'TF':tf,'cataloged_TF':tf in TFnames,'positive_sampled_raw':tf in measured,'in_model':tf in modelpos,'measured_model_prior_outedges':int((ps[tf]>0).sum()) if tf in ps else 0})
            pd.DataFrame(tfcov).to_csv(d/'candidate_TF_pre_effect_coverage.tsv',sep='\t',index=False)
            builds.append({'study':st,'n_original_final_FB':len(fb),'n_known_FB':len(known),'n_unknown_FB':len(fb)-len(known),'n_sampled_cells':len(obs),'n_original_states':len(states),'expected_patients':c['patients'],'expected_states':states,'HVG_requested':2000,'HVG_actual':int(hvg.sum()),'n_model_genes':len(modelgenes),'embedding_key':embedding_key,'input_sources':inp,'files':files_receipt(d),'selected_TFs':tfs,'TF_coverage':tfcov});del raw,linear,hv,ainput
        atom(out/'COMPLETE.json',{'status':'COMPLETE_FROZEN_CELLORACLE_INPUTS_BEFORE_FIT_OR_EFFECTS','source':rec(__file__),'common':rec(Path(__file__).with_name('finite_common_v1.py')),'candidate_plan':rec(a.candidate_plan),'S2_complete':p['S2_complete'],'selected_TFs':tfs,'seed_sampling_and_shuffle':SEED,'official_internal_Ridge_bagging_seed':123,'cohorts':builds,'normalization':'Integer original counts/full original library UMI *10000 -> linear model input; log1p exactly once inside official import','fit_or_simulate_executed':False,'UTC':utc()})
    except BaseException as ex:atom(out/'FAILED.json',{'status':'FAILED_PRESERVED_NO_FIT','error':repr(ex),'traceback':traceback.format_exc(),'UTC':utc()});raise
if __name__=='__main__':main()
