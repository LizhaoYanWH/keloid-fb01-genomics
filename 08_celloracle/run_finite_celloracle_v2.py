"""One cohort / one fixed CellOracle fit / <=3 TF, numeric propagation only. No transitions."""
from finite_common_v1 import *
import argparse,random,traceback,warnings,gzip

def numeric_outputs(base,pred,genes,obs,programs,tf,condition,out,expected_patients,expected_states,failure_status=None,state_outedge_counts=None):
    import numpy as np,pandas as pd
    delta=pred-base if pred is not None else np.full(base.shape,np.nan);sd=base.std(axis=0,ddof=1);gidx={x:i for i,x in enumerate(genes)};rows=[];gene_rows=[]
    groups=[((pat,state),obs[obs.independent_patient_id.eq(pat)&obs.FB_state.astype(str).eq(state)]) for pat in expected_patients for state in expected_states]
    for (pat,state),fr in groups:
        ci=obs.index.get_indexer(fr.index);d=delta[ci];b=base[ci];means=d.mean(axis=0) if len(ci) and pred is not None else np.full(len(genes),np.nan)
        for j,gene in enumerate(genes):gene_rows.append({'patient':pat,'state':str(state),'TF':tf,'condition':condition,'gene_symbol':gene,'n_cells':len(ci),'WT_mean_log_normalized':float(b[:,j].mean()) if len(ci) else None,'WT_SD_all_sampled_cells':float(sd[j]),'mean_delta_log_normalized':float(means[j]) if np.isfinite(means[j]) else None,'mean_delta_WT_SD_standardized':float(means[j]/sd[j]) if sd[j]>0 and np.isfinite(means[j]) else None,'n_positive_cells':int((d[:,j]>0).sum()) if pred is not None else None,'n_negative_cells':int((d[:,j]<0).sum()) if pred is not None else None,'n_zero_cells':int((d[:,j]==0).sum()) if pred is not None else None,'perturbed_gene_itself':gene==tf,'status':failure_status or ('NA_NO_CELLS_IN_ORIGINAL_PATIENT_STATE' if not len(ci) else 'DESCRIPTIVE_PREDICTION')})
        for name,original in programs.items():
            original=list(dict.fromkeys(original))
            for view in ['full_model_covered','without_perturbed_TF']:
                ix=[gidx[x] for x in original if x in gidx and (view=='full_model_covered' or x!=tf)];valid=[j for j in ix if sd[j]>0]
                rows.append({'patient':pat,'state':str(state),'TF':tf,'condition':condition,'program':name,'view':view,'n_cells':len(ci),'original_mapped_symbols':len(original),'model_covered_members':len(ix),'positive_WT_SD_members':len(valid),'mean_delta_log_normalized':float(d[:,ix].mean()) if ix and len(ci) and pred is not None else None,'mean_delta_WT_SD_standardized':float((d[:,valid]/sd[valid]).mean()) if valid and len(ci) and pred is not None else None,'coverage_label':'FULL_MAPPED_MEMBERS' if len(ix)==len(original) else 'MODEL_COVERED_SUBSET_NOT_FULL_SCORE','status':failure_status or ('NA_NO_CELLS_IN_ORIGINAL_PATIENT_STATE' if not len(ci) else ('DESCRIPTIVE_PREDICTION' if ix else 'NA_NO_MODEL_GENE_COVERAGE'))})
    frame=pd.DataFrame(rows);gene_frame=pd.DataFrame(gene_rows)
    # Preserve raw model aggregates, including mechanical TF self-decrease and structural zero.
    # Main program/gene summaries exclude states with no non-self regulatory outedge for this TF.
    for table,value_columns in [(frame,['mean_delta_log_normalized','mean_delta_WT_SD_standardized']),(gene_frame,['mean_delta_log_normalized','mean_delta_WT_SD_standardized','n_positive_cells','n_negative_cells','n_zero_cells'])]:
        table['actual_TF_nonself_outedges_for_state']=table.state.map(state_outedge_counts or {}) if tf else np.nan
        structural=table.state.map(lambda state:bool(tf) and state_outedge_counts is not None and state_outedge_counts.get(str(state),0)==0)&table.n_cells.gt(0)
        table['STRUCTURAL_NO_OUTEDGE']=structural
        raw=table.copy()
        if not failure_status:raw.loc[structural,'status']='STRUCTURAL_NO_OUTEDGE_RAW_MODEL_VALUES_NOT_EFFECT_EVIDENCE'
        raw_name='RAW_model_patient_state_program_deltas.tsv' if table is frame else 'RAW_model_all_gene_patient_state_deltas.tsv.gz'
        raw.to_csv(out/raw_name,sep='\t',index=False,na_rep='NA')
        table.loc[structural,value_columns]=np.nan
        if not failure_status:table.loc[structural,'status']='STRUCTURAL_NO_OUTEDGE'
    frame.to_csv(out/'patient_state_program_deltas.tsv',sep='\t',index=False,na_rep='NA');gene_frame.to_csv(out/'all_gene_patient_state_deltas.tsv.gz',sep='\t',index=False,na_rep='NA')
    paired=[]
    for key,fr in frame.groupby(['patient','TF','condition','program','view'],sort=True,dropna=False):
        fb01=fr[fr.state.eq('FB01')];other=fr[~fr.state.eq('FB01')]
        for field in ['mean_delta_log_normalized','mean_delta_WT_SD_standardized']:
            x=float(fb01[field].iloc[0]) if len(fb01) and pd.notna(fb01[field].iloc[0]) else None;y=float(other[field].mean()) if other[field].notna().any() else None
            paired.append(dict(zip(['patient','TF','condition','program','view'],key),readout=field,FB01=x,otherFB_state_equal_mean=y,FB01_minus_otherFB_state_equal_mean=x-y if x is not None and y is not None else None,otherFB_states_with_value=int(other[field].notna().sum()),otherFB_states_structurally_unsupported=int(other.STRUCTURAL_NO_OUTEDGE.sum()),status='DESCRIPTIVE_NOT_HEALTHY_RESTORATION'))
    pd.DataFrame(paired).to_csv(out/'paired_FB01_otherFB_state_equal_deltas.tsv',sep='\t',index=False,na_rep='NA')
    return frame

def export_coef(coefficients,out,name):
    import numpy as np,pandas as pd
    summaries=[]
    for state,c in coefficients.items():
        need(c.index.is_unique and c.columns.is_unique and c.index.tolist()==c.columns.tolist(),'Coefficient gene axes mismatch');a=c.to_numpy();need(np.isfinite(a).all(),'Nonfinite fit coefficients preserved as failed fit');ii,jj=np.nonzero(a);pd.DataFrame({'source_TF':np.array(c.index)[ii],'target_gene':np.array(c.columns)[jj],'coefficient':a[ii,jj]}).to_csv(out/f'{name}_{state}_nonzero_coefficients.tsv.gz',sep='\t',index=False);summaries.append({'state':str(state),'n_rows':a.shape[0],'n_columns':a.shape[1],'nonzero_n':len(ii),'implicit_zero_n':a.size-len(ii)})
    atom(out/(name+'_coefficient_axis.json'),{'states':summaries,'gene_axis':next(iter(coefficients.values())).columns.tolist(),'orientation':'rows source regulator; columns target. Unlisted entries exactly zero.'})

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--manifest',required=True);ap.add_argument('--root-grant',required=True);ap.add_argument('--study',choices=STUDIES,required=True);ap.add_argument('--out',required=True);a=ap.parse_args();guard(own_GiB=16,output_reserve_MiB=256)
    g=grant(a.root_grant,'APPROVED_FINITE_CELLORACLE_ONE_COHORT',__file__,{'manifest_sha256':sha(a.manifest),'study':a.study});m=load(a.manifest);need(m['status']=='COMPLETE_FROZEN_CELLORACLE_INPUTS_BEFORE_FIT_OR_EFFECTS','Actual input manifest incomplete');need(g['candidate_plan_sha256']==m['candidate_plan']['sha256'],'Candidate plan identity mismatch');p,f,reg=check_plan(checkref(m['candidate_plan']));cohort=next(x for x in m['cohorts'] if x['study']==a.study)
    for r in cohort['files']:checkref(r)
    need(set(g['installed_celloracle_source_sha256'])=={'trajectory/oracle_core.py','trajectory/modified_VelocytoLoom_class.py','trajectory/oracle_GRN.py','network/regression_models.py','network_analysis/links_object.py'},'Exact reviewed source-family binding required')
    for rel,expected in g['installed_celloracle_source_sha256'].items():
        import importlib.util
        spec=importlib.util.find_spec('celloracle');need(spec is not None,'CellOracle missing');need(sha(Path(next(iter(spec.submodule_search_locations)))/rel)==expected,'Installed reviewed CellOracle source changed: '+rel)
    out=Path(a.out);need(not out.exists(),'Run output exists; never duplicate fit');out.mkdir(parents=True);atom(out/'RUNNING.json',{'PID':os.getpid(),'UTC':utc(),'source':rec(__file__),'common':rec(Path(__file__).with_name('finite_common_v1.py')),'manifest':rec(a.manifest),'study':a.study,'host':platform.node()})
    try:
        import numpy as np,pandas as pd,scanpy as sc,anndata as ad,celloracle as co
        need(co.__version__=='0.20.0' and sc.__version__=='1.9.8','Reviewed runtime versions changed');random.seed(SEED);np.random.seed(SEED);src=Path(cohort['files'][0]['path']).parent
        aa=ad.read_h5ad(src/'sampled_linear_CP10K.h5ad');prior=pd.read_parquet(src/'model_promoter_prior.parquet');programs=load(src/'programs.json')['programs'];genes=aa.var_names.tolist();tfs=m['selected_TFs'];need(1<=len(tfs)<=3 and aa.n_obs==cohort['n_sampled_cells'],'Frozen shape/TF mismatch');need(aa.obs.independent_patient_id.isin(f['cohorts'][a.study]['patients']).all(),'Unknown patient in fit');need(aa.obs.FB_state.astype(str).str.startswith('FB').all(),'Original FB states required');need(np.isfinite(aa.X).all() and (aa.X>=0).all(),'Linear input invalid')
        expected_log=np.log1p(aa.X);oracle=co.Oracle();fit_warning=[]
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always');oracle.import_anndata_as_raw_count(aa,cluster_column_name='FB_state',embedding_name=cohort['embedding_key'],transform='natural_log');need(np.array_equal(np.asarray(oracle.adata.layers['normalized_count']),expected_log),'Official natural_log transform is not exact log1p linear CP10K');oracle.import_TF_data(TF_info_matrix=prior);n_pc=min(30,aa.n_obs-1,aa.n_vars-1);oracle.perform_PCA(n_components=n_pc);k=max(2,int(.025*aa.n_obs));oracle.knn_imputation(n_pca_dims=n_pc,k=k,balanced=True,b_sight=8*k,b_maxl=4*k,n_jobs=1);guard(own_GiB=16,start=False);links=oracle.get_links(cluster_name_for_GRN_unit='FB_state',alpha=10,bagging_number=20,model_method='bagging_ridge',verbose_level=1,n_jobs=1);links.filter_links(p=.001,weight='coef_abs',threshold_number=2000);oracle.get_cluster_specific_TFdict_from_Links(links_object=links);oracle.fit_GRN_for_simulation(alpha=10,use_cluster_specific_TFdict=True,verbose_level=1);fit_warning=[str(w.message) for w in caught]
        atom(out/'FIT_DIAGNOSTICS.json',{'warnings':fit_warning,'n_components':n_pc,'KNN_k':k,'balanced':True,'b_sight':8*k,'b_maxl':4*k,'alpha':10,'bagging_number':20,'n_jobs':1,'internal_Ridge_bagging_seed':123,'global_sampling_and_shuffle_seed':SEED,'no_transition_or_quiver':True,'resources':guard(own_GiB=16,start=False)})
        for state,z in links.links_dict.items():z.to_csv(out/f'full_unfiltered_links_{state}.tsv.gz',sep='\t',index=False)
        for state,z in links.filtered_links.items():z.to_csv(out/f'fixed_recipe_filtered_links_{state}.tsv.gz',sep='\t',index=False)
        export_coef(oracle.coef_matrix_per_cluster,out,'real_GRN');base=np.asarray(oracle.adata.layers['imputed_count']);need(base.shape==(aa.n_obs,aa.n_vars) and np.isfinite(base).all(),'Baseline invalid');np.savez_compressed(out/'WT_imputed_log_normalized.npz',expression=base,gene_symbols=np.array(genes),cell_ids=np.array(oracle.adata.obs_names));oracle.adata.obs.to_csv(out/'modeled_cell_metadata.tsv.gz',sep='\t',index_label='cell_id')
        sd=base.std(axis=0,ddof=1);coverage=[];prior_targets=set(prior.gene_short_name);prior_cols=set(prior.columns)-{'peak_id','gene_short_name'}
        for name,original in programs.items():
            for gene in original:
                j=genes.index(gene) if gene in genes else None;coverage.append({'program':name,'gene_symbol':gene,'in_model':j is not None,'in_prior_targets':gene in prior_targets,'in_prior_regulators':gene in prior_cols,'WT_SD_positive':bool(sd[j]>0) if j is not None else False,'actual_active_regulator':gene in oracle.active_regulatory_genes})
        pd.DataFrame(coverage).to_csv(out/'actual_model_program_coverage.tsv',sep='\t',index=False)
        conditions=[('',False,'UNPERTURBED_ENGINEERING_CONTROL')]+[(tf,r,label) for tf in tfs for label,r in [('TF_KO',False),('SHUFFLED_GRN_SAME_TF_KO',True)]];oracle.calculate_randomized_coef_table(random_seed=SEED);export_coef(oracle.coef_matrix_per_cluster_randomized,out,'shuffled_GRN');condition_receipts=[];all_programs=[]
        state_outedges={};outedge_coverage=[]
        for randomized,coefficients in [(False,oracle.coef_matrix_per_cluster),(True,oracle.coef_matrix_per_cluster_randomized)]:
            for tf in tfs:
                counts={}
                for state in cohort['expected_states']:
                    c=coefficients.get(state);n=int(np.count_nonzero(c.loc[tf,c.columns!=tf].to_numpy())) if c is not None and tf in c.index else 0;counts[str(state)]=n
                    n_cells=int(oracle.adata.obs.FB_state.astype(str).eq(str(state)).sum());status='ACTIVE_NONSELF_OUTEDGES' if n else ('STRUCTURAL_NO_OUTEDGE' if n_cells else 'NA_NO_CELLS_IN_ORIGINAL_STATE')
                    outedge_coverage.append({'TF':tf,'state':str(state),'randomized_GRN':randomized,'nonself_outedge_n':n,'state_sampled_cells':n_cells,'status':status})
                state_outedges[(tf,randomized)]=counts
        pd.DataFrame(outedge_coverage).to_csv(out/'TF_state_actual_outedge_coverage.tsv',sep='\t',index=False)
        for number,(tf,randomized,label) in enumerate(conditions):
            guard(own_GiB=16,start=False);dest=out/f'condition_{number:02d}_{tf or "EMPTY"}_{label}';dest.mkdir();diag={'TF':tf,'condition':label,'randomized_GRN':randomized,'n_propagation':3,'clip_delta_X':False,'ignore_warning':False};state_support=state_outedges.get((tf,randomized)) if tf else None
            if tf and (tf not in genes or tf not in oracle.active_regulatory_genes):
                diag['status']='NA_NOT_MEASURED_OR_NO_ACTIVE_REGULATORY_OUTEDGE';all_programs.append(numeric_outputs(base,None,genes,oracle.adata.obs,programs,tf,label,dest,cohort['expected_patients'],cohort['expected_states'],failure_status=diag['status'],state_outedge_counts=state_support));atom(dest/'COMPLETE.json',{'status':diag['status'],'diagnostics':diag,'outputs':files_receipt(dest)});condition_receipts.append(rec(dest/'COMPLETE.json'));continue
            try:
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter('always');oracle.simulate_shift(perturb_condition={tf:0.0} if tf else {},n_propagation=3,use_randomized_GRN=randomized,ignore_warning=False,clip_delta_X=False)
                pred=np.asarray(oracle.adata.layers['simulated_count']);need(pred.shape==base.shape and np.isfinite(pred).all(),'Nonfinite/invalid predicted expression; retained as failed condition');diag['warning_messages']=[str(w.message) for w in caught];diag['max_abs_delta']=float(np.abs(pred-base).max());diag['status']='EXPLORATORY_NUMERIC_PREDICTION_NOT_REAL_KO'
                if not tf:need(np.array_equal(pred,base),'Empty intervention must be exactly zero')
                oracle.evaluate_simulated_gene_distribution_range().to_csv(dest/'official_all_gene_OOD_distribution.tsv.gz',sep='\t',index=True);all_programs.append(numeric_outputs(base,pred,genes,oracle.adata.obs,programs,tf,label,dest,cohort['expected_patients'],cohort['expected_states'],state_outedge_counts=state_support));atom(dest/'DIAGNOSTICS.json',diag);atom(dest/'COMPLETE.json',{'status':'COMPLETE_NUMERIC_CONDITION','diagnostics':diag,'outputs':files_receipt(dest)})
            except ValueError as ex:
                if not tf:raise
                diag.update(status='NA_OFFICIAL_OR_NUMERIC_CONDITION_REJECTED_NOT_BIOLOGICAL_NEGATIVE',error=repr(ex),traceback=traceback.format_exc());all_programs.append(numeric_outputs(base,None,genes,oracle.adata.obs,programs,tf,label,dest,cohort['expected_patients'],cohort['expected_states'],failure_status=diag['status'],state_outedge_counts=state_support));atom(dest/'COMPLETE.json',{'status':diag['status'],'diagnostics':diag,'outputs':files_receipt(dest)})
            condition_receipts.append(rec(dest/'COMPLETE.json'))
        if all_programs:
            allframes=pd.concat(all_programs,ignore_index=True);allframes.to_csv(out/'all_conditions_patient_state_program_deltas.tsv.gz',sep='\t',index=False,na_rep='NA');summaries=[]
            for key,fr in allframes.groupby(['TF','condition','state','program','view'],sort=True,dropna=False):
                for metric in ['mean_delta_log_normalized','mean_delta_WT_SD_standardized']:
                    vals=fr[metric].dropna();summaries.append(dict(zip(['TF','condition','state','program','view'],key),readout=metric,patient_denominator=len(cohort['expected_patients']),patients_with_value=len(vals),patients_NA=len(fr)-len(vals),positive=int((vals>0).sum()),negative=int((vals<0).sum()),zero=int((vals==0).sum()),equal_patient_mean=float(vals.mean()) if len(vals) else None,median=float(vals.median()) if len(vals) else None,status='DESCRIPTIVE_NO_PATIENT_P_NO_MODEL_VOTING'))
            pd.DataFrame(summaries).to_csv(out/'equal_patient_summary_all_signs_NA.tsv',sep='\t',index=False,na_rep='NA')
        atom(out/'COMPLETE.json',{'status':'COMPLETE_FINITE_CELLORACLE_ATTEMPTS_WITH_ALL_NA_AND_OPPOSITE_RESULTS_RETAINED','study':a.study,'manifest':rec(a.manifest),'source':rec(__file__),'common':rec(Path(__file__).with_name('finite_common_v1.py')),'candidate_plan':m['candidate_plan'],'condition_receipts':condition_receipts,'outputs':files_receipt(out),'boundary':'Model-based expression propagation, not real KO, clinical efficacy, healthy restoration, causal genetic mediation or independent patient validation','UTC':utc()})
    except BaseException as ex:atom(out/'FAILED.json',{'status':'FAILED_PRESERVED','error':repr(ex),'traceback':traceback.format_exc(),'UTC':utc()});raise
if __name__=='__main__':main()
