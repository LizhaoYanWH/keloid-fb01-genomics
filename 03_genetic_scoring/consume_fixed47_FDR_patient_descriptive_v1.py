"""Root-approved existing47 MC tests -> one BH; all known-patient/context descriptions, no new model."""
from pathlib import Path
import os,sys,json,hashlib,time,argparse,fcntl,shutil,resource
P=Path('/root/five_new_studies_20261020/projects/N02_keloid_genetic_cell_state')
R=P/'reviews/scDRS_fixed_family_consumption_v1'
RAW=R/'complete_raw_audit_v1'
PREP=P/'reviews/scDRS_production_preparation_v1'
OUT=P/'results/scDRS_fixed47_family_patient_descriptive_v1'
APP=R/'ROOT_APPROVAL_CONSUMER_v1.json'
COHORTS=['GSE241132_D0','GSE220300','GSE335482']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def need(ok,reason):
    if not ok:raise RuntimeError(reason)
def atom(p,j):
    q=p.with_suffix(p.suffix+'.tmp')
    q.write_text(json.dumps(j,indent=2,allow_nan=False)+'\n')
    q.replace(p)
def resource_cut():
    limit=int(Path('/sys/fs/cgroup/memory.max').read_text())
    used=int(Path('/sys/fs/cgroup/memory.current').read_text())
    mem={x.split(':')[0]:int(x.split()[1])*1024 for x in Path('/proc/meminfo').read_text().splitlines()}
    return {'effective_free_bytes':min(limit-used,mem['MemAvailable']),'cgroup_current_bytes':used,'system_free_bytes':shutil.disk_usage('/root').free,'data_free_bytes':shutil.disk_usage('/root/autodl-tmp').free}
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run',action='store_true');ap.add_argument('--approval-sha256');a=ap.parse_args()
    if not a.run:
        print(json.dumps({'status':'PREPARATION_ONLY_NO_FDR','source_SHA256':sha(__file__),'fixed_family':47}));return
    need(sys.platform=='linux','AutoDL only')
    need(sha(APP)==a.approval_sha256,'Root approval SHA differs')
    approval=json.loads(APP.read_text())
    need(approval['approved'] is True and approval['source_SHA256']==sha(__file__),'Exact root source approval missing')
    need(sha(R/'FAMILY_FROZEN_BEFORE_FDR.json')==approval['family_freeze_SHA256'],'Group family changed')
    need(sha(RAW/'COMPLETE.json')==approval['complete_raw_audit_SHA256'],'Full-cohort raw audit changed')
    for c in COHORTS:
        need(sha(P/f'results/scDRS_whole_cohort_v1/{c}/COMPLETE.json')==approval['cohort_COMPLETE_SHA256'][c],'Cohort changed '+c)
    os.sched_setaffinity(0,sorted(os.sched_getaffinity(0))[-2:])
    for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
    lock=(P.parents[1]/'metadata/assess01_CPU_worker.lock').open('a+')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    cut=resource_cut()
    need(min(cut['system_free_bytes'],cut['data_free_bytes'])>=5*2**30,'Both disks5GiB guard')
    need(cut['effective_free_bytes']>=approval['job_RAM_cap_bytes']+approval['effective_memory_floor_bytes']+approval['other_unmaterialized_bytes'],'Actual memory plus root peer reservation guard')
    need(not OUT.exists(),'Never overwrite a prior partial or complete consumer')
    claim=R/'consumer_once_claim_v1.json'
    with claim.open('x') as f:json.dump({'PID':os.getpid(),'approval_SHA256':a.approval_sha256,'source_SHA256':sha(__file__)},f)
    OUT.mkdir()
    atom(OUT/'ROOT_APPROVAL_USED.json',approval)
    start=time.monotonic()
    try:
        import numpy as np,pandas as pd,importlib.metadata as im
        from statsmodels.stats.multitest import multipletests
        raw_receipt=json.loads((RAW/'COMPLETE.json').read_text())
        raw_file=RAW/'all47_unadjusted_official_MC_with_frozen_patient_support.tsv'
        need(sha(raw_file)==raw_receipt['files'][raw_file.name]['SHA256'],'Combined47 MC table changed')
        family=pd.read_csv(R/'all_frozen_group_family_v1.tsv',sep='\t',dtype=str,keep_default_na=False)
        units=pd.read_csv(R/'all_patient_and_context_units_v1.tsv',sep='\t',dtype=str,keep_default_na=False)
        z=pd.read_csv(raw_file,sep='\t',dtype={'group':str})
        need(len(z)==47 and z.family_key.is_unique and list(z.family_key)==list(family.family_key),'Full frozen47/order mismatch')
        p=pd.to_numeric(z.assoc_mcp,errors='coerce').to_numpy()
        finite=np.isfinite(p)
        need(((p[finite]>=0)&(p[finite]<=1)).all(),'Invalid finite MC P')
        # Every fixed row stays in multiplicity accounting. Missing input P/q stay NA externally.
        q=multipletests(np.where(finite,p,1.),method='fdr_bh')[1]
        q[~finite]=np.nan
        z['FDR_BH_all47_fixed_MC']=q
        z['MC_family_size']=47
        z['MC_P_type']='control-gene-set association, not independent-patient disease test'
        descriptions=[]
        marker_rows=[]
        references=[]
        for c in COHORTS:
            src=P/f'results/scDRS_whole_cohort_v1/{c}'
            producer=json.loads((src/'COMPLETE.json').read_text())
            fn='patient_capture_lineage_FBstate_P2_views_descriptive.tsv'
            need(sha(src/fn)==producer['outputs'][fn]['sha256'],'Patient cache-derived summary changed '+c)
            t=pd.read_csv(src/fn,sep='\t',dtype={'independent_patient_id':str,'capture_id':str,'FB_state':str},keep_default_na=False)
            t=t.loc[t.view.eq('original_GEO')].copy()
            t['n_cells']=pd.to_numeric(t.n_cells)
            t['mean_norm_score']=pd.to_numeric(t.mean_norm_score)
            def describe(part):
                n=int(part.n_cells.sum())
                return n,float(np.average(part.mean_norm_score,weights=part.n_cells)) if n else None
            for row in units.loc[units.cohort.eq(c)].to_dict('records'):
                base=t.loc[t.independent_patient_id.eq(row['independent_patient_id'])]
                if not row['independent_patient_id']:
                    base=base.loc[base.capture_id.eq(row.get('context_capture_id',''))]
                group=base.loc[base.lineage.eq(row['group'])] if row['group_family']=='lineage' else base.loc[base.lineage.eq('fibroblast')&base.FB_state.eq(row['group'])]
                n,mean=describe(group)
                need(n==int(row['n_cells']),'Frozen patient/context cell count changed '+row['family_key'])
                other_n=0;other_mean=None;delta=None
                if row['group_family']=='FB_state':
                    other_n,other_mean=describe(base.loc[base.lineage.eq('fibroblast')&base.FB_state.ne(row['group'])])
                    if mean is not None and other_mean is not None:delta=mean-other_mean
                descriptions.append({**row,'mean_norm_score':mean,'other_FB_cells_same_patient_or_context':other_n,'other_FB_mean_norm_score':other_mean,'delta_vs_other_FB_same_patient_or_context':delta,'description_only_no_new_P_test':True,'P2_existing_views_preserved_in_original_producer':True})
            if c!='GSE241132_D0':
                ref=Path('/root/autodl-tmp/Five_Studies_Execution_20261020/projects/01_regional/versions/lineage_v2/results')/f'expansion_readouts/{c}/FB_state_top100_annotation_markers.tsv'
                m=pd.read_csv(ref,sep='\t',dtype={'group':str},keep_default_na=False)
                references.append({'use':'Existing state annotation only, never new gene-set scoring','path':str(ref),'SHA256':sha(ref)})
                for state,g in m.groupby('group',sort=True):
                    marker_rows.append({'cohort':c,'group_family':'FB_state','group':state,'existing_annotation_top10':','.join(g.gene_symbol.astype(str).iloc[:10]),'source_path':str(ref),'source_SHA256':sha(ref)})
        desc=pd.DataFrame(descriptions)
        desc.to_csv(OUT/'all_frozen_patient_and_context_descriptions_no_new_tests.tsv',sep='\t',index=False)
        markers=pd.DataFrame(marker_rows)
        z=z.merge(markers.drop(columns=['source_path','source_SHA256']),on=['cohort','group_family','group'],how='left',validate='1:1')
        z.to_csv(OUT/'all47_official_MC_one_BH_family_with_existing_annotation.tsv',sep='\t',index=False)
        markers.to_csv(OUT/'existing_all15_state_annotation_reference_only.tsv',sep='\t',index=False)
        state=[]
        known=desc.loc[desc.group_family.eq('FB_state')&desc.independent_patient_id.ne('')]
        for key,t in known.groupby('family_key',sort=False):
            d=pd.to_numeric(t.delta_vs_other_FB_same_patient_or_context,errors='coerce')
            state.append({'family_key':key,'cohort':t.cohort.iloc[0],'FB_state':t.group.iloc[0],'known_patient_slots':len(t),'patients_with_state_and_other_FB':int(d.notna().sum()),'patients_positive_delta':int(d.gt(0).sum()),'patients_negative_delta':int(d.lt(0).sum()),'patient_ids':'|'.join(t.independent_patient_id),'deltas_in_same_patient_order':'|'.join('NA' if pd.isna(v) else format(v,'.9g') for v in d),'descriptive_only_not_population_genetic_replication':True})
        pd.DataFrame(state).to_csv(OUT/'all15_FB_states_known_patient_consistency_descriptive.tsv',sep='\t',index=False)
        cross=Path('/root/autodl-tmp/Five_Studies_Execution_20261020/projects/01_regional/versions/lineage_v2/results/state_correspondence_v1/all_cross_study_state_correspondence.tsv')
        references.append({'use':'Reuse all54 prior RNA state correspondence, state IDs not equated across studies','path':str(cross),'SHA256':sha(cross),'rows':len(pd.read_csv(cross,sep='\t'))})
        summaries=[]
        for (c,kind),g in z.groupby(['cohort','group_family'],sort=True):
            summaries.append({'cohort':c,'family':kind,'fixed_groups':len(g),'raw_MC_lt05':int(g.assoc_mcp.lt(.05).sum()),'one47family_FDR_lt05':int(g.FDR_BH_all47_fixed_MC.lt(.05).sum()),'missing_P':int(g.assoc_mcp.isna().sum()),'missing_q':int(g.FDR_BH_all47_fixed_MC.isna().sum())})
        pd.DataFrame(summaries).to_csv(OUT/'complete_family_summary.tsv',sep='\t',index=False)
        peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024
        need(peak<=approval['job_RAM_cap_bytes'],'Consumer owned RSS exceeded approval')
        r={'status':'COMPLETE_EXISTING47_MC_ONE_BH_AND_ALL_PATIENT_CONTEXT_DESCRIPTIONS_NO_NEW_MODEL','source_SHA256':sha(__file__),'approval_SHA256':a.approval_sha256,'family_freeze_SHA256':approval['family_freeze_SHA256'],'complete_raw_audit_SHA256':approval['complete_raw_audit_SHA256'],'family_size':47,'tested_P_count':int(finite.sum()),'missing_P_rows_preserved':int((~finite).sum()),'NA_multiplicity_convention':'Internal NA-as1 only for conservative47 BH accounting; external inputP/outputq stayNA','no_recalculation_of_scores_or_MC':True,'no_new_signature_cluster_or_population_P_test':True,'all15_state_patient_comparators_predefined_other_FB_within_same_patient':'Descriptive difference of cache-derived normalized scores; not genetic causal effect','full_family_summary':summaries,'existing_annotation_reference_manifest':references,'versions':{k:im.version(k) for k in ['numpy','pandas','statsmodels']},'elapsed_seconds':time.monotonic()-start,'peak_RSS_bytes':peak,'resources_before':cut,'resources_after':resource_cut(),'files':{p.name:{'bytes':p.stat().st_size,'SHA256':sha(p)} for p in OUT.iterdir() if p.is_file()}}
        atom(OUT/'COMPLETE.json',r)
        print(json.dumps({'COMPLETE_SHA256':sha(OUT/'COMPLETE.json'),'summary':summaries,'elapsed_seconds':r['elapsed_seconds'],'peak_RSS_bytes':peak}))
    except BaseException as e:
        atom(OUT/'FAILED_PRESERVED_NO_AUTOMATIC_RETRY.json',{'error':repr(e),'source_SHA256':sha(__file__)})
        raise
if __name__=='__main__':main()
