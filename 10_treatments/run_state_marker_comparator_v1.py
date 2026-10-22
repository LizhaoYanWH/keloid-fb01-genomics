"""Cached native-state marker panels matched without looking at drug effects.

No new RNA fit/state/DE/P; no drug-effect/low-expression-P selection of members.
"""
from pathlib import Path
import argparse, importlib.util, json, os, signal, subprocess, sys, time
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[key] = '2'
BASE = Path('/root/five_new_studies_20261020/projects/N02_keloid_genetic_cell_state/mechanisms_20261021/regulatory_axis_v1')
HERE = BASE / 'state_marker_comparator_v1'; UP = BASE.parent / 'literature_intervention_v1'
R = HERE / 'score_fixed_state_comparators_v1.R'; RS = Path('/root/five_studies_envs/r/bin/Rscript')
BASE_SHA = 'ed261586a4ba3d9a882c6ae1b534691584c9dcd9cd51d34cf87216495f325012'
BASE_FREEZE = 'a40845ba6d763f53191cf417149c6d94e0b515c3b3b826f784be2302a6cca446'
UP_COMPLETE = '889f7ccfe8d31a16d39a3ffd75909921cd58ae2b89744db34f78584987961dcd'
ORIGINAL_R_SHA = 'a29c2289881da9bd7be959d7e0dff53d0e6fb0d625898e92b5551bf249f9992c'
TABLE = UP / 'results_v1/all7_contrasts_full_HGNC_t_logFC_for_fixed_prior_footprints.tsv.gz'
CAP = 16*2**20
spec = importlib.util.spec_from_file_location('frozen_regulatory_context',BASE/'run_regulatory_axis_v1.py')
b = importlib.util.module_from_spec(spec); spec.loader.exec_module(b)
b.need(b.sha(BASE/'run_regulatory_axis_v1.py')==BASE_SHA and b.sha(BASE/'FREEZE.json')==BASE_FREEZE,'Original regulatory source/freeze changed')


def freeze():
    import pandas as pd, importlib.metadata
    b.need(not (HERE/'FREEZE.json').exists(),'Existing freeze; no overwrite'); HERE.mkdir(parents=True,exist_ok=True)
    old=json.loads((BASE/'FREEZE.json').read_text());up=json.loads((UP/'FREEZE.json').read_text())
    b.need(b.sha(UP/'results_v1/COMPLETE.json')==UP_COMPLETE,'Upstream real intervention changed')
    complete=json.loads((UP/'results_v1/COMPLETE.json').read_text())
    b.need(b.sha(TABLE)==complete['all_outputs'][TABLE.name]['sha256'],'Need complete original allHGNC gene admission axis')
    b.need(b.sha(UP/'run_paired_intervention_analysis_v1.R')==ORIGINAL_R_SHA,'Original input loaders changed')
    ext=b.P/'extensions_20261021/external_RNA_v1/FREEZE.json';ef=json.loads(ext.read_text())
    original25=ef['fallback_shared_FB01_nonGS_module'];b.need(len(original25)==25,'Original shared25 required')
    h,e2s,_,_=b.network_and_genes();ban=set(old['GS_Ensembl'])|set(old['candidates'].values())
    axis=pd.read_csv(TABLE,sep='\t',dtype=str,usecols=['dataset','gene_id','HGNC_symbol'],keep_default_na=False)
    assay={ds:sorted(set(z.gene_id)) for ds,z in axis.groupby('dataset')};b.need(set(assay)=={'GSE232079','GSE232435','GSE7890'},'Original three assay axes differ')
    common=set.intersection(*(set(v) for v in assay.values()))
    cohorts={};pools={};ledger=[];files=[R,RS,BASE/'run_regulatory_axis_v1.py',BASE/'FREEZE.json',ext,TABLE,UP/'FREEZE.json',UP/'results_v1/COMPLETE.json',UP/'results_v1/allPS_complete.tsv.gz',UP/'run_paired_intervention_analysis_v1.R',b.HGNC]
    # Freeze all already-frozen original source/annotation/sample inputs, no new source.
    for r in up['inputs']:
        if '/sources/' in r['path'] or Path(r['path']).name in ('all56_GSM_sample_metadata.tsv','all_unique_approved_HGNC_gene_axis.tsv','GPL570_fixed_all_probe_to_HGNC_Ensembl.tsv.gz'):
            files.append(Path(r['path']))
    for st in ('GSE220300','GSE335482'):
        root=b.O/'projects/01_regional/versions/lineage_v2/results';pb=root/'pseudobulk'/st
        marker=root/'expansion_readouts'/st/'FB_state_top100_annotation_markers.tsv'
        m=pd.read_csv(marker,sep='\t',dtype=str,keep_default_na=False)
        genes=pd.read_csv(pb/'genes.tsv',sep='\t',dtype=str,keep_default_na=False);geneaxis=set(genes.source_gene_id)
        states=sorted(m.group.unique());b.need(len(states)==(9 if st=='GSE220300' else 6) and all(m.groupby('group').size()==100),'Original complete native top100 family changed')
        cohorts[st]={'marker':str(marker),'counts':str(pb/'all_FB_counts.tsv.gz'),'samples':str(pb/'all_FB_samples.tsv'),'genes':str(pb/'genes.tsv'),'patients':old['cohorts'][st]['patients'],'states':states}
        for state in states:
            rows=m[m.group.eq(state)];key=st+'::'+state;chosen=[]
            for rank,row in enumerate(rows.itertuples(),1):
                import re
                g=row.names;symbol=e2s.get(g,'')
                reason='HOUSEKEEPING_EXCLUDED' if re.match(r'^(MT-|RPL|RPS|HB[ABDEGMQZ])',row.gene_symbol) else 'ORIGINAL_GS9_EXCLUDED' if g in ban else 'HGNC_NONUNIQUE_OR_UNRESOLVED_NA' if not symbol else 'DISCOVERY_ASSAY_ABSENT_NA' if g not in geneaxis else 'NOT_ASSAYED_ALL3_DRUG_AXES_NA' if g not in common else 'ELIGIBLE_FIXED_TOP100_POOL'
                ledger.append({'study':st,'state':state,'program':key,'original_marker_rank':rank,'Ensembl':g,'HGNC_symbol':symbol,'status':reason})
                if reason=='ELIGIBLE_FIXED_TOP100_POOL':chosen.append(g)
            pools[key]=sorted(set(chosen))
        files += [marker,pb/'all_FB_counts.tsv.gz',pb/'all_FB_samples.tsv',pb/'genes.tsv',pb/'pseudobulk_receipt.json']
    files += [b.O/'projects/01_regional/scripts/15_expansion_readouts.py',b.O/'projects/01_regional/scripts/05_freeze_and_pseudobulk.py',b.O/'projects/01_regional/scripts/04_atlas_annotation.py']
    f={'status':'STATIC_NATIVE_MARKER_AND_ASSAY_AXES_FROZEN_NO_MATCHING_OR_DRUG_SCORES','source_sha256':b.sha(__file__),'R_source_sha256':b.sha(R),
       'original_shared25':original25,'GS9_excluded_Ensembl':sorted(ban),'original_ECM321':old['fixed_modules']['ECM'],'Ensembl_to_HGNC':e2s,
       'assay_gene_axes':assay,'common_all3_assayed_Ensembl':sorted(common),'candidate_pools':pools,'cohorts':cohorts,'all1500_marker_qualification':ledger,
       'matching':'For each of15 native states independently: exactly cached top100, remove originalGS9/housekeeping/nonuniqueHGNC/outside original3 measured axes. Do not inspect t/P/sign or expand to top101. Each documented patient sums all own capture raw allFB pseudobulks, normalizes by its FULL original library toCP10K; average four patients equally, then log2(meanCP10K+.1). Match the25 original sharedFB01 anchor expression levels to25 distinct candidate genes by scipy linear_sum_assignment minimum total absolute log-expression distance. Lexically sorted anchor/candidate inputs, deterministic installed solver. No caliper/late dropping; report all distances. If<25 or missing anchor/background, retain whole programNA. FB01 itself included; identifiers never aligned acrossstudies.',
       'estimator':'Same original unfiltered measured-HGNC background TMMlogCPM/prior2 or log2 provider normalizedHC probe-average, perassay per-sample percentile rank mean. Exact frozen original input loader functions reused via AST whitelist; original25 scores must reproduce<1e-12. Preserve all16 programs/all7legal contrasts and all22 source deltas; t/Welch95interval descriptive only, no newP/limma/CAMERA/FRY.',
       'boundary':'Post-discovery descriptive comparator. Drug effects were previously seen, but membership algorithm never reads them. Original25 inHC may be24assayed while new panels25: retain unequal coverage. Background expression matching does not match network correlation, ECMmembership, specificity or biological importance. Rare/native states and shared gene panels are not independent replicated states; cultured bulk has no in-tissue FB01 labels. No clinical/state-selective efficacy inference.',
       'versions':{k:importlib.metadata.version(k) for k in ('numpy','scipy','pandas')},'resource':{'RSS_GiB':4,'CPU':[6,7],'output_cap_bytes':CAP},
       'inputs':[b.receipt(p) for p in sorted(set(files))]}
    b.atom(HERE/'FREEZE.json',f)
    print(json.dumps({'source_sha256':b.sha(__file__),'R_source_sha256':b.sha(R),'freeze_sha256':b.sha(HERE/'FREEZE.json'),'eligible_pool_sizes':{k:len(v) for k,v in pools.items()}}),flush=True)


def worker(f):
    import numpy as np,pandas as pd
    from scipy.optimize import linear_sum_assignment
    os.sched_setaffinity(0,{6,7});out=HERE/'results';b.need(not out.exists(),'Existing output; do not overwrite');out.mkdir()
    for r in f['inputs']:
        p=Path(r['path']);b.need(p.stat().st_size==r['bytes'] and b.sha(p)==r['sha256'],'Frozen source changed: '+str(p))
    members=[{'program':'ORIGINAL_SHARED_FB01_25','ensembl_gene_id':g,'weight':1.} for g in f['original_shared25']]
    catalog=[];pairings=[];baseline=[];support=[];panels={'ORIGINAL_SHARED_FB01_25':set(f['original_shared25'])}
    for st,c in f['cohorts'].items():
        x=pd.read_csv(c['counts'],sep='\t',index_col=0);m=pd.read_csv(c['samples'],sep='\t',dtype=str,keep_default_na=False)
        b.need(x.index.is_unique and x.columns.is_unique and set(x.columns)==set(m.pseudobulk_id),'Original aggregate axis mismatch')
        b.need(np.isfinite(x.to_numpy()).all() and (x.to_numpy()>=0).all() and (x.to_numpy()==np.rint(x.to_numpy())).all(),'Original integer aggregate counts required')
        mm=m.set_index('pseudobulk_id').loc[x.columns];b.need(np.array_equal(x.sum(0).to_numpy(np.int64),mm.total_gene_UMI.to_numpy(np.int64)),'Full original library differs')
        patient=[]
        for p in c['patients']:
            selected=mm.index[mm.patient_id.eq(p)];b.need(len(selected)>0,'Documented patient lacks pseudobulk')
            raw=x.loc[:,selected].sum(1);lib=int(raw.sum());b.need(lib>0,'Zero patient library')
            patient.append(raw*1e4/lib);support.append({'study':st,'patient':p,'captures':len(selected),'full_gene_library':lib,'FB_cells':int(mm.loc[selected,'n_FB_cells'].astype(int).sum()),'capture_ids':'|'.join(selected)})
        b.need(len(patient)==4,'Exactly four independent patients required')
        mean=pd.concat(patient,axis=1).mean(1);metric=np.log2(mean+.1);anchors=sorted(f['original_shared25']);b.need(set(anchors).issubset(metric.index),'Original25 anchor absent; preserve before matching')
        for g in sorted(set(anchors)|set().union(*(set(f['candidate_pools'][st+'::'+s]) for s in c['states']))):
            baseline.append({'study':st,'Ensembl':g,'HGNC_symbol':f['Ensembl_to_HGNC'].get(g,''),'mean_patient_CP10K':float(mean[g]),'matching_log2_meanCP10K_plus0p1':float(metric[g])})
        for state in c['states']:
            key=st+'::'+state;candidate=f['candidate_pools'][key];status='MATCHED_FIXED25' if len(candidate)>=25 else 'INSUFFICIENT_ORIGINAL_TOP100_CANDIDATES_NA_NO_BACKFILL'
            record={'study':st,'state':state,'program':key,'eligible_pool_genes':len(candidate),'status':status,'fixed_members':0,'median_abs_log2_expression_distance':np.nan,'max_abs_log2_expression_distance':np.nan}
            if len(candidate)>=25:
                cost=np.abs(metric.loc[anchors].to_numpy()[:,None]-metric.loc[candidate].to_numpy()[None,:]);row,col=linear_sum_assignment(cost)
                chosen=[candidate[j] for j in col];b.need(len(row)==25 and len(set(chosen))==25,'One-to-one match failure')
                panels[key]=set(chosen);record.update({'fixed_members':25,'median_abs_log2_expression_distance':float(np.median(cost[row,col])),'max_abs_log2_expression_distance':float(np.max(cost[row,col]))})
                for i,j in zip(row,col):
                    target=anchors[i];gene=candidate[j];pairings.append({'study':st,'state':state,'program':key,'original25_anchor':target,'selected_marker':gene,'anchor_metric':float(metric[target]),'selected_metric':float(metric[gene]),'abs_log2_expression_distance':float(cost[i,j])})
                    members.append({'program':key,'ensembl_gene_id':gene,'weight':1.})
            else:panels[key]=set()
            catalog.append(record)
    pd.DataFrame(catalog).to_csv(out/'all15_native_state_matching_catalog.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(members).to_csv(out/'fixed16_program_members.tsv',sep='\t',index=False)
    pd.DataFrame(pairings).to_csv(out/'all15_program_expression_matching_pairs.tsv',sep='\t',index=False)
    pd.DataFrame(baseline).to_csv(out/'original_patient_equal_FB_expression_matching_background.tsv.gz',sep='\t',index=False)
    pd.DataFrame(support).to_csv(out/'all8patient_expression_background_units.tsv',sep='\t',index=False)
    pd.DataFrame(f['all1500_marker_qualification']).to_csv(out/'all1500_cached_marker_qualification.tsv',sep='\t',index=False)
    h,_,_,_=b.network_and_genes();ecm=set(h.loc[h.entrez_id.isin(set(map(str,f['original_ECM321']['canonical_entrez']))),'ensembl_gene_id'])-{''}
    overlap=[]
    for name,a in panels.items():
        for other,z in panels.items():overlap.append({'program':name,'other_program':other,'n_program':len(a),'n_other':len(z),'intersection':len(a&z),'union':len(a|z),'Jaccard':len(a&z)/len(a|z) if a|z else np.nan,'program_ECM321_overlap':len(a&ecm),'program_original25_overlap':len(a&panels['ORIGINAL_SHARED_FB01_25'])})
    pd.DataFrame(overlap).to_csv(out/'all16_program_overlap_ECM_original25_ledger.tsv',sep='\t',index=False,na_rep='NA')
    with (out/'R_descriptive.log').open('x') as log:
        child=subprocess.Popen([str(RS),'--vanilla',str(R),str(HERE)],stdout=log,stderr=subprocess.STDOUT)
        b.atom(out/'R_RUNNING.json',{'worker':b.identity(os.getpid()),'R_child':b.identity(child.pid)});rc=child.wait()
    b.need(rc==0 and (out/'R_DESCRIPTIVE_COMPLETE.json').exists(),'Descriptive R failed; preserve no retry')
    outputs=[b.receipt(p) for p in sorted(out.iterdir()) if p.is_file()];b.need(sum(r['bytes'] for r in outputs)<=CAP,'16MiB output cap exceeded')
    b.atom(out/'COMPLETE.json',{'status':'COMPLETE_ALL15_NATIVE_STATE_MARKER_DESCRIPTIVE_COMPARATORS_ROOT_REVIEW_PENDING','source_sha256':b.sha(__file__),'R_source_sha256':b.sha(R),'freeze_sha256':b.sha(HERE/'FREEZE.json'),'outputs':outputs,'native_states':15,'matched_states':sum(v['fixed_members']==25 for v in catalog),'inference_boundary':f['boundary']})


def tree_identities(root):
    found={root};identities=[]
    for _ in range(5):
        for p in Path('/proc').glob('[0-9]*/stat'):
            try:
                t=p.read_text();v=t[t.rfind(')')+2:].split()
                if int(v[1]) in found:found.add(int(p.parent.name))
            except (FileNotFoundError,ProcessLookupError):pass
    for pid in found:
        try:identities.append(b.identity(pid))
        except FileNotFoundError:pass
    return identities


def run(approval):
    f=json.loads((HERE/'FREEZE.json').read_text());g=json.loads(Path(approval).read_text())
    b.need(g['source_sha256']==f['source_sha256']==b.sha(__file__) and g['R_source_sha256']==b.sha(R)==f['R_source_sha256'] and g['freeze_sha256']==b.sha(HERE/'FREEZE.json'),'Root/source/freeze disagreement')
    b.need(not (HERE/'CLAIM.json').exists() and not (HERE/'results').exists(),'Existing claim/output; no duplicate')
    r=b.resources();peer=float(g.get('peer_unmaterialized_GiB',20));b.need(r['effective_free_bytes']>=(4+peer+8)*2**30 and r['system_free_bytes']>=5*2**30+CAP and r['data_free_bytes']>=5*2**30,'Actual resource HOLD')
    fd=os.open(HERE/'CLAIM.json',os.O_WRONLY|os.O_CREAT|os.O_EXCL);os.write(fd,json.dumps({'supervisor':b.identity(os.getpid()),'grant':b.receipt(approval),'resources':r}).encode());os.close(fd)
    begin=time.monotonic();peak=0
    with (HERE/'worker.log').open('x') as log:
        child=subprocess.Popen([sys.executable,str(Path(__file__).resolve()),'_worker'],stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        b.atom(HERE/'RUNNING.json',{'supervisor':b.identity(os.getpid()),'child':b.identity(child.pid)})
        while child.poll() is None:
            own=tree_identities(child.pid);rss=0
            for ident in own:
                try:
                    kv={v.split(':',1)[0]:v.split(':',1)[1].strip() for v in Path(f'/proc/{ident["pid"]}/status').read_text().splitlines() if ':' in v};rss+=int(kv.get('VmRSS','0 kB').split()[0])*1024
                except FileNotFoundError:pass
            peak=max(peak,rss);q=b.resources()
            if rss>4*2**30 or q['effective_free_bytes']<8*2**30 or min(q['system_free_bytes'],q['data_free_bytes'])<5*2**30 or time.monotonic()-begin>g.get('timeout_seconds',600):
                os.killpg(child.pid,signal.SIGTERM);child.wait();b.atom(HERE/'RESOURCE_STOP.json',{'owned_tree':own,'RSS':rss,'resources':q})
            time.sleep(2)
        rc=child.wait()
    b.atom(HERE/'RETURN.json',{'child_pid':child.pid,'returncode':rc,'seconds':time.monotonic()-begin,'peak_owned_tree_RSS_bytes':peak,'COMPLETE_exists':(HERE/'results/COMPLETE.json').exists(),'resources':b.resources()})
    b.need(rc==0 and (HERE/'results/COMPLETE.json').exists(),'Failure/incomplete; preserve no retry')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['freeze','run','_worker']);ap.add_argument('--approval');a=ap.parse_args()
    if a.stage=='freeze':freeze()
    elif a.stage=='run':run(a.approval)
    else:worker(json.loads((HERE/'FREEZE.json').read_text()))
