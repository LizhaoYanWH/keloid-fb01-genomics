"""One bounded descriptive temporal check; existing pseudobulks and fixed programmes.

Execute only on the existing AutoDL host. No scDRS refit, genetic association test,
new signature, cell clustering, differential expression, or model selection.
"""
from pathlib import Path
import argparse, csv, hashlib, json, os, platform, resource, shutil, time

OLD = Path('/root/autodl-tmp/Five_Studies_Execution_20261020')
ROOT = Path('/root/five_new_studies_20261020/projects/N02_keloid_genetic_cell_state')
OUT = ROOT / 'normal_repair_timecourse_20261021'
PB = OLD / 'shared_inputs/reference_profiles/GSE241132_author_labels_v1'
DRUG = ROOT / 'mechanisms_20261021/literature_intervention_v1'
PROGRAMS = ('fixed_shared_nonGS_FB01_25', 'original_MAGMA1000_Z_RNA')
EXPECTED = {
    'counts_profiles_by_genes.npz': 'fdb26d3434054bf8fd1753a7494379bd2b25fe6bc0241ef0663cd8cee3aac265',
    'profiles.tsv': '8aa2d1510a295a05cd6b5e2b29e562feabf8eb012dffc94714364578d90b94f1',
    'genes.tsv': 'fb0e69ea0524266e33bab961ac68109786a7d46e9f437820115c67e7e79a2370',
    'reference_profile_receipt.json': 'debd43ad2504c7c3e939f80805a8ffb0b1f42d6d4c239ebeb071d4440c308782',
}
CELLS = {'PWH26': [724,255,920,843], 'PWH27': [521,497,1992,2096], 'PWH28': [907,260,909,2335]}

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def record(p):
    return dict(path=str(p), bytes=p.stat().st_size, sha256=sha(p))

def need(ok, message):
    if not ok:
        raise RuntimeError(message)

def write_json(p, obj):
    temp = p.with_suffix(p.suffix + '.tmp')
    temp.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
    temp.replace(p)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--stage', choices=['freeze','run'], required=True)
    stage = ap.parse_args().stage
    need(platform.system() == 'Linux' and platform.node().startswith('autodl-container-gh1ez37ml9'), 'Existing AutoDL host only')
    for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
        os.environ[name] = '1'
    os.sched_setaffinity(0, [min(os.sched_getaffinity(0))])
    need(min(shutil.disk_usage('/root').free, shutil.disk_usage('/root/autodl-tmp').free) >= 5*2**30 + 32*2**20, 'Both disks require 5GiB plus32MiB headroom')
    cg = Path('/sys/fs/cgroup')
    limit = int((cg/'memory.max').read_text())
    current = int((cg/'memory.current').read_text())
    need(limit-current >= 1.5*2**30, 'Require 1GiB task plan plus0.5GiB free reserve before imports')
    import numpy as np
    import pandas as pd
    from scipy import sparse
    from scipy.stats import rankdata
    OUT.mkdir(parents=True, exist_ok=True)
    sources = [PB/n for n in EXPECTED]
    for p in sources:
        need(sha(p) == EXPECTED[p.name], 'Frozen pseudobulk changed: '+p.name)
    sources += [DRUG/'fixed4_program_members_and_weights.tsv', DRUG/'all_unique_approved_HGNC_gene_axis.tsv', ROOT/'extensions_20261021/external_RNA_v1/FREEZE.json']
    need(sha(sources[4]) == 'a701950c5133a1615c2520b525d7529072ac3e58a7e11289dc838a5170f499d8', 'Original programme membership/weights changed')
    programme = pd.read_csv(sources[4], sep='\t', keep_default_na=False)
    programme = programme.loc[programme.program.isin(PROGRAMS)].copy()
    need(programme.groupby('program').size().to_dict() == {PROGRAMS[0]:25, PROGRAMS[1]:1000}, 'Exactly original25 and1000 rows required')
    ext = json.loads(sources[-1].read_text())
    need(set(programme.loc[programme.program.eq(PROGRAMS[0]),'ensembl_gene_id']) == set(ext['fallback_shared_FB01_nonGS_module']), 'Shared25 identities changed')
    g = programme.loc[programme.program.eq(PROGRAMS[1])]
    need({str(r.source_id):float(r.weight) for r in g.itertuples()} == {str(k):float(v) for k,v in ext['original_GS_weights'].items()}, 'Original1000 Z weights changed')
    meta = pd.read_csv(PB/'profiles.tsv', sep='\t')
    genes = pd.read_csv(PB/'genes.tsv', sep='\t', keep_default_na=False)
    need(genes.source_gene_id.is_unique and len(genes)==36600, 'Original complete36600 gene axis required')
    fb = meta.loc[meta.author_cell_type.eq('Fibroblast')].sort_values(['patient_id','day']).copy()
    need(len(fb)==12 and int(fb.n_cells.sum())==12259, 'All12 author-label FB pseudobulks required')
    for donor, counts in CELLS.items():
        sub = fb.loc[fb.patient_id.eq('GSE241132::'+donor)]
        need(sub.day.tolist()==[0,1,7,30] and sub.n_cells.tolist()==counts, 'Frozen donor/time/cell count mismatch')
    hg = pd.read_csv(sources[5], sep='\t', dtype=str, keep_default_na=False)
    need(hg.ensembl_gene_id.is_unique, 'HGNC background must be unambiguous')
    bg = genes.source_gene_id.isin(hg.ensembl_gene_id).to_numpy()
    ids = genes.loc[bg,'source_gene_id'].tolist()
    need(len(ids)>0 and len(ids)==len(set(ids)), 'Measured HGNC background invalid')
    frozen = dict(status='FROZEN_BEFORE_TEMPORAL_SCORES', source_sha256=sha(__file__), inputs=[record(p) for p in sources],
        original_profiles=fb.to_dict('records'), n_independent_donors=3, normal_times_days=[0,1,7,30], total_FB_cells=12259,
        D0_FB_cells=2152, raw_genes=36600, HGNC_background_genes=len(ids),
        score='For each full pseudobulk, average-tie rank(counts)/G over every measured one-to-one Approved HGNC gene, including zeros; programme mean or original-Z weighted mean over measured fixed members. Monotone library normalization leaves rank unchanged.',
        endpoint='All24 programme-profile scores and every donor difference D1-D0,D7-D0,D30-D0,D30-D7; equal-donor mean/min/max, no P-values/FDR/threshold fitting',
        interpretation='RNA expression proxy; not scDRS, individualPRS, change in genotype/risk, FB01 abundance, causal trajectory, or KL-specific persistence. Author labels/exclusion of author Doublet only, no new QC claims.',
        no_new_algorithms_or_genes=True, ancestry_extension=False)
    if stage=='freeze':
        need(not (OUT/'FREEZE.json').exists(), 'Existing freeze; do not regenerate')
        programme.to_csv(OUT/'frozen_two_programmes_1025_members.tsv', sep='\t', index=False)
        write_json(OUT/'FREEZE.json', frozen)
        print(json.dumps({'status':frozen['status'],'freeze_sha256':sha(OUT/'FREEZE.json'),'background':len(ids)}))
        return
    original = json.loads((OUT/'FREEZE.json').read_text())
    need(original == frozen, 'Current sources/inputs/design differ from frozen record')
    if (OUT/'results/COMPLETE.json').exists():
        print('Existing complete result retained; no rerun'); return
    with (OUT/'ONCE_CLAIM.json').open('x') as handle:
        json.dump(dict(pid=os.getpid(),boot_id=(Path('/proc/sys/kernel/random/boot_id').read_text().strip()),start_ticks=Path('/proc/self/stat').read_text().split()[21],source_sha256=sha(__file__),freeze_sha256=sha(OUT/'FREEZE.json')),handle)
    started=time.monotonic()
    X = sparse.load_npz(PB/'counts_profiles_by_genes.npz').tocsr()
    need(X.shape == (len(meta),len(genes)) and np.isfinite(X.data).all() and (X.data>=0).all() and np.array_equal(X.data,np.rint(X.data)), 'Original integer UMI matrix invalid')
    raw = X[fb.index.to_numpy(),:].toarray()
    need(np.array_equal(raw.sum(axis=1),fb.total_captured_gene_UMI.to_numpy()), 'Full original library UMI sums mismatch')
    arr=raw[:,bg]; ranks=rankdata(arr,axis=1,method='average')/len(ids)
    index={name:i for i,name in enumerate(ids)}
    coverage=[]; records=[]; member_records=[]
    for pr in PROGRAMS:
        members=programme.loc[programme.program.eq(pr)].copy()
        measured=members.ensembl_gene_id.isin(index)
        use=members.loc[measured]
        weights=use.weight.to_numpy(float)
        need(len(use)>0 and np.isfinite(weights).all() and (weights>0).all(), 'Programme has no valid measured weights')
        ix=[index[z] for z in use.ensembl_gene_id]
        scores=(ranks[:,ix]*weights).sum(axis=1)/weights.sum()
        for member in members.itertuples():
            coverage.append(dict(program=pr,source_id=member.source_id,ensembl_gene_id=member.ensembl_gene_id,weight=float(member.weight),assayed=member.ensembl_gene_id in index))
        for j, row in enumerate(fb.to_dict('records')):
            records.append(dict(**row,program=pr,score=float(scores[j]),requested_members=len(members),assayed_members=len(use),original_weight_sum=float(members.weight.sum()),assayed_weight_sum=float(weights.sum()),nonzero_members=int((arr[j,ix]>0).sum()),HGNC_background_genes=len(ids)))
            for member in use.itertuples():
                k=index[member.ensembl_gene_id]
                member_records.append(dict(program=pr,profile_id=row['profile_id'],patient_id=row['patient_id'],day=row['day'],ensembl_gene_id=member.ensembl_gene_id,weight=float(member.weight),UMI=int(arr[j,k]),rank_fraction=float(ranks[j,k])))
    score_df=pd.DataFrame(records); contrasts=[]; summaries=[]
    for pr in PROGRAMS:
        for donor in CELLS:
            q=score_df.loc[score_df.program.eq(pr)&score_df.patient_id.eq('GSE241132::'+donor)].set_index('day')
            for new,old in [(1,0),(7,0),(30,0),(30,7)]:
                contrasts.append(dict(program=pr,patient_id='GSE241132::'+donor,contrast=f'D{new}_minus_D{old}',delta=float(q.loc[new,'score']-q.loc[old,'score'])))
    cd=pd.DataFrame(contrasts)
    for (pr,c), group in cd.groupby(['program','contrast'],sort=False):
        values=group.delta.to_numpy();need(len(values)==3,'Three complete paired donors required')
        summaries.append(dict(program=pr,contrast=c,n_donors=3,mean_delta=float(values.mean()),min_delta=float(values.min()),max_delta=float(values.max()),positive=int((values>0).sum()),negative=int((values<0).sum()),zero=int((values==0).sum()),inference='DESCRIPTIVE_NO_P_NO_FDR'))
    result=OUT/'results';result.mkdir(exist_ok=False)
    for name,df in [('all24_profile_program_scores.tsv',score_df),('all24_donor_program_contrasts.tsv',cd),('all8_equal_donor_descriptive_summaries.tsv',pd.DataFrame(summaries)),('all1025_requested_member_assay_coverage.tsv',pd.DataFrame(coverage)),('all_measured_member_UMI_rank_values.tsv.gz',pd.DataFrame(member_records))]:
        df.to_csv(result/name,sep='\t',index=False)
    # A portable source/UMI audit for this small descriptive check; no matrices downloaded locally.
    pd.DataFrame(dict(profile_id=fb.profile_id,full_raw_UMI=raw.sum(axis=1).astype('int64'),HGNC_UMI=arr.sum(axis=1).astype('int64'),nonzero_HGNC_genes=(arr>0).sum(axis=1))).to_csv(result/'all12_full_library_UMI_QA.tsv',sep='\t',index=False)
    complete=dict(status='COMPLETE_FIXED_TWO_PROGRAMME_NORMAL_TIMECOURSE_DESCRIPTIVE',source_sha256=sha(__file__),freeze_sha256=sha(OUT/'FREEZE.json'),
        independent_donors=3,profiles=12,scores=len(records),paired_donor_differences=len(contrasts),summaries=len(summaries),
        wall_seconds=time.monotonic()-started,VmHWM_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
        files=[record(p) for p in sorted(result.iterdir())],versions=dict(numpy=np.__version__,pandas=pd.__version__),
        caveats=frozen['interpretation'], no_tests_or_new_signatures=True)
    write_json(result/'COMPLETE.json',complete)
    print(json.dumps(complete))

if __name__=='__main__':
    main()
