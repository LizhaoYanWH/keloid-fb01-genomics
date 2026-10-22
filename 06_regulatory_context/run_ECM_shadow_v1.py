"""One post-discovery ECM-footprint exclusion view; no raw RNA or limma rerun.

Uses complete original patient gene deltas and seven complete moderated-t vectors.
Static freeze must precede root-authorized once-only AutoDL execution.
"""
from pathlib import Path
import argparse, importlib.util, json, os, resource, subprocess, sys, time, traceback
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[key] = '2'
BASE = Path('/root/five_new_studies_20261020/projects/N02_keloid_genetic_cell_state/mechanisms_20261021/regulatory_axis_v1')
HERE = BASE / 'ECM_shadow_v1'
UP = BASE.parent / 'literature_intervention_v1'
TABLE = UP / 'results_v1/all7_contrasts_full_HGNC_t_logFC_for_fixed_prior_footprints.tsv.gz'
BASE_SHA = 'ed261586a4ba3d9a882c6ae1b534691584c9dcd9cd51d34cf87216495f325012'
BASE_FREEZE = 'a40845ba6d763f53191cf417149c6d94e0b515c3b3b826f784be2302a6cca446'
BASE_COMPLETE = 'b33fb138d962b46089a498e20004703fcd61d22d3d0eb51fbb84914404e59db1'
UP_COMPLETE = '889f7ccfe8d31a16d39a3ffd75909921cd58ae2b89744db34f78584987961dcd'
TABLE_SHA = 'e3086073af2c85b1da7c34bc1fbcbe0f258a8d009f1cf4c03f0f7f84b7354794'
COHORTS = ('GSE220300', 'GSE335482', 'GSE163973')
CONTRASTS = ['Corin_minus_DMSO_primary4', 'TGFbeta_minus_Control', 'HS276_TGFbeta_minus_TGFbeta',
             'HS276_TGFbeta_minus_Control', 'HC_minus_NoHC_keloid5', 'HC_minus_NoHC_normalScar4',
             'HC_response_keloid_minus_normalScar']
VIEW = 'no_GS9_mapping_OR_original_ECM321_targets'
CAP = 48 * 2**20
spec = importlib.util.spec_from_file_location('frozen_regulatory_context', BASE / 'run_regulatory_axis_v1.py')
b = importlib.util.module_from_spec(spec); spec.loader.exec_module(b)
b.need(b.sha(BASE / 'run_regulatory_axis_v1.py') == BASE_SHA, 'Original regulatory source changed')
b.need(b.sha(BASE / 'FREEZE.json') == BASE_FREEZE, 'Original regulatory freeze changed')


def freeze():
    import pandas as pd
    b.need(not (HERE / 'FREEZE.json').exists(), 'Existing freeze; no overwrite')
    HERE.mkdir(parents=True, exist_ok=True)
    old = json.loads((BASE / 'FREEZE.json').read_text())
    complete = BASE / 'results/COMPLETE.json'; upstream = UP / 'results_v1/COMPLETE.json'
    b.need(b.sha(complete) == BASE_COMPLETE and b.sha(upstream) == UP_COMPLETE, 'Original complete producer changed')
    cp = json.loads(complete.read_text()); up = json.loads(upstream.read_text())
    b.need(up['returncode'] == 0 and b.sha(TABLE) == TABLE_SHA == up['all_outputs'][TABLE.name]['sha256'], 'Seven real intervention results are not complete and bound')
    b.need(json.loads((UP / 'FREEZE.json').read_text())['contrasts7'] == CONTRASTS, 'Original seven contrasts differ')
    h, _, net, _ = b.network_and_genes()
    ecm = old['fixed_modules']['ECM']
    b.need(ecm['id'] == 'R-HSA-1474244' and len(set(ecm['canonical_entrez'])) == 321, 'Original fixed ECM321 changed')
    ids = set(map(str, ecm['canonical_entrez']))
    ledger = []
    for eid in sorted(ids):
        z = h[h.entrez_id.eq(eid)]
        if z.empty:
            ledger.append({'original_ECM_Entrez': eid, 'Approved_HGNC': '', 'status': 'NO_APPROVED_HGNC_MAPPING_NA'})
        else:
            for r in z.itertuples():
                ledger.append({'original_ECM_Entrez': eid, 'Approved_HGNC': r.symbol,
                               'status': 'APPROVED_EXACT_ENTREZ' if len(z) == 1 else 'MULTIPLE_APPROVED_SYMBOLS_RETAIN_ALL'})
    mapped = sorted({v['Approved_HGNC'] for v in ledger if v['Approved_HGNC']})
    original = set(old['excluded_symbols_no_direct_overlap']); union = sorted(original | set(mapped))
    metadata = BASE / 'results/all14_patient_group_support.tsv'
    receipts = {r['path']: r for r in cp['outputs']}; files = [complete, upstream, TABLE, metadata, UP / 'FREEZE.json', BASE / 'run_regulatory_axis_v1.py', BASE / 'FREEZE.json']
    delta_paths = {}
    for st in COHORTS:
        p = BASE / f'results/{st}_patient_FB01_otherFB_log2CP10K_delta.tsv.gz'
        b.need(str(p) in receipts and b.sha(p) == receipts[str(p)]['sha256'], 'Patient delta not complete-producer bound')
        columns = pd.read_csv(p, sep='\t', nrows=0).columns.tolist()
        b.need(len(columns[1:]) == len(set(columns[1:])), 'Patient delta gene header duplicate')
        delta_paths[st] = str(p); files.append(p)
    files += [Path(r['path']) for r in old['inputs'] if 'CollecTRI_human_interactions' in r['path'] or 'hgnc_complete_set' in r['path'] or 'site-packages/decoupler/' in r['path'] or 'specificity_locus_v1/FREEZE.json' in r['path']]
    b.need(len(net.source.unique()) == 1185 and len(original) == 1432, 'Original full family or exclusions differ')
    f = {'status': 'STATIC_INPUTS_MEMBERSHIP_FROZEN_NO_NEW_TF_EFFECTS', 'source_sha256': b.sha(__file__),
         'base_source_sha256': BASE_SHA, 'base_freeze_sha256': BASE_FREEZE, 'view': VIEW,
         'network_sources': old['network_sources'], 'axes': old['axes'], 'tmin': old['tmin'],
         'cohorts': delta_paths, 'patient_support': str(metadata), 'contrasts': CONTRASTS,
         'original_excluded_symbols': sorted(original), 'original_ECM321': ecm,
         'ECM_Entrez_Approved_HGNC_ledger': ledger, 'ECM_symbols': mapped, 'excluded_symbols_union': union,
         'design': 'One additional post-discovery sensitivity. Same cached patient FB01-otherFB full-library log2(CP10K+.1) and complete seven moderated-t vectors. Exclude original GS1000/9/mapping targets UNION exact-HGNC-mapped original ECM321; never backfill or rerun RNA/limma. All original1185 TFs, signs, missingness and native statistic retained.',
         'inference': 'ULM native is within-vector allTF BHq; MLM native rawP is explicitly BH-corrected over all finite returned TFs. These are gene-network enrichment statistics, not biological-patient P. No new patient-level association inference or pooled family. Retention after removing direct ECM targets does not establish TF causality, cell-state uniqueness or drug efficacy; loss may reflect target coverage/power. Different views have different footprint backgrounds and t scales: no effect-retention percentage.',
         'resource': {'own_RSS_GiB': 4, 'CPU_affinity': [6, 7], 'output_cap_bytes': CAP},
         'inputs': [b.receipt(p) for p in sorted(set(files))]}
    b.atom(HERE / 'FREEZE.json', f)
    print(json.dumps({'source_sha256': b.sha(__file__), 'freeze_sha256': b.sha(HERE / 'FREEZE.json'),
                      'original_exclusions': len(original), 'original_ECM321_mapped_symbols': len(mapped), 'union_exclusions': len(union)}), flush=True)


def worker(f):
    import numpy as np, pandas as pd, decoupler as dc
    from scipy.stats import false_discovery_control
    os.sched_setaffinity(0, {6, 7})
    out = HERE / 'results'; b.need(not out.exists(), 'Existing output; no duplicate')
    out.mkdir()
    for r in f['inputs']:
        b.need(Path(r['path']).stat().st_size == r['bytes'] and b.sha(r['path']) == r['sha256'], 'Frozen input/source changed: '+r['path'])
    h, _, net, _ = b.network_and_genes(); approved = set(h.symbol)
    excluded = set(f['excluded_symbols_union']); ecm = set(f['ECM_symbols']); old_excluded = set(f['original_excluded_symbols'])
    axis = {tf: name for name, tfs in f['axes'].items() for tf in tfs}
    support = pd.read_csv(f['patient_support'], sep='\t', dtype={'patient': str})
    contrasts = []; admission = []
    for st, p in f['cohorts'].items():
        x = pd.read_csv(p, sep='\t', index_col=0)
        b.need(x.index.is_unique and x.columns.is_unique and set(x.columns).issubset(approved), 'Cached delta needs original unique Approved HGNC')
        b.need(np.isfinite(x.to_numpy(float)).all(), 'Cached patient gene delta nonfinite')
        meta = support[support.study.eq(st)].set_index('patient')
        b.need(set(x.index) == set(meta.index[(meta.FB01_n.gt(0) & meta.otherFB_n.gt(0))]), 'Complete patient support mismatch')
        contrasts.append(('state', st, x, meta))
    tab = pd.read_csv(TABLE, sep='\t', keep_default_na=False)
    b.need(set(tab.contrast) == set(CONTRASTS), 'Exactly original seven complete contrasts required')
    for contrast in CONTRASTS:
        z = tab[tab.contrast.eq(contrast)].copy()
        b.need(z.HGNC_symbol.is_unique and set(z.HGNC_symbol).issubset(approved), 'No t-vector remapping or averaging')
        v = pd.to_numeric(z.t, errors='coerce'); valid = np.isfinite(v.to_numpy(float))
        x = pd.DataFrame([dict(zip(z.loc[valid, 'HGNC_symbol'], v[valid]))], index=[contrast])
        contrasts.append(('drug', str(z.dataset.unique().item()), x, None))
        for row, vt in zip(z.itertuples(), v):
            admission.append({'contrast': contrast, 'gene': row.HGNC_symbol, 'source_status': row.test_status,
                              'finite_t': bool(np.isfinite(vt)), 'original_exclusion': row.HGNC_symbol in old_excluded,
                              'original_ECM321': row.HGNC_symbol in ecm, 'excluded_union': row.HGNC_symbol in excluded})
    rows = []; cov = []; errors = []
    for kind, context, x, meta in contrasts:
        use = x.loc[:, ~x.columns.isin(excluded)]
        n = net[net.target.isin(use.columns)].copy(); counts = n.groupby('source').target.nunique()
        for tf in f['network_sources']:
            targets = set(net.loc[net.source.eq(tf), 'target']); observed = set(x.columns)
            cov.append({'kind': kind, 'context': context, 'view': VIEW, 'TF': tf, 'prior_axis': axis.get(tf, ''),
                        'finite_gene_background_before': x.shape[1], 'finite_gene_background_after': use.shape[1],
                        'full_regulon_targets': len(targets), 'assayed_targets_before': len(targets & observed),
                        'assayed_targets_after': int(counts.get(tf, 0)),
                        'original_exclusion_assayed_targets': len(targets & observed & old_excluded),
                        'additional_ECM_assayed_targets_removed': len((targets & observed & ecm) - old_excluded),
                        'status': 'ADMITTED_TMIN10' if counts.get(tf, 0) >= f['tmin'] else 'INSUFFICIENT_TARGETS_NA'})
        for method in ('ulm', 'mlm'):
            try:
                score, native = getattr(dc.mt, method)(use, net=n, tmin=f['tmin'], empty=False, verbose=True)
                b.need(score.index.equals(x.index), 'Official sample order differs')
                q = native.copy()
                if method == 'mlm':
                    for sample in q.index:
                        v = native.loc[sample].to_numpy(float); valid = np.isfinite(v); vals = np.full(v.shape, np.nan)
                        vals[valid] = false_discovery_control(v[valid], method='bh'); q.loc[sample] = vals
                error = None
            except (np.linalg.LinAlgError, ValueError) as e:
                error = str(e); errors.append({'kind': kind, 'context': context, 'method': method, 'error': error, 'traceback': traceback.format_exc()})
            for sample in x.index:
                m = meta.loc[sample] if kind == 'state' else None
                for tf in f['network_sources']:
                    exists = error is None and tf in score.columns
                    rows.append({'kind': kind, 'context': context, 'sample': sample, 'view': VIEW, 'method': method, 'TF': tf,
                                 'prior_axis': axis.get(tf, ''), 'condition': m.condition if m is not None else 'public_culture_intervention',
                                 'FB01_n': int(m.FB01_n) if m is not None else np.nan, 'otherFB_n': int(m.otherFB_n) if m is not None else np.nan,
                                 'enrichment_t': score.loc[sample, tf] if exists else np.nan,
                                 'official_returned_value': native.loc[sample, tf] if exists else np.nan,
                                 'official_return_kind': 'BH_q' if method == 'ulm' else 'RAW_P',
                                 'network_BH_q_NOT_patient_P': q.loc[sample, tf] if exists else np.nan,
                                 'status': 'OFFICIAL_COMPLETE' if exists else ('OFFICIAL_TECHNICAL_NA_NO_SELECTIVE_TF_REMOVAL' if error else 'INSUFFICIENT_TARGETS_NA')})
            print(json.dumps({'kind': kind, 'context': context, 'method': method, 'samples': len(x), 'target_genes': use.shape[1]}), flush=True)
    a = pd.DataFrame(rows)
    b.need(len(a) == (13+7)*1185*2, 'Complete family row count mismatch')
    a.to_csv(out / 'all3cohorts7contrasts_all1185TF_new_ECM_exclusion_view.tsv.gz', sep='\t', index=False, na_rep='NA')
    a[a.prior_axis.ne('')].to_csv(out / 'fixed15_all3cohorts7contrasts_new_ECM_exclusion_view.tsv', sep='\t', index=False, na_rep='NA')
    pd.DataFrame(cov).to_csv(out / 'allTF_background_coverage_ECM_exclusion.tsv.gz', sep='\t', index=False)
    pd.DataFrame(admission).to_csv(out / 'all7_complete_input_gene_admission.tsv.gz', sep='\t', index=False)
    pd.DataFrame(f['ECM_Entrez_Approved_HGNC_ledger']).to_csv(out / 'original321_ECM_entrez_HGNC_ledger.tsv', sep='\t', index=False)
    summary = []
    for key, z in a[a.kind.eq('state') & a.condition.eq('keloid')].groupby(['context', 'method', 'TF']):
        v = z.enrichment_t.to_numpy(float); finite = np.isfinite(v)
        summary.append(dict(zip(['study', 'method', 'TF'], key), view=VIEW, prior_axis=axis.get(key[2], ''),
                            n_patients=len(v), finite_patients=int(finite.sum()), positive=int((v[finite]>0).sum()), negative=int((v[finite]<0).sum()),
                            mean_enrichment_t=float(v[finite].mean()) if finite.any() else np.nan, patient_inference_P=np.nan))
    pd.DataFrame(summary).to_csv(out / 'allTF_keloid_patient_direction_summary.tsv.gz', sep='\t', index=False, na_rep='NA')
    b.atom(out / 'OFFICIAL_METHOD_TECHNICAL_FAILURES.json', errors)
    outputs = [b.receipt(p) for p in sorted(out.iterdir()) if p.is_file()]
    b.need(sum(r['bytes'] for r in outputs) <= CAP, '48MiB declared output exceeded')
    b.atom(out / 'COMPLETE.json', {'status': 'COMPLETE_ONE_ADDITIONAL_ECM_EXCLUSION_FOOTPRINT_VIEW_ROOT_REVIEW_PENDING',
           'source_sha256': b.sha(__file__), 'freeze_sha256': b.sha(HERE / 'FREEZE.json'), 'outputs': outputs,
           'all_TF_rows': len(a), 'technical_failures': errors, 'peak_RSS_bytes': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
           'inference_boundary': f['inference']})


def run(approval):
    grant = json.loads(Path(approval).read_text()); f = json.loads((HERE / 'FREEZE.json').read_text())
    b.need(grant['source_sha256'] == b.sha(__file__) == f['source_sha256'] and grant['freeze_sha256'] == b.sha(HERE / 'FREEZE.json'), 'Root/source/freeze mismatch')
    b.need(not (HERE / 'CLAIM.json').exists() and not (HERE / 'results').exists(), 'Existing claim/output; no relaunch')
    r = b.resources(); own = 4; peer = float(grant.get('peer_unmaterialized_GiB', 22))
    b.need(r['effective_free_bytes'] >= (own+peer+8)*2**30 and r['system_free_bytes'] >= 5*2**30+CAP and r['data_free_bytes'] >= 5*2**30, 'Actual RAM/disk HOLD')
    fd = os.open(HERE / 'CLAIM.json', os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    os.write(fd, json.dumps({'supervisor': b.identity(os.getpid()), 'grant': b.receipt(approval), 'resources': r}).encode()); os.close(fd)
    begin = time.monotonic()
    with (HERE / 'worker.log').open('x') as log:
        child = subprocess.Popen([sys.executable, str(Path(__file__).resolve()), '_worker'], stdout=log, stderr=subprocess.STDOUT)
        b.atom(HERE / 'RUNNING.json', {'supervisor': b.identity(os.getpid()), 'child': b.identity(child.pid), 'UTC': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())})
        while child.poll() is None:
            try:
                kv = {v.split(':',1)[0]: v.split(':',1)[1].strip() for v in Path(f'/proc/{child.pid}/status').read_text().splitlines() if ':' in v}
                rss = int(kv.get('VmRSS','0 kB').split()[0])*1024; q = b.resources()
                if rss > own*2**30 or q['effective_free_bytes'] < 8*2**30 or min(q['system_free_bytes'],q['data_free_bytes']) < 5*2**30 or time.monotonic()-begin > grant.get('timeout_seconds',600):
                    child.terminate(); child.wait(); b.atom(HERE / 'RESOURCE_STOP.json', {'own_child_pid': child.pid, 'RSS': rss, 'resources': q})
            except FileNotFoundError:
                pass
            time.sleep(3)
        rc = child.wait()
    b.atom(HERE / 'RETURN.json', {'child_pid': child.pid, 'returncode': rc, 'seconds': time.monotonic()-begin, 'COMPLETE_exists': (HERE / 'results/COMPLETE.json').exists(), 'resources': b.resources()})
    b.need(rc == 0 and (HERE / 'results/COMPLETE.json').exists(), 'Incomplete; preserve all source/partial outputs')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('stage', choices=['freeze','run','_worker']); ap.add_argument('--approval'); args = ap.parse_args()
    if args.stage == 'freeze': freeze()
    elif args.stage == 'run': run(args.approval)
    else: worker(json.loads((HERE / 'FREEZE.json').read_text()))
