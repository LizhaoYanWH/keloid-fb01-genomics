"""One fixed union-of-published-windows sensitivity; reuse immutable official runner helpers."""
from pathlib import Path
import argparse, fcntl, gc, hashlib, importlib.util, json, multiprocessing as mp
import os, resource, sys, time, traceback

ROOT = Path('/root/five_new_studies_20261020/projects/N02_keloid_genetic_cell_state')
EXT = ROOT / 'extensions_20261021'
HERE = EXT / 'regional_exclusion_v1'
BASE_DIR = EXT / 'specificity_locus_v1'
BASE_PATH = BASE_DIR / 'run_specificity_locus_v1.py'
BASE_SHA = '2a4b3a9b119820b3abf21363f5096846f4742214bfd722c2d251733eeafbb570'
BASE_FREEZE_SHA = 'd7d8d74a641078f35baa3b0cfbdee3a23ffff8daf9f9c11bd82d6c1e7c1d83ab'
LOCUS = EXT / 'locus_membership_v1'
LOUT = LOCUS / 'results_v2'
SPEC = LOCUS / 'all22_union_v1/FIXED_ALL22_UNION_DELETION_SPEC.json'
PYTHON = ROOT / 'tools/scdrs104_env/bin/python'
COHORTS = ('GSE220300', 'GSE335482')
CAP = 128 * 2**20


def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda: f.read(2**20), b''):
            h.update(block)
    return h.hexdigest()


def atom(p, j):
    q = p.with_suffix(p.suffix + '.tmp')
    q.write_text(json.dumps(j, indent=2, allow_nan=False) + '\n')
    q.replace(p)


def need(b, s):
    if not b:
        raise RuntimeError(s)


def base_module():
    need(sha(BASE_PATH) == BASE_SHA, 'Pinned original scientific helper source changed')
    spec = importlib.util.spec_from_file_location('n02_pinned_scientific_helpers', BASE_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    # Only the provenance/output-location constant changes in this imported
    # process-local namespace. No physical base source or running B is changed.
    module.HERE = HERE
    return module


def freeze():
    import pandas as pd
    b = base_module()
    need(not (HERE / 'FREEZE.json').exists(), 'No overwrite of static C freeze')
    need(sha(BASE_DIR / 'FREEZE.json') == BASE_FREEZE_SHA, 'Original frozen inputs changed')
    original = json.loads((BASE_DIR / 'FREEZE.json').read_text())
    lc = json.loads((LOUT / 'COMPLETE.json').read_text())
    need(lc['status'] == 'COMPLETE_STATIC_LOCUS_MEMBERSHIP_NO_RNA_EFFECT_OR_ASSOCIATION',
         'Exact successful same-build static membership required')
    need(sha(LOCUS / 'STATIC_FREEZE_v2.json') == lc['freeze_sha256'], 'Static locus freeze changed')
    for r in lc['outputs'].values():
        b.verify(r)
    membership = pd.read_csv(LOUT / 'all1000_GS_gene_location_and_membership.tsv',
                             sep='\t', dtype={'canonical_entrez': str})
    need(len(membership) == 1000 and set(membership.canonical_entrez) ==
         set(original['original_gene_weights']), 'Locus table is not the complete original1000GS')
    need((membership.coordinate_status == 'MAPPED_NCBI37').all(), 'Missing same-build gene body')
    need(all(float(r.original_Z) == original['original_gene_weights'][r.canonical_entrez]
             for _, r in membership.iterrows()), 'Original pure Z weights changed')
    removed = sorted(membership.loc[membership.n_regions > 0, 'canonical_entrez'])
    spec = json.loads(SPEC.read_text())
    need(len(removed) == 61 and lc['n_outside_all22'] == 939 and lc['n_missing_gene_locations'] == 0,
         'Root-reviewed61/939 full mapping differs')
    need(spec['base_COMPLETE_sha256'] == sha(LOUT / 'COMPLETE.json') and
         removed == sorted(spec['remove_original_GS_entrez']) and
         spec['n_remove'] == 61 and spec['n_remain'] == 939,
         'Exact root-declared static all22 union spec differs')
    f = {
        'status': 'FROZEN_ALL22_PUBLISHED_WINDOW_UNION_EXCLUSION_BEFORE_RNA_EFFECTS',
        'UTC': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        'driver_sha256': sha(__file__), 'base_scientific_helper_sha256': BASE_SHA,
        'base_input_freeze_sha256': BASE_FREEZE_SHA,
        'original_gene_weights': original['original_gene_weights'],
        'variants_removed_genes': {'remove_all22_published_window_union': removed},
        'cohorts': list(COHORTS), 'n_original': 1000, 'n_removed': 61, 'n_retained': 939,
        'mapping': 'All22 fixed source GWAS lead-SNP windows; actual GRCh37 same-build '
                   'original NCBI37 gene body +/-10kb overlap. Static full membership v2 '
                   'was independently root-read; no RNA/score/P used for selection.',
        'estimand': 'Remaining original genes outside these22 published windows; '
                    'not deletion of all true causal loci, not a non-risk gene set.',
        'retained_weights': 'Original remaining MAGMA Z weights; never rerank/backfill from RNA',
        'official_recipe': 'Unchanged pinned helpers: original whole-cohort loader, covariates, '
                           'preprocess, mean_var matched1000controls, vs weights, seed20261020; '
                           'raw baseline exact float32 identity check; no cell batch score',
        'family': 'Exactly15 frozen FB-state tests, one BH15 family; all60 state-patient '
                  'contrasts descriptive. Also report combined B30+C15 BH45 as declared '
                  'before allB/C effects, with original B30/C15 q separately retained. '
                  'Original47 primary tests unchanged.',
        'retention': 'Only six-score gzip/control membership/coverage/whole-family descriptions; '
                     'no new Nx1006 disk cache, system-only new directory <=128MiB',
        'interpretation': 'Loss of significance does not establish dependence on these windows; '
                          'membership/weight/coverage, MC magnitude and patient contrasts retained.',
        'inputs': original['inputs'] + [b.receipt(BASE_PATH), b.receipt(BASE_DIR / 'FREEZE.json'),
                  b.receipt(LOCUS / 'STATIC_FREEZE_v2.json'), b.receipt(LOUT / 'COMPLETE.json'),
                  b.receipt(SPEC)] +
                  list(lc['outputs'].values()),
    }
    HERE.mkdir(parents=True, exist_ok=True)
    atom(HERE / 'FREEZE.json', f)
    print(json.dumps({'freeze_sha256': sha(HERE / 'FREEZE.json'), 'removed': removed,
                      'n_removed': 61, 'no_scDRS_fit_performed': True}))


def worker(f, out):
    try:
        import numpy as np
        import pandas as pd
        from statsmodels.stats.multitest import multipletests
        b = base_module()
        rows, patients, coverage = [], [], []
        for st in COHORTS:
            atom(out / 'CURRENT_STAGE.json', {'cohort': st, 'stage': 'official_all22_exclusion'})
            a, cov, original = b.load_rna(st, f)
            r, p, c = b.score_sensitivity(st, a, cov, original, f, out)
            rows += r
            patients += p
            coverage += c
            del a, cov, original
            gc.collect()
        frame = pd.DataFrame(rows)
        need(len(frame) == 15 and len(patients) == 60 and len(coverage) == 2,
             'Complete fixed15/60/two-case family lost')
        finite = frame.assoc_mcp.notna()
        q = multipletests(frame.assoc_mcp.fillna(1).to_numpy(dtype=float), method='fdr_bh')[1]
        frame['sensitivity_family_BH15_q'] = np.where(finite, q, np.nan)
        frame.to_csv(out / 'all15_state_union_exclusion_tests_BH15.tsv', sep='\t',
                     index=False, na_rep='NA')
        b_internal = json.loads((BASE_DIR / 'results_B_v1/INTERNAL_COMPLETE.json').read_text())
        b_key = 'all30_state_sensitivity_tests_BH30.tsv'
        b.verify(b_internal['outputs'][b_key])
        b_table = pd.read_csv(BASE_DIR / 'results_B_v1' / b_key, sep='\t',
                              float_precision='round_trip')
        need(len(b_table) == 30, 'The original complete B30 family changed')
        joint = pd.concat([b_table.assign(source_family='B30'),
                           frame.assign(source_family='C15')], ignore_index=True, sort=False)
        need(len(joint) == 45, 'Combined planned45-family lost slots')
        joint_finite = joint.assoc_mcp.notna()
        joint_q = multipletests(joint.assoc_mcp.fillna(1).to_numpy(dtype=float),
                                method='fdr_bh')[1]
        joint['combined_B30_C15_BH45_q'] = np.where(joint_finite, joint_q, np.nan)
        joint.to_csv(out / 'all45_B30_C15_sensitivity_tests_BH45.tsv', sep='\t',
                     index=False, na_rep='NA')
        pd.DataFrame(patients).to_csv(out / 'all60_state_patient_union_exclusion_contrasts.tsv',
                                     sep='\t', index=False, na_rep='NA')
        atom(out / 'BOTH_OFFICIAL_SCORE_DIAGNOSTICS.json', coverage)
        need(b.written_bytes(HERE) <= CAP, 'Combined C128MiB output cap exceeded')
        atom(out / 'INTERNAL_COMPLETE.json', {
            'status': 'COMPLETE_TWO_ALL22_WINDOW_EXCLUSIONS_NO_PRIMARY_REPLACEMENT',
            'driver_sha256': sha(__file__), 'base_scientific_helper_sha256': BASE_SHA,
            'freeze_sha256': sha(HERE / 'FREEZE.json'),
            'n_score_runs': 2, 'n_state_tests': 15, 'n_patient_descriptions': 60,
            'n_joint_B30_C15_tests': 45,
            'B30_complete_sha256': sha(BASE_DIR / 'results_B_v1/COMPLETE.json'),
            'B30_source_table_sha256': sha(BASE_DIR / 'results_B_v1' / b_key),
            'peak_RSS_bytes': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
            'outputs': {str(p.relative_to(out)): b.receipt(p) for p in out.rglob('*')
                        if p.is_file() and p.name not in ['INTERNAL_COMPLETE.json', 'RUNNING.json']}})
    except BaseException as exc:
        atom(out / 'WORKER_FAILURE_PRESERVED.json', {'error': repr(exc),
                                                   'traceback': traceback.format_exc()})
        raise


def run(args):
    b = base_module()
    approval = Path(args.approval)
    need(sha(approval) == args.approval_sha256, 'Root approval SHA differs')
    j = json.loads(approval.read_text())
    need(j['action'] == 'RUN_N02_ALL22_WINDOW_EXCLUSION_ONCE' and
         j['driver_sha256'] == sha(__file__) and j['base_scientific_helper_sha256'] == BASE_SHA and
         j['freeze_sha256'] == sha(HERE / 'FREEZE.json'), 'Root driver/base/freeze not bound')
    need(sys.executable == str(PYTHON), 'Exact existing official scDRS environment required')
    f = json.loads((HERE / 'FREEZE.json').read_text())
    need(f['driver_sha256'] == sha(__file__), 'Freeze driver differs')
    for r in f['inputs']:
        b.verify(r)
    bc = json.loads((BASE_DIR / 'results_B_v1/COMPLETE.json').read_text())
    br = json.loads((BASE_DIR / 'B_ACTUAL_RETURN.json').read_text())
    need(bc['returncode'] == 0 and bc['stage'] == 'sensitivity' and br['returncode'] == 0 and
         bc['source_sha256'] == BASE_SHA and bc['freeze_sha256'] == BASE_FREEZE_SHA,
         'B must truly complete/return0 before C; never compete/restart it')
    need(j['threads'] == 4 and j['RAM_cap_bytes'] <= 24 * 2**30 and
         j['memory_floor_bytes'] >= 8 * 2**30 and j['max_wall_seconds'] <= 3600 and
         j['output_reserve_bytes'] == CAP, 'Finite C resource scope changed')
    for name in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
        os.environ[name] = '4'
    allowed = sorted(os.sched_getaffinity(0))
    quota, period = Path('/sys/fs/cgroup/cpu.max').read_text().split()
    need(len(allowed) >= 4 and (quota == 'max' or int(quota) / int(period) >= 4), 'CPU quota')
    os.sched_setaffinity(0, allowed[-4:])
    out = HERE / 'results_v1'
    need(not out.exists(), 'Never overwrite partial/complete C results')
    with (HERE / 'ONCE_CLAIM.json').open('x') as w:
        json.dump({**b.pid_identity(os.getpid()), 'approval_sha256': args.approval_sha256}, w)
    with (HERE / 'owner.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        r = b.resources()
        need(r['data_free_bytes'] >= 5 * 2**30 and r['system_free_bytes'] >= 5 * 2**30 + CAP,
             'Both5GiB and system128MiB additional output reservation required')
        need(r['effective_free_bytes'] >= j['RAM_cap_bytes'] + j['memory_floor_bytes'] +
             j['peer_reserved_bytes'], 'Actual memory admission failed; do not auto-retry')
        out.mkdir()
        atom(out / 'RUNNING.json', {**b.pid_identity(os.getpid()), 'approval': j,
                                  'resources': r, 'B_complete_sha256': sha(BASE_DIR / 'results_B_v1/COMPLETE.json'),
                                  'B_return_sha256': sha(BASE_DIR / 'B_ACTUAL_RETURN.json')})
        started = time.monotonic()
        proc = mp.get_context('spawn').Process(target=worker, args=(f, out))
        proc.start()
        atom(out / 'ACTUAL_CHILD_LAUNCH.json', b.pid_identity(proc.pid))
        peak, failed = 0, None
        while proc.is_alive():
            proc.join(5)
            try:
                r = b.resources()
                need(min(r['system_free_bytes'], r['data_free_bytes']) >= 5 * 2**30, 'Disk floor')
                need(r['effective_free_bytes'] >= j['memory_floor_bytes'], 'Memory floor')
                need(b.written_bytes(HERE) <= CAP, 'Combined output cap')
                need(time.monotonic() - started <= j['max_wall_seconds'], 'Finite C wall cap')
                text = Path(f'/proc/{proc.pid}/status').read_text()
                rss = next((int(z.split()[1]) * 1024 for z in text.splitlines()
                            if z.startswith('VmRSS:')), 0)
                peak = max(peak, rss)
                need(rss <= j['RAM_cap_bytes'], 'Owned C RSS cap')
            except FileNotFoundError:
                pass
            except BaseException as exc:
                failed = repr(exc)
                proc.terminate()
                proc.join(30)
                break
        if failed or proc.exitcode != 0:
            atom(out / 'FAILED_PRESERVED_NO_AUTO_RETRY.json', {'error': failed,
                 'child_returncode': proc.exitcode, 'peak_RSS_bytes': peak})
            raise RuntimeError('C failure retained; no automatic retry')
        c = json.loads((out / 'INTERNAL_COMPLETE.json').read_text())
        for r in c['outputs'].values():
            b.verify(r)
        need(sha(BASE_PATH) == BASE_SHA and sha(__file__) == j['driver_sha256'], 'Source mutated')
        atom(out / 'COMPLETE.json', {
            'status': 'COMPLETE_N02_ALL22_WINDOW_EXCLUSION_ROOT_REVIEW_PENDING',
            'driver_sha256': sha(__file__), 'base_scientific_helper_sha256': BASE_SHA,
            'freeze_sha256': sha(HERE / 'FREEZE.json'), 'approval_sha256': sha(approval),
            'internal_complete_sha256': sha(out / 'INTERNAL_COMPLETE.json'),
            'returncode': proc.exitcode, 'elapsed_seconds': time.monotonic() - started,
            'peak_RSS_bytes': max(peak, c['peak_RSS_bytes']),
            'combined_directory_bytes': b.written_bytes(HERE), 'final_resources': b.resources()})
        print(json.dumps({'COMPLETE': str(out / 'COMPLETE.json'),
                          'sha256': sha(out / 'COMPLETE.json')}))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('mode', choices=['freeze', 'run'])
    p.add_argument('--approval')
    p.add_argument('--approval-sha256')
    args = p.parse_args()
    need(ROOT.exists() and Path('/proc').exists(), 'Existing AutoDL only')
    if args.mode == 'freeze':
        freeze()
    else:
        need(args.approval and args.approval_sha256, 'Root scientific grant required')
        run(args)


if __name__ == '__main__':
    main()
