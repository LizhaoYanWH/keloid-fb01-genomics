"""Finite N02 specificity/region diagnostics; official scDRS sensitivity, never a new primary GS.

Run only on the existing AutoDL host. `freeze` consumes static memberships and
complete source receipts, not RNA effects. `run` needs a root approval binding
the source and FREEZE SHA. Old results and manuscripts are strictly read-only.
"""
from pathlib import Path
import argparse, fcntl, gc, gzip, hashlib, json, multiprocessing as mp, os
import resource, shutil, sys, time, traceback

ROOT = Path('/root/five_new_studies_20261020/projects/N02_keloid_genetic_cell_state')
OLD = Path('/root/autodl-tmp/Five_Studies_Execution_20261020')
HERE = ROOT / 'extensions_20261021/specificity_locus_v1'
PREP = ROOT / 'reviews/scDRS_production_preparation_v1'
RNA = PREP / 'RNA_freeze_v1'
GS = PREP / 'official_MAGMA_GS_v1'
PRIMARY = ROOT / 'results/scDRS_whole_cohort_v1'
TABLE = ROOT / 'manuscript/full_candidate_v1_20261020/tables'
PYTHON = ROOT / 'tools/scdrs104_env/bin/python'
SITE = ROOT / 'tools/scdrs104_env/lib/python3.12/site-packages'
COHORTS = ('GSE220300', 'GSE335482')
MODULES = {
    'ECM': ('R-HSA-1474244', 'Extracellular matrix organization'),
    'TGFbeta': ('R-HSA-170834', 'Signaling by TGF-beta Receptor Complex'),
    'CellCycle': ('R-HSA-1640170', 'Cell Cycle'),
}
OUTPUT_CAP = 256 * 2**20


def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda: f.read(2**20), b''):
            h.update(b)
    return h.hexdigest()


def atom(p, j):
    p = Path(p)
    q = p.with_suffix(p.suffix + '.tmp')
    q.write_text(json.dumps(j, indent=2, ensure_ascii=False, allow_nan=False) + '\n')
    q.replace(p)


def need(b, why):
    if not b:
        raise RuntimeError(why)


def receipt(p):
    p = Path(p)
    return {'path': str(p), 'bytes': p.stat().st_size, 'sha256': sha(p)}


def verify(r):
    p = Path(r['path'])
    need(p.is_file() and p.stat().st_size == r['bytes'] and sha(p) == r['sha256'],
         'Input changed: ' + str(p))


def resources():
    m = {x.split(':')[0]: int(x.split()[1]) * 1024
         for x in Path('/proc/meminfo').read_text().splitlines()}
    limit = int(Path('/sys/fs/cgroup/memory.max').read_text())
    used = int(Path('/sys/fs/cgroup/memory.current').read_text())
    return {'effective_free_bytes': min(m['MemAvailable'], limit - used),
            'cgroup_current_bytes': used, 'cgroup_limit_bytes': limit,
            'system_free_bytes': shutil.disk_usage('/root').free,
            'data_free_bytes': shutil.disk_usage('/root/autodl-tmp').free}


def pid_identity(pid):
    return {'pid': pid, 'start_ticks': int(Path(f'/proc/{pid}/stat').read_text().split()[21]),
            'boot_id': Path('/proc/sys/kernel/random/boot_id').read_text().strip()}


def written_bytes(out):
    return sum(p.stat().st_size for p in out.rglob('*') if p.is_file())


def read_gs():
    import pandas as pd
    s = pd.read_csv(GS / 'Keloid_EUR.top1000_Z.gs', sep='\t', dtype=str)
    members = [x.rsplit(':', 1) for x in s.GENESET.iloc[0].split(',')]
    return {g: float(w) for g, w in members}


def freeze():
    import numpy as np
    import pandas as pd
    need(not (HERE / 'FREEZE.json').exists(), 'Freeze exists; do not overwrite')
    need(HERE.resolve().is_relative_to(Path('/root/five_new_studies_20261020')),
         'Outputs must use system volume')
    HERE.mkdir(parents=True, exist_ok=True)
    rf = json.loads((RNA / 'COMPLETE.json').read_text())
    gf = json.loads((GS / 'COMPLETE.json').read_text())
    need(rf['status'] == 'COMPLETE_INPUT_FREEZE_NO_SCDRS_FIT', 'Incomplete RNA freeze')
    need(gf['status'] == 'COMPLETE_OFFICIAL_TOP1000_GS_NO_RNA_SELECTION_NO_SCORE',
         'Incomplete original gene set')
    shared = OLD / 'shared_inputs/programs_01_03'
    program = json.loads((shared / 'program_registry.json').read_text())
    reactome = shared / 'sources/NCBI2Reactome_All_Levels.txt'
    need(sha(reactome) == program['sources'][reactome.name]['sha256'],
         'Reactome differs from already frozen official source')
    pathways = {name: set() for name in MODULES}
    by_id = {v[0]: k for k, v in MODULES.items()}
    with reactome.open() as f:
        for line in f:
            a = line.rstrip('\n').split('\t')
            if len(a) >= 6 and a[1] in by_id and a[5] == 'Homo sapiens':
                pathways[by_id[a[1]]].add(a[0])
    need(all(pathways.values()), 'A fixed official pathway has no members')
    gs = read_gs()
    need(len(gs) == 1000 and all(np.isfinite(w) and w > 0 for w in gs.values()),
         'Original 1000 weighted genes changed')
    member = pd.read_csv(TABLE / 'S2a_fixed1000_gene_members_and_weights.tsv', sep='\t',
                         dtype={'canonical_entrez': str, 'CHR': str})
    need(len(member) == 1000 and set(member.canonical_entrez) == set(gs),
         'Genetic membership table differs from official gene set')
    need(all(float(r.official_GS_weight_string) == gs[r.canonical_entrez]
             for _, r in member.iterrows()), 'Official GS weights differ from membership table')
    regions = json.loads((TABLE / 'S9b_published_region_freeze_manifest.json').read_text())
    need(regions['n_merged_windows'] == 22, 'Original 22 regions changed')
    chromosome_of = {}
    for _, row in member.iterrows():
        chrom = str(row.CHR).replace('23', 'X')
        need(chrom in [str(i) for i in range(1, 23)] + ['X'], 'Unexpected original GS chromosome')
        chromosome_of[row.canonical_entrez] = chrom
    # The existing published-region manifest labels its intervals GRCh37, but the
    # actual source table/variant harmonization may use GRCh38. Do NOT overlap
    # coordinates until independently closed. The bounded approved fallback is
    # whole-chromosome concentration, which is invariant to hg19/hg38 position.
    ids = [str(i) for i in range(1, 23)] + ['X']
    pure_mass = {r: sum(w for g, w in gs.items() if chromosome_of[g] == r) for r in ids}
    largest = sorted(ids, key=lambda r: (-pure_mass[r], ids.index(r)))[0]
    need(pure_mass[largest] > 0, 'No original gene has an admitted chromosome')
    variants = {
        'remove_ECM': sorted(set(gs) & pathways['ECM']),
        'remove_largest_RNA_blind_chromosome': sorted(g for g in gs if chromosome_of[g] == largest),
    }
    inputs = [RNA / 'COMPLETE.json', GS / 'COMPLETE.json',
              GS / 'Keloid_EUR.top1000_Z.gs', TABLE / 'S2a_fixed1000_gene_members_and_weights.tsv',
              TABLE / 'S9b_published_region_freeze_manifest.json', shared / 'program_registry.json',
              reactome, PREP / 'environment_install_receipt_v1.json']
    for st in COHORTS:
        src = json.loads((PRIMARY / st / 'COMPLETE.json').read_text())
        need(src['status'] == 'COMPLETE_ONE_WHOLE_COHORT_OFFICIAL_API_NO_BATCH_SPLITTING',
             'Primary whole-cohort score not complete')
        inputs += [PRIMARY / st / 'COMPLETE.json',
                   PRIMARY / st / 'six_cell_scores.tsv.gz',
                   PRIMARY / st / 'actual_gene_mean_variance_bins.tsv.gz',
                   PRIMARY / st / 'cell_axis.tsv']
        for name, r in rf['cohorts'][st]['files'].items():
            p = RNA / st / name
            need(p.stat().st_size == r['bytes'] and sha(p) == r['sha256'],
                 'Original frozen RNA artifact changed')
            inputs.append(p)
        p = Path(rf['cohorts'][st]['input']['counts_path'])
        need(sha(p) == rf['cohorts'][st]['input']['counts_sha256'], 'Counts changed')
        inputs.append(p)
        for name, r in src['outputs'].items():
            if name in ['six_cell_scores.tsv.gz', 'actual_gene_mean_variance_bins.tsv.gz',
                        'cell_axis.tsv']:
                p = PRIMARY / st / name
                need(p.stat().st_size == r['bytes'] and sha(p) == r['sha256'],
                     'Original complete score artifact changed')
    env = json.loads((PREP / 'environment_install_receipt_v1.json').read_text())
    for r in env['verified_source_files']:
        if r['path'].startswith('scdrs/') and len(Path(r['path']).parts) == 2:
            p = SITE / r['path']
            need(sha(p) == r['sha256'], 'Installed official scDRS source changed')
            inputs.append(p)
    j = {
        'status': 'FROZEN_FINITE_SPECIFICITY_LOCUS_INPUTS_BEFORE_EXTENSION_EFFECTS',
        'UTC': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        'source_sha256': sha(__file__), 'cohorts': list(COHORTS),
        'fixed_modules': {k: {'id': MODULES[k][0], 'name': MODULES[k][1],
                              'canonical_entrez': sorted(v)} for k, v in pathways.items()},
        'original_gene_weights': gs, 'original_gene_chromosome': chromosome_of,
        'all_chromosomes': ids,
        'original_22_region_build_status': 'Coordinate build not independently closed; '
                                           'no coordinate overlap or region deletion performed. '
                                           'The original22 manifest is preserved as a source '
                                           'receipt, not used as a valid gene mapping.',
        'largest_chromosome_choice': {'id': largest, 'rule': 'Largest sum of original MAGMA GS Z '
                                     'weights among complete1-22/X chromosomes, no RNA/statistic '
                                     'used; ties in natural chromosome order. Native MAGMA CHR '
                                     'only, no uncertain coordinate overlap.',
                                     'pure_Z_mass_by_all23': pure_mass},
        'variants_removed_genes': variants,
        'sensitivity_gene_rule': 'Keep original remaining genes/weights; never backfill from RNA '
                                 'or rerank MAGMA; measured intersection is official API only',
        'full_cohort_background': True, 'n_ctrl_each': 1000, 'seed': 20261020,
        'sensitivity_FDR_family': 'Exactly30 FB-state x2variant tests across two cohorts; '
                                  'NA retained; original47 tests unchanged',
        'diagnostic_inference': 'No new P/FDR for module correlations or residual contrasts. '
                                'These use cells as descriptive measurements, not patients.',
        'region_inference': 'Chromosome concentration only, not single-locus exclusion. Pure Z '
                            'and official raw-kernel coefficient masses are separate; neither '
                            'decomposes normalized scDRS score/causal variance. Loss of '
                            'significance after deletion is not proof of exclusive mediation; '
                            'coverage, weight mass, MC magnitude and patient contrasts coexist.',
        'output_rule': 'system volume only; never persist Nx1006 new caches; max256MiB',
        'inputs': [receipt(p) for p in dict.fromkeys(inputs)],
    }
    atom(HERE / 'FREEZE.json', j)
    print(json.dumps({'freeze': str(HERE / 'FREEZE.json'),
                      'sha256': sha(HERE / 'FREEZE.json'), 'largest_chromosome': largest,
                      'module_sizes': {k: len(v) for k, v in pathways.items()},
                      'removed_counts': {k: len(v) for k, v in variants.items()}}))


def load_rna(st, frozen):
    import numpy as np
    import pandas as pd
    import h5py
    import scipy.sparse as sp
    import scanpy as sc
    import anndata as ad
    rf = json.loads((RNA / 'COMPLETE.json').read_text())['cohorts'][st]
    obs = pd.read_csv(RNA / st / 'all_final_cell_metadata.tsv.gz', sep='\t', dtype=str,
                      keep_default_na=False).set_index('cell_id')
    genes = pd.read_csv(RNA / st / 'admitted_gene_axis.tsv', sep='\t',
                        dtype={'canonical_entrez': str})
    cov = pd.read_csv(RNA / st / 'technical_covariates.tsv.gz', sep='\t', index_col=0,
                      float_precision='round_trip')
    scores = pd.read_csv(PRIMARY / st / 'six_cell_scores.tsv.gz', sep='\t', index_col=0,
                         float_precision='round_trip')
    need(obs.index.equals(cov.index) and obs.index.equals(scores.index),
         'Frozen metadata/cov/primary score order differs')
    with h5py.File(rf['input']['counts_path']) as h:
        z = h['layers/counts']
        x = sp.csr_matrix((z['data'][:], z['indices'][:], z['indptr'][:]),
                          shape=tuple(z.attrs['shape']))
    x.sum_duplicates()
    x.eliminate_zeros()
    need(x.shape == (rf['n_cells'], rf['input_genes']), 'Raw count shape changed')
    need(np.issubdtype(x.dtype, np.integer) and (x.data >= 0).all(), 'Counts no longer raw integer')
    umi = np.asarray(x.sum(1)).ravel()
    det = np.diff(x.indptr)
    need(np.array_equal(umi, obs.full_input_UMI.astype('int64')) and
         np.array_equal(det, obs.full_input_nGenes.astype('int64')), 'Input QC totals differ')
    need(np.allclose(cov.iloc[:, 0], np.log1p(umi), rtol=0, atol=1e-14) and
         np.allclose(cov.iloc[:, 1], np.log1p(det), rtol=0, atol=1e-14), 'Cov values differ')
    x = x[:, genes.source_index.to_numpy()].astype(np.float64).tocsr()
    a = ad.AnnData(X=x, obs=obs,
                   var=pd.DataFrame(index=genes.canonical_entrez))
    del x
    gc.collect()
    before = a.obs_names.copy()
    sc.pp.normalize_per_cell(a, counts_per_cell_after=1e4)
    sc.pp.log1p(a)
    need(a.obs_names.equals(before), 'Official normalization dropped a frozen cell')
    return a, cov, scores


def diagnostics(st, a, cov, original, frozen, out):
    import numpy as np
    import pandas as pd
    from scipy.stats import rankdata
    obs = a.obs.copy()
    fb = obs['FB_state'].isin(json.loads((RNA / 'COMPLETE.json').read_text())
                              ['cohorts'][st]['FB_state_family'])
    x = a.X[fb.to_numpy()]
    d = obs.loc[fb, ['independent_patient_id', 'capture_id', 'FB_state', 'region_original',
                    'region_P2_swap', 'include_P2_omit']].copy()
    d['original_norm_score'] = original.loc[d.index, 'norm_score'].astype(np.float64)
    for c in cov:
        d[c] = cov.loc[d.index, c]
    coverage = []
    gs = set(frozen['original_gene_weights'])
    module_cols = []
    for name, definition in frozen['fixed_modules'].items():
        all_members = set(definition['canonical_entrez'])
        for subset, request in [('full', all_members), ('nonGS', all_members - gs)]:
            members = sorted(request & set(a.var_names))
            col = name + '_' + subset + '_mean_logCP10k'
            module_cols.append(col)
            d[col] = np.asarray(x[:, a.var_names.get_indexer(members)].mean(1)).ravel() \
                if members else np.nan
            coverage.append({'cohort': st, 'module': name, 'fixed_pathway': definition['id'],
                             'subset': subset, 'requested_genes': len(request),
                             'measurable_genes': len(members),
                             'GS_overlap_in_full_module': len(all_members & gs),
                             'measured_members': ','.join(members),
                             'proxy_not_signaling_activity_or_proliferation_rate': True})
    pd.DataFrame(coverage).to_csv(out / f'{st}_fixed_module_coverage.tsv', sep='\t', index=False)
    d.to_csv(out / f'{st}_FB_cell_diagnostic_values.tsv.gz', sep='\t')

    def correlation(v, w):
        keep = np.isfinite(v) & np.isfinite(w)
        if keep.sum() < 3:
            return None, int(keep.sum())
        v, w = rankdata(np.asarray(v)[keep]), rankdata(np.asarray(w)[keep])
        if v.std() == 0 or w.std() == 0:
            return None, int(keep.sum())
        return float(np.corrcoef(v, w)[0, 1]), int(keep.sum())

    correlations = []
    for group_kind, keys in [('patient', ['independent_patient_id']),
                             ('patient_capture', ['independent_patient_id', 'capture_id']),
                             ('patient_state', ['independent_patient_id', 'FB_state']),
                             ('patient_capture_state', ['independent_patient_id', 'capture_id',
                                                        'FB_state'])]:
        for key, sub in d.groupby(keys, sort=True, dropna=False):
            key = key if isinstance(key, tuple) else (key,)
            meta = dict(zip(keys, key))
            for col in module_cols + list(cov):
                r, n = correlation(sub.original_norm_score.to_numpy(), sub[col].to_numpy())
                correlations.append({'cohort': st, 'view': group_kind, **meta,
                                     'variable': col, 'spearman_r': r, 'finite_cells': n,
                                     'n_cells': len(sub), 'descriptive_no_cell_P': True})
    pd.DataFrame(correlations).to_csv(out / f'{st}_within_patient_state_correlations.tsv',
                                    sep='\t', index=False, na_rep='NA')
    summaries = []
    states = json.loads((RNA / 'COMPLETE.json').read_text())['cohorts'][st]['FB_state_family']
    known = json.loads((RNA / 'COMPLETE.json').read_text())['cohorts'][st]['independent_patients']
    for patient in known:
        s = d.loc[d.independent_patient_id == patient].copy()
        n = len(s)
        base = np.column_stack([np.ones(n), s[list(cov)].to_numpy(dtype=float),
                                pd.get_dummies(s.capture_id, drop_first=True,
                                               dtype=float).to_numpy()])
        non = [k + '_nonGS_mean_logCP10k' for k in MODULES]
        designs = {'technical_capture': base,
                   'technical_capture_plus3_nonGS_modules':
                       np.column_stack([base, s[non].to_numpy(dtype=float)])}
        # Rank/residual contrasts are descriptive. No residual receives an old or new scDRS P.
        for design_name, design in designs.items():
            finite = np.isfinite(design).all(1) & np.isfinite(s.original_norm_score.to_numpy())
            y = s.original_norm_score.to_numpy(dtype=float)
            residual = np.full(n, np.nan)
            rank = None
            if finite.all() and n:
                fit, _, rank, _ = np.linalg.lstsq(design, y, rcond=None)
                residual = y - design @ fit
            for state in states:
                inside = s.FB_state.eq(state).to_numpy()
                other = ~inside
                ok = inside.any() and other.any() and np.isfinite(residual).all()
                summaries.append({'cohort': st, 'patient': patient, 'state': state,
                                  'design': design_name, 'n_state': int(inside.sum()),
                                  'n_other_FB': int(other.sum()), 'design_columns': design.shape[1],
                                  'design_rank': int(rank) if rank is not None else None,
                                  'original_mean_difference': float(y[inside].mean() - y[other].mean())
                                  if inside.any() and other.any() else None,
                                  'residual_mean_difference': float(residual[inside].mean() -
                                                                    residual[other].mean()) if ok else None,
                                  'status': 'DESCRIPTIVE' if ok else 'NA_MISSING_SUPPORT',
                                  'not_inherited_risk_not_conditional_scDRS_test': True})
    pd.DataFrame(summaries).to_csv(out / f'{st}_all_state_patient_residual_contrasts.tsv',
                                 sep='\t', index=False, na_rep='NA')
    stats = pd.read_csv(PRIMARY / st / 'actual_gene_mean_variance_bins.tsv.gz', sep='\t',
                        index_col=0, float_precision='round_trip')
    stats.index = stats.index.astype(str)
    measured = sorted(set(frozen['original_gene_weights']) & set(stats.index))
    raw = {g: frozen['original_gene_weights'][g] /
              float(np.sqrt(stats.loc[g, 'var_tech'] + 1e-2)) for g in measured}
    denominator = sum(raw.values())
    pure_all = sum(frozen['original_gene_weights'].values())
    pure_measured = sum(frozen['original_gene_weights'][g] for g in measured)
    regions = frozen['all_chromosomes']
    rows = []
    for region in regions:
        all_g = [g for g in frozen['original_gene_weights']
                 if frozen['original_gene_chromosome'][g] == region]
        m = [g for g in all_g if g in raw]
        z_all = sum(frozen['original_gene_weights'][g] for g in all_g)
        z_measured = sum(frozen['original_gene_weights'][g] for g in m)
        rows.append({'cohort': st, 'chromosome': region, 'original_GS_genes': len(all_g),
                     'measured_GS_genes': len(m), 'all_GS_Z_sum_in_region': z_all,
                     'all_GS_Z_sum_denominator': pure_all,
                     'all_GS_Z_fraction': z_all / pure_all,
                     'measured_Z_sum_in_region': z_measured,
                     'measured_Z_sum_denominator': pure_measured,
                     'measured_Z_fraction': z_measured / pure_measured,
                     'official_raw_kernel_coefficient_fraction':
                         sum(raw[g] for g in m) / denominator,
                     'raw_coefficient_denominator': denominator,
                     'not_normalized_score_decomposition': True})
    pd.DataFrame(rows).to_csv(out / f'{st}_all23_chromosome_weight_quality.tsv',
                             sep='\t', index=False)


def score_sensitivity(st, a, cov, original, frozen, out):
    import numpy as np
    import pandas as pd
    import scdrs
    scdrs.preprocess(a, cov=cov, adj_prop=None, n_mean_bin=20, n_var_bin=20, n_chunk=128)
    primary_genes = sorted(set(frozen['original_gene_weights']) & set(a.var_names))
    primary_weights = [frozen['original_gene_weights'][g] for g in primary_genes]
    # Kernel-only identity check: no baseline1000-control recomputation.
    a.X = a.X.tocsc()
    raw, _ = scdrs.method._compute_raw_score(a, primary_genes, primary_weights, 'vs')
    need(np.array_equal(np.asarray(raw, dtype=np.float32),
                        original.raw_score.to_numpy(dtype=np.float32)),
         'Recreated original preprocessing/raw kernel is not exact at official float32 output')
    rows, patients, coverage = [], [], []
    rf = json.loads((RNA / 'COMPLETE.json').read_text())['cohorts'][st]
    for variant, removed in frozen['variants_removed_genes'].items():
        dest = out / f'{st}__{variant}'
        dest.mkdir()
        row_start, patient_start = len(rows), len(patients)
        keep = {g: w for g, w in frozen['original_gene_weights'].items() if g not in set(removed)}
        measurable = sorted(set(keep) & set(a.var_names))
        need(measurable, 'All remaining genes unmeasurable; retain NA via explicit failure review')
        started = time.monotonic()
        original_selector = scdrs.method._select_ctrl_geneset
        calls = []
        def capture(*args, **kwargs):
            lists, weights = original_selector(*args, **kwargs)
            with gzip.open(dest / 'actual1000_control_members_and_input_weights.tsv.gz', 'wt') as f:
                f.write('control_id\tposition\tcanonical_entrez\tinput_weight\n')
                for k in sorted(lists):
                    for position, (g, w) in enumerate(zip(lists[k], weights[k])):
                        f.write(f'{k}\t{position}\t{g}\t{w:.17g}\n')
            calls.append({'n_ctrl': len(lists), 'measured': len(args[1]), 'seed': args[6]})
            return lists, weights
        scdrs.method._select_ctrl_geneset = capture
        try:
            y = scdrs.score_cell(a, list(keep), gene_weight=list(keep.values()),
                                ctrl_match_key='mean_var', n_ctrl=1000, weight_opt='vs',
                                return_ctrl_raw_score=False, return_ctrl_norm_score=True,
                                random_seed=20261020)
        finally:
            scdrs.method._select_ctrl_geneset = original_selector
        need(len(calls) == 1 and calls[0]['n_ctrl'] == 1000 and
             calls[0]['seed'] == 20261020, 'Official selection coverage changed')
        need(y.index.equals(a.obs_names) and y.shape == (a.n_obs, 1006),
             'Whole-cohort sensitivity cell/score coverage changed')
        need(np.isfinite(y.to_numpy()).all(), 'Nonfinite scores; preserve failure, do not delete cells')
        y.iloc[:, :6].to_csv(dest / 'six_cell_scores.tsv.gz', sep='\t')
        controls = [f'ctrl_norm_score_{k}' for k in range(1000)]
        for state in rf['FB_state_family']:
            mask = a.obs.FB_state.eq(state).to_numpy()
            r = {'cohort': st, 'variant': variant, 'state': state, 'n_cells': int(mask.sum()),
                 'requested_genes_after_removal': len(keep), 'measurable_genes': len(measurable),
                 'n_ctrl': 1000}
            if mask.any():
                q = float(np.quantile(y.loc[mask, 'norm_score'], .95))
                z = np.quantile(y.loc[mask, controls], .95, axis=0)
                r.update(status='COMPLETE', score_q95=q,
                         assoc_mcp=float(((z >= q).sum() + 1) / 1001),
                         assoc_mcz=float((q - z.mean()) / z.std()) if z.std() else None)
            else:
                r.update(status='NA_EMPTY_FROZEN_STATE', score_q95=None,
                         assoc_mcp=None, assoc_mcz=None)
            rows.append(r)
        for patient in rf['independent_patients']:
            pmask = a.obs.independent_patient_id.eq(patient).to_numpy()
            fbmask = pmask & a.obs.FB_state.isin(rf['FB_state_family']).to_numpy()
            for state in rf['FB_state_family']:
                smask = fbmask & a.obs.FB_state.eq(state).to_numpy()
                omask = fbmask & ~a.obs.FB_state.eq(state).to_numpy()
                patients.append({'cohort': st, 'variant': variant, 'patient': patient,
                                 'state': state, 'n_state': int(smask.sum()),
                                 'n_other_FB': int(omask.sum()),
                                 'mean_difference_vs_other_FB':
                                     float(y.loc[smask, 'norm_score'].mean() -
                                           y.loc[omask, 'norm_score'].mean())
                                     if smask.any() and omask.any() else None,
                                 'descriptive_not_genotype_replication': True})
        c = {'cohort': st, 'variant': variant, 'removed_original_genes': removed,
             'requested_genes_after_removal': len(keep), 'measurable_genes': measurable,
             'original_pure_Z_sum': sum(frozen['original_gene_weights'].values()),
             'retained_pure_Z_sum': sum(keep.values()),
             'retained_pure_Z_fraction': sum(keep.values()) /
                                         sum(frozen['original_gene_weights'].values()),
             'original_measurable_genes': len(primary_genes),
             'retained_measurable_gene_fraction': len(measurable) / len(primary_genes),
             'all_final_cells': a.n_obs, 'n_ctrl': 1000, 'seed': 20261020,
             'elapsed_seconds': time.monotonic() - started, 'selection_calls': calls,
             'no_Nx1006_cache_saved': True, 'nonfinite': 0,
             'original_raw_kernel_exact_float32_check': True,
             'whole_cohort_background_preserved': True,
             'peak_process_RSS_bytes': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024}
        atom(dest / 'SCORE_DIAGNOSTICS.json', c)
        pd.DataFrame(rows[row_start:]).to_csv(dest / 'all_frozen_state_MC_tests_unadjusted.tsv',
                                             sep='\t', index=False, na_rep='NA')
        pd.DataFrame(patients[patient_start:]).to_csv(dest / 'all_state_patient_contrasts.tsv',
                                                    sep='\t', index=False, na_rep='NA')
        atom(dest / 'COMPLETE.json', {
            'status': 'COMPLETE_ONE_FIXED_SENSITIVITY_WHOLE_COHORT_1000CONTROLS',
            'cohort': st, 'variant': variant, 'source_sha256': sha(__file__),
            'freeze_sha256': sha(HERE / 'FREEZE.json'),
            'outputs': {p.name: receipt(p) for p in dest.iterdir() if p.is_file()},
            'score_diagnostics': c, 'FDR_pending_all30': True,
            'non_significance_not_proof_of_ECM_or_chromosome_mediation': True})
        coverage.append(c)
        del y
        gc.collect()
    return rows, patients, coverage


def worker(frozen, out, stage):
    try:
        import numpy as np
        import pandas as pd
        from statsmodels.stats.multitest import multipletests
        all_rows, all_patients, all_coverage = [], [], []
        for st in COHORTS:
            atom(out / 'CURRENT_STAGE.json', {'cohort': st, 'stage': 'load_original_sparse_RNA'})
            a, cov, original = load_rna(st, frozen)
            if stage == 'diagnostics':
                diagnostics(st, a, cov, original, frozen, out)
                atom(out / f'{st}_DIAGNOSTICS_COMPLETE.json',
                     {'status': 'COMPLETE_DESCRIPTIVE_ONLY_NO_NEW_P', 'cohort': st,
                      'source_sha256': sha(__file__), 'freeze_sha256': sha(HERE / 'FREEZE.json'),
                      'outputs': {p.name: receipt(p) for p in out.glob(st + '_*')
                                  if p.is_file()}})
            else:
                atom(out / 'CURRENT_STAGE.json', {'cohort': st, 'stage': 'two_fixed_official_sensitivities'})
                rows, patients, coverage = score_sensitivity(st, a, cov, original, frozen, out)
                all_rows += rows
                all_patients += patients
                all_coverage += coverage
            del a, cov, original
            gc.collect()
        if stage == 'sensitivity':
            frame = pd.DataFrame(all_rows)
            need(len(frame) == 30, 'Must retain all30 sensitivity rows')
            finite = frame.assoc_mcp.notna()
            # Fixed planned family size30 includes NA slots. P=1 placeholders only
            # reserve multiplicity; displayed P/Q for unavailable states stay NA.
            p = frame.assoc_mcp.fillna(1).to_numpy(dtype=float)
            q = multipletests(p, method='fdr_bh')[1]
            frame['sensitivity_family_BH30_q'] = np.where(finite, q, np.nan)
            frame.to_csv(out / 'all30_state_sensitivity_tests_BH30.tsv', sep='\t',
                         index=False, na_rep='NA')
            pd.DataFrame(all_patients).to_csv(out / 'all120_state_patient_sensitivity_contrasts.tsv',
                                             sep='\t', index=False, na_rep='NA')
            atom(out / 'ALL4_OFFICIAL_SCORE_DIAGNOSTICS.json', all_coverage)
        need(written_bytes(HERE) <= OUTPUT_CAP, 'Combined A+B256MiB cap exceeded')
        atom(out / 'INTERNAL_COMPLETE.json', {
            'status': 'COMPLETE_FINITE_' + stage.upper(), 'stage': stage,
            'original47_tests_and1000_GS_unchanged': True,
            'n_sensitivity_rows': len(all_rows), 'n_sensitivity_score_runs': len(all_coverage),
            'n_patient_sensitivity_rows': len(all_patients), 'n_ctrl_each': 1000,
            'not_causal_not_clinical_validation': True,
            'peak_RSS_bytes': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
            'outputs': {str(p.relative_to(out)): receipt(p)
                        for p in out.rglob('*') if p.is_file() and
                        p.name not in ['INTERNAL_COMPLETE.json', 'RUNNING.json']}})
    except BaseException as e:
        atom(out / 'WORKER_FAILURE_PRESERVED.json', {'error': repr(e),
                                                    'traceback': traceback.format_exc()})
        raise


def run(args):
    approved = Path(args.approval)
    need(sha(approved) == args.approval_sha256, 'Root approval changed')
    approval = json.loads(approved.read_text())
    stage = args.stage
    need(stage in ['diagnostics', 'sensitivity'], 'Explicit A/B stage required')
    need(approval['action'] == 'RUN_N02_FINITE_SPECIFICITY_LOCUS_' + stage.upper() + '_ONCE',
         'Wrong action')
    need(approval['source_sha256'] == sha(__file__), 'Source not root-approved')
    need(approval['freeze_sha256'] == sha(HERE / 'FREEZE.json'), 'Freeze not root-approved')
    need(sys.executable == str(PYTHON), 'Use existing exact scDRS environment')
    frozen = json.loads((HERE / 'FREEZE.json').read_text())
    need(frozen['source_sha256'] == sha(__file__), 'Freeze source differs')
    for r in frozen['inputs']:
        verify(r)
    nthreads = approval['threads']
    need(nthreads == 4 and approval['output_reserve_bytes'] == OUTPUT_CAP,
         'Fixed4CPU/256MiB output design changed')
    for name in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
        os.environ[name] = str(nthreads)
    allowed = sorted(os.sched_getaffinity(0))
    quota, period = Path('/sys/fs/cgroup/cpu.max').read_text().split()
    need(len(allowed) >= nthreads and (quota == 'max' or int(quota) / int(period) >= nthreads),
         'CPU quota insufficient')
    os.sched_setaffinity(0, allowed[-nthreads:])
    out = HERE / ('results_A_v1' if stage == 'diagnostics' else 'results_B_v1')
    need(not out.exists(), 'No overwrite/retry of partial or complete results')
    with (HERE / (stage + '_ONCE_CLAIM.json')).open('x') as f:
        json.dump({**pid_identity(os.getpid()), 'approval_sha256': args.approval_sha256}, f)
    # Isolated owner lock; existing workers are not stopped or modified.
    with (HERE / 'owner.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        r = resources()
        need(r['data_free_bytes'] >= 5 * 2**30, 'Data disk5GiB floor')
        need(r['system_free_bytes'] >= 5 * 2**30 + OUTPUT_CAP, 'System output reserve unavailable')
        need(r['effective_free_bytes'] >= approval['RAM_cap_bytes'] +
             approval['memory_floor_bytes'] + approval['peer_reserved_bytes'],
             'Actual cgroup/host memory admission failed; do not auto-retry')
        need(approval['RAM_cap_bytes'] <= 24 * 2**30 and
             approval['memory_floor_bytes'] >= 8 * 2**30 and
             approval['max_wall_seconds'] <= 5400, 'Approved finite resource design exceeded')
        out.mkdir()
        atom(out / 'RUNNING.json', {**pid_identity(os.getpid()), 'resources': r,
                                  'approval': approval, 'UTC': time.strftime('%Y-%m-%dT%H:%M:%SZ',
                                                                           time.gmtime())})
        child = mp.get_context('spawn').Process(target=worker, args=(frozen, out, stage))
        started = time.monotonic()
        child.start()
        peak = 0
        error = None
        while child.is_alive():
            child.join(5)
            try:
                r = resources()
                need(min(r['system_free_bytes'], r['data_free_bytes']) >= 5 * 2**30,
                     'A disk violated5GiB floor')
                need(r['effective_free_bytes'] >= approval['memory_floor_bytes'],
                     'Effective memory floor violated')
                need(written_bytes(HERE) <= OUTPUT_CAP, 'Combined A+B output cap exceeded')
                need(time.monotonic() - started <= approval['max_wall_seconds'],
                     'Finite wall budget exceeded')
                status = Path(f'/proc/{child.pid}/status').read_text()
                rss = next((int(x.split()[1]) * 1024 for x in status.splitlines()
                            if x.startswith('VmRSS:')), 0)
                peak = max(peak, rss)
                need(rss <= approval['RAM_cap_bytes'], 'Owned RSS cap exceeded')
            except FileNotFoundError:
                pass
            except BaseException as e:
                error = repr(e)
                child.terminate()
                child.join(30)
                break
        if error or child.exitcode != 0:
            atom(out / 'FAILED_PRESERVED_NO_AUTO_RETRY.json',
                 {'error': error, 'child_returncode': child.exitcode,
                  'peak_RSS_bytes': peak, 'resources': resources()})
            raise RuntimeError('Finite extension failed; all partial evidence retained')
        complete = json.loads((out / 'INTERNAL_COMPLETE.json').read_text())
        for r in complete['outputs'].values():
            verify(r)
        need(sha(__file__) == approval['source_sha256'], 'Source changed while running')
        atom(out / 'COMPLETE.json', {
            'status': 'COMPLETE_N02_FINITE_SPECIFICITY_LOCUS_ROOT_REVIEW_PENDING',
            'stage': stage,
            'source_sha256': sha(__file__), 'freeze_sha256': sha(HERE / 'FREEZE.json'),
            'approval_sha256': args.approval_sha256,
            'internal_complete_sha256': sha(out / 'INTERNAL_COMPLETE.json'),
            'returncode': child.exitcode, 'elapsed_seconds': time.monotonic() - started,
            'peak_RSS_bytes': max(peak, complete['peak_RSS_bytes']),
            'outputs_bytes': written_bytes(out), 'combined_directory_bytes': written_bytes(HERE),
            'final_resources': resources()})
        print(json.dumps({'COMPLETE': str(out / 'COMPLETE.json'),
                          'sha256': sha(out / 'COMPLETE.json')}))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('mode', choices=['freeze', 'run'])
    p.add_argument('--stage', choices=['diagnostics', 'sensitivity'])
    p.add_argument('--approval')
    p.add_argument('--approval-sha256')
    args = p.parse_args()
    need(Path('/proc').exists() and ROOT.exists(), 'Scientific steps require existing AutoDL')
    if args.mode == 'freeze':
        freeze()
    else:
        need(args.approval and args.approval_sha256, 'Root scientific approval required')
        run(args)


if __name__ == '__main__':
    main()
