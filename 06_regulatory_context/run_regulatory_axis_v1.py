"""Bounded patient-pseudobulk CollecTRI regulatory context, not causal TF tests.

freeze reads static memberships and metadata only. run needs root source/freeze
authorization. No old input/result/manuscript is modified. Scientific run: AutoDL.
"""
from pathlib import Path
import argparse, collections, hashlib, importlib.metadata, inspect, json, os
import resource, shutil, subprocess, sys, time, traceback
for _k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[_k] = '2'
P = Path('/root/five_new_studies_20261020/projects/N02_keloid_genetic_cell_state')
O = Path('/root/autodl-tmp/Five_Studies_Execution_20261020')
HERE = P / 'mechanisms_20261021/regulatory_axis_v1'
RNA = P / 'reviews/scDRS_production_preparation_v1/RNA_freeze_v1'
GS = P / 'reviews/scDRS_production_preparation_v1/official_MAGMA_GS_v1'
EXT = P / 'extensions_20261021/external_RNA_v1'
NET = O / 'shared/programs/curated_regulatory_LR_v1/CollecTRI_human_interactions.tsv.gz'
HGNC = O / 'shared_inputs/programs_01_03/sources/hgnc_complete_set.txt'
AXES = {'NFkB': ['RELA', 'NFKB1', 'NFKB2', 'REL', 'RELB'],
        'AP1_JNK_downstream': ['JUN', 'JUNB', 'JUND', 'FOS', 'FOSB', 'FOSL1', 'FOSL2'],
        'TGFbeta_SMAD': ['SMAD2', 'SMAD3', 'SMAD4']}
CANDIDATES = {'TAB2': 'ENSG00000055208', 'NEDD4': 'ENSG00000069869',
              'ASIP': 'ENSG00000101440', 'LSP1': 'ENSG00000130592',
              'TNNT3': 'ENSG00000130595', 'SETBP1': 'ENSG00000152217',
              'CARS1-AS1': 'ENSG00000247473', 'SETBP1-DT': 'ENSG00000267414',
              'PRR33': 'ENSG00000283787'}
CAP = 64 * 2**20


def need(ok, why):
    if not ok:
        raise RuntimeError(why)


def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda: f.read(2**20), b''):
            h.update(b)
    return h.hexdigest()


def receipt(p):
    p = Path(p)
    return {'path': str(p), 'bytes': p.stat().st_size, 'sha256': sha(p)}


def atom(p, j):
    p = Path(p); q = p.with_suffix(p.suffix + '.tmp')
    q.write_text(json.dumps(j, ensure_ascii=False, indent=2, allow_nan=False) + '\n')
    q.replace(p)


def clean(x):
    import numpy as np
    if isinstance(x, dict):
        return {str(k): clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, np.generic):
        x = x.item()
    if isinstance(x, float) and not np.isfinite(x):
        return None
    return x


def resources():
    m = {x.split(':')[0]: int(x.split()[1]) * 1024
         for x in Path('/proc/meminfo').read_text().splitlines()}
    lim = Path('/sys/fs/cgroup/memory.max').read_text().strip()
    free = m['MemAvailable']
    if lim != 'max':
        free = min(free, int(lim) - int(Path('/sys/fs/cgroup/memory.current').read_text()))
    return {'effective_free_bytes': free, 'system_free_bytes': shutil.disk_usage('/root').free,
            'data_free_bytes': shutil.disk_usage('/root/autodl-tmp').free}


def identity(pid):
    s = Path(f'/proc/{pid}/stat').read_text()
    return {'pid': pid, 'start_ticks': int(s[s.rfind(')')+2:].split()[19]),
            'boot_id': Path('/proc/sys/kernel/random/boot_id').read_text().strip()}


def harray(n):
    import h5py, numpy as np
    if isinstance(n, h5py.Dataset):
        a = n[:]
        return np.asarray([v.decode() if isinstance(v, bytes) else str(v) for v in a]) if a.dtype.kind in 'OSU' else a
    a = harray(n['categories']); c = n['codes'][:]
    return np.asarray([a[i] if i >= 0 else '' for i in c])


def network_and_genes():
    import pandas as pd
    h = pd.read_csv(HGNC, sep='\t', dtype=str, keep_default_na=False)
    h = h[h.status.eq('Approved')].copy()
    unique = h.ensembl_gene_id.ne('') & ~h.ensembl_gene_id.duplicated(False) & ~h.symbol.duplicated(False)
    ens_to_symbol = dict(zip(h.loc[unique, 'ensembl_gene_id'], h.loc[unique, 'symbol']))
    approved = set(h.symbol); aliases = collections.defaultdict(set)
    for r in h.itertuples():
        for col in ('alias_symbol', 'prev_symbol'):
            for a in str(getattr(r, col)).split('|'):
                if a:
                    aliases[a].add(r.symbol)
    def mapper(x):
        if x in approved:
            return x
        return next(iter(aliases[x])) if len(aliases[x]) == 1 else x
    n = pd.read_csv(NET, sep='\t')
    n['source_original'] = n.source; n['target_original'] = n.target
    n.source = n.source.map(mapper); n.target = n.target.map(mapper)
    # The established 04 recipe: signed duplicate edge weights average; zero cancels.
    net = n.groupby(['source', 'target'], as_index=False).weight.mean()
    net = net[net.weight.ne(0)].copy()
    return h, ens_to_symbol, net, n


def freeze():
    import decoupler as dc, pandas as pd
    need(not (HERE / 'FREEZE.json').exists(), 'Existing freeze; do not overwrite')
    HERE.mkdir(parents=True, exist_ok=True)
    rf = json.loads((RNA / 'COMPLETE.json').read_text())
    ef = json.loads((EXT / 'FREEZE.json').read_text())
    ec = json.loads((EXT / 'COMPLETE.json').read_text())
    need('COMPLETE' in ec['status'], 'External pipeline incomplete')
    h, e2s, net, original = network_and_genes()
    need(set(CANDIDATES.values()).issubset(e2s), 'Candidate gene identifiers not HGNC-resolved')
    need(all(e2s[v] == k for k, v in CANDIDATES.items()), 'Candidate symbol/Ensembl disagreement')
    g = pd.read_csv(GS / 'Keloid_EUR.top1000_Z.gs', sep='\t', dtype=str)
    gsids = {x.rsplit(':', 1)[0] for x in g.GENESET.iloc[0].split(',')}
    need(len(gsids) == 1000, 'Original 1000 genetic members differ')
    gsens = set(h.loc[h.entrez_id.isin(gsids), 'ensembl_gene_id']) - {''}
    mapping = set().union(*[set(v) for v in ef['mapping_features'].values()])
    excluded_ens = gsens | set(CANDIDATES.values()) | mapping
    excluded_symbols = sorted(e2s[x] for x in excluded_ens if x in e2s)
    needed = [NET, HGNC, RNA / 'COMPLETE.json', GS / 'COMPLETE.json', GS / 'Keloid_EUR.top1000_Z.gs',
              EXT / 'FREEZE.json', EXT / 'COMPLETE.json', EXT / 'frozen41449_cell_manifest.tsv.gz',
              EXT / 'mapping_v1/all12249_two_reference_state_transport.tsv.gz',
              EXT / 'mapping_v1/MAPPING_TRANSPORT_GATES.json',
              P / 'extensions_20261021/gene_state_bridge_v1/results/all_primary_PP4_ge08_unique_genes.tsv',
              P / 'extensions_20261021/specificity_locus_v1/FREEZE.json']
    cohorts = {}
    for st in ('GSE220300', 'GSE335482'):
        p = Path(rf['cohorts'][st]['input']['FB_state_path'])
        need(sha(p) == rf['cohorts'][st]['input']['FB_state_sha256'], 'Discovery FB raw source changed')
        m = RNA / st / 'all_final_cell_metadata.tsv.gz'
        needed += [p, m]
        cohorts[st] = {'raw': str(p), 'node': 'layers/counts', 'metadata': str(m),
                       'patients': rf['cohorts'][st]['independent_patients']}
    for c in ef['captures']:
        needed.append(Path(c['raw_path']))
    for method in ('ulm', 'mlm'):
        needed.append(Path(inspect.getfile(getattr(dc.mt, method).func)))
    import decoupler._Method as official_method
    needed += [Path(inspect.getfile(official_method._run)), Path(inspect.getfile(type(dc.mt.ulm)))]
    sf = json.loads((P / 'extensions_20261021/specificity_locus_v1/FREEZE.json').read_text())
    programs = sf['fixed_modules']
    need(programs is not None, 'Existing fixed module freeze schema missing')
    f = {'status': 'STATIC_FROZEN_NO_NEW_RNA_EFFECTS', 'source_sha256': sha(__file__),
         'decoupler_version': importlib.metadata.version('decoupler'),
         'cohorts': cohorts, 'external_captures': ef['captures'], 'candidates': CANDIDATES,
         'axes': AXES, 'network_sources': sorted(net.source.unique()),
         'ensembl_to_symbol': e2s, 'original_GS_Entrez': sorted(gsids),
         'GS_Ensembl': sorted(gsens), 'external_mapping_Ensembl': sorted(mapping),
         'excluded_symbols_no_direct_overlap': excluded_symbols, 'fixed_modules': programs,
         'views': ['all_targets', 'no_GS9_mapping_targets'], 'tmin': 10,
         'contrast': 'log2(patient FB01 full-library CP10K +0.1) minus log2(patient otherFB full-library CP10K +0.1); cohort-wide any-FB-count>0 gene background',
         'external_definition': 'Frozen two-reference FB01_consensus vs all other original own-QC FB including Unknown. Normal scars 0/1/1 FB01 preserved as no/tiny-support diagnostic, not replicated patient mechanism.',
         'inference_boundary': 'Gene-network ULM/MLM q is not patient association P. Cohort/patient directions and small-n correlations descriptive. All original candidates and all network TFs preserved; no chosen favorable axis.',
         'correlation_family': 'Each9candidate x15priorTF x3cohort x2views x2methods, Spearman of patient FB01-other contrasts; only documented keloid patients; no inferential P at n4/4/3.',
         'resource': {'own_RSS_GiB': 4, 'threads': 2, 'output_cap_bytes': CAP},
         'inputs': [receipt(p) for p in sorted(set(needed))]}
    original.to_csv(HERE / 'canonical_network_edge_ledger.tsv.gz', sep='\t', index=False)
    atom(HERE / 'FREEZE.json', f)
    print(json.dumps({'freeze_sha256': sha(HERE / 'FREEZE.json'), 'source_sha256': sha(__file__),
                      'TFs': len(f['network_sources']), 'excluded_symbols': len(excluded_symbols),
                      'axis_TFs': sum(map(len, AXES.values()))}))


def aggregate(path, node, metadata, groups, e2s, sums):
    import h5py, numpy as np
    from scipy.sparse import csr_matrix
    with h5py.File(path) as h:
        obs = harray(h['obs'][h['obs'].attrs['_index']]); var = harray(h['var'][h['var'].attrs['_index']])
        need(len(set(var)) == len(var), 'Duplicate raw Ensembl axis')
        meta = metadata.reindex(obs)
        z = h[node]; need(z.attrs['encoding-type'] == 'csr_matrix', 'Raw CSR required')
        ix = [j for j, e in enumerate(var) if e in e2s]
        syms = [e2s[var[j]] for j in ix]
        need(len(syms) == len(set(syms)), 'Ambiguous canonical count axis')
        ptr = z['indptr']
        for start in range(0, len(obs), 512):
            end = min(start + 512, len(obs)); p = ptr[start:end+1]; lo, hi = int(p[0]), int(p[-1])
            d = z['data'][lo:hi]; ids = z['indices'][lo:hi]
            need(np.isfinite(d).all() and (d >= 0).all() and np.equal(d, np.rint(d)).all(), 'Original integer raw counts required')
            x = csr_matrix((d, ids, p-lo), shape=(end-start, len(var)))
            lib = np.asarray(x.sum(1)).ravel().astype(np.int64)
            for key, mask in groups(meta.iloc[start:end]).items():
                rows = np.flatnonzero(mask)
                if not len(rows):
                    continue
                v = np.asarray(x[rows][:, ix].sum(0)).ravel().astype(np.int64)
                cur = sums.setdefault(key, {'n': 0, 'lib': 0, 'raw': collections.Counter()})
                cur['n'] += len(rows); cur['lib'] += int(lib[rows].sum())
                cur['raw'].update(dict(zip(syms, v)))
            # Verify observed original full libraries against previously frozen units.
            column = 'full_input_UMI' if 'full_input_UMI' in meta else 'total_counts'
            have = meta.iloc[start:end][column].notna().to_numpy()
            need(np.array_equal(lib[have], meta.iloc[start:end].loc[have, column].to_numpy(dtype=np.int64)), 'Full-library count source differs')
    return sums


def cohort_groups(meta):
    import numpy as np
    result = {}
    for p in meta['patient'].dropna().unique():
        for group in ('FB01', 'otherFB'):
            result[(str(p), group)] = (meta.patient.eq(p) & meta.group.eq(group)).to_numpy()
    return result


def worker(f):
    import decoupler as dc, numpy as np, pandas as pd
    from scipy.stats import spearmanr, false_discovery_control
    os.sched_setaffinity(0, {6, 7})
    out = HERE / 'results'; need(not out.exists(), 'Existing output; no overwrite/retry')
    out.mkdir()
    for r in f['inputs']:
        need(Path(r['path']).stat().st_size == r['bytes'] and sha(r['path']) == r['sha256'], 'Frozen input/source changed: '+r['path'])
    _, e2s, net, _ = network_and_genes()
    axis_of = {tf: axis for axis, tfs in AXES.items() for tf in tfs}
    groups_all = []; activity = []; coverage = []; gene_rows = []; errors = []; correlations = []
    inputs = []
    for st, definition in f['cohorts'].items():
        m = pd.read_csv(definition['metadata'], sep='\t', dtype=str, keep_default_na=False).set_index('cell_id')
        m = m[m.lineage.eq('fibroblast')].copy()
        m['patient'] = m.independent_patient_id.where(m.independent_patient_id.isin(definition['patients']))
        m['group'] = np.where(m.FB_state.eq('FB01'), 'FB01', 'otherFB'); m['condition'] = 'keloid'
        inputs.append((st, [(definition['raw'], definition['node'])], m))
    m = pd.read_csv(EXT / 'frozen41449_cell_manifest.tsv.gz', sep='\t', dtype=str, keep_default_na=False).set_index('cell_id')
    m = m[m.is_own_QC_FB.eq('True')].copy()
    transfer = pd.read_csv(EXT / 'mapping_v1/all12249_two_reference_state_transport.tsv.gz', sep='\t', dtype=str, keep_default_na=False).set_index('cell_id')
    need(set(m.index) == set(transfer.index), 'External FB and frozen mapping identity mismatch')
    m['patient'] = m.patient_id; m['group'] = np.where(transfer.loc[m.index].transport_label.eq('FB01_consensus'), 'FB01', 'otherFB')
    inputs.append(('GSE163973', [(c['raw_path'], 'X') for c in f['external_captures']], m))
    for st, sources, m in inputs:
        sums = {}
        for path, node in sources:
            aggregate(path, node, m, cohort_groups, e2s, sums)
        patients = sorted(m.patient.dropna().unique()); genes = sorted(set().union(*(set(x['raw']) for x in sums.values())))
        raw = pd.DataFrame({k: [v['raw'][g] for g in genes] for k, v in sums.items()}, index=genes).T
        measured = raw.sum(0).gt(0)
        genes = list(raw.columns[measured]); cp = raw[genes].astype(float)
        for k, v in sums.items():
            cp.loc[k] *= 1e4 / v['lib']
        delta = []; sample_metadata = []
        for p in patients:
            condition = str(m.loc[m.patient.eq(p), 'condition'].iloc[0])
            a, b = (p, 'FB01'), (p, 'otherFB')
            na = sums.get(a, {}).get('n', 0); nb = sums.get(b, {}).get('n', 0)
            groups_all.append({'study': st, 'patient': p, 'condition': condition, 'FB01_n': na, 'otherFB_n': nb,
                               'FB01_library': sums.get(a, {}).get('lib', 0), 'otherFB_library': sums.get(b, {}).get('lib', 0),
                               'support': 'ADEQUATE_DESCRIPTIVE' if min(na, nb) >= 20 else 'TINY_OR_ABSENT_FB01_DIAGNOSTIC'})
            if not na or not nb:
                continue
            d = np.log2(cp.loc[a].to_numpy()+.1) - np.log2(cp.loc[b].to_numpy()+.1)
            delta.append(d); sample_metadata.append((p, condition, na, nb))
            for symbol in CANDIDATES:
                j = genes.index(symbol) if symbol in genes else -1
                gene_rows.append({'study': st, 'patient': p, 'condition': condition, 'gene': symbol,
                                  'assay_status': 'ASSAYED_DETECTED' if j >= 0 else 'UNDETECTED_OR_UNMAPPED_NA',
                                  'log2_CP10K_plus0p1_delta': d[j] if j >= 0 else np.nan})
        X = pd.DataFrame(delta, index=[v[0] for v in sample_metadata], columns=genes)
        X.to_csv(out / f'{st}_patient_FB01_otherFB_log2CP10K_delta.tsv.gz', sep='\t')
        for view in f['views']:
            excluded = set(f['excluded_symbols_no_direct_overlap']) if view != 'all_targets' else set()
            features = [g for g in genes if g not in excluded]
            small_net = net[net.target.isin(features)].copy()
            counts = small_net.groupby('source').target.nunique().to_dict()
            for tf in f['network_sources']:
                targets = set(net.loc[net.source.eq(tf), 'target'])
                coverage.append({'study': st, 'view': view, 'TF': tf, 'prior_axis': axis_of.get(tf, ''),
                                 'full_regulon_targets': len(targets), 'assayed_targets': counts.get(tf, 0),
                                 'excluded_direct_overlap_targets': len(targets & excluded),
                                 'status': 'ADMITTED_TMIN10' if counts.get(tf, 0) >= f['tmin'] else 'INSUFFICIENT_TARGETS_NA'})
            for method in ('ulm', 'mlm'):
                begin = time.monotonic()
                try:
                    scores, official_values = getattr(dc.mt, method)(X[features], net=small_net, tmin=f['tmin'], empty=False, verbose=True)
                    need(scores.index.equals(X.index), 'Official score patient order differs')
                    # decoupler2.2.0 _run: ULM already BH, MLM explicitly skipped.
                    # Keep the returned quantity and add BH only for raw MLM P.
                    qvalues = official_values.copy()
                    if method == 'mlm':
                        for p in qvalues.index:
                            values = official_values.loc[p].to_numpy(float); valid = np.isfinite(values)
                            corrected = np.full(values.shape, np.nan)
                            corrected[valid] = false_discovery_control(values[valid], method='bh')
                            qvalues.loc[p] = corrected
                    for p, condition, na, nb in sample_metadata:
                        for tf in f['network_sources']:
                            exists = tf in scores.columns
                            activity.append({'study': st, 'view': view, 'method': method, 'patient': p,
                                             'condition': condition, 'FB01_n': na, 'otherFB_n': nb,
                                             'TF': tf, 'prior_axis': axis_of.get(tf, ''),
                                             'activity_delta_enrichment_t': scores.loc[p, tf] if exists else np.nan,
                                             'official_returned_network_statistic': official_values.loc[p, tf] if exists else np.nan,
                                             'official_return_kind': 'BH_q' if method == 'ulm' else 'RAW_P',
                                             'network_BH_q_within_patient_NOT_patient_P': qvalues.loc[p, tf] if exists else np.nan,
                                             'status': 'OFFICIAL_COMPLETE' if exists else 'INSUFFICIENT_TARGETS_NA'})
                except (np.linalg.LinAlgError, ValueError) as e:
                    errors.append({'study': st, 'view': view, 'method': method, 'status': 'OFFICIAL_TECHNICAL_NA_NO_SELECTIVE_TF_REMOVAL',
                                   'error_type': type(e).__name__, 'error': str(e), 'traceback': traceback.format_exc()})
                    for p, condition, na, nb in sample_metadata:
                        for tf in f['network_sources']:
                            activity.append({'study': st, 'view': view, 'method': method, 'patient': p, 'condition': condition,
                                             'FB01_n': na, 'otherFB_n': nb, 'TF': tf, 'prior_axis': axis_of.get(tf, ''),
                                             'activity_delta_enrichment_t': np.nan, 'official_returned_network_statistic': np.nan,
                                             'official_return_kind': 'BH_q' if method == 'ulm' else 'RAW_P',
                                             'network_BH_q_within_patient_NOT_patient_P': np.nan,
                                             'status': 'OFFICIAL_TECHNICAL_NA_NO_SELECTIVE_TF_REMOVAL'})
                print(json.dumps({'cohort': st, 'view': view, 'method': method, 'seconds': time.monotonic()-begin}), flush=True)
        del sums, raw, cp, X
    a = pd.DataFrame(activity); g = pd.DataFrame(gene_rows)
    # All 9 x15 axes, two views/methods, three cohorts. No n3/n4 association P.
    for (st, view, method), z in a[a.condition.eq('keloid')].groupby(['study', 'view', 'method']):
        for gene in CANDIDATES:
            gg = g[g.study.eq(st) & g.condition.eq('keloid') & g.gene.eq(gene)].set_index('patient')
            for tf in axis_of:
                zz = z[z.TF.eq(tf)].set_index('patient')
                both = gg[['log2_CP10K_plus0p1_delta']].join(zz[['activity_delta_enrichment_t']], how='inner').dropna()
                rho = spearmanr(both.iloc[:, 0], both.iloc[:, 1]).statistic if len(both) >= 3 and (both.nunique() > 1).all() else np.nan
                correlations.append({'study': st, 'view': view, 'method': method, 'gene': gene, 'TF': tf,
                                     'prior_axis': axis_of[tf], 'patients': len(both), 'Spearman_descriptive_only': rho,
                                     'association_P': np.nan, 'reason_no_P': 'n3_or4_shared_state_and_possible_coexpression_not_causality'})
    a.to_csv(out / 'all_TF_all_patient_all_view_method_activity.tsv.gz', sep='\t', index=False, na_rep='NA')
    pd.DataFrame(coverage).to_csv(out / 'all_TF_coverage_and_direct_overlap.tsv.gz', sep='\t', index=False, na_rep='NA')
    g.to_csv(out / 'all9_candidate_patient_expression_deltas.tsv', sep='\t', index=False, na_rep='NA')
    pd.DataFrame(groups_all).to_csv(out / 'all14_patient_group_support.tsv', sep='\t', index=False)
    pd.DataFrame(correlations).to_csv(out / 'all9candidate_15axisTF_correlations_DESCRIPTIVE.tsv', sep='\t', index=False, na_rep='NA')
    prior = net[net.target.isin(CANDIDATES)].copy()
    prior.to_csv(out / 'prior_TF_to_all9_candidate_signed_edges.tsv', sep='\t', index=False)
    summaries = []
    for key, z in a[a.condition.eq('keloid')].groupby(['study', 'view', 'method', 'TF'], sort=True):
        v = z.activity_delta_enrichment_t.to_numpy(float); finite = np.isfinite(v)
        summaries.append(dict(zip(['study', 'view', 'method', 'TF'], key),
                              prior_axis=axis_of.get(key[-1], ''), n_patients=len(v), finite_patients=int(finite.sum()),
                              mean_enrichment_t=float(np.mean(v[finite])) if finite.any() else np.nan,
                              median_enrichment_t=float(np.median(v[finite])) if finite.any() else np.nan,
                              positive=int((v[finite] > 0).sum()), negative=int((v[finite] < 0).sum()),
                              patient_inference_P=np.nan))
    pd.DataFrame(summaries).to_csv(out / 'all_TF_cohort_direction_summary.tsv.gz', sep='\t', index=False, na_rep='NA')
    atom(out / 'OFFICIAL_METHOD_TECHNICAL_FAILURES.json', errors)
    outputs = [receipt(p) for p in sorted(out.iterdir()) if p.is_file()]
    need(sum(r['bytes'] for r in outputs) <= CAP, 'Output exceeded declared64MiB')
    atom(out / 'COMPLETE.json', clean({'status': 'COMPLETE_REGULATORY_CONTEXT_ROOT_SCIENTIFIC_REVIEW_PENDING',
         'source_sha256': sha(__file__), 'freeze_sha256': sha(HERE / 'FREEZE.json'), 'outputs': outputs,
         'TF_activity_rows': len(a), 'all_candidates': list(CANDIDATES), 'method_technical_failures': errors,
         'inference_boundary': f['inference_boundary'], 'peak_RSS_bytes': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024}))


def run(approval):
    need(HERE.resolve().is_relative_to(Path('/root/five_new_studies_20261020')), 'System output only')
    grant = json.loads(Path(approval).read_text()); f = json.loads((HERE / 'FREEZE.json').read_text())
    need(grant['source_sha256'] == sha(__file__) == f['source_sha256'], 'Root/source disagreement')
    need(grant['freeze_sha256'] == sha(HERE / 'FREEZE.json'), 'Root/freeze disagreement')
    need(not (HERE / 'results').exists() and not (HERE / 'CLAIM.json').exists(), 'Existing claim/output; do not relaunch')
    r = resources(); own = float(grant.get('own_RSS_GiB', 4)); peer = float(grant.get('peer_unmaterialized_GiB', 30)); floor = 8
    need(r['effective_free_bytes'] >= (own+peer+floor)*2**30, 'Initial RAM reserve HOLD')
    need(r['system_free_bytes'] >= 5*2**30+CAP and r['data_free_bytes'] >= 5*2**30, 'Dual disk floor HOLD')
    fd = os.open(HERE / 'CLAIM.json', os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    os.write(fd, json.dumps({'supervisor': identity(os.getpid()), 'approval': receipt(approval), 'resources': r}).encode()); os.close(fd)
    start = time.monotonic()
    with (HERE / 'worker.log').open('x') as log:
        child = subprocess.Popen([sys.executable, str(Path(__file__).resolve()), '_worker'], stdout=log, stderr=subprocess.STDOUT)
        atom(HERE / 'RUNNING.json', {'supervisor': identity(os.getpid()), 'child': identity(child.pid), 'start_UTC': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())})
        while child.poll() is None:
            # Own child only. Never signal any peer; root fixes resource contention.
            try:
                fields = dict((l.split(':', 1)[0], l.split(':', 1)[1].strip()) for l in Path(f'/proc/{child.pid}/status').read_text().splitlines() if ':' in l)
                rss = int(fields.get('VmRSS', '0 kB').split()[0])*1024
                rr = resources()
                violation = rss > own*2**30 or rr['effective_free_bytes'] < floor*2**30 or rr['system_free_bytes'] < 5*2**30 or rr['data_free_bytes'] < 5*2**30
                violation |= time.monotonic()-start > float(grant.get('timeout_seconds', 3600))
                if violation:
                    child.terminate(); child.wait()
                    atom(HERE / 'RESOURCE_STOP.json', {'own_child': identity(child.pid) if Path(f'/proc/{child.pid}').exists() else {'pid': child.pid}, 'RSS': rss, 'resources': rr})
                    break
            except FileNotFoundError:
                pass
            time.sleep(3)
        rc = child.wait()
    atom(HERE / 'RETURN.json', {'child_pid': child.pid, 'returncode': rc, 'elapsed_seconds': time.monotonic()-start,
                               'COMPLETE_exists': (HERE / 'results/COMPLETE.json').exists(), 'resources': resources()})
    need(rc == 0 and (HERE / 'results/COMPLETE.json').exists(), 'Worker failed/incomplete; preserve whole scene')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('stage', choices=['freeze', 'run', '_worker']); ap.add_argument('--approval')
    args = ap.parse_args()
    if args.stage == 'freeze':
        freeze()
    elif args.stage == 'run':
        need(args.approval, 'Root approval required'); run(args.approval)
    else:
        worker(json.loads((HERE / 'FREEZE.json').read_text()))
