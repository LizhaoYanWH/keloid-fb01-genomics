"""Fixed external RNA transport: input-only freeze, then root-approved mapping->score.

No new QC/clustering. The six prior-QC captures and discovery state labels stay fixed.
"""
from pathlib import Path
import argparse, collections, csv, gc, gzip, hashlib, json, math, multiprocessing as mp, os, re, resource, shutil, sys, time, traceback, warnings

P = Path('/root/five_new_studies_20261020/projects/N02_keloid_genetic_cell_state')
OLD = Path('/root/autodl-tmp/Five_Studies_Execution_20261020')
HERE = P / 'extensions_20261021/external_RNA_v1'
PREP = P / 'reviews/scDRS_production_preparation_v1'
RF = PREP / 'RNA_freeze_v1'
GS = PREP / 'official_MAGMA_GS_v1'
COHORTS = ('GSE220300', 'GSE335482')
CAPTURES = ('GSM4994379_KL1', 'GSM4994380_KL2', 'GSM4994381_KL3', 'GSM4994382_NS1', 'GSM4994383_NS2', 'GSM4994384_NS3')
EXPECTED_SINGLETS = (7726, 7265, 6896, 6715, 7434, 5413)
PYTHON = P / 'tools/scdrs104_env/bin/python'
CAP = 128 * 2**20


def need(ok, why):
    if not ok:
        raise RuntimeError(why)


def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda: f.read(2**20), b''):
            h.update(b)
    return h.hexdigest()


def atom(p, obj):
    q = Path(str(p) + '.tmp')
    q.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False) + '\n')
    q.replace(p)


def receipt(p):
    p = Path(p)
    return {'path': str(p), 'bytes': p.stat().st_size, 'sha256': sha(p)}


def verify(r):
    p = Path(r['path'])
    need(p.stat().st_size == r['bytes'] and sha(p) == r['sha256'], 'Input changed: ' + str(p))


def resources():
    mem = {r.split(':')[0]: int(r.split()[1]) * 1024 for r in Path('/proc/meminfo').read_text().splitlines()}
    used = int(Path('/sys/fs/cgroup/memory.current').read_text())
    maximum = int(Path('/sys/fs/cgroup/memory.max').read_text())
    return {'effective_free_bytes': min(mem['MemAvailable'], maximum - used),
            'system_free_bytes': shutil.disk_usage('/root').free, 'data_free_bytes': shutil.disk_usage('/root/autodl-tmp').free,
            'owned_peak_RSS_bytes': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024}


def guard(job=24, start=False, peer=12):
    r = resources()
    need(r['effective_free_bytes'] >= ((job if start else 0) + 8 + peer) * 2**30, 'Job+8floor+peer RAM guard')
    need(r['owned_peak_RSS_bytes'] <= job * 2**30, 'Owned RSS cap')
    need(r['data_free_bytes'] >= 5 * 2**30 and r['system_free_bytes'] >= 5 * 2**30 + CAP, 'Both5GiB; system128MiB reserve')
    return r


def strings(node):
    import numpy as np
    import h5py
    if isinstance(node, h5py.Dataset):
        return np.asarray([v.decode() if isinstance(v, bytes) else str(v) for v in node[:]])
    categories = strings(node['categories']); codes = node['codes'][:]
    return np.asarray([categories[k] if k >= 0 else '' for k in codes])


def hindex(group):
    return strings(group[group.attrs['_index']])


def write_complete(out, status, inputs, extra=None):
    atom(out / 'COMPLETE.json', {'status': status, 'source_sha256': sha(__file__),
                                'freeze_sha256': sha(HERE / 'FREEZE.json'), 'input_bindings': inputs,
                                **(extra or {}), 'outputs': {p.name: receipt(p) for p in out.iterdir()
                                                         if p.is_file() and p.name != 'COMPLETE.json'}})


def freeze():
    import numpy as np
    import pandas as pd
    import h5py
    need(not (HERE / 'FREEZE.json').exists(), 'Freeze exists; no overwrite')
    guard(job=2, start=True, peer=32)
    HERE.mkdir(parents=True, exist_ok=True)
    own = OLD / 'projects/03_susceptibility'
    fb_source = own / 'results/mature_scar_state_bridge_v1/mature_FB_ownQC_state_mapping.h5ad'
    need(sha(fb_source) == '8f256ccf91f1646a6967ccd5b20aa5c481b8a3eaae9f753b9afb2a8654268de9', 'Prior FB identity changed')
    with h5py.File(fb_source) as h:
        original_fb = set(hindex(h['obs']))
    need(len(original_fb) == 12249, 'Prior12249 FB required')
    hgnc = OLD / 'shared_inputs/programs_01_03/sources/hgnc_complete_set.txt'
    need(sha(hgnc) == '12363916be36f8f6f91879b0faade4eeaac973722de7136ec61329b77f2c2f91', 'Original HGNC changed')
    h = pd.read_csv(hgnc, sep='\t', dtype=str, keep_default_na=False)
    h = h.loc[h.status.eq('Approved') & h.ensembl_gene_id.ne('')]
    unique_ens = h.ensembl_gene_id.value_counts()
    gene_rows = h.loc[h.ensembl_gene_id.map(unique_ens).eq(1)]
    canonical = dict(zip(gene_rows.ensembl_gene_id, gene_rows.entrez_id))
    gs_complete = json.loads((GS / 'COMPLETE.json').read_text())
    gp = GS / 'Keloid_EUR.top1000_Z.gs'
    need(sha(gp) == gs_complete['files'][gp.name]['sha256'], 'Original GS changed')
    original_gs = dict(v.rsplit(':', 1) for v in pd.read_csv(gp, sep='\t').GENESET.iloc[0].split(','))
    need(len(original_gs) == 1000, 'Original1000GS required')
    trait_ens = set(h.loc[h.entrez_id.isin(original_gs), 'ensembl_gene_id'])
    candidate_path = P / 'extensions_20261021/gene_state_bridge_v1/results/all_primary_PP4_ge08_unique_genes.tsv'
    candidate_complete = P / 'extensions_20261021/gene_state_bridge_v1/results/COMPLETE.json'
    cc = json.loads(candidate_complete.read_text())
    need(sha(candidate_path) == cc['files'][candidate_path.name]['sha256'], 'Original complete9 candidates changed')
    candidates = pd.read_csv(candidate_path, sep='\t', dtype=str, keep_default_na=False).to_dict('records')
    candidate_genes = {r['gene_id'] for r in candidates}
    need(len(candidate_genes) == 9, 'Keep complete original9 candidates')
    manifests = []; captures = []; inputs = [fb_source, hgnc, gp, GS / 'COMPLETE.json', RF / 'COMPLETE.json', own / 'results/mature_scar_QC_v1/COMPLETE.json', candidate_path, candidate_complete]
    axis = None
    for capture, expected in zip(CAPTURES, EXPECTED_SINGLETS):
        raw = own / f'inputs/processed/GSE163973/{capture}/raw.h5ad'
        qp = own / f'results/scRNA_qc/{capture}/cell_qc_audit.tsv.gz'
        dp = own / f'results/scRNA_qc/{capture}/scDblFinder_calls.tsv'
        q = pd.read_csv(qp, sep='\t', dtype=str, keep_default_na=False).set_index('cell_id')
        d = pd.read_csv(dp, sep='\t', dtype=str, keep_default_na=False).set_index('cell_id')
        with h5py.File(raw) as f:
            cells = hindex(f['obs']); genes = hindex(f['var'])
            need(len(cells) == len(set(cells)), 'Raw cell IDs duplicated')
            need(f['X'].attrs['encoding-type'] == 'csr_matrix', 'Raw CSR required')
            need(len(genes) == 33538 and (axis is None or np.array_equal(axis, genes)), 'Original complete33538 axis differs')
            axis = genes
            raw_dtype = str(f['X/data'].dtype)
            raw_conditions = sorted(set(strings(f['obs/condition'])))
            expected_condition = 'keloid' if '_KL' in capture else 'normal_scar'
            QC_conditions = sorted(set(q.condition))
            need(raw_conditions == [expected_condition] and QC_conditions == raw_conditions, 'Frozen KL/ordinary-scar condition coding differs')
            need(set(cells) == set(q.index) and set(d.index) == set(q.index[q.preliminary_pass.eq('True')]), 'Prior full QC/preliminary-only doublet identities differ')
            keep = q.reindex(cells).preliminary_pass.eq('True') & d.reindex(cells).scDblFinder_class.eq('singlet')
            need(int(keep.sum()) == expected, 'Prior QC singlet count differs')
            qm = q.loc[cells[keep.to_numpy()]].copy()
            qm['is_own_QC_FB'] = qm.index.isin(original_fb)
            need(np.array_equal(qm.is_own_QC_FB.to_numpy(), qm.author_cell_type.eq('Fibroblast').to_numpy()), 'Frozen author-FB/ownQC identity differs')
            qm['source_capture_directory'] = capture
            manifests.append(qm.reset_index())
        captures.append({'capture': capture, 'raw_path': str(raw), 'QC_path': str(qp), 'doublet_path': str(dp), 'singlets': expected, 'raw_X_data_dtype': raw_dtype, 'raw_condition_values': raw_conditions, 'QC_condition_values': QC_conditions})
        inputs += [raw, qp, dp]
    manifest = pd.concat(manifests, ignore_index=True)
    need(len(manifest) == 41449 and manifest.cell_id.is_unique and int(manifest.is_own_QC_FB.sum()) == 12249, 'Full41449/FB12249 changed')
    need(set(manifest.loc[manifest.is_own_QC_FB, 'cell_id']) == original_fb, 'Frozen FB cells differ')
    manifest.to_csv(HERE / 'frozen41449_cell_manifest.tsv.gz', sep='\t', index=False)
    rn = json.loads((RF / 'COMPLETE.json').read_text())
    features, ledger, references = {}, [], {}
    for st in COHORTS:
        marker = OLD / f'projects/01_regional/versions/lineage_v2/results/expansion_readouts/{st}/FB_state_top100_annotation_markers.tsv'
        raw = Path(rn['cohorts'][st]['input']['FB_state_path'])
        need(sha(raw) == rn['cohorts'][st]['input']['FB_state_sha256'], 'Frozen discovery FB counts changed')
        with h5py.File(raw) as f:
            ref_axis = set(hindex(f['var']))
        m = pd.read_csv(marker, sep='\t', dtype={'names': str, 'gene_symbol': str})
        states = rn['cohorts'][st]['FB_state_family']; chosen = set()
        for state in states:
            subset = m.loc[m.group.eq(state)]
            old_top50 = [r for r in subset.itertuples() if not re.match(r'^(MT-|RPL|RPS|HB[ABDEGMQZ])', str(r.gene_symbol))][:50]
            for r in old_top50:
                g = str(r.names)
                reason = 'ORIGINAL_GS_EXCLUDED' if g in trait_ens else 'ORIGINAL9_CANDIDATE_EXCLUDED' if g in candidate_genes else 'HGNC_UNRESOLVED_EXCLUDED' if g not in canonical else 'EXTERNAL_ASSAY_MISSING' if g not in set(axis) else 'REFERENCE_ASSAY_MISSING' if g not in ref_axis else 'MAPPING_FEATURE'
                ledger.append({'reference': st, 'state': state, 'ensembl_gene_id': g, 'symbol': str(r.gene_symbol), 'status': reason})
                if reason == 'MAPPING_FEATURE':
                    chosen.add(g)
        need(len(chosen) >= 30 and not(chosen & trait_ens), 'Original mapper technical>=30 nonGS features unavailable')
        features[st] = sorted(chosen)
        references[st] = {'raw_FB_h5ad': str(raw), 'metadata_path': str(RF / st / 'all_final_cell_metadata.tsv.gz'),
                          'states': states, 'documented_patients': rn['cohorts'][st]['independent_patients'], 'marker_path': str(marker)}
        inputs += [raw, marker, RF / st / 'all_final_cell_metadata.tsv.gz']
    pd.DataFrame(ledger).to_csv(HERE / 'all_reference_top50_nonGS_feature_ledger.tsv', sep='\t', index=False)
    fallback = set(r['ensembl_gene_id'] for r in ledger if r['reference'] == COHORTS[0] and r['state'] == 'FB01' and r['status'] == 'MAPPING_FEATURE') & set(r['ensembl_gene_id'] for r in ledger if r['reference'] == COHORTS[1] and r['state'] == 'FB01' and r['status'] == 'MAPPING_FEATURE')
    identity = P / 'extensions_20261021/gene_state_bridge_v1/SOURCE_IDENTITY_PREFLIGHT_v1.json'
    inputs += [identity, PREP / 'environment_install_receipt_v1.json', HERE / 'frozen41449_cell_manifest.tsv.gz', HERE / 'all_reference_top50_nonGS_feature_ledger.tsv']
    env = json.loads((PREP / 'environment_install_receipt_v1.json').read_text())
    site = P / 'tools/scdrs104_env/lib/python3.12/site-packages'
    for r in env['verified_source_files']:
        if r['path'].startswith('scdrs/') and len(Path(r['path']).parts) == 2:
            need(sha(site / r['path']) == r['sha256'], 'Pinned official scDRS source changed')
            inputs.append(site / r['path'])
    import sklearn, sklearn.linear_model._logistic as sklearn_logistic, sklearn.utils.class_weight as sklearn_class_weight
    inputs += [Path(sklearn_logistic.__file__), Path(sklearn_class_weight.__file__)]
    frozen = {'status': 'INPUT_ONLY_FROZEN_EXTERNAL_RNA_NO_MAPPING_OR_SCORE', 'source_sha256': sha(__file__),
              'captures': captures, 'all_cells': 41449, 'FB_cells': 12249, 'references': references,
              'mapping_features': features, 'fallback_shared_FB01_nonGS_module': sorted(fallback),
              'fallback_not_automatically_run': True, 'canonical_Ensembl_to_Entrez': canonical,
              'original_GS_weights': {g: float(w) for g, w in original_gs.items()},
              'original9_candidates': candidates,
              'mapping_recipe': 'Original all-state nonMT/RP/HB top50 then remove original1000GS and all9 described candidates before external trait scoring; logistic C.5 balanced, donor-equal base sample weights; per-reference probability>=.5 and two-reference FB01 consensus; no retraining from external labels',
              'LODO_stop_rule': 'Both discovery references must recover FB01 with donor-equal aggregate precision greater than donor-equal FB01 prevalence and recall>0; report every heldout donor including all rare labels. This is label-transport calibration, not independent truth.',
              'LODO_feature_selection_limitation': 'Features/states were originally selected using all discovery donors. This fixed-panel conditional LODO is not nested feature selection, leakage-free accuracy or independent validation.',
              'external_stop_rule': 'No consensusFB01-like cells in any KL patient with a nonempty otherFB comparator -> stop state replication before score; keep Unknown/all6patient support rows',
              'score_recipe': 'One whole41449 official scDRS, all original1000 weights, fixed min50+uniqueEntrez, full raw UMI/nGene cov before filtering, admitted-count normalization1e4+log1p,1000controls,seed20261020,vs,adj_propNone',
              'family': 'Exactly two MC association groups(all ownQC authorFB; consensusFB01-like), one BH2. Six patient score differences descriptive; KL3 main context, NS3 separately. No clinical3vs3 P.',
              'not_blind': 'GSE163973 previously used in old03; ordinary scars are not healthy skin; current discovery source input identity checked separately',
              'sklearn_version': sklearn.__version__,
              'class_weight_diagnostic': 'Existing sklearn1.9.1 balanced logistic passes base sample_weight into compute_class_weight; report official weighted-count formula as a diagnostic, not bitwise per-iteration optimizer weights.',
              'output_cap_bytes': CAP, 'inputs': [receipt(p) for p in dict.fromkeys(inputs)]}
    atom(HERE / 'FREEZE.json', frozen)
    print(json.dumps({'status': frozen['status'], 'FREEZE_sha256': sha(HERE / 'FREEZE.json'),
                      'features': {k: len(v) for k, v in features.items()}, 'fallback_features': len(fallback)}))


def feature_matrix(path, ids, features, node):
    import h5py
    import numpy as np
    from scipy.sparse import csr_matrix
    wanted = {v: i for i, v in enumerate(ids)}
    result = np.empty((len(ids), len(features)), dtype=np.float32); seen = set()
    with h5py.File(path) as h:
        cells = hindex(h['obs']); genes = hindex(h['var']); lookup = {g: i for i, g in enumerate(genes)}
        ix = [lookup[g] for g in features]; z = h[node]; ptr = z['indptr']
        for start in range(0, len(cells), 512):
            end = min(start + 512, len(cells)); p = ptr[start:end + 1]
            block = csr_matrix((z['data'][int(p[0]):int(p[-1])], z['indices'][int(p[0]):int(p[-1])], p - p[0]), shape=(end - start, len(genes)))
            lib = np.asarray(block.sum(1)).ravel(); sub = block[:, ix].toarray().astype(np.float64)
            for j, cell in enumerate(cells[start:end]):
                if cell in wanted:
                    need(lib[j] > 0, 'Zero full library')
                    result[wanted[cell]] = np.log1p(sub[j] / lib[j] * 1e4).astype(np.float32); seen.add(cell)
            guard()
    need(seen == set(ids), 'Mapping cell coverage incomplete')
    return result


def mapping(frozen, out):
    import numpy as np
    import pandas as pd
    from sklearn.linear_model import LogisticRegression
    from sklearn.exceptions import ConvergenceWarning
    from sklearn.utils.class_weight import compute_class_weight
    import joblib
    manifest = pd.read_csv(HERE / 'frozen41449_cell_manifest.tsv.gz', sep='\t', dtype=str, keep_default_na=False).set_index('cell_id')
    fb = manifest.loc[manifest.is_own_QC_FB.eq('True')].copy(); predicted = fb[['patient_id', 'condition']].copy()
    diagnostic = []; all_weights = []; gates = []
    for st in COHORTS:
        r = frozen['references'][st]; features = frozen['mapping_features'][st]
        meta = pd.read_csv(r['metadata_path'], sep='\t', dtype=str, keep_default_na=False).set_index('cell_id')
        meta = meta.loc[meta.independent_patient_id.isin(r['documented_patients']) & meta.FB_state.isin(r['states'])]
        A = feature_matrix(r['raw_FB_h5ad'], list(meta.index), features, 'layers/counts')
        labels = meta.FB_state.to_numpy(); donor = meta.independent_patient_id.to_numpy()
        calibrations = []
        for heldout in [*r['documented_patients'], 'FULL_REFERENCE']:
            train = donor != heldout if heldout != 'FULL_REFERENCE' else np.ones(len(donor), dtype=bool)
            known = sorted(set(donor[train])); n = int(train.sum())
            donor_counts = collections.Counter(donor[train])
            sw = np.asarray([n / (len(known) * donor_counts[p]) for p in donor[train]], dtype=float)
            model = LogisticRegression(C=.5, class_weight='balanced', max_iter=1000, random_state=20261020, n_jobs=1)
            with warnings.catch_warnings(record=True) as captured:
                warnings.simplefilter('always'); model.fit(A[train], labels[train], sample_weight=sw)
            need(not any(issubclass(w.category, ConvergenceWarning) for w in captured), 'Original finite logistic fit did not converge; no solver switch')
            cw = dict(zip(model.classes_, compute_class_weight('balanced', classes=model.classes_, y=labels[train], sample_weight=sw)))
            for p in known:
                ix = donor[train] == p
                all_weights.append({'reference': st, 'heldout': heldout, 'train_patient': p,
                                    'base_weight_sum': float(sw[ix].sum()),
                                    'weighted_count_formula_class_adjusted_weight_sum_diagnostic': float(sum(w * cw[s] for w, s in zip(sw[ix], labels[train][ix]))),
                                    'patient_equal_before_additional_class_balancing': True})
            if heldout != 'FULL_REFERENCE':
                test = ~train; prob = model.predict_proba(A[test]); best = prob.argmax(1)
                accepted = np.where(prob.max(1) >= .5, model.classes_[best], 'Unknown')
                truth = labels[test] == 'FB01'; guess = accepted == 'FB01'; tp = int((truth & guess).sum())
                rec = {'reference': st, 'heldout_patient': heldout, 'n_cells': int(test.sum()),
                       'true_FB01': int(truth.sum()), 'accepted_FB01': int(guess.sum()), 'true_positive': tp,
                       'FB01_precision': tp / int(guess.sum()) if guess.any() else None,
                       'FB01_recall': tp / int(truth.sum()) if truth.any() else None,
                       'accepted_any_state_fraction': float((accepted != 'Unknown').mean()),
                       'classes_absent_from_training': ','.join(sorted(set(r['states']) - set(model.classes_))),
                       'LODO_label_reproduction_not_external_truth': True,
                       'fixed_panel_selected_using_all_discovery_donors_not_nested': True}
                diagnostic.append(rec); calibrations.append(rec)
            else:
                joblib.dump(model, out / f'{st}_fixed_nonGS_mapper.joblib')
                B = np.empty((len(fb), len(features)), dtype=np.float32)
                for c in frozen['captures']:
                    ids = list(fb.loc[fb.source_capture_directory.eq(c['capture'])].index)
                    ix = fb.index.get_indexer(ids); B[ix] = feature_matrix(c['raw_path'], ids, features, 'X')
                prob = model.predict_proba(B); best = prob.argmax(1)
                predicted[st + '_raw_prediction'] = model.classes_[best]
                predicted[st + '_max_probability'] = prob.max(1)
                predicted[st + '_accepted_prediction'] = np.where(prob.max(1) >= .5, model.classes_[best], 'Unknown')
                del B
        tp = sum(v['true_positive'] / v['n_cells'] for v in calibrations)
        guess = sum(v['accepted_FB01'] / v['n_cells'] for v in calibrations)
        true = sum(v['true_FB01'] / v['n_cells'] for v in calibrations)
        precision = tp / guess if guess else 0.; recall = tp / true if true else 0.; prevalence = true / len(calibrations)
        gates.append({'reference': st, 'donor_equal_precision': precision, 'donor_equal_recall': recall,
                      'donor_equal_prevalence': prevalence, 'gate_PASS': precision > prevalence and recall > 0,
                      'rule': 'Predeclared target-state transport gate; no macroF1 cutoff or external trait score used'})
        del A; gc.collect()
    a = predicted[COHORTS[0] + '_accepted_prediction']; b = predicted[COHORTS[1] + '_accepted_prediction']
    predicted['transport_label'] = np.where(a.eq('FB01') & b.eq('FB01'), 'FB01_consensus',
                                          np.where(a.ne('Unknown') & b.ne('Unknown') & a.ne('FB01') & b.ne('FB01'), 'other_mapped_FB', 'Unknown'))
    predicted.to_csv(out / 'all12249_two_reference_state_transport.tsv.gz', sep='\t')
    pd.DataFrame(diagnostic).to_csv(out / 'all8_discovery_LODO_FB01_diagnostics.tsv', sep='\t', index=False, na_rep='NA')
    pd.DataFrame(all_weights).to_csv(out / 'all_reference_training_patient_weights.tsv', sep='\t', index=False)
    support = []
    for patient in sorted(fb.patient_id.unique()):
        sub = predicted.loc[predicted.patient_id.eq(patient)]; n = int(sub.transport_label.eq('FB01_consensus').sum())
        support.append({'patient': patient, 'condition': str(sub.condition.iloc[0]), 'all_FB': len(sub),
                        'FB01_consensus': n, 'other_FB_including_Unknown': len(sub) - n,
                        'Unknown': int(sub.transport_label.eq('Unknown').sum())})
    pd.DataFrame(support).to_csv(out / 'all6_patient_mapping_support.tsv', sep='\t', index=False)
    supported_KL = sum(v['condition'] == 'keloid' and v['FB01_consensus'] > 0 and v['other_FB_including_Unknown'] > 0 for v in support)
    passed = all(v['gate_PASS'] for v in gates) and supported_KL > 0
    atom(out / 'MAPPING_TRANSPORT_GATES.json', {'status': 'PASS_TO_FIXED_SCORE' if passed else 'STOP_STATE_REPLICATION_NO_SCORE',
                                             'reference_gates': gates, 'supported_KL_patients': supported_KL,
                                             'all6_patient_support': support, 'no_threshold_tuning_or_reclustering': True})
    return passed


def score(frozen, mapping_dir, out):
    import numpy as np
    import pandas as pd
    import h5py, scipy.sparse as sp
    import anndata as ad, scanpy as sc, scdrs
    manifest = pd.read_csv(HERE / 'frozen41449_cell_manifest.tsv.gz', sep='\t', dtype=str, keep_default_na=False).set_index('cell_id')
    pieces = []; axes = None; metadata = []
    for capture in frozen['captures']:
        with h5py.File(capture['raw_path']) as h:
            ids = hindex(h['obs']); genes = hindex(h['var']); z = h['X']
            need(axes is None or np.array_equal(axes, genes), 'Complete raw gene axis changed')
            axes = genes
            X = sp.csr_matrix((z['data'][:], z['indices'][:], z['indptr'][:]), shape=tuple(z.attrs['shape']))
            use = np.asarray([cell in manifest.index for cell in ids]); ids = ids[use]
            pieces.append(X[use].copy()); metadata.extend(ids)
    X = sp.vstack(pieces, format='csr'); del pieces; gc.collect(); X.sum_duplicates(); X.eliminate_zeros()
    need(len(metadata) == len(manifest) == 41449 and list(manifest.index) == metadata, 'All41449 cells/order changed')
    need(np.issubdtype(X.dtype, np.integer) and (X.data >= 0).all(), 'Original raw integer counts required')
    U = np.asarray(X.sum(1)).ravel(); N = np.diff(X.indptr)
    old_u = pd.to_numeric(manifest.total_counts).to_numpy(dtype=np.float64)
    old_n = pd.to_numeric(manifest.n_genes_by_counts).to_numpy(dtype=np.float64)
    need(np.array_equal(old_u, np.rint(old_u)) and np.array_equal(old_n, np.rint(old_n)) and
         np.array_equal(U, old_u.astype(np.int64)) and np.array_equal(N, old_n.astype(np.int64)), 'Original full raw QC totals differ')
    transfer = pd.read_csv(mapping_dir / 'all12249_two_reference_state_transport.tsv.gz', sep='\t', index_col=0, dtype=str, keep_default_na=False)
    labels = transfer.transport_label.reindex(manifest.index).fillna('NOT_FB')
    ownfb = manifest.is_own_QC_FB.eq('True').to_numpy(); consensus = labels.eq('FB01_consensus').to_numpy()
    # Original9 fixed genetic candidates, including NEDD4/TAB2: raw full-library
    # measurements, separate from later normalized scDRS and no new patient P.
    lookup = {g: i for i, g in enumerate(axes)}; descriptions = []; differences = []
    for pid in sorted(manifest.patient_id.unique()):
        pm = manifest.patient_id.eq(pid).to_numpy()
        for candidate in frozen['original9_candidates']:
            gene = candidate['gene_id']; pair = {}
            for name, mask in [('FB01_consensus', pm & consensus), ('other_FB_including_Unknown', pm & ownfb & ~consensus)]:
                n = int(mask.sum()); lib = int(U[mask].sum()); measured = gene in lookup
                if measured and n:
                    values = np.asarray(X[mask, lookup[gene]].toarray()).ravel()
                    raw = int(values.sum()); detected = int((values > 0).sum())
                    cp10k = raw / lib * 1e4 if lib else None
                    mean = float((values / U[mask] * 1e4).mean()); fraction = detected / n
                else:
                    raw = detected = cp10k = mean = fraction = None
                rec = {'patient': pid, 'condition': str(manifest.loc[pm, 'condition'].iloc[0]),
                       'group': name, 'gene_id': gene, 'approved_symbol': candidate['approved_symbol'],
                       'n_cells': n, 'full_raw_library_UMI': lib, 'gene_raw_UMI': raw,
                       'pseudobulk_CP10K': cp10k, 'mean_cell_CP10K': mean, 'detected_cells': detected,
                       'detection_fraction': fraction,
                       'status': 'ASSAYED' if measured and n else 'NO_GROUP_CELLS_NA' if measured else 'NOT_ASSAYED_NA'}
                descriptions.append(rec); pair[name] = rec
            a, b = pair['FB01_consensus'], pair['other_FB_including_Unknown']
            d = {'patient': pid, 'condition': a['condition'], 'gene_id': gene,
                 'approved_symbol': candidate['approved_symbol'], 'n_FB01': a['n_cells'], 'n_other_FB': b['n_cells'],
                 'descriptive_no_patient_P': True}
            for metric in ['pseudobulk_CP10K', 'mean_cell_CP10K', 'detection_fraction']:
                d[metric + '_difference_FB01_minus_otherFB'] = a[metric] - b[metric] if a[metric] is not None and b[metric] is not None else None
            differences.append(d)
    pd.DataFrame(descriptions).to_csv(out / 'all9_candidates_all6patients_two_group_raw_expression.tsv', sep='\t', index=False, na_rep='NA')
    pd.DataFrame(differences).to_csv(out / 'all9_candidates_all6patients_FB01_expression_contrasts.tsv', sep='\t', index=False, na_rep='NA')
    canonical = [frozen['canonical_Ensembl_to_Entrez'].get(g, '') for g in axes]
    frequency = collections.Counter(g for g in canonical if g); detection = np.asarray((X > 0).sum(0)).ravel()
    use = np.asarray([bool(g) and frequency[g] == 1 and detection[i] >= 50 for i, g in enumerate(canonical)])
    admitted = [g for g, keep in zip(canonical, use) if keep]
    cov = pd.DataFrame({'log1p_full_input_UMI': np.log1p(U), 'log1p_full_input_nGenes': np.log1p(N)}, index=manifest.index)
    a = ad.AnnData(X=X[:, use].astype(np.float64), obs=manifest.copy(), var=pd.DataFrame(index=admitted)); del X; gc.collect()
    need((np.asarray(a.X.sum(1)).ravel() > 0).all(), 'Zero admitted library; do not delete cells')
    sc.pp.normalize_per_cell(a, counts_per_cell_after=1e4); sc.pp.log1p(a)
    scdrs.preprocess(a, cov=cov, adj_prop=None, n_mean_bin=20, n_var_bin=20, n_chunk=128)
    original = scdrs.method._select_ctrl_geneset; calls = []
    def capture(*args, **kwargs):
        lists, weights = original(*args, **kwargs)
        with gzip.open(out / 'actual1000_control_members_and_weights.tsv.gz', 'wt') as f:
            f.write('control_id\tposition\tcanonical_entrez\tinput_weight\n')
            for k in sorted(lists):
                for pos, (g, w) in enumerate(zip(lists[k], weights[k])):
                    f.write(f'{k}\t{pos}\t{g}\t{w:.17g}\n')
        calls.append({'n_ctrl': len(lists), 'measured': len(args[1]), 'seed': args[6]}); return lists, weights
    scdrs.method._select_ctrl_geneset = capture
    try:
        y = scdrs.score_cell(a, list(frozen['original_GS_weights']), gene_weight=list(frozen['original_GS_weights'].values()), ctrl_match_key='mean_var', n_ctrl=1000, weight_opt='vs', return_ctrl_raw_score=False, return_ctrl_norm_score=True, random_seed=20261020)
    finally:
        scdrs.method._select_ctrl_geneset = original
    need(y.shape == (41449, 1006) and y.index.equals(manifest.index) and np.isfinite(y.to_numpy()).all(), 'Wholecohort score finite/coverage contract failed')
    need(len(calls) == 1 and calls[0]['n_ctrl'] == 1000 and calls[0]['seed'] == 20261020, 'Official1000control single-call contract failed')
    y.iloc[:, :6].to_csv(out / 'all41449_six_cell_scores.tsv.gz', sep='\t')
    pd.DataFrame({'source_ensembl': axes, 'canonical_entrez': canonical, 'detected_cells': detection, 'admitted': use}).to_csv(out / 'all33538_gene_admission.tsv.gz', sep='\t', index=False)
    rows = []; controls = [f'ctrl_norm_score_{k}' for k in range(1000)]
    for name, mask in [('all_ownQC_author_FB', ownfb), ('FB01_consensus', consensus)]:
        if mask.any():
            q = float(np.quantile(y.loc[mask, 'norm_score'], .95)); null = np.quantile(y.loc[mask, controls], .95, axis=0)
            p = float(((null >= q).sum() + 1) / 1001)
            rows.append({'group': name, 'n_cells': int(mask.sum()), 'MC_P': p, 'MC_Z': float((q - null.mean()) / null.std()) if null.std() else None, 'not_patient_risk_P': True})
        else:
            rows.append({'group': name, 'n_cells': 0, 'MC_P': None, 'MC_Z': None, 'not_patient_risk_P': True})
    from statsmodels.stats.multitest import multipletests
    q = multipletests([r['MC_P'] if r['MC_P'] is not None else 1 for r in rows], method='fdr_bh')[1]
    for row, value in zip(rows, q):
        row['BH2_q'] = float(value) if row['MC_P'] is not None else None
    pd.DataFrame(rows).to_csv(out / 'fixed2_external_group_MC_BH2.tsv', sep='\t', index=False, na_rep='NA')
    patient = []
    for pid in sorted(manifest.patient_id.unique()):
        pm = manifest.patient_id.eq(pid).to_numpy(); sm = pm & consensus; other = pm & ownfb & ~consensus
        patient.append({'patient': pid, 'condition': str(manifest.loc[pm, 'condition'].iloc[0]), 'n_FB01': int(sm.sum()), 'n_other_FB_including_Unknown': int(other.sum()),
                        'mean_difference_FB01_minus_otherFB': float(y.loc[sm, 'norm_score'].mean() - y.loc[other, 'norm_score'].mean()) if sm.any() and other.any() else None,
                        'descriptive_only_not_new_genotype_risk_or_clinical_P': True})
    pd.DataFrame(patient).to_csv(out / 'all6_patient_FB01_context_score_descriptions.tsv', sep='\t', index=False, na_rep='NA')
    atom(out / 'OFFICIAL_SCORE_DIAGNOSTICS.json', {'all_final_cells': 41449, 'original_requested_genes': 1000, 'measurable_original_genes': calls[0]['measured'], 'n_ctrl': 1000, 'seed': 20261020, 'selection_calls': calls, 'nonfinite': 0, 'no_whole_control_cache_persisted': True, 'ordinary_scar_not_healthy': True, 'not_blind_external_material': True})


def pipeline(frozen, j, approval_sha):
    start = time.monotonic(); mapping_dir = HERE / 'mapping_v1'; score_dir = HERE / j.get('score_directory', 'score_v1')
    try:
        need(re.fullmatch(r'score_v[0-9]+', score_dir.name) is not None and not score_dir.exists(), 'Explicit new score phase directory required')
        reuse = j.get('reuse_mapping_COMPLETE_sha256')
        if reuse:
            need(sha(mapping_dir / 'COMPLETE.json') == reuse, 'Root reuse mapping SHA differs')
            m = json.loads((mapping_dir / 'COMPLETE.json').read_text())
            need(m['status'] == 'COMPLETE_MAPPING_READY' and m['source_sha256'] == sha(__file__) and m['freeze_sha256'] == sha(HERE / 'FREEZE.json'), 'Cannot reuse partial or differently frozen mapping')
            for r in m['outputs'].values():
                verify(r)
            passed = True
        else:
            need(not mapping_dir.exists(), 'No mapping overwrite/refit')
            mapping_dir.mkdir(); passed = mapping(frozen, mapping_dir)
            write_complete(mapping_dir, 'COMPLETE_MAPPING_READY' if passed else 'COMPLETE_MAPPING_SCIENTIFIC_STOP', frozen['inputs'])
        if not passed:
            atom(HERE / 'SCIENTIFIC_STOP_NO_EXTERNAL_SCORE.json', {'reason': 'Fixed target-state LODO/consensus mapping gate failed', 'mapping_COMPLETE_sha256': sha(mapping_dir / 'COMPLETE.json'), 'no_signature_or_threshold_change': True}); return
        guard(); need(time.monotonic() - start <= j['max_wall_seconds'], 'Mapping wall budget exhausted')
        score_dir.mkdir(); score(frozen, mapping_dir, score_dir)
        need(sum(p.stat().st_size for p in HERE.rglob('*') if p.is_file()) <= CAP, '128MiB complete output cap')
        guard(); need(time.monotonic() - start <= j['max_wall_seconds'], 'Finite pipeline wall cap')
        write_complete(score_dir, 'COMPLETE_ONE_WHOLE41449_EXTERNAL_SCORE_FIXED2GROUPS', frozen['inputs'], {'mapping_COMPLETE_sha256': sha(mapping_dir / 'COMPLETE.json')})
        atom(HERE / 'COMPLETE.json', {'status': 'COMPLETE_FIXED_EXTERNAL_RNA_MAPPING_AND_SCORE_ROOT_REVIEW_PENDING', 'source_sha256': sha(__file__), 'freeze_sha256': sha(HERE / 'FREEZE.json'), 'approval_sha256': approval_sha, 'mapping_COMPLETE_sha256': sha(mapping_dir / 'COMPLETE.json'), 'score_COMPLETE_sha256': sha(score_dir / 'COMPLETE.json'), 'elapsed_seconds': time.monotonic() - start, 'resources': resources()})
    except BaseException as e:
        atom(HERE / ('FAILED_PRESERVED_' + j.get('attempt_tag', 'v1') + '.json'), {'error': repr(e), 'traceback': traceback.format_exc(), 'source_sha256': sha(__file__)}); raise


def run(args):
    need(sys.executable == str(PYTHON), 'Exact existing scDRS environment required')
    ap = Path(args.approval); j = json.loads(ap.read_text())
    need(sha(ap) == args.approval_sha256 and j['action'] == 'RUN_FIXED_EXTERNAL_RNA_MAPPING_THEN_ONE_SCORE', 'Root action/SHA binding required')
    need(j['source_sha256'] == sha(__file__) and j['freeze_sha256'] == sha(HERE / 'FREEZE.json') and j['C_scientific_support_reviewed'] is True, 'Source/freeze/C scientific gate not root-approved')
    c = P / 'extensions_20261021/regional_exclusion_v1'
    need(sha(c / 'results_v1/COMPLETE.json') == j['C_COMPLETE_sha256'] and json.loads((c / 'C_ACTUAL_RETURN.json').read_text())['returncode'] == 0 and json.loads((c / 'results_v1/COMPLETE.json').read_text())['returncode'] == 0, 'C not complete/approved')
    frozen = json.loads((HERE / 'FREEZE.json').read_text())
    need(frozen['source_sha256'] == sha(__file__), 'Frozen executor changed')
    for r in frozen['inputs']:
        verify(r)
    need(j['threads'] == 4 and j['RAM_cap_GiB'] == 24 and j['max_wall_seconds'] <= 7200, 'Root finite resource recipe changed')
    for name in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
        os.environ[name] = '4'
    quota, period = Path('/sys/fs/cgroup/cpu.max').read_text().split()
    need(len(os.sched_getaffinity(0)) >= 4 and (quota == 'max' or int(quota) / int(period) >= 4), 'Insufficient4CPU quota')
    os.sched_setaffinity(0, sorted(os.sched_getaffinity(0))[-4:]); guard(start=True)
    tag = j.get('attempt_tag', 'v1'); need(re.fullmatch(r'v[0-9]+', tag) is not None, 'Explicit root attempt tag required')
    with (HERE / ('ONCE_CLAIM_' + tag + '.json')).open('x') as f:
        json.dump({'PID': os.getpid(), 'source_sha256': sha(__file__), 'approval_sha256': sha(ap)}, f)
    proc = mp.get_context('spawn').Process(target=pipeline, args=(frozen, j, sha(ap)))
    proc.start(); start = time.monotonic(); problem = None
    atom(HERE / ('ACTUAL_CHILD_' + tag + '.json'), {'PID': proc.pid, 'parent_PID': os.getpid()})
    while proc.is_alive():
        proc.join(5)
        try:
            rr = resources()
            need(rr['effective_free_bytes'] >= 8 * 2**30, 'Runtime8GiB floor')
            need(min(rr['system_free_bytes'], rr['data_free_bytes']) >= 5 * 2**30, 'Runtime disk5GiB floor')
            need(time.monotonic() - start <= j['max_wall_seconds'], 'Runtime wall cap')
            text = Path(f'/proc/{proc.pid}/status').read_text()
            rss = next((int(z.split()[1]) * 1024 for z in text.splitlines() if z.startswith('VmRSS:')), 0)
            need(rss <= 24 * 2**30, 'Runtime owned24GiB RSS cap')
            need(sum(p.stat().st_size for p in HERE.rglob('*') if p.is_file()) <= CAP, 'Runtime128MiB output cap')
        except FileNotFoundError:
            pass
        except BaseException as e:
            problem = repr(e); proc.terminate(); proc.join(30); break
    atom(HERE / ('ACTUAL_RETURN_' + tag + '.json'), {'returncode': proc.exitcode, 'error': problem, 'source_sha256': sha(__file__), 'approval_sha256': sha(ap)})
    need(problem is None and proc.exitcode == 0, 'Owned external run failed; partial phase retained, no automatic retry')


def main():
    p = argparse.ArgumentParser(); p.add_argument('mode', choices=['freeze', 'run']); p.add_argument('--expected-source-sha256'); p.add_argument('--approval'); p.add_argument('--approval-sha256'); args = p.parse_args()
    need(P.exists() and Path('/proc').exists(), 'Existing AutoDL only')
    if args.mode == 'freeze':
        for name in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
            os.environ[name] = '1'
        os.sched_setaffinity(0, [min(os.sched_getaffinity(0))])
        need(args.expected_source_sha256 == sha(__file__), 'Expected freeze source SHA required'); freeze()
    else:
        need(args.approval and args.approval_sha256, 'Root scientific grant required'); run(args)


if __name__ == '__main__':
    main()
