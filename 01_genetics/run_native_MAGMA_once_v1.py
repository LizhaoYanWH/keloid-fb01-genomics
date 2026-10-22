#!/usr/bin/env python3
"""Root-gated, once-only native MAGMA annotation and full-family association.

No scDRS, new GWAS join, native-method replacement, or scientific result filter.
Default --plan is read-only. --run requires a separate, exact root approval.
"""
import argparse
import collections
import csv
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time

sys.dont_write_bytecode = True
ROOT = Path('/root/five_new_studies_20261020')
P = ROOT / 'projects/N02_keloid_genetic_cell_state'
S = P / 'engineering/MAGMA_inputs_v1'
M = S / 'whole_EUR_mapping_v1'
R = P / 'reviews/MAGMA_native_association_v1'
AUDIT = R / 'input_audit_v1/frozen_inputs_and_audit_v1.json'
AUDIT_SHA = '530cf483464db12bbeb89eff66905d9ac7e167b4c34054cab709466a8befb173'
LEDGER = R / 'input_audit_v1/all19427_gene_admission_ledger.tsv'
OUT = P / 'results/MAGMA_complete_EUR_native_v1'
APPROVAL = R / 'ROOT_APPROVAL_TO_RUN.json'
RUN_ID = 'N02_MAGMA_complete_EUR_native_v1'
GIB = 1024**3
RAM_CAP = 32 * GIB
MEMORY_FLOOR = 8 * GIB
OUTPUT_RESERVE = GIB
WALL_SECONDS = 6 * 3600
BINARY = S / 'MAGMALinux1.10_verified_unpacked/magma'


def need(ok, reason):
    if not ok:
        raise RuntimeError(reason)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(4 * 1024**2), b''):
            h.update(b)
    return h.hexdigest()


def canonical_sha(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def atomic_new(path, value):
    """No replacement of any prior marker, including a failed partial attempt."""
    raw = (json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n').encode()
    tmp = path.with_name(path.name + '.tmp')
    with tmp.open('xb') as f:
        f.write(raw)
        f.flush()
        os.fsync(f.fileno())
    os.link(tmp, path)
    tmp.unlink()
    fd = os.open(path.parent, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def commands():
    return {
        'annotate': [str(BINARY), '--annotate', 'window=10,10',
                     '--snp-loc', str(M / 'complete_mapped_EUR_SNP_location_GRCh37.tsv'),
                     '--gene-loc', str(S / 'gene37_verified_unpacked/NCBI37.3.gene.loc'),
                     '--out', str(OUT / 'native')],
        'association': [str(BINARY), '--bfile', str(S / 'EUR_LD_verified_unpacked/g1000_eur'),
                        '--pval', str(M / 'complete_mapped_EUR_SNP_P.tsv'),
                        'use=SNP,P', 'N=1282582', 'duplicate=drop',
                        '--gene-annot', str(OUT / 'native.genes.annot'),
                        '--gene-model', 'snp-wise=mean', '--seed', '20261020',
                        '--genes-only', '--out', str(OUT / 'complete_EUR')]
    }


def plan():
    return dict(run_id=RUN_ID, frozen_inputs_SHA256=AUDIT_SHA,
                argv=commands(), argv_canonical_SHA256=canonical_sha(commands()),
                threads=2, affinity='last two CPUs of inherited allowed set',
                RAM_cap_bytes=RAM_CAP, effective_memory_floor_bytes=MEMORY_FLOOR,
                required_start_available_bytes=RAM_CAP + MEMORY_FLOOR,
                both_disks_floor_bytes=5 * GIB, system_output_reserve_bytes=OUTPUT_RESERVE,
                max_wall_seconds=WALL_SECONDS, official_gene_family=19427,
                admitted_SNPs=8862187, constant_total_N=1282582,
                scDRS=False, no_automatic_retry=True,
                approval_path=str(APPROVAL), output=str(OUT))


def proc_stat(pid):
    text = Path(f'/proc/{pid}/stat').read_text()
    z = text[text.rfind(')') + 2:].split()
    return dict(state=z[0], pgrp=int(z[2]), start_ticks=z[19])


def resource_cut(pgrp=None):
    disk = {str(p): shutil.disk_usage(p).free for p in ('/', '/root/autodl-tmp')}
    mem = int(next(x.split()[1] for x in Path('/proc/meminfo').read_text().splitlines()
                   if x.startswith('MemAvailable:'))) * 1024
    cg = Path('/sys/fs/cgroup')
    limit = (cg / 'memory.max').read_text().strip()
    current = int((cg / 'memory.current').read_text())
    available = mem if limit == 'max' else min(mem, int(limit) - current)
    rss = 0
    identities = []
    for d in Path('/proc').iterdir():
        if not d.name.isdigit():
            continue
        try:
            pid = int(d.name)
            st = proc_stat(pid)
            if pid != os.getpid() and (pgrp is None or st['pgrp'] != pgrp):
                continue
            r = int(next(x.split()[1] for x in (d / 'status').read_text().splitlines()
                         if x.startswith('VmRSS:'))) * 1024
            rss += r
            identities.append(dict(PID=pid, start_ticks=st['start_ticks'], RSS_bytes=r))
        except (FileNotFoundError, ProcessLookupError, StopIteration, PermissionError):
            continue
    return dict(disk_free_bytes=disk, host_available_bytes=mem,
                cgroup_max_bytes=limit, cgroup_current_bytes=current,
                effective_available_bytes=available, owned_RSS_bytes=rss, processes=identities)


def guard(cut, start=False):
    need(min(cut['disk_free_bytes'].values()) >= 5 * GIB, 'Both disks must retain 5GiB')
    need(cut['disk_free_bytes']['/'] >= 5 * GIB + (OUTPUT_RESERVE if start else 0),
         'System disk output reserve not available')
    need(cut['effective_available_bytes'] >= (RAM_CAP + MEMORY_FLOOR if start else MEMORY_FLOOR),
         'Effective min(host,cgroup) memory guard')
    need(cut['owned_RSS_bytes'] < RAM_CAP, 'Owned process group RSS reached approved cap')


def frozen_inputs(audit, rehash=True):
    need(sha(AUDIT) == AUDIT_SHA, 'Frozen input receipt SHA differs')
    need(audit['status'] == 'COMPLETE_ENGINEERING_AUDIT_FROZEN_INPUTS_NO_ASSOCIATION',
         'Engineering input receipt is not complete')
    need(audit['status_counts'] == {'ADMITTED_UNDIRECTED_P': 8862187,
         'NO_REFERENCE_COORDINATE_ALLELE_MATCH': 896968,
         'AMBIGUOUS_SOURCE_DUPLICATE_ESTIMATES': 720}, 'Original family differs')
    checked = []
    for item in audit['files']:
        p = Path(item['path'])
        need(p.is_file() and not p.is_symlink(), 'Missing or substituted frozen input: ' + str(p))
        st = p.stat()
        need(st.st_size == item['bytes'] and st.st_mtime_ns == item['mtime_ns'],
             'Frozen input bytes/mtime differ: ' + str(p))
        if rehash:
            need(sha(p) == item['SHA256'], 'Frozen input SHA differs: ' + str(p))
            guard(resource_cut())
        checked.append(dict(path=str(p), bytes=st.st_size, SHA256=item['SHA256']))
    return checked


def terminate_owned(child):
    """Only this executor's still-live, new-session child group may be signalled."""
    if child is None or child.poll() is not None:
        return
    os.killpg(child.pid, signal.SIGTERM)
    try:
        child.wait(timeout=10)
    except subprocess.TimeoutExpired:
        os.killpg(child.pid, signal.SIGKILL)
        child.wait(timeout=10)


def native_phase(name, argv, run_start, peak):
    child = None
    log = OUT / (name + '.stdout_stderr.txt')
    started = time.monotonic()
    guard(resource_cut(), start=True)
    with log.open('xb') as f:
        try:
            child = subprocess.Popen(argv, cwd=OUT, stdin=subprocess.DEVNULL,
                                     stdout=f, stderr=subprocess.STDOUT, start_new_session=True)
            child_identity = dict(PID=child.pid, **proc_stat(child.pid))
            atomic_new(OUT / (name + '.launch.json'), dict(argv=argv, **child_identity))
            with (OUT / (name + '.resources.jsonl')).open('x') as monitor:
                while child.poll() is None:
                    cut = resource_cut(child.pid)
                    peak[0] = max(peak[0], cut['owned_RSS_bytes'])
                    monitor.write(json.dumps(dict(elapsed_seconds=time.monotonic() - run_start, **cut)) + '\n')
                    monitor.flush()
                    guard(cut)
                    need(time.monotonic() - run_start < WALL_SECONDS, 'Approved finite wall time reached')
                    time.sleep(1)
            need(child.returncode == 0, name + ' native returncode=' + str(child.returncode))
            result = dict(phase=name, argv=argv, returncode=child.returncode,
                          elapsed_seconds=time.monotonic() - started, identity=child_identity)
            atomic_new(OUT / (name + '.return.json'), result)
            return result
        except BaseException:
            terminate_owned(child)
            raise


def check_annotation(genes):
    found = {}
    with (OUT / 'native.genes.annot').open() as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            z = line.split()
            need(len(z) >= 3 and z[0] in genes and z[0] not in found, 'Native annotation gene identity/schema')
            g = genes[z[0]]
            chr_, left, right = z[1].split(':')
            expected_chr = '23' if g['chromosome'] == 'X' else g['chromosome']
            actual_chr = '23' if chr_ == 'X' else chr_
            need(actual_chr == expected_chr and max(1, int(left)) == int(g['window_start_1based'])
                 and int(right) == int(g['window_stop_1based']), 'Native annotation coordinate differs')
            need(len(z) - 2 == int(g['n_unique_admitted_SNP']) and len(set(z[2:])) == len(z) - 2,
                 'Native annotation coverage differs from frozen full family')
            found[z[0]] = len(z) - 2
    expected = {g for g, v in genes.items() if int(v['n_unique_admitted_SNP']) > 0}
    need(set(found) == expected, 'Native annotation silently omitted or added a covered gene')
    need(sum(found.values()) == 4719623, 'Complete annotation assignment total differs')
    return dict(annotated_genes=len(found), gene_SNP_assignments=sum(found.values()))


def parse_native_gene_output(path, known):
    header = None
    rows = {}
    with path.open() as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            z = line.split()
            if header is None:
                header = z
                need({'GENE', 'CHR', 'NSNPS', 'ZSTAT', 'P'} <= set(header), 'Native gene output schema differs')
                continue
            need(len(z) == len(header), 'Truncated native gene output row')
            row = dict(zip(header, z))
            gene = row['GENE']
            need(gene in known and gene not in rows, 'Unknown or duplicate native gene ID')
            expected_chr = '23' if known[gene]['chromosome'] == 'X' else known[gene]['chromosome']
            actual_chr = '23' if row['CHR'] == 'X' else row['CHR']
            need(actual_chr == expected_chr, 'Native output gene chromosome differs')
            # Preserve exact original values; nonfinite is a diagnostic, never a success filter.
            for key in ('P', 'ZSTAT'):
                if row[key].upper() not in ('NA', 'NAN', '.'):
                    float(row[key])
            int(row['NSNPS'])
            rows[gene] = row
    need(header is not None and rows, 'Native output has no association rows')
    return header, rows


def complete_ledger(genes):
    header, native = parse_native_gene_output(OUT / 'complete_EUR.genes.out', genes)
    counts = collections.Counter()
    diagnostic = collections.Counter()
    with (OUT / 'all19427_gene_results_and_NA_ledger.tsv').open('x', newline='') as f:
        cols = list(next(iter(genes.values()))) + ['native_output_status'] + ['native_' + x for x in header]
        writer = csv.DictWriter(f, fieldnames=cols, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for gene, g in genes.items():
            n = native.get(gene)
            if n is not None:
                status = 'NATIVE_RESULT_PRESENT_RAW_VALUES_RETAINED'
                for key in ('P', 'ZSTAT'):
                    try:
                        value = float(n[key])
                        diagnostic[key + ('_nonfinite' if not math.isfinite(value) else '_zero' if value == 0 else '_finite_nonzero')] += 1
                    except ValueError:
                        diagnostic[key + '_NA'] += 1
            else:
                status = g['native_preanalysis_status'] if g['native_preanalysis_status'] != 'ELIGIBLE_FOR_NATIVE_MODEL_INTERNAL_QC_PENDING' else 'NATIVE_OMITTED_AFTER_INTERNAL_PROCESSING_SEE_FULL_LOGS'
            counts[status] += 1
            writer.writerow({**g, 'native_output_status': status,
                             **{'native_' + key: n[key] if n else 'NA' for key in header}})
    return dict(official_gene_family=len(genes), native_output_genes=len(native),
                full_family_status_counts=dict(counts), raw_numeric_diagnostics=dict(diagnostic),
                statistical_or_project_acceptance=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', action='store_true')
    parser.add_argument('--plan', action='store_true')
    args = parser.parse_args()
    if not args.run:
        print(json.dumps(plan(), indent=2))
        return
    need(sys.platform == 'linux' and os.getenv('FIVE_STUDIES_COMPUTE') == 'autodl', 'Existing AutoDL only')
    need(Path(__file__).resolve() == R / 'run_native_MAGMA_once_v1.py', 'Exact approved remote source path only')
    approved = json.loads(APPROVAL.read_text())
    this_sha = sha(__file__)
    need(approved.get('approved') is True and approved.get('run_id') == RUN_ID
         and approved.get('executor_SHA256') == this_sha
         and approved.get('frozen_inputs_SHA256') == AUDIT_SHA
         and approved.get('plan_canonical_SHA256') == canonical_sha(plan()), 'Root approval exact binding is absent or differs')
    os.umask(0o077)
    allowed = sorted(os.sched_getaffinity(0))
    need(len(allowed) >= 2, 'Two allowed CPUs required')
    os.sched_setaffinity(0, set(allowed[-2:]))
    os.nice(10)
    for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        os.environ[key] = '2'
    lock = (ROOT / 'metadata/assess01_CPU_worker.lock').open('a')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    guard(resource_cut(), start=True)
    audit = json.loads(AUDIT.read_text())
    need(not OUT.exists(), 'Prior output, including any partial attempt, must be preserved')
    need(sha(LEDGER) == approved.get('gene_ledger_SHA256'), 'Full gene ledger SHA approval differs')
    checked = frozen_inputs(audit)
    genes = {g['official_gene_id']: g for g in csv.DictReader(LEDGER.open(), delimiter='\t')}
    need(len(genes) == 19427 and os.access(BINARY, os.X_OK), 'Full gene ledger or official executable missing')
    # O_EXCL claim remains even if the later run fails; there is no automatic second attempt.
    claim = R / 'native_once_launchclaim.json'
    atomic_new(claim, dict(run_id=RUN_ID, PID=os.getpid(), start_ticks=proc_stat(os.getpid())['start_ticks'],
                          source_SHA256=this_sha, frozen_inputs_SHA256=AUDIT_SHA, approval_SHA256=sha(APPROVAL)))
    OUT.mkdir()
    started = time.monotonic()
    peak = [0]
    atomic_new(OUT / 'RUNNING.json', dict(**plan(), PID=os.getpid(), identity=proc_stat(os.getpid()),
                                       source_SHA256=this_sha, effective_affinity=sorted(os.sched_getaffinity(0)),
                                       UTC=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), resource_cut=resource_cut()))
    shutil.copyfile(APPROVAL, OUT / 'ROOT_APPROVAL_USED.json')
    try:
        phases = [native_phase('annotate', commands()['annotate'], started, peak)]
        annotation = check_annotation(genes)
        atomic_new(OUT / 'annotation_structure_check.json', annotation)
        phases.append(native_phase('association', commands()['association'], started, peak))
        ledger = complete_ledger(genes)
        # Same full input SHA verification after both phases, before any COMPLETE marker.
        frozen_inputs(audit)
        need(sha(__file__) == this_sha and sha(LEDGER) == approved['gene_ledger_SHA256'], 'Source or gene ledger changed during execution')
        outputs = [dict(path=str(f), bytes=f.stat().st_size, SHA256=sha(f)) for f in sorted(OUT.iterdir()) if f.is_file()]
        atomic_new(OUT / 'COMPLETE.json', dict(status='COMPLETE_NATIVE_MAGMA_FULL_FAMILY_TECHNICAL_ONLY',
                   source_SHA256=this_sha, frozen_inputs_SHA256=AUDIT_SHA, approval_SHA256=sha(OUT / 'ROOT_APPROVAL_USED.json'),
                   plan=plan(), phases=phases, annotation=annotation, gene_ledger=ledger,
                   frozen_inputs_checked=checked, outputs=outputs, max_observed_owned_RSS_bytes=peak[0],
                   elapsed_seconds=time.monotonic() - started, final_resources=resource_cut(),
                   scDRS_run=False, original_GWAS_MAGMA_already_reported=True,
                   new_genetic_discovery_claim=False, root_scientific_acceptance_pending=True))
    except BaseException as error:
        atomic_new(OUT / 'FAILED_PRESERVED.json', dict(status='FAILED_PRESERVED_NO_AUTOMATIC_RETRY',
                   error=type(error).__name__, reason=str(error), elapsed_seconds=time.monotonic() - started,
                   source_SHA256=this_sha, frozen_inputs_SHA256=AUDIT_SHA, max_observed_owned_RSS_bytes=peak[0]))
        raise


if __name__ == '__main__':
    main()
