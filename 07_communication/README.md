# Paired receiver responses and NicheNet

**Input:** Frozen discovery counts and labels; official ligand-target prior.

**Analysis/output:** Full paired gene responses; all 1,226-ligand rankings and receptor support.

**Paper connection:** Table 2; communication supporting results.

## Source files

- [`01_freeze_inputs.py`](01_freeze_inputs.py)
- [`02_aggregate_profiles.py`](02_aggregate_profiles.py)
- [`03_paired_response_and_NicheNet_resume_v5.R`](03_paired_response_and_NicheNet_resume_v5.R)
- [`audit_candidate_complex_TF_edges_v1.py`](audit_candidate_complex_TF_edges_v1.py)
- [`common.py`](common.py)
- [`run_TGFbeta_existing_response_calibration_v1.R`](run_TGFbeta_existing_response_calibration_v1.R)

## Reproduction

Read [the execution guide](../docs/REPRODUCTION.md) before running these sources. They preserve the scientific algorithms and the original input-integrity checks. Historical absolute paths and local execution grants are documented execution scaffolding, not portable user commands. Configure an isolated working copy against downloaded inputs; do not remove scientific checks or overwrite completed results.

Completed supporting tables are in [`results/`](../results/). Complete families retain opposite directions and unavailable entries.
