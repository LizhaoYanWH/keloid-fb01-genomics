# Genetic gene-set localization

**Input:** Fixed MAGMA1000-Z; complete RNA cohorts.

**Analysis/output:** scDRS scores, 1,000 controls, full 47-test family and patient descriptions.

**Paper connection:** Figure 1E-F; Figure 2E-H.

## Source files

- [`consume_fixed47_FDR_patient_descriptive_v1.py`](consume_fixed47_FDR_patient_descriptive_v1.py)
- [`freeze_RNA_inputs_v2.py`](freeze_RNA_inputs_v2.py)
- [`scDRS_official_API_once_v2.py`](scDRS_official_API_once_v2.py)

## Reproduction

Read [the execution guide](../docs/REPRODUCTION.md) before running these sources. They preserve the scientific algorithms and the original input-integrity checks. Historical absolute paths and local execution grants are documented execution scaffolding, not portable user commands. Configure an isolated working copy against downloaded inputs; do not remove scientific checks or overwrite completed results.

Completed supporting tables are in [`results/`](../results/). Complete families retain opposite directions and unavailable entries.
