# Whole-GWAS gene prioritization and colocalization

**Input:** GCST90652488, European LD, three GTEx contexts.

**Analysis/output:** 19,427-gene ledger; fixed 1,000 genes; complete three-prior coloc family.

**Paper connection:** Figure 1B; Figure 4A.

## Source files

- [`audit_and_freeze_inputs_v1.py`](audit_and_freeze_inputs_v1.py)
- [`coloc_complete_family_driver_v1.R`](coloc_complete_family_driver_v1.R)
- [`freeze_official_MAGMA_GS_v1.py`](freeze_official_MAGMA_GS_v1.py)
- [`run_native_MAGMA_once_v1.py`](run_native_MAGMA_once_v1.py)

## Reproduction

Read [the execution guide](../docs/REPRODUCTION.md) before running these sources. They preserve the scientific algorithms and the original input-integrity checks. Historical absolute paths and local execution grants are documented execution scaffolding, not portable user commands. Configure an isolated working copy against downloaded inputs; do not remove scientific checks or overwrite completed results.

Completed supporting tables are in [`results/`](../results/). Complete families retain opposite directions and unavailable entries.
