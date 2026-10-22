# RNA quality control and fibroblast-state correspondence

**Input:** Raw counts and original capture/lineage metadata.

**Analysis/output:** Frozen cell axes, pseudobulk profiles and 54 state correspondences.

**Paper connection:** Figure 1C-D; Figure 2A-D.

## Source files

- [`02_qc_prepare.py`](02_qc_prepare.py)
- [`05_freeze_and_pseudobulk.py`](05_freeze_and_pseudobulk.py)
- [`18_state_correspondence_and_purity.py`](18_state_correspondence_and_purity.py)

## Reproduction

Read [the execution guide](../docs/REPRODUCTION.md) before running these sources. They preserve the scientific algorithms and the original input-integrity checks. Historical absolute paths and local execution grants are documented execution scaffolding, not portable user commands. Configure an isolated working copy against downloaded inputs; do not remove scientific checks or overwrite completed results.

Completed supporting tables are in [`results/`](../results/). Complete families retain opposite directions and unavailable entries.
