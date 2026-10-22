# Controlled public-treatment transcriptomes

**Input:** GSE232435, GSE7890, GSE232079; fixed gene programs.

**Analysis/output:** Seven treatment contrasts, matched state comparators and full gene-set tests.

**Paper connection:** Figure 5.

## Source files

- [`prepare_intervention_inputs_v1.py`](prepare_intervention_inputs_v1.py)
- [`run_intervention_pipeline_v1.py`](run_intervention_pipeline_v1.py)
- [`run_paired_intervention_analysis_v1.R`](run_paired_intervention_analysis_v1.R)
- [`run_state_marker_comparator_v1.py`](run_state_marker_comparator_v1.py)
- [`score_fixed_state_comparators_v1.R`](score_fixed_state_comparators_v1.R)

## Reproduction

Read [the execution guide](../docs/REPRODUCTION.md) before running these sources. They preserve the scientific algorithms and the original input-integrity checks. Historical absolute paths and local execution grants are documented execution scaffolding, not portable user commands. Configure an isolated working copy against downloaded inputs; do not remove scientific checks or overwrite completed results.

Completed supporting tables are in [`results/`](../results/). Complete families retain opposite directions and unavailable entries.
