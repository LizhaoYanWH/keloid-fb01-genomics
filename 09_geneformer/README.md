# Embedding perturbation sensitivity

**Input:** Frozen candidates, cells, reference axes and official Geneformer.

**Analysis/output:** All candidate/operator/raw/centered/NA responses.

**Paper connection:** Figure 4E-F.

## Source files

- [`tf_models/prepare_TF_followup_v1.py`](tf_models/prepare_TF_followup_v1.py)
- [`tf_models/run_TF_geneformer_v1.py`](tf_models/run_TF_geneformer_v1.py)
- [`tf_models/runtime_common_v1.py`](tf_models/runtime_common_v1.py)
- [`candidate_models/prepare_multimodel_v1.py`](candidate_models/prepare_multimodel_v1.py)
- [`candidate_models/run_geneformer_v1.py`](candidate_models/run_geneformer_v1.py)
- [`candidate_models/runtime_common_v1.py`](candidate_models/runtime_common_v1.py)
- [`rank_controls/prepare_genetic_rank4_v1.py`](rank_controls/prepare_genetic_rank4_v1.py)
- [`rank_controls/run_genetic_rankcontrol_GF_v1.py`](rank_controls/run_genetic_rankcontrol_GF_v1.py)
- [`rank_controls/runtime_common_v1.py`](rank_controls/runtime_common_v1.py)

## Reproduction

Read [the execution guide](../docs/REPRODUCTION.md) before running these sources. They preserve the scientific algorithms and the original input-integrity checks. Historical absolute paths and local execution grants are documented execution scaffolding, not portable user commands. Configure an isolated working copy against downloaded inputs; do not remove scientific checks or overwrite completed results.

Completed supporting tables are in [`results/`](../results/). Complete families retain opposite directions and unavailable entries.
