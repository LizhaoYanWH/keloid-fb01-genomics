# Keloid FB01 genomics

A reader-oriented code and result companion for **Scar-associated genetic signals map to a repair-linked fibroblast state in keloid single-cell transcriptomes**.

The workflow asks where a gene set selected from population genetics is expressed within keloid fibroblast heterogeneity, then places the associated FB01 state in normal-repair, regulatory and controlled-treatment contexts.

## Analysis roadmap

| Step | Folder | What it answers |
|---|---|---|
| 1 | [01_genetics](01_genetics/) | Which genes are prioritized independently of lesion RNA? |
| 2 | [02_cell_states](02_cell_states/) | How are fibroblasts and comparable RNA states defined? |
| 3 | [03_genetic_scoring](03_genetic_scoring/) | Where does the fixed genetic gene-set signal localize? |
| 4 | [04_validation](04_validation/) | Does the association persist after fixed exclusions and external RNA mapping? |
| 5 | [05_normal_repair](05_normal_repair/) | How does the same program behave during normal healing? |
| 6 | [06_regulatory_context](06_regulatory_context/) | Which candidates and TF footprints characterize the state? |
| 7 | [07_communication](07_communication/) | Which extracellular signals are compatible with its expression response? |
| 8 | [08_celloracle](08_celloracle/) | What do finite TF simulations predict for fixed response programs? |
| 9 | [09_geneformer](09_geneformer/) | How sensitive are embedding predictions to candidates, operators and references? |
| 10 | [10_treatments](10_treatments/) | How do fixed programs respond in public controlled-treatment datasets? |

Each stage README lists **inputs → analysis → outputs → paper figures**. Start with the first seven stages to follow the biological argument; the model folders are explicitly computational follow-up.

## Find what you need

- **Read the analysis:** [methods](docs/METHODS.md) and [execution guide](docs/REPRODUCTION.md).
- **Locate data:** [public dataset list](config/datasets.tsv).
- **Inspect findings:** [result index](docs/RESULTS_INDEX.md) and [completed numeric tables](results/).
- **Generate figures:** [plotting/composition source only](figures/).
- **Check the package:** `python validate_repository.py` (integrity and syntax; no analysis).

The repository intentionally contains **no image files, manuscript drafts, model-review conversations, credentials, cloud orchestration or queue/retry scripts**. Download original matrices, priors and model weights from their official sources. Complete source snapshots retain the algorithms used; production paths and frozen-input checks require configuration before a new run. Full clean-environment re-execution has not been certified.

## Scope

The primary result is genetic gene-set localization to a reproducible fibroblast expression context. Normal repair, communication and virtual perturbation describe biological context and testable hypotheses. They do not establish an inherited genotype in RNA donors, a validated causal axis or clinical treatment efficacy. All supplied complete result families retain opposite directions and NA entries.

Maintainer: Lizhao Yan. [Resource and reuse terms](RIGHTS_AND_RESOURCES.md). The earlier full numerical archive is linked from the execution guide for traceability; this LizhaoYanWH repository is the organized code companion.
