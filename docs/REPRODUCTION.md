# Reproducing the analysis

## Start with the completed tables

1. Open `config/workflow_steps.tsv` for the study order and the paper figure served by each stage.
2. `results/` contains complete primary tables and selected completed intermediate families. Their file contents are unchanged from the production analyses. Consult `config/source_files.tsv` for provenance and `FILE_MANIFEST.json` for checksums.
3. Run `python validate_repository.py` to check integrity, Python syntax, forbidden image files and table headers. This check does not run scientific analyses.

## Obtain inputs and software

Download inputs from the official accessions listed in `config/datasets.tsv`. Raw counts, full nominal QTL rows, the European LD reference, HGNC mapping, NicheNet priors, CollecTRI, the CellOracle promoter prior and Geneformer weights must be obtained under their original access terms. Do not substitute a differently filtered input or reconstruct values from rounded manuscript figures.

Use separate R/Python environments for scDRS, communication and the perturbation models. `config/environment_versions.tsv` reports versions confirmed in the study record. It is not a complete dependency lock. `docs/METHODS.md` specifies the scientific preprocessing, candidate families, estimators and controls.

## Follow the numbered stages

The numbered folders preserve executed scientific drivers. Some drivers retain production paths such as `/root/five_new_studies_20261020`, input-hash freezes and local approval/resource checks. They are source snapshots, not a claim that a fresh clone runs without configuration. Prepare a separate working copy with an explicit input-path map and freshly verified source hashes. Obtain original frozen input/annotation records before repeating an input-dependent step; do not generate dummy approval files or weaken identity checks to make it run.

The native MAGMA commands are recorded in `01_genetics/run_native_MAGMA_once_v1.py` (`commands()`). The coloc source backend is recorded by official commit in `config/analysis_parameters.json`; install/obtain that backend under its licence. Subsequent scDRS runs score each entire cohort once, then reuse the cache across the fixed state and patient views. Do not score each state separately and merge normalized outputs.

The full earlier numeric reproduction record remains in the [archived supporting-data release](https://github.com/Jaqen00/keloid-fb01-genomics/releases/tag/v1.0.0). This is a historical data-source link, not the maintained code directory. It provides full output families beyond the compact tables included here. A clean raw-input-to-final-figure rerun has not been certified; source inspection, table integrity and end-to-end computational validation are distinct checks.

## Figures

`figures/` contains plotting and Adobe composition code only. No generated PNG, PDF, SVG, AI, JPEG or other image files are stored in the repository. The scientific plotting modules consume completed numeric tables; they do not refit models. Adobe conceptual-panel sources require artwork provided separately by the author; that artwork is not a measured dataset. No image-service credential or hosted image API helper is included.

## Scientific interpretation

Cell-level gene-set null tests, patient descriptions, reference-dependent embedding outputs and expression-network simulations address different questions. Simulations are hypotheses about defined readouts, not measured treatment efficacy. Keep the complete control families, opposite directions, zeros and NA entries when adapting the workflow.
