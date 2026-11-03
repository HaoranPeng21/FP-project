# FP-project

R code for studying the genomic diversity, functions and ecology of
*Faecalibacterium*. The final study catalogue contains **17,103 genomes and
83 species**. Individual analyses use the relevant prespecified subsets.

Each analysis is one independent R script in `analyses/`. Open a script in
RStudio, set `data_dir` and `out_dir` at the top, and run it from top to bottom.
There is no shared function file, configuration file or master runner.

```r
data_dir <- "/path/to/private/data/01_catalogue_tree"
out_dir <- "/path/to/private/results/01_catalogue_tree"
```

Alternatively, run a script directly with `Rscript`; `FP_DATA_DIR` and
`FP_OUTPUT_DIR` can supply these two paths. Package installation is described in
[DEPENDENCIES.md](DEPENDENCIES.md); input filenames and formats are in
[INPUTS.md](INPUTS.md).

| Figure | Independent R script |
|---|---|
| 1A | [01_catalogue_tree.R](analyses/01_catalogue_tree.R) |
| 1B-D; S2D-E | [02_benchmark.R](analyses/02_benchmark.R) |
| 2A | [03_pathway_prevalence.R](analyses/03_pathway_prevalence.R) |
| 2B | [04_host_distribution.R](analyses/04_host_distribution.R) |
| 2C; S3 | [05_animal_trees.R](analyses/05_animal_trees.R) |
| 3A | [06_human_distribution.R](analyses/06_human_distribution.R) |
| 3B-C | [07_infant_trajectories.R](analyses/07_infant_trajectories.R) |
| 3D | [08_infant_associations.R](analyses/08_infant_associations.R) |
| 4 | [09_coabundance.R](analyses/09_coabundance.R) |
| 5A | [10_phenotype_meta.R](analyses/10_phenotype_meta.R) |
| 5B-F | [11_abo.R](analyses/11_abo.R) |
| 6 | [12_within_species.R](analyses/12_within_species.R) |
| S2A-C | [13_catalogue_qc.R](analyses/13_catalogue_qc.R) |
| S5 | [14_geography.R](analyses/14_geography.R) |
| S6A-B | [15_gene_sharing.R](analyses/15_gene_sharing.R) |
| S6C-F | [16_gene_pair_prevalence.R](analyses/16_gene_pair_prevalence.R) |
| S7 | [17_maast_tree_permanova.R](analyses/17_maast_tree_permanova.R) |
| S4 | [18_host_pathways.R](analyses/18_host_pathways.R) |
| S1 | [19_taxonomy_history.R](analyses/19_taxonomy_history.R) |

This repository distributes **code and documentation only**. Private data,
participant/sample records, trees, statistical tables and generated figures
are not included. Input filenames are neutral; the original directory layout
is not required.

The scripts start from the input tables, matrices and trees. For analyses that
use frozen model results, those results are explicit inputs; the scripts do
not repeat upstream sequence processing or model fitting. Scientific scopes
and validation details are summarized in [VALIDATION.md](VALIDATION.md).
