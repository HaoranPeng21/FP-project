# FP-project

**Faecalibacterium comparative genomics and microbiome ecology**

FP-project provides R analysis and visualization code for investigating the genomic
diversity, host-associated functions and ecological distribution of
*Faecalibacterium*. It connects genome phylogeny and metagenomic profiling
benchmarks with host and geographic distributions, early-life colonization,
co-abundance patterns, host phenotype associations, ABO/FUT2 comparisons,
gene-repertoire sharing and within-species genetic associations.

The repository contains **19 independent R analysis modules**, covering the study's
**six main figures and seven supplementary figures**. Each file in `analyses/`
corresponds to one analysis, with shared functions in `R/` and a single runner,
`run.R`. All 19 modules were executed successfully with the existing private
inputs; the checked outputs comprised 27 PDF panels or panel groups. See
[VALIDATION.md](VALIDATION.md) for the verification scope and source-version notes.

This is a **code-only repository**. Participant records, abundance matrices,
genome metadata, statistical result tables, trees, generated figures and private
configurations are not distributed. Public species names, pathway names and
schematic labels remain in the code where they explain an analysis. Input schemas
and private-data configuration are described in [INPUTS.md](INPUTS.md).

This is a **data-to-figure release**. It starts from the private catalogue,
annotation matrices, benchmark predictions and frozen statistical results listed
in [INPUTS.md](INPUTS.md). Genome assembly, sequence alignment, profiling,
pangenome construction and the upstream genome-level association scans are not
rerun. No Python scripts, shell programs or historical project scripts are called.
The code produces analysis panels and natural panel groups; the manuscript's
manual page arrangement is not reproduced pixel for pixel.

## Run

Use R 4.5.1 or a compatible environment; see [DEPENDENCIES.md](DEPENDENCIES.md).

```bash
Rscript install_dependencies.R
Rscript validate_code.R
Rscript run.R --list
```

Keep the private data outside this repository. Set the project root and a separate
private output directory. `R/inputs.R` gives paths relative to the data root.

```bash
export FP_DATA_ROOT=/path/to/private/project
export FP_OUTPUT_DIR=/path/to/private/figure_outputs
export FP_CONFIG_FILE=/path/to/private/config.local.R
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
Rscript run.R all
```

Copy `config.example.R` outside the repository and replace its example paths.
Three inputs located outside the original project tree require explicit overrides:
the ABO abundance matrix, ABO/FUT2 table and pairwise gene-prevalence table.
Any other input key can also be overridden, so the original directory structure
is optional. Missing inputs cause an explicit error; data are never downloaded or
replaced with synthetic values.

Run one analysis or a subset by its numeric prefix:

```bash
Rscript run.R 02
Rscript run.R 03 04 05
```

Optional settings: `FP_PNG=0` disables optional PNGs, `FP_FONT` selects a local
font, and `FP_R_LIBS` supplies additional R library directories separated by the
platform path separator. The default font is the portable `sans` family.

Each run writes PDFs, result summaries, `run_receipt.tsv`, `sessionInfo.txt` and a
private input-path/checksum receipt under `FP_OUTPUT_DIR`. The output directory
must be outside the repository. These outputs can contain private scientific
results or paths and are not part of this code-only release.

## Figure-to-code map

| Figure/panels | R script | Starting point |
|---|---|---|
| 1A | `01_catalogue_tree.R` | bac120 tree and frozen 17,103-tip annotations |
| 1B–D; S2D–E | `02_benchmark.R` | strict83/shared24 predictions and frozen 901-run selection |
| 2A | `03_pathway_prevalence.R` | genome counts for all 5,256 pathway comparisons |
| 2B | `04_host_distribution.R` | 19-species × 12-host genome counts |
| 2C; S3 | `05_animal_trees.R` | four species subtrees and full 998-tip selected overview |
| 3A | `06_human_distribution.R` | audited host-stratum summaries and 14-species order |
| 3B–C | `07_infant_trajectories.R` | 10-species timepoint summaries, infant-bootstrap intervals and frozen two-part GEE tests |
| 3D | `08_infant_associations.R` | selected adjusted longitudinal-model coefficients |
| 4 | `09_coabundance.R` | bac120 partner tree, cohort direction counts and frozen REML meta-results |
| 5A | `10_phenotype_meta.R` | frozen phenotype meta-analysis, full prevalence and species order |
| 5B–F | `11_abo.R` | abundance, ABO/FUT2, covariates and GalNAc counts |
| 6 | `12_within_species.R` | frozen genome scans and PC5 logistic/Firth gene-panel tables |
| S2A–C | `13_catalogue_qc.R` | catalogue quality metrics, rarefaction estimates and SGB mapping |
| S5 | `14_geography.R` | country counts, boundaries and frozen country prevalence |
| S6A–B | `15_gene_sharing.R` | 25-species repertoire Jaccard, representatives and bac120 tree |
| S6C–F | `16_gene_pair_prevalence.R` | shared-family prevalence for highlighted pairs 1, 3, 4, 7 |
| S7 | `17_maast_tree_permanova.R` | Maast trees and genome-tip continent annotations |
| S4 | `18_host_pathways.R` | retained nominal host-associated pathway screen |
| S1 | `19_taxonomy_history.R` | public timeline labels and catalogue release crosswalk |

Script filenames are relative to `analyses/`. The supplementary numbering above
follows the supplied PDFs; S5 is a country distribution, S6 gene sharing and S7
Maast phylogeography. S1 is an illustrative R reconstruction: its original drawing
script was not found. Its release crosswalk uses the modal catalogue label tuple
for each displayed lineage and does not assert that every genome has that tuple.

## Preserved analysis definitions

- **Benchmark:** all83, same21 and shared24 are separate scopes. The supplied PDF's
  S2D–E titles say same21, but its TP/FN/FP/TN and Figure 1C metrics identify
  **shared24** (1,300 targets; 21,624 cells per tool). The R code follows those
  numerical inputs and labels the comparison shared24. The 901-run selection is read unchanged; redundant original
  group-23 simulations are excluded upstream. Zero predictions are displayed at
  a labelled positive proxy on logarithmic axes and remain zero in scoring.
- **Pathways:** Figure 2A recomputes Fisher tests and BH across all 5,256 comparisons.
  The retained display is 14 nonredundant profiles plus two shared backgrounds;
  stars additionally require an absolute prevalence difference of at least 25
  percentage points. Figure S4 retains its nominal screen and inherited global
  q-values; no correction is performed on the selected rows alone.
- **Distribution:** human bubbles use prevalence for area and mean abundance
  among detected samples for colour, retaining all 14 eligible species in 3A.
  The supplied 3B–C uses the earlier ten-species two-part analysis: C shows the
  arithmetic mean among detected samples, with bootstrap CIs; stars retain the
  binomial-logit detection and log-abundance Gaussian GEE tests, with infant
  clustering and endpoint-wise ten-species BH correction. The later 14-species
  zero-inclusive trajectory revision is not the supplied figure. S5 uses countries, at least 50 samples per country, and
  overall prevalence strictly greater than 5%.
- **Co-abundance:** colour is the positive direction fraction among cohorts
  with nominal P<0.05. It is not a proportion significant in all cohorts. Thin
  and thick frames retain nominal meta-P and global BH-q meanings respectively.
- **Phenotypes:** Figure 5A retains upstream phenome-wide q-values, at least
  two contributing studies and the full-sample prevalence filter. Display
  filtering never defines a new FDR family.
- **ABO:** CLR uses all 33 named Faecalibacterium species and pseudocount 1e-6;
  tests and plots are carriers-only. Four-group pairwise Wilcoxon tests use one
  84-test BH family. HC3 models adjust age, sex, BMI and log10 read depth; the
  three-comparator and pooled contrasts retain 42- and 14-test families.
- **Within species:** Figure 6 preserves the three supplied cases and 26 shown
  gene families. Family q-values are inherited from PC5 logistic/Firth tests
  across the complete eligible family set per comparison. A physical-neighbourhood
  carriage bracket is descriptive and is not a cluster-level q-value.
- **Gene sharing:** this is the 25-species pangenome subset, not all 83 species.
  Native R recomputes representative distances, Jaccard checks and a 99,999-label
  QAP permutation test. The original frozen QAP P is retained in the figure;
  both that P and the R recomputation are exported. R and NumPy use different
  random-number generators, so finite-permutation P-values can differ slightly.
- **Maast:** PERMANOVA, PERMDISP and Fitch tests are recomputed in R with 1,999
  permutations and the original seeds. Known-continent groups require at least
  five genomes. BH is across four species separately for each test. Significant
  PERMDISP remains visible alongside PERMANOVA.

The implementation and local verification record are described in
[VALIDATION.md](VALIDATION.md). Raw data access remains governed by the owners of
the original datasets; this repository does not provide access to those data.
