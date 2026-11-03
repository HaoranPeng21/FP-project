# Copy this file OUTSIDE the repository and set FP_CONFIG_FILE to that copy.
# FP_DATA_ROOT contains the private project tree described in R/inputs.R.
# Use absolute paths only in your private copy; do not commit that copy.
input_overrides <- list(
  abo_abundance = "/path/to/private/abo_abundance.tsv",
  abo_phenotype = "/path/to/private/abo_phenotype.tsv",
  gene_pair_prevalence = "/path/to/private/species_pair_shared_prevalence.tsv"
)
# Any key in R/inputs.R can be overridden in the same way.
