# Validation and analysis scope

All 19 standalone R scripts were executed with the existing private inputs;
the run produced 27 PDF panels or panel groups. Each script can be run directly,
without sourcing another project script. Local logs, input aliases, checksums
and generated figures remain outside this repository.

The final inclusion list and main-tree annotation table independently contain
**17,103 unique genomes and 83 species**. Runtime checks enforce these counts.
Catalogue QC and taxonomic-history summaries use that same inclusion list.
QC shows completeness/contamination for the final included genomes.

The following analysis definitions are retained:

- Benchmark: the frozen 901-run selection; the full 83-species and shared24
  comparison scopes remain separate. Shared24 follows the supplied figure's
  numerical values, correcting its stale same21 title.
- Pathways: Fisher tests and BH correction over the complete 5,256-test family;
  the displayed subset does not define the correction family.
- Infant trajectories: the supplied ten-species two-part GEE analysis;
  abundance means use detected samples, with the original bootstrap intervals
  and endpoint-wise correction families.
- Co-abundance: direction ratios count only nominally significant cohorts.
  Frames retain the frozen REML P-values and global BH q-values.
- Phenotype meta-analysis: retain the original phenome-wide q-values before
  filtering the display; the script reads the frozen model results.
- ABO/FUT2: carriers-only CLR over all 33 named genus rows, pseudocount 1e-6;
  all six pairwise contrasts per displayed species share one BH family.
  HC3 models adjust age, sex, BMI and log10 read depth, without a cohort term.
  Pooled and three-comparator models keep their own correction families.
- Within-species associations: the supplied PC5 BMI, tricyclic-antidepressant
  and PPI cases. Brackets describe physical-neighbourhood carriage; gene stars
  use the original member-family q-values.
- Gene sharing: the 25-species pangenome subset and 300 dependent pairs.
  QAP permutes species labels jointly. The native R and original NumPy runs
  have different random generators; both P-values are exported, and the
  supplied figure's verified P-value is retained in its label.
- Maast: 1,999 permutations, original R seeds, known continents with at least
  five genomes, and BH correction across the four species per test.
- Host pathways: the retained nominal Fisher screen is descriptive; stars
  use upstream global q-values.
- Taxonomic history: S1 is an illustrative R reconstruction because the
  original drawing script was unavailable.

Only two named statistical functions remain: HC3 model calculations in the
ABO script and the repeatedly evaluated Fitch score in the Maast script.
Plotting-library callbacks remain where required for heatmap cell labels.
All other analysis and plotting steps run at the top level or in ordinary loops.

Statistical result tables were checked against the previously validated release.
The benchmark, pathway, ABO and gene-sharing results are unchanged; Maast
statistics agree within numerical rounding, with unchanged P/q-values.
These checks cover the data-to-figure scripts. Upstream assembly, profiling,
pangenome construction and model estimation remain explicit private inputs.
Fonts, spacing and manual manuscript-page assembly may differ from the PDFs.
