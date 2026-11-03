# Private inputs

Use a separate input folder for each analysis and set `data_dir` in its script.
The files below remain private and are not distributed with this repository.
Filenames are neutral aliases of the existing inputs; no source data are changed.

The final catalogue inclusion list and tree annotations must contain **17,103
unique genomes and 83 species**. QC and taxonomy scripts select this inclusion
list explicitly. Other figure inputs cover the subsets stated in their scripts.

| Analysis | Files in `data_dir` |
|---|---|
| 01_catalogue_tree | `bac120_tree.nwk`, `catalogue_tips.csv` |
| 02_benchmark | `benchmark_all83.tsv`, `benchmark_shared24.tsv`, `benchmark_runs.tsv` |
| 03_pathway_prevalence | `species_scope.tsv`, `pathway_prevalence.tsv`, `pathway_tests.tsv`, `pathway_selection.tsv` |
| 04_host_distribution | `species_scope.tsv`, `host_counts.tsv` |
| 05_animal_trees | `bac120_tree.nwk`, `animal_tips.csv`, `animal_species_tips.csv` |
| 06_human_distribution | `human_strata.tsv`, `human_species_order.tsv` |
| 07_infant_trajectories | `infant_summary.tsv`, `infant_tests.tsv`, `infant_palette.tsv` |
| 08_infant_associations | `infant_models.tsv` |
| 09_coabundance | `coabundance_tree.nwk`, `coabundance_ratios.tsv`, `coabundance_taxonomy.txt`, `coabundance_order.txt`, `coabundance_meta.tsv` |
| 10_phenotype_meta | `phenotype_meta.tsv`, `phenotype_scope.tsv`, `phenotype_info.tsv`, `species_phylogeny_order.tsv` |
| 11_abo | `abo_abundance.tsv`, `abo_phenotype.txt`, `phenotype_master.csv`, `phenotype_scope.tsv`, `abo_species.tsv`, `galnac_prevalence.tsv` |
| 12_within_species | `within_main`, `within_supplementary`, `within_presence` |
| 13_catalogue_qc | `catalogue.csv`, `rarefaction.tsv`, `benchmark_mapping.tsv`, `analysis_set.tsv` |
| 14_geography | `world_counts.tsv`, `world_boundaries.geojson`, `geography_country_plot.tsv` |
| 15_gene_sharing | `gene_label_map.tsv`, `gene_representatives.tsv`, `bac120_tree.nwk`, `gene_ratio.tsv`, `gene_qap.tsv` |
| 16_gene_pair_prevalence | `gene_highlights.tsv`, `gene_label_map.tsv`, `gene_pair_prevalence.tsv` |
| 17_maast_tree_permanova | `maast_tips.csv`, `maast_trees` |
| 18_host_pathways | `host_pathways.tsv` |
| 19_taxonomy_history | `catalogue.csv`, `analysis_set.tsv` |

Tabular inputs use a header row. TSV files are tab-separated; CSV files are
comma-separated; `.nwk` files are Newick trees. Input directories for the
within-species and Maast analyses retain the filenames read explicitly in the
corresponding R script. Genome labels must match the tree tips after the
suffix normalization shown in the script.

| File | Required columns or structure |
|---|---|
| `bac120_tree.nwk` | GTDB-Tk bac120 Newick; finite branch lengths; labels must match the private annotation tables |
| `catalogue_tips.csv` | CSV; label, ID, ID_original, haoran_species, haoran_species_base, haoran_species_group, haoran_sp, Continent, Host1, Country, Host, Host2, Type, Type_isolate |
| `benchmark_all83.tsv` | TSV; run_id, tool, species_id, standardized_name, metaphlan_unique_map, selected_genome_id, is_database_representative, expected_percent, expected_reads, observed_percent, truth_present, detected, error_percent |
| `benchmark_shared24.tsv` | TSV; run_id, tool, species_id, standardized_name, metaphlan_unique_map, selected_genome_id, is_database_representative, expected_percent, expected_reads, observed_percent, truth_present, detected, error_percent, sgb, gtdb_species |
| `benchmark_runs.tsv` | TSV; run_id, run_type, selection_seed |
| `species_scope.tsv` | TSV; species_id, display_name, n_genomes, catalogue_genomes, known_host_genomes, unknown_host_genomes, display_name_updated, short_name, figure_order, is_reference |
| `pathway_prevalence.tsv` | TSV; species_id, feature_kind, feature_id, n_present, n_genomes, prevalence, wilson_lo, wilson_hi, display_name, pathway_name, species_name, figure_order |
| `pathway_tests.tsv` | TSV; species_id, species, species_name, reference, feature_id, pathway_name, n_species, k_species, n_reference, k_reference, prevalence_species, prevalence_reference, delta, p_value, q_value, abs_delta, direction, eligible (additional result fields permitted) |
| `pathway_selection.tsv` | TSV; feature_id, pathway_name, n_equivalent_pathways, equivalent_feature_ids, equivalent_pathway_names, selected_by |
| `host_counts.tsv` | TSV; host, s__Faecalibacterium longum, s__Faecalibacterium duncaniae, s__Faecalibacterium taiwanense, s__Faecalibacterium prausnitzii, s__Faecalibacterium prausnitzii_E, s__Faecalibacterium prausnitzii_I, s__Faecalibacterium intestinale, s__Faecalibacterium sp900539945, s__Faecalibacterium hattorii, s__Faecalibacterium butyricigenerans, s__Faecalibacterium tardum, s__Faecalibacterium sp900540455, s__Faecalibacterium wellingii, s__Faecalibacterium prausnitzii_F, s__Faecalibacterium gallistercoris, s__Faecalibacterium sp003449675, s__Faecalibacterium prausnitzii_M (additional result fields permitted) |
| `animal_tips.csv` | CSV; label, ID, ID_original, haoran_species, haoran_species_base, haoran_species_group, haoran_sp, Continent, Host1, Country, Host, Host2, Type, Type_isolate |
| `animal_species_tips.csv` | CSV; label, ID, ID_original, haoran_species, haoran_species_base, haoran_species_group, haoran_sp, Continent, Host1, Country, Host, Host2, Type, Type_isolate |
| `human_strata.tsv` | TSV; panel, panel_order, group_order, group_key, group, n_samples, n_subjects, species, species_rank, n_positive, prevalence_pct, mean_positive_abundance_pct, log2_mean_positive_abundance_vs_reference |
| `human_species_order.tsv` | TSV; display_rank, prevalence_rank, species, display_species, n_samples, n_positive, prevalence_pct |
| `infant_summary.tsv` | TSV; timepoint, species, n_samples, detected_samples, prevalence, prevalence_ci_low, prevalence_ci_high, mean_abundance_if_detected, mean_abundance_if_detected_ci_low, mean_abundance_if_detected_ci_high; 10 species, 7 timepoints |
| `infant_tests.tsv` | TSV; species, endpoint, model, working_correlation, cluster_unit, wald_chi2, df, p_value, q_value, significance; frozen 10-species two-part GEE tests |
| `infant_palette.tsv` | TSV; species, display_species, source_color, display_color |
| `infant_models.tsv` | TSV; species, trait, endpoint, estimand, adjustment_model, coefficient, display_model_class, heatmap_row, heatmap_column, adjustment_covariates, n_rows, n_subjects, omnibus_p_value, q_value, status, effect_measure, effect_ratio, ci_low (additional result fields permitted) |
| `coabundance_tree.nwk` | Newick partner phylogeny; labels map through coabundance_taxonomy |
| `coabundance_ratios.tsv` | TSV; otu_a, otu_b, n_cohorts_total, n_tested, n_sig, n_sig_pos, n_sig_neg, pct_all_sig_pos, pct_all_sig_neg, pct_tested_sig_pos, pct_tested_sig_neg, ratio_pos_among_sig, edge_class |
| `coabundance_taxonomy.txt` | Two-column headerless TSV: leaf, GTDB taxonomy |
| `coabundance_order.txt` | One species identifier per line; frozen 13-species phylogenetic order |
| `coabundance_meta.tsv` | TSV; otu_a, otu_b, species_a, species_b, edge_key, n_cohorts_used, n_total_cohorts, n_positive_cohorts, n_negative_cohorts, positive_direction_fraction, tau2_reml, meta_z_reml, meta_se_reml, meta_r_reml, meta_r_ci_low, meta_r_ci_high, meta_p_reml, cochran_q (additional result fields permitted) |
| `phenotype_meta.tsv` | TSV; association_source, phenotype_raw, phenotype_canonical, phenotype, feature, k_study, beta, se, z, p_meta, p_meta_fdr |
| `phenotype_scope.tsv` | TSV; feature, species_label, full_nonmissing, full_positive, full_prevalence, full_strict_gt_0.05, full_exactly_eq_0.05, eligible_nonmissing, eligible_positive, eligible_prevalence, eligible_strict_gt_0.05, eligible_exactly_eq_0.05, eligible_minus_full_prevalence, threshold_classification |
| `phenotype_info.tsv` | TSV; SuperClass, Class, SubClass, UnifiedName, Description, Unit, ValueType, Total_N |
| `species_phylogeny_order.tsv` | TSV; phylogeny_order, species, node, tree_y, tree_x, node_label |
| `abo_abundance.tsv` | TSV; clade_name and one column per private sample; species-level relative abundance |
| `abo_phenotype.txt` | TSV; sample, Bloodtype, Blood_genotype, FUT2_secretor_status |
| `phenotype_master.csv` | CSV; sample_id, study, age, gender, BMI, number_reads (only these six columns are read) |
| `abo_species.tsv` | TSV; phenotype_heatmap_order, species_id, species_label, n_genomes, phylogeny_order, row_label |
| `galnac_prevalence.tsv` | TSV; species_gtdb_r226, n_genomes_total, n_pathway_positive, pct_pathway_positive, ci95_low_pct, ci95_high_pct, n_core_cluster_before_step4, n_step4_genomewide, figure_role |
| `within_main` | Directory of frozen compressed genome scans, gene tracks, family statistics and long gene-presence tables |
| `within_supplementary` | Directory of frozen PPI-case genome scans, gene tracks, family statistics and long gene-presence tables |
| `within_presence` | Directory of frozen phenotype-class counts and BMI-case PC5 gene-panel tables |
| `catalogue.csv` | CSV; ID, Host1, Haoran_species, CheckM2.completeness, CheckM2.contamination, Species, Species_v220, Species_v226 |
| `rarefaction.tsv` | TSV; genomes_sampled, mean_species, q025_species, q975_species |
| `benchmark_mapping.tsv` | TSV; sgb, species_id, standardized_name, project_species, representative_id, mapping_basis, ani_percent, ani_reference_release, current_sgb_database, gtdb_species, gtdb_name_support, n_hq_genomes_supporting_gtdb_name, native_species, native_genus, mapping_note |
| `world_counts.tsv` | TSV; iso3, country, count, share_of_all_samples, dataset_total, country_resolved_total |
| `world_boundaries.geojson` | Natural Earth Admin-0 GeoJSON (polygons, ISO country identifiers); external input, not bundled |
| `geography_country_plot.tsv` | TSV; continent, country, country_order, country_total_samples, species, species_display, species_phylogeny_order, species_overall_prevalence, detected_samples, prevalence, wilson95_low, wilson95_high |
| `gene_label_map.tsv` | TSV; source_label, species_id, skani_95_average, species_name, gff_genomes, mapped_genomes, unmapped_genomes, one_to_one_label_mapping, plot_label |
| `gene_representatives.tsv` | TSV; species_id, representative_ID; additional upstream columns may be present |
| `gene_ratio.tsv` | TSV numeric square Jaccard matrix; first column is the source species key, followed by one column per source species |
| `gene_qap.tsv` | TSV; n_species, n_pairs, spearman_rho, ordinary_spearman_p, qap_permutations, qap_seed, qap_lower_extreme, qap_two_sided_extreme, qap_one_sided_p, qap_two_sided_p, qap_null_mean, qap_null_sd, qap_null_min, qap_null_max |
| `gene_highlights.tsv` | TSV; key_id, highlight_group, species_i_short, species_j_short, phylo_distance, shared_jaccard |
| `gene_pair_prevalence.tsv` | TSV; species_i, species_j, gene, prevalence_i, prevalence_j |
| `maast_tips.csv` | CSV; label, ID, ID_original, Continent, Host1, Country, Type, Type_isolate, Species, tree_file |
| `maast_trees` | Directory of Newick species trees; filenames are resolved through maast_tips.tree_file |
| `host_pathways.tsv` | TSV; species_id, feature_kind, feature_id, n_animal, n_human, present_animal, present_human, prevalence_animal, prevalence_human, prevalence_diff, fisher_odds_ratio, fisher_p, adj_score_z_animal, adj_p, fisher_q_global, adj_q_global, fisher_q_within_species, adj_q_within_species (additional result fields permitted) |
| `analysis_set.tsv` | TSV; ID and catalogue QC fields; frozen 17,103-genome inclusion set |

For the gene-sharing label map, `species_name` is the display-name field.
An existing column ending in `standardized_name` is also accepted and renamed
inside the script. The representative table is matched only to the
25-species pangenome subset; it does not define the full catalogue scope.

Keep private input files and output directories outside this code repository.
