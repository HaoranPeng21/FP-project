# Private input contract

No input files or example participant rows are supplied. The table below describes
column names and formats only. Paths relative to `FP_DATA_ROOT` are defined once
in `R/inputs.R`; private absolute-path overrides belong in an external
`FP_CONFIG_FILE`. Existing source links can be used privately, but none are
distributed in this repository.

| Input key | Format / required fields |
|---|---|
| `bac120_tree` | GTDB-Tk bac120 Newick; finite branch lengths; labels must match the private annotation tables |
| `catalogue_tips` | CSV; label, ID, ID_original, haoran_species, haoran_species_base, haoran_species_group, haoran_sp, Continent, Host1, Country, Host, Host2, Type, Type_isolate |
| `animal_tips` | CSV; label, ID, ID_original, haoran_species, haoran_species_base, haoran_species_group, haoran_sp, Continent, Host1, Country, Host, Host2, Type, Type_isolate |
| `animal_species_tips` | CSV; label, ID, ID_original, haoran_species, haoran_species_base, haoran_species_group, haoran_sp, Continent, Host1, Country, Host, Host2, Type, Type_isolate |
| `catalogue` | CSV; ID, Host1, Haoran_species, CheckM2.completeness, CheckM2.contamination, Species, Species_v220, Species_v226 |
| `analysis_set` | TSV; ID and catalogue QC fields; frozen 17,103-genome inclusion set |
| `rarefaction` | TSV; genomes_sampled, mean_species, q025_species, q975_species |
| `benchmark_all83` | TSV; run_id, tool, species_id, standardized_name, metaphlan_unique_map, selected_genome_id, is_database_representative, expected_percent, expected_reads, observed_percent, truth_present, detected, error_percent |
| `benchmark_shared24` | TSV; run_id, tool, species_id, standardized_name, metaphlan_unique_map, selected_genome_id, is_database_representative, expected_percent, expected_reads, observed_percent, truth_present, detected, error_percent, sgb, gtdb_species |
| `benchmark_runs` | TSV; run_id, run_type, selection_seed |
| `benchmark_mapping` | TSV; sgb, species_id, standardized_name, project_species, representative_id, mapping_basis, ani_percent, ani_reference_release, current_sgb_database, gtdb_species, gtdb_name_support, n_hq_genomes_supporting_gtdb_name, native_species, native_genus, mapping_note |
| `species_scope` | TSV; species_id, display_name, n_genomes, catalogue_genomes, known_host_genomes, unknown_host_genomes, display_name_updated, short_name, figure_order, is_reference |
| `pathway_prevalence` | TSV; species_id, feature_kind, feature_id, n_present, n_genomes, prevalence, wilson_lo, wilson_hi, display_name, pathway_name, species_name, figure_order |
| `pathway_tests` | TSV; species_id, species, species_name, reference, feature_id, pathway_name, n_species, k_species, n_reference, k_reference, prevalence_species, prevalence_reference, delta, p_value, q_value, abs_delta, direction, eligible (additional result fields permitted) |
| `pathway_selection` | TSV; feature_id, pathway_name, n_equivalent_pathways, equivalent_feature_ids, equivalent_pathway_names, selected_by |
| `host_counts` | TSV; host, s__Faecalibacterium longum, s__Faecalibacterium duncaniae, s__Faecalibacterium taiwanense, s__Faecalibacterium prausnitzii, s__Faecalibacterium prausnitzii_E, s__Faecalibacterium prausnitzii_I, s__Faecalibacterium intestinale, s__Faecalibacterium sp900539945, s__Faecalibacterium hattorii, s__Faecalibacterium butyricigenerans, s__Faecalibacterium tardum, s__Faecalibacterium sp900540455, s__Faecalibacterium wellingii, s__Faecalibacterium prausnitzii_F, s__Faecalibacterium gallistercoris, s__Faecalibacterium sp003449675, s__Faecalibacterium prausnitzii_M (additional result fields permitted) |
| `human_strata` | TSV; panel, panel_order, group_order, group_key, group, n_samples, n_subjects, species, species_rank, n_positive, prevalence_pct, mean_positive_abundance_pct, log2_mean_positive_abundance_vs_reference |
| `human_species_order` | TSV; display_rank, prevalence_rank, species, display_species, n_samples, n_positive, prevalence_pct |
| `infant_summary` | TSV; timepoint, species, n_samples, detected_samples, prevalence, prevalence_ci_low, prevalence_ci_high, mean_abundance_if_detected, mean_abundance_if_detected_ci_low, mean_abundance_if_detected_ci_high; 10 species, 7 timepoints |
| `infant_tests` | TSV; species, endpoint, model, working_correlation, cluster_unit, wald_chi2, df, p_value, q_value, significance; frozen 10-species two-part GEE tests |
| `infant_palette` | TSV; species, display_species, source_color, display_color |
| `infant_models` | TSV; species, trait, endpoint, estimand, adjustment_model, coefficient, display_model_class, heatmap_row, heatmap_column, adjustment_covariates, n_rows, n_subjects, omnibus_p_value, q_value, status, effect_measure, effect_ratio, ci_low (additional result fields permitted) |
| `coabundance_tree` | Newick partner phylogeny; labels map through coabundance_taxonomy |
| `coabundance_taxonomy` | Two-column headerless TSV: leaf, GTDB taxonomy |
| `coabundance_order` | One species identifier per line; frozen 13-species phylogenetic order |
| `coabundance_ratios` | TSV; otu_a, otu_b, n_cohorts_total, n_tested, n_sig, n_sig_pos, n_sig_neg, pct_all_sig_pos, pct_all_sig_neg, pct_tested_sig_pos, pct_tested_sig_neg, ratio_pos_among_sig, edge_class |
| `coabundance_meta` | TSV; otu_a, otu_b, species_a, species_b, edge_key, n_cohorts_used, n_total_cohorts, n_positive_cohorts, n_negative_cohorts, positive_direction_fraction, tau2_reml, meta_z_reml, meta_se_reml, meta_r_reml, meta_r_ci_low, meta_r_ci_high, meta_p_reml, cochran_q (additional result fields permitted) |
| `phenotype_meta` | TSV; association_source, phenotype_raw, phenotype_canonical, phenotype, feature, k_study, beta, se, z, p_meta, p_meta_fdr |
| `phenotype_info` | TSV; SuperClass, Class, SubClass, UnifiedName, Description, Unit, ValueType, Total_N |
| `phenotype_scope` | TSV; feature, species_label, full_nonmissing, full_positive, full_prevalence, full_strict_gt_0.05, full_exactly_eq_0.05, eligible_nonmissing, eligible_positive, eligible_prevalence, eligible_strict_gt_0.05, eligible_exactly_eq_0.05, eligible_minus_full_prevalence, threshold_classification |
| `species_phylogeny_order` | TSV; phylogeny_order, species, node, tree_y, tree_x, node_label |
| `abo_abundance` | TSV; clade_name and one column per private sample; species-level relative abundance |
| `abo_phenotype` | TSV; sample, Bloodtype, Blood_genotype, FUT2_secretor_status |
| `phenotype_master` | CSV; sample_id, study, age, gender, BMI, number_reads (only these six columns are read) |
| `galnac_prevalence` | TSV; species_gtdb_r226, n_genomes_total, n_pathway_positive, pct_pathway_positive, ci95_low_pct, ci95_high_pct, n_core_cluster_before_step4, n_step4_genomewide, figure_role |
| `abo_species` | TSV; phenotype_heatmap_order, species_id, species_label, n_genomes, phylogeny_order, row_label |
| `within_main` | Directory of frozen compressed genome scans, gene tracks, family statistics and long gene-presence tables |
| `within_supplementary` | Directory of frozen PPI-case genome scans, gene tracks, family statistics and long gene-presence tables |
| `within_presence` | Directory of frozen phenotype-class counts and BMI-case PC5 gene-panel tables |
| `geography_country_plot` | TSV; continent, country, country_order, country_total_samples, species, species_display, species_phylogeny_order, species_overall_prevalence, detected_samples, prevalence, wilson95_low, wilson95_high |
| `world_boundaries` | Natural Earth Admin-0 GeoJSON (polygons, ISO country identifiers); external input, not bundled |
| `world_counts` | TSV; iso3, country, count, share_of_all_samples, dataset_total, country_resolved_total |
| `gene_representatives` | TSV; species_id, representative_ID; additional upstream columns may be present |
| `gene_label_map` | TSV; source_label, species_id, skani_95_average, final84_standardized_name, gff_genomes, genomes_with_final84_group, genomes_missing_final84_group, one_to_one_label_mapping, plot_label |
| `gene_ratio` | TSV numeric square Jaccard matrix; first column is the source species key, followed by one column per source species |
| `gene_pair_prevalence` | TSV; species_i, species_j, gene, prevalence_i, prevalence_j |
| `gene_qap` | TSV; n_species, n_pairs, spearman_rho, ordinary_spearman_p, qap_permutations, qap_seed, qap_lower_extreme, qap_two_sided_extreme, qap_one_sided_p, qap_two_sided_p, qap_null_mean, qap_null_sd, qap_null_min, qap_null_max |
| `gene_highlights` | TSV; key_id, highlight_group, species_i_short, species_j_short, phylo_distance, shared_jaccard |
| `maast_trees` | Directory of Newick species trees; filenames are resolved through maast_tips.tree_file |
| `maast_tips` | CSV; label, ID, ID_original, Continent, Host1, Country, Type, Type_isolate, Species, tree_file |
| `host_pathways` | TSV; species_id, feature_kind, feature_id, n_animal, n_human, present_animal, present_human, prevalence_animal, prevalence_human, prevalence_diff, fisher_odds_ratio, fisher_p, adj_score_z_animal, adj_p, fisher_q_global, adj_q_global, fisher_q_within_species, adj_q_within_species (additional result fields permitted) |

Identifiers are join keys inside private inputs. Use consistent identifiers across
matrices, metadata and trees; independently anonymizing each table breaks these
joins. The public release does not export identifiers or anonymized records.

Frozen-result inputs retain the original model fits and multiple-testing universe.
They are required for the corresponding data-to-figure entry points; this code
does not infer missing coefficients from a PDF or refit a different model.

For within-species panels the three input directories must contain the filenames
explicitly read in `analyses/12_within_species.R`. They are the PC5 release, not
earlier Fisher-only or exploratory module releases. For gene sharing, the label map
and Jaccard matrix refer to the same 25-species Roary i100 subset.

Tree tips containing quotes or FASTA suffixes are normalized before joining.
Duplicate normalized tips, missing required joins and incompatible fixed scopes
are errors. Public display-name substitutions occur after the scientific joins.
