# Retained analysis sources

The supplied PDFs determine the figure numbering. The paths below identify the
scientific source implementations inside the private project, not runtime code
dependencies. This release includes its own R implementation and shared R helpers.
Historical render variants, audit experiments and manuscript-writing scripts are
not included.

| Analysis | Scientific source in the private project | Refactoring boundary |
|---|---|---|
| 01 | `04_animal_phylogeny/code/annotate_tree_with_skani_rerun.R` | Frozen full-tree annotations; native ape fan and rings |
| 02 | `12_Manuscript/SFigure_v0.3/code/build_ani_priority24.py`; `code/render_ani_priority24.R` | Native R scoring of frozen shared24/all83 predictions; PDF counts resolve stale same21 titles |
| 03–04 | `13_Manuscript/Figure2/build_Figure2_composite.R` | Retained 19-species heatmaps; full Fisher/BH recomputation |
| 05 | `13_Manuscript/Figure2/tree_helpers_from_reference.R` | Same selected tips and four focal species; native R fan geometry |
| 06 | `08_geography/scripts/summarize_host_strata_bubbles.py`; `13_Manuscript/Figure3/build_Figure3.py` | R rendering from the frozen stratum summaries |
| 07 | `12_Manuscript/Figure_V0.3/build_figure_v03.py`; `06a_early_life_mother_infant/code/plot_main2_composite.R` | Supplied ten-species detection/positive-abundance curves, bootstrap CIs and two-part GEE results |
| 08 | `13_Manuscript/Figure3/data/main2_selected_models.tsv`; `06a_early_life_mother_infant/code/plot_main2_composite.R` | Selected model estimates and their original full-family q-values |
| 09 | `12_Manuscript/Figure_V0.3/panel_d_render.R` | Retained fan/ring geometry with explicit meta markers and legends |
| 10 | `09_phenotype_meta_analysis/scripts/plot_fullphen_heatmap_from_meta.R` | Frozen meta-results, global q-values, prevalence filter and tree order |
| 11 | `10_blood_group/scripts/plot_four_group_and_galnac.R` | Native R matrix extraction, CLR, all 84 Wilcoxon tests and HC3 contrasts |
| 12 | `12_Manuscript/Figure_V0.5/code/plot_main_figure_presence_portrait.R` | PC5 final rendering; no model refit, legacy caption or font-installation code |
| 13 | `12_Manuscript/Figure_V0.0/make_figure_genome_overview.R`; SGB mapping table | Native R QC counts, frozen rarefaction, mapping display |
| 14 | `12_Manuscript/SFigure_v0.3.1/build_sfigure_v031_country_detection.R`; `Figure_world_sample_distribution/ggplot_world_map.R` | Native R country map and prevalence bubbles |
| 15 | `05_gene_sharing/code/plot_final84_phylogeny_gene_sharing.R`; `species_pair_qap_test.py` | Native R tree distances/QAP; preserve original figure QAP label |
| 16 | `05_gene_sharing/code/plot_highlight_pair_prevalence_heatmaps.R` | Native R 100×100 prevalence bins and marginal histograms |
| 17 | `04_animal_phylogeny/code/analyze_human_continent_phylogenetic_association.R` | Same R tests and seeds; precomputed Fitch traversal; native R fan display |
| 18 | `12_Manuscript/SFigure_v0.2/Figure_S3/code/plot_Figure_S3.R` | Same nominal screen and inherited complete-family q-values |
| 19 | Supplied `Figure_S1_New.pdf` and catalogue release labels | Illustrative reconstruction; original drawing code not located |

Rendering has been simplified where possible. Scientific input/result tables are
not embedded in scripts. The two larger renderers, co-abundance and within-species
gene architecture, retain their necessary tree/coordinate logic rather than
depending on external legacy scripts or pre-rendered panels.
