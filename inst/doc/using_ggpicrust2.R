## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----installation, eval = FALSE-----------------------------------------------
# # install.packages("devtools")
# devtools::install_github("cafferychen777/ggpicrust2")
# 
# library(ggpicrust2)
# library(tibble)

## ----one-command, eval = FALSE------------------------------------------------
# data("ko_abundance")
# data("metadata")
# 
# results <- ggpicrust2(
#   data = ko_abundance,
#   metadata = metadata,
#   group = "Environment",
#   pathway = "KO",
#   daa_method = "LinDA",
#   ko_to_kegg = TRUE,
#   order = "pathway_class",
#   p_values_bar = TRUE,
#   x_lab = "pathway_name"
# )
# 
# # Access the main outputs
# results[[1]]$plot
# head(results[[1]]$results)

## ----ko-to-kegg, eval = FALSE-------------------------------------------------
# kegg_pathway_abundance <- ko2kegg_abundance(data = ko_abundance)
# head(kegg_pathway_abundance[, 1:3])

## ----daa, eval = FALSE--------------------------------------------------------
# daa_results <- pathway_daa(
#   abundance = kegg_pathway_abundance,
#   metadata = metadata,
#   group = "Environment",
#   daa_method = "ALDEx2"
# )
# 
# head(daa_results)

## ----annotation, eval = FALSE-------------------------------------------------
# annotated_daa <- pathway_annotation(
#   pathway = "KO",
#   daa_results_df = daa_results,
#   ko_to_kegg = TRUE
# )
# 
# head(annotated_daa)

## ----errorbar, eval = FALSE---------------------------------------------------
# pathway_errorbar(
#   abundance = kegg_pathway_abundance,
#   daa_results_df = annotated_daa,
#   Group = "Environment"
# )

## ----heatmap-pca, eval = FALSE------------------------------------------------
# sig_pathways <- annotated_daa$feature[annotated_daa$p_adjust < 0.05]
# 
# if (length(sig_pathways) > 0) {
#   pathway_heatmap(
#     abundance = kegg_pathway_abundance[sig_pathways, , drop = FALSE],
#     metadata = metadata,
#     group = "Environment"
#   )
# }
# 
# pathway_pca(
#   abundance = kegg_pathway_abundance,
#   metadata = metadata,
#   group = "Environment"
# )

## ----contrib-readers, eval = FALSE--------------------------------------------
# # For pred_metagenome_contrib.tsv
# contrib_data <- read_contrib_file("pred_metagenome_contrib.tsv")
# 
# # For pred_metagenome_strat.tsv
# strat_data <- read_strat_file("pred_metagenome_strat.tsv")

## ----contrib-aggregate, eval = FALSE------------------------------------------
# taxa_contrib <- aggregate_taxa_contributions(
#   contrib_data = contrib_data,
#   taxonomy = your_taxonomy_table,
#   tax_level = "Genus",
#   top_n = 10,
#   daa_results_df = daa_results
# )
# 
# head(taxa_contrib)

## ----contrib-plots, eval = FALSE----------------------------------------------
# taxa_contribution_bar(
#   contrib_agg = taxa_contrib,
#   metadata = metadata,
#   group = "Environment",
#   facet_by = "function"
# )
# 
# taxa_contribution_heatmap(
#   contrib_agg = taxa_contrib,
#   n_functions = 20
# )

## ----gsea, eval = FALSE-------------------------------------------------------
# gsea_results <- pathway_gsea(
#   abundance = ko_abundance %>% column_to_rownames("#NAME"),
#   metadata = metadata,
#   group = "Environment",
#   pathway_type = "KEGG",
#   method = "camera"
# )
# 
# annotated_gsea <- gsea_pathway_annotation(
#   gsea_results = gsea_results,
#   pathway_type = "KEGG"
# )
# 
# visualize_gsea(
#   gsea_results = annotated_gsea,
#   plot_type = "barplot",
#   n_pathways = 15
# )

