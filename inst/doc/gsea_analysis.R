## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----installation, eval=FALSE-------------------------------------------------
# install.packages(c("ggpicrust2", "MicrobiomeStat", "ggridges", "ggVennDiagram",
#                    "circlize", "igraph", "BiocManager"))
# BiocManager::install(c("limma", "fgsea", "ComplexHeatmap"))

## ----setup, eval=FALSE--------------------------------------------------------
# library(ggpicrust2)

## ----basic-gsea, eval=FALSE---------------------------------------------------
# # Load example data
# data(ko_abundance)
# data(metadata)
# metadata$Environment <- factor(
#   metadata$Environment, levels = c("Pro-inflammatory", "Pro-survival")
# )
# 
# # Prepare abundance data
# abundance_data <- as.data.frame(ko_abundance)
# rownames(abundance_data) <- abundance_data[, "#NAME"]
# abundance_data <- abundance_data[, -1]
# 
# # Run the competitive camera test
# gsea_results <- pathway_gsea(
#   abundance = abundance_data,
#   metadata = metadata,
#   group = "Environment",
#   pathway_type = "KEGG",
#   method = "camera",
#   inter.gene.cor = NA_real_,
#   min_size = 5,
#   max_size = 500,
#   p_adjust_method = "BH"
# )
# 
# # View the top results
# head(gsea_results)

## ----inspect-results, eval=FALSE----------------------------------------------
# table(gsea_results$direction)
# sum(gsea_results$p.adjust < 0.05, na.rm = TRUE)
# unique(gsea_results[, c("method", "score_type", "score_label")])

## ----covariate-gsea, eval=FALSE-----------------------------------------------
# # Mouse_Sex is observed in the bundled metadata and varies within both groups.
# 
# gsea_results_adjusted <- pathway_gsea(
#   abundance = abundance_data,
#   metadata = metadata,
#   group = "Environment",
#   covariates = "Mouse_Sex",
#   pathway_type = "KEGG",
#   method = "camera",
#   inter.gene.cor = NA_real_
# )
# 
# # The results now reflect the group effect after adjusting for confounders
# head(gsea_results_adjusted)

## ----fry-gsea, eval=FALSE-----------------------------------------------------
# # Fast rotation gene set test
# gsea_results_fry <- pathway_gsea(
#   abundance = abundance_data,
#   metadata = metadata,
#   group = "Environment",
#   pathway_type = "KEGG",
#   method = "fry",
#   min_size = 5,
#   max_size = 500
# )
# 
# head(gsea_results_fry)

## ----fgsea, eval=FALSE--------------------------------------------------------
# # Preranked testing uses a different null from camera/fry.
# gsea_results_fgsea <- pathway_gsea(
#   abundance = abundance_data,
#   metadata = metadata,
#   group = "Environment",
#   pathway_type = "KEGG",
#   method = "fgsea",
#   rank_method = "signal2noise",
#   comparison = c("Pro-survival", "Pro-inflammatory"),
#   min_size = 10,
#   max_size = 500,
#   p_adjust_method = "BH",
#   seed = 42
# )
# 
# # View the top results
# head(gsea_results_fgsea)

## ----annotate-gsea, eval=FALSE------------------------------------------------
# # Annotate GSEA results
# annotated_results <- gsea_pathway_annotation(
#   gsea_results = gsea_results,
#   pathway_type = "KEGG"
# )
# 
# # View the annotated results
# head(annotated_results)

## ----pathway-labels, eval=FALSE-----------------------------------------------
# # Option 1: Use raw GSEA results (shows pathway IDs)
# plot_with_ids <- visualize_gsea(
#   gsea_results = gsea_results,
#   plot_type = "barplot",
#   n_pathways = 10
# )
# 
# # Option 2: Use annotated results (automatically shows pathway names)
# plot_with_names <- visualize_gsea(
#   gsea_results = annotated_results,
#   plot_type = "barplot",
#   n_pathways = 10
# )
# 
# # Option 3: Explicitly specify which column to use for labels
# plot_custom_labels <- visualize_gsea(
#   gsea_results = annotated_results,
#   plot_type = "barplot",
#   pathway_label_column = "pathway_name",
#   n_pathways = 10
# )
# 
# # Compare the plots
# plot_with_ids
# plot_with_names
# plot_custom_labels

## ----barplot, eval=FALSE------------------------------------------------------
# # Create a barplot of the top-ranked pathways
# barplot <- visualize_gsea(
#   gsea_results = annotated_results,
#   plot_type = "barplot",
#   n_pathways = 20,
#   sort_by = "p.adjust"
# )
# 
# # Display the plot
# barplot

## ----dotplot, eval=FALSE------------------------------------------------------
# # Create a dotplot of the top-ranked pathways
# dotplot <- visualize_gsea(
#   gsea_results = annotated_results,
#   plot_type = "dotplot",
#   n_pathways = 20,
#   sort_by = "p.adjust"
# )
# 
# # Display the plot
# dotplot

## ----enrichment-plot, eval=FALSE----------------------------------------------
# # This is a score-summary bar chart, not a running enrichment curve.
# enrichment_plot <- visualize_gsea(
#   gsea_results = annotated_results,
#   plot_type = "enrichment_plot",
#   n_pathways = 10,
#   sort_by = "p.adjust"
# )
# 
# # Display the plot
# enrichment_plot

## ----ridge-plot, eval=FALSE---------------------------------------------------
# # Create a ridge plot for GSEA results
# # Note: Requires ggridges package to be installed
# ridge_plot <- pathway_ridgeplot(
#   gsea_results = gsea_results,
#   abundance = abundance_data,
#   metadata = metadata,
#   group = "Environment",
#   pathway_type = "KEGG",
#   comparison = c("Pro-inflammatory", "Pro-survival"),
#   n_pathways = 10,
#   sort_by = "p.adjust",
#   show_direction = TRUE,
#   colors = c("Down" = "#3182bd", "Up" = "#de2d26")
# )
# 
# # Display the plot
# ridge_plot

## ----leading-edge-plots, eval=FALSE-------------------------------------------
# annotated_fgsea <- gsea_pathway_annotation(gsea_results_fgsea, pathway_type = "KEGG")
# leading_results <- annotated_fgsea[
#   !is.na(annotated_fgsea$leading_edge) & nzchar(annotated_fgsea$leading_edge), , drop = FALSE
# ]
# if (nrow(leading_results) > 0) {
#   network_plot <- visualize_gsea(
#     leading_results, plot_type = "network", n_pathways = 10,
#     network_params = list(similarity_measure = "jaccard", similarity_cutoff = 0.2)
#   )
#   print(network_plot)
#   leading_heatmap <- visualize_gsea(
#     leading_results, plot_type = "heatmap", n_pathways = 10,
#     abundance = abundance_data, metadata = metadata, group = "Environment",
#     heatmap_params = list(cluster_rows = TRUE, cluster_columns = TRUE,
#                           show_rownames = TRUE)
#   )
#   ComplexHeatmap::draw(leading_heatmap)
# }

## ----compare-gsea-daa, eval=FALSE---------------------------------------------
# # Compare KEGG pathways to KEGG pathways, not individual KO identifiers.
# kegg_pathway_abundance <- ko2kegg_abundance(data = ko_abundance)
# daa_results <- pathway_daa(
#   abundance = kegg_pathway_abundance,
#   metadata = metadata,
#   group = "Environment",
#   daa_method = "LinDA"
# )
# 
# # Compare only pathways that both procedures actually tested.
# # Keep each analysis's original multiple-testing adjustment.
# common_pathways <- intersect(annotated_results$pathway_id, daa_results$feature)
# gsea_common <- annotated_results[annotated_results$pathway_id %in% common_pathways, , drop = FALSE]
# daa_common <- daa_results[daa_results$feature %in% common_pathways, , drop = FALSE]
# comparison <- compare_gsea_daa(
#   gsea_results = gsea_common,
#   daa_results = daa_common,
#   plot_type = "venn",
#   p_threshold = 0.05
# )
# 
# # Display the comparison plot
# comparison$plot
# 
# # View the comparison results
# comparison$results

