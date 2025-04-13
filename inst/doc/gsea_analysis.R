## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----setup, eval=FALSE--------------------------------------------------------
# # Install ggpicrust2
# if (!requireNamespace("ggpicrust2", quietly = TRUE)) {
#   devtools::install_github("cafferychen777/ggpicrust2")
# }
# 
# # Install required Bioconductor packages
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# 
# BiocManager::install(c("fgsea", "clusterProfiler", "enrichplot", "DOSE", "pathview"))
# 
# # Load the package
# library(ggpicrust2)
# library(dplyr)
# library(ggplot2)

## ----basic-gsea, eval=FALSE---------------------------------------------------
# # Load example data
# data(ko_abundance)
# data(metadata)
# 
# # Prepare abundance data
# abundance_data <- as.data.frame(ko_abundance)
# rownames(abundance_data) <- abundance_data[, "#NAME"]
# abundance_data <- abundance_data[, -1]
# 
# # Run GSEA analysis
# gsea_results <- pathway_gsea(
#   abundance = abundance_data,
#   metadata = metadata,
#   group = "Environment",
#   pathway_type = "KEGG",
#   method = "fgsea",
#   rank_method = "signal2noise",
#   nperm = 1000,
#   min_size = 10,
#   max_size = 500,
#   p.adjust = "BH",
#   seed = 42
# )
# 
# # View the top results
# head(gsea_results)

## ----annotate-gsea, eval=FALSE------------------------------------------------
# # Annotate GSEA results
# annotated_results <- gsea_pathway_annotation(
#   gsea_results = gsea_results,
#   pathway_type = "KEGG"
# )
# 
# # View the annotated results
# head(annotated_results)

## ----barplot, eval=FALSE------------------------------------------------------
# # Create a barplot of the top enriched pathways
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
# # Create a dotplot of the top enriched pathways
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
# # Create an enrichment plot for a specific pathway
# enrichment_plot <- visualize_gsea(
#   gsea_results = annotated_results,
#   plot_type = "enrichment_plot",
#   n_pathways = 10,
#   sort_by = "NES"
# )
# 
# # Display the plot
# enrichment_plot

## ----compare-gsea-daa, eval=FALSE---------------------------------------------
# # Run DAA analysis
# daa_results <- pathway_daa(
#   abundance = abundance_data,
#   metadata = metadata,
#   group = "Environment",
#   daa_method = "ALDEx2"
# )
# 
# # Annotate DAA results
# annotated_daa_results <- pathway_annotation(
#   pathway = "KO",
#   daa_results_df = daa_results,
#   ko_to_kegg = TRUE
# )
# 
# # Compare GSEA and DAA results
# comparison <- compare_gsea_daa(
#   gsea_results = annotated_results,
#   daa_results = annotated_daa_results,
#   plot_type = "venn",
#   p_threshold = 0.05
# )
# 
# # Display the comparison plot
# comparison$plot
# 
# # View the comparison results
# comparison$results

## ----integrated-analysis, eval=FALSE------------------------------------------
# # Run integrated analysis
# integrated_results <- ggpicrust2_extended(
#   data = ko_abundance,
#   metadata = metadata,
#   group = "Environment",
#   pathway = "KO",
#   daa_method = "LinDA",
#   ko_to_kegg = TRUE,
#   run_gsea = TRUE,
#   gsea_params = list(
#     method = "fgsea",
#     rank_method = "signal2noise",
#     nperm = 1000
#   )
# )
# 
# # Access DAA results
# daa_results <- integrated_results$daa_results
# 
# # Access GSEA results
# gsea_results <- integrated_results$gsea_results
# 
# # Access plots
# daa_plot <- integrated_results$daa_plot
# gsea_plot <- integrated_results$gsea_plot

