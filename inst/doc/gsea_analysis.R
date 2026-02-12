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
# BiocManager::install(c("limma", "fgsea", "clusterProfiler", "enrichplot", "DOSE", "pathview"))
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
# # Run GSEA analysis with camera (recommended)
# gsea_results <- pathway_gsea(
#   abundance = abundance_data,
#   metadata = metadata,
#   group = "Environment",
#   pathway_type = "KEGG",
#   method = "camera",
#   min_size = 5,
#   max_size = 500,
#   p.adjust = "BH"
# )
# 
# # View the top results
# head(gsea_results)

## ----covariate-gsea, eval=FALSE-----------------------------------------------
# # Example with covariate adjustment
# # Assuming metadata has columns: group, age, sex, BMI
# 
# gsea_results_adjusted <- pathway_gsea(
#   abundance = abundance_data,
#   metadata = metadata,
#   group = "Disease",
#   covariates = c("age", "sex", "BMI"),  # Adjust for these confounders
#   pathway_type = "KEGG",
#   method = "camera"
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
# # Note: Preranked methods don't account for inter-gene correlations
# # P-values may be less reliable (Wu et al., 2012)
# gsea_results_fgsea <- pathway_gsea(
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

## ----ridge-plot, eval=FALSE---------------------------------------------------
# # Create a ridge plot for GSEA results
# # Note: Requires ggridges package to be installed
# ridge_plot <- pathway_ridgeplot(
#   gsea_results = gsea_results,
#   abundance = abundance_data,
#   metadata = metadata,
#   group = "Environment",
#   pathway_type = "KEGG",
#   n_pathways = 10,
#   sort_by = "p.adjust",
#   show_direction = TRUE,
#   colors = c("Down" = "#3182bd", "Up" = "#de2d26")
# )
# 
# # Display the plot
# ridge_plot

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

