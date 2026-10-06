## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----installation, eval = FALSE-----------------------------------------------
# install.packages(c("ggpicrust2", "MicrobiomeStat", "BiocManager"))
# BiocManager::install(c("KEGGREST", "limma"))

## ----example-data, eval = FALSE-----------------------------------------------
# library(ggpicrust2)
# library(tibble)
# data("ko_abundance")
# data("metadata")
# alpha <- 0.05

## ----one-command, eval = FALSE------------------------------------------------
# results <- ggpicrust2(
#   data = ko_abundance,
#   metadata = metadata,
#   group = "Environment",
#   pathway = "KO",
#   daa_method = "LinDA",
#   ko_to_kegg = TRUE,
#   order = "pathway_class",
#   p_values_bar = TRUE,
#   p_values_threshold = alpha,
#   x_lab = "pathway_name"
# )
# 
# # A method's plot is NULL when no pathways can be plotted.
# results[[1]]$plot
# head(results[[1]]$results)

## ----ko-to-kegg, eval = FALSE-------------------------------------------------
# kegg_pathway_abundance <- ko2kegg_abundance(data = ko_abundance)
# head(kegg_pathway_abundance[, 1:3])

## ----sample-groups, eval = FALSE----------------------------------------------
# stopifnot(
#   !anyNA(metadata$sample_name),
#   !anyDuplicated(metadata$sample_name),
#   setequal(colnames(kegg_pathway_abundance), metadata$sample_name)
# )
# sample_groups <- setNames(metadata$Environment, metadata$sample_name)

## ----daa, eval = FALSE--------------------------------------------------------
# daa_results <- pathway_daa(
#   abundance = kegg_pathway_abundance,
#   metadata = metadata,
#   group = "Environment",
#   daa_method = "LinDA"
# )
# 
# head(daa_results)

## ----annotation, eval = FALSE-------------------------------------------------
# annotated_daa <- pathway_annotation(
#   pathway = "KO",
#   daa_results_df = daa_results,
#   ko_to_kegg = TRUE,
#   p_adjust_threshold = alpha
# )
# 
# head(annotated_daa)

## ----errorbar, eval = FALSE---------------------------------------------------
# sig_pathways <- unique(annotated_daa$feature[
#   !is.na(annotated_daa$p_adjust) & annotated_daa$p_adjust < alpha
# ])
# 
# p <- NULL
# if (length(sig_pathways) > 0) {
#   p <- pathway_errorbar(
#     abundance = kegg_pathway_abundance,
#     daa_results_df = annotated_daa,
#     Group = sample_groups,
#     ko_to_kegg = TRUE,
#     p_values_threshold = alpha,
#     order = "pathway_class",
#     x_lab = "pathway_name"
#   )
# } else {
#   message("No pathways pass the adjusted p-value threshold; skipping the error bar plot.")
# }
# p

## ----heatmap-pca, eval = FALSE------------------------------------------------
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

## ----aldex2-alternative, eval = FALSE-----------------------------------------
# # Install once with BiocManager::install("ALDEx2").
# set.seed(207)
# aldex_results <- pathway_daa(
#   abundance = kegg_pathway_abundance,
#   metadata = metadata,
#   group = "Environment",
#   daa_method = "ALDEx2"
# )
# daa_results <- aldex_results[
#   aldex_results$method == "ALDEx2_Welch's t test", , drop = FALSE
# ]

## ----contrib-example, eval = FALSE--------------------------------------------
# contrib_input <- expand.grid(
#   sample = paste0("S", 1:4),
#   function_id = c("K00001", "K00002"),
#   taxon = c("ASV1", "ASV2"),
#   stringsAsFactors = FALSE
# )
# contrib_input$taxon_function_abun <- seq_len(nrow(contrib_input))
# contrib_data <- read_contrib_file(data = contrib_input)
# contrib_metadata <- data.frame(
#   sample_name = paste0("S", 1:4),
#   Environment = rep(c("Control", "Treatment"), each = 2)
# )
# taxonomy <- data.frame(
#   ASV = c("ASV1", "ASV2"),
#   Genus = c("ExampleGenusA", "ExampleGenusB")
# )
# taxa_contrib <- aggregate_taxa_contributions(
#   contrib_data, taxonomy = taxonomy, tax_level = "Genus", top_n = 2
# )
# head(taxa_contrib)

## ----contrib-plots, eval = FALSE----------------------------------------------
# taxa_contribution_bar(
#   contrib_agg = taxa_contrib,
#   metadata = contrib_metadata,
#   group = "Environment",
#   facet_by = "function"
# )
# taxa_contribution_heatmap(contrib_agg = taxa_contrib, n_functions = 2)

## ----pathway-contrib-example, eval = FALSE------------------------------------
# path_input <- contrib_input
# path_input$function_id <- ifelse(path_input$function_id == "K00001",
#                                   "GLYCOLYSIS", "PWY-5484")
# path_data <- read_pathway_contrib_file(data = path_input)
# path_taxa_contrib <- aggregate_taxa_contributions(
#   path_data, taxonomy = taxonomy, tax_level = "Genus", top_n = 2
# )
# pathway_annotation_df <- pathway_annotation(
#   data = data.frame(function_id = unique(path_taxa_contrib$function_id)),
#   pathway = "MetaCyc"
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

