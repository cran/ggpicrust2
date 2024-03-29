% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ko2kegg_abundance.R
\name{ko2kegg_abundance}
\alias{ko2kegg_abundance}
\title{Convert KO abundance in picrust2 export files to KEGG pathway abundance}
\usage{
ko2kegg_abundance(file = NULL, data = NULL)
}
\arguments{
\item{file}{A character string representing the file path of the input file containing KO abundance data in picrust2 export format. The input file should have KO identifiers in the first column and sample identifiers in the first row. The remaining cells should contain the abundance values for each KO-sample pair.}

\item{data}{An optional data.frame containing KO abundance data in the same format as the input file. If provided, the function will use this data instead of reading from the file. By default, this parameter is set to NULL.}
}
\value{
A data frame with KEGG pathway abundance values. Rows represent KEGG pathways, identified by their KEGG pathway IDs. Columns represent samples, identified by their sample IDs from the input file. Each cell contains the abundance of a specific KEGG pathway in a given sample, calculated by summing the abundances of the corresponding KOs in the input file.
}
\description{
This function takes a file containing KO (KEGG Orthology) abundance data in picrust2 export format and converts it to KEGG pathway abundance data.
The input file should be in .tsv, .txt, or .csv format.
}
\examples{
\dontrun{
library(ggpicrust2)
library(readr)

# Example 1: Demonstration with a hypothetical input file

# Prepare an input file path
input_file <- "path/to/your/picrust2/results/pred_metagenome_unstrat.tsv"

# Run ko2kegg_abundance function
kegg_abundance <- ko2kegg_abundance(file = input_file)

# Alternatively, read the data from a file and use the data argument
file_path <- "path/to/your/picrust2/results/pred_metagenome_unstrat.tsv"
ko_abundance <- read_delim(file_path, delim = "\t")
kegg_abundance <- ko2kegg_abundance(data = ko_abundance)

# Example 2: Working with real data
# In this case, we're using an existing dataset from the ggpicrust2 package.

# Load the data
data(ko_abundance)

# Apply the ko2kegg_abundance function to our real dataset
kegg_abundance <- ko2kegg_abundance(data = ko_abundance)
}
}
