% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pathway_annotation.R
\name{safe_extract}
\alias{safe_extract}
\title{Safely Extract Elements from a List}
\usage{
safe_extract(list, field, index = 1)
}
\arguments{
\item{list}{A list object from which to extract elements}

\item{field}{The name of the field to extract from the list}

\item{index}{The index position to extract from the field. Default is 1}
}
\value{
The extracted element if successful, NA if extraction fails
}
\description{
Safely extracts elements from a list, returning NA if the extraction fails
}
\examples{
# Create a sample list
my_list <- list(
  a = list(x = 1:3),
  b = list(y = 4:6)
)

# Extract existing element
safe_extract(my_list, "a", 1)

# Extract non-existing element (returns NA)
safe_extract(my_list, "c", 1)
}
