% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bcClassify.R, R/bcClassify1.R
\name{bcClassify}
\alias{bcClassify}
\title{Bladder Cancer Immune Subtype Classification}
\usage{
bcClassify(expr_data)

bcClassify(expr_data)
}
\arguments{
\item{expr_data}{Expression matrix with genes as rows and samples as columns}
}
\value{
A character vector containing classification results for each sample

A character vector containing classification results for each sample
}
\description{
Classifies bladder cancer samples into immune subtypes using gene expression signatures

Classifies bladder cancer samples into immune subtypes using gene expression signatures
}
\examples{
\dontrun{
# Load expression data
expr_data <- read.csv("expression_data.csv", row.names = 1)
# Perform classification
results <- bcClassify(expr_data)
print(table(results))
}
\dontrun{
# Load expression data
expr_data <- read.csv("expression_data.csv", row.names = 1)
# Perform classification
results <- bcClassify(expr_data)
print(table(results))
}
}
