% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotPCA.R
\name{plotPCA}
\alias{plotPCA}
\title{Plot PCA Data}
\usage{
plotPCA(object, intgroup = "condition", ntop = 500, x = "PC1", y = "PC2")
}
\arguments{
\item{object}{a ‘DESeqDataSet’ object with results metadata after running \code{DESeq}.}

\item{intgroup}{interesting groups: a character vector of names in ‘colData(x)’ to use for grouping}

\item{ntop}{number of top genes to use for principal components, selected by highest row variance}

\item{x}{x-axis principal component}

\item{y}{y-axis principal component}
}
\value{
A ggplot object of the PCA data
}
\description{
This function takes a ‘DESeqDataSet’ object object, performs a variance stabilizing transform, and returns a PCA plot. This code is modified from the DESeq2 function 'plotPCA'.
}
