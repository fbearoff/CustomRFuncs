% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotMA.R
\name{plotMA}
\alias{plotMA}
\title{Plot MA}
\usage{
plotMA(condition1, condition2, top_n_fc = 12)
}
\arguments{
\item{condition1}{Reference condition}

\item{condition2}{Experimental condition}

\item{top_n_fc}{Number of DE genes to be highlighted per direction}
}
\value{
A ggplot of the MA plot
}
\description{
This function returns a MA plot of the supplied conditions from a DESeq2 object.
}
