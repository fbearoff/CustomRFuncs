#' Return PCA Data
#'
#' This function takes a 'DESeqTransform' object, performs PCA, and returns the PCA data for PC1-4 as a dataframe. This code is modified from the DESeq2 function 'plotPCA'.
#'
#' @param object a ‘DESeqTransform’ object, with data in ‘assay(x)’, produced for example by either ‘rlog’ or ‘varianceStabilizingTransformation’.
#' @param intgroup interesting groups: a character vector of names in ‘colData(x)’ to use for grouping
#' @param ntop number of top genes to use for principal components, selected by highest row variance
#' @return A dataframe of the PCA data
#' @export
returnPCA <- function(object, intgroup = "condition", ntop = 500) {
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(
    ntop,
    length(rv)
  ))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2 / sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup,
    drop = FALSE
  ])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = ":"))
  } else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    PC3 = pca$x[, 3],
    PC4 = pca$x[, 4],
    group = group,
    intgroup.df,
    name = colnames(object)
  )
  attr(d, "percentVar") <- percentVar[1:4]
  return(d)
}
