#' Plot PCA Data in 3D
#'
#' This function takes a ‘DESeqDataSet’ object and returns a 3D PCA plot.
#'
#' @param object a ‘DESeqDataSet’ object with results metadata after running `DESeq`.
#' @param intgroup interesting groups: a character vector of names in ‘colData(x)’ to use for grouping
#' @param ntop number of top genes to use for principal components, selected by highest row variance
#' @return A plotly object of the PCA data
#' @export
plot3DPCA <- function(object, intgroup = "condition", ntop = 500) {
  pca_data <- CustomRFuncs::returnPCA(object, intgroup, ntop)
  plotly::plot_ly(
    data = pca_data,
    x = ~PC1,
    y = ~PC2,
    z = ~PC3,
    color = ~condition,
    type = "scatter3d",
    text = ~name,
    mode = "markers+text"
  )
}
