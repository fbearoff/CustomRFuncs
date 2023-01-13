#' Plot PCA Data
#'
#' This function takes a ‘DESeqDataSet’ object object, performs a variance stabilizing transform, and returns a PCA plot. This code is modified from the DESeq2 function 'plotPCA'.
#'
#' @param object a ‘DESeqDataSet’ object with results metadata after running `DESeq`.
#' @param intgroup interesting groups: a character vector of names in ‘colData(x)’ to use for grouping
#' @param ntop number of top genes to use for principal components, selected by highest row variance
#' @param x x-axis principal component
#' @param y y-axis principal component
#' @return A ggplot object of the PCA data
#' @export
plotPCA <- function(object = dds, intgroup = "condition", ntop = 500, x = "PC1", y = "PC2") {
  pca_data <- CustomRFuncs::returnPCA(object, intgroup, ntop)

  percent_var <- round(100 * attr(pca_data, "percentVar"))
  pca <- ggplot2::ggplot(pca_data, ggplot2::aes(!!rlang::sym(x), !!rlang::sym(y),
    color = condition,
    fill = condition,
    label = condition
  )) +
    ggplot2::geom_point() +
    ggrepel::geom_text_repel(ggplot2::aes(label = name)) +
    ggforce::geom_mark_ellipse(
      show.legend = FALSE,
      con.cap = 0,
      expand = grid::unit(1.5, "mm")
    ) +
  ggplot2::xlab(paste0(x, ": ", percent_var[as.numeric(stringr::str_sub(x, -1))], "% variance")) +
  ggplot2::ylab(paste0(y, ": ", percent_var[as.numeric(stringr::str_sub(y, -1))], "% variance")) +
    ggplot2::coord_fixed() +
    ggplot2::theme_classic(base_size = 20) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
  return(pca)
}
