#' Plot Gene
#'
#' This function returns a data plot of the queried gene.
#'
#' @param gene_name the gene to be plotted
#' @export
plotGene <- function(gene_name) {
  gene_id <- tx2gene[grep(gene_name, gene_symbol, ignore.case = TRUE) & chr %in% c(1:22, "X", "Y", "MT"), gene_id][1]
  plot_data <- DESeq2::plotCounts(dds,
    gene = as.character(gene_id),
    returnData = TRUE
  )
  ggplot2::ggplot(plot_data, ggplot2::aes(
    x = condition,
    y = count,
    color = condition
  )) +
    ggplot2::geom_boxplot(show.legend = TRUE) +
    ggplot2::scale_color_manual(values = viridis::plasma(5)) +
    ggplot2::ggtitle(paste0(gene_name, "\n", gene_id)) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = .5))
}
