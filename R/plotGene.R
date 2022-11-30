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
  ggplot(plot_data, aes(
    x = condition,
    y = count,
    color = condition
  )) +
    geom_boxplot(show.legend = TRUE) +
    scale_color_manual(values = plasma(5)) +
    ggtitle(paste0(gene_name, "\n", gene_id)) +
    theme(plot.title = element_text(hjust = .5))
}
