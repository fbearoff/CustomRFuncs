#' Plot Gene
#'
#' This function returns a data plot of the queried gene.
#'
#' @param gene_name the gene to be plotted
#' @importFrom dplyr "%>%"
#' @importFrom rlang .data
#' @export
plotGene <- function(gene_name) {
  chr <- count <- condition <- gene_symbol <- NULL

  if (!exists("tx2gene")) {
    stop("`tx2gene` variable not defined in global environment")
  }
  tx2gene <- tx2gene
  if (!exists("dds")) {
    stop("`dds` variable not defined in global environment")
  }
  dds <- dds
  gene_id <- tx2gene %>%
    dplyr::filter(.data$chr %in% c(1:22, "X", "Y", "MT")) %>%
    dplyr::filter(.data$gene_symbol != "") %>%
    dplyr::filter(grepl("IL1B", .data$gene_symbol, ignore.case = TRUE)) %>%
    dplyr::select(.data$gene_id) %>%
    dplyr::slice_head(n = 1) %>%
    as.character()
  plot_data <- DESeq2::plotCounts(dds,
    gene = gene_id,
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
