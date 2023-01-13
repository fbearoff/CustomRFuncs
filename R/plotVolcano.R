#' Plot Volcano
#'
#' This function returns a volcano plot of the supplied conditions from a DESeq2 object.
#'
#' @param condition1 Reference condition
#' @param condition2 Experimental condition
#' @returns A ggplot of the volcano plot
#' @export
plotVolcano <- function(condition1, condition2) {

  z <- CustomRFuncs::compDESeq2(condition1, condition2)
  comparison <- paste(condition1, "vs", condition2)
  EnhancedVolcano::EnhancedVolcano(z,
    lab = z$gene_symbol,
    title = bquote(italic(.(comparison))),
    subtitle = "Volcano plot",
    caption = NULL,
    x = "log2FoldChange",
    y = "padj",
    ylab = bquote(~ -Log[10] ~ italic("Padj")),
    pCutoff = 0.05,
    FCcutoff = 1,
    legendLabels = c(
      "NS",
      bquote(~ Log[2] ~ "FC"),
      "p-adj",
      bquote("p-adj &" ~ Log[2] ~ "FC")
    ),
    drawConnectors = FALSE,
    pointSize = c(ifelse(abs(z$log2FoldChange) > 2 & z$padj < 0.05, 8, 1)),
  ) +
    ggplot2::coord_cartesian(xlim = c(-6, 6)) +
    ggplot2::scale_x_continuous(
      breaks = seq(-6, 6, 1)
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5)
    )
}
