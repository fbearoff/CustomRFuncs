#' Plot MA
#'
#' This function returns a MA plot of the supplied conditions from a DESeq2 object.
#'
#' @param condition1 Reference condition
#' @param condition2 Experimental condition
#' @param top_n_fc Number of DE genes to be highlighted per direction
#' @returns A ggplot of the MA plot
#' @importFrom dplyr "%>%"
#' @export
plotMA <- function(condition1, condition2, top_n_fc = 12) {
  z <- CustomRFuncs::compDESeq2(condition1, condition2)
  top_fc <- z %>%
    dplyr::filter(padj <= 0.05, baseMean > 10) %>%
    dplyr::arrange(dplyr::desc(log2FoldChange)) %>%
    dplyr::slice_head(n = 12)
  bottom_fc <- z %>%
    dplyr::filter(padj <= 0.05, baseMean > 10) %>%
    dplyr::arrange(dplyr::desc(log2FoldChange)) %>%
    dplyr::slice_tail(n = 12)

  top_fc <- dplyr::bind_rows(top_fc, bottom_fc)
  p <- ggplot2::ggplot(
    z,
    ggplot2::aes(
      x = baseMean + 0.01,
      y = log2FoldChange,
      color = ifelse(is.na(padj),
        padj > 0.05,
        padj < 0.05 &
          abs(log2FoldChange) > 1 &
          baseMean > 10
      ),
      shape = ifelse(is.na(padj),
        padj > 0.05,
        padj < 0.05 &
          abs(log2FoldChange) > 1 &
          baseMean > 10
      ),
      size = ifelse(is.na(padj), .1,
        ifelse(padj > 0.05 | abs(log2FoldChange) < 1 | baseMean < 10,
          0.1,
          log10(padj) * -1
        )
      )
    )
  )

  p <- p + ggplot2::geom_point(na.rm = TRUE) +
    ggrepel::geom_text_repel(
      data = top_fc,
      ggplot2::aes(
        label = paste0(
          gene_symbol,
          " (",
          round(log2FoldChange,
            digits = 1
          ),
          ")"
        ),
        size = 12
      ),
      show.legend = FALSE,
      color = "black",
      fontface = "bold"
    ) +
    ggplot2::scale_x_log10(limit = c(0.1, 1e6)) +
    ggplot2::scale_y_continuous(trans = "pseudo_log") +
    ggplot2::scale_color_manual(
      labels = c("padj>0.05", "DE Transcripts", "padj=NA"),
      values = viridis::plasma(2, begin = .2, end = .7),
      ggplot2::guide_legend(title = "")
    ) +
    ggplot2::scale_shape_manual(
      guide = "none",
      na.value = 2,
      values = c(0, 1, 2)
    ) +
    ggplot2::scale_size(
      guide = "none"
    ) +
    ggplot2::labs(
      x = "Mean Expression",
      y = "Log2 Fold Change"
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(
      shape = c(0, 1, 2),
      size = 4
    ))) +
    ggplot2::annotate(
      geom = ggtext::GeomRichText,
      y = max(z$log2FoldChange, na.rm = TRUE) / 2,
      x = log(mean(z$baseMean, na.rm = TRUE)) / 10,
      label = paste0(
        nrow(dplyr::filter(
          z,
          padj <= 0.05,
          padj != "NA",
          log2FoldChange > 1,
          baseMean > 10
        )),
        " Genes <i>Upregulated</i> in ", "<b>",
        condition1, "</b>"
      ),
      hjust = "middle",
      size = 6
    ) +
    ggplot2::annotate(
      geom = ggtext::GeomRichText,
      y = min(z$log2FoldChange, na.rm = TRUE) / 2,
      x = log(mean(z$baseMean, na.rm = TRUE)) / 10,
      label = paste0(
        nrow(dplyr::filter(
          z,
          padj <= 0.05,
          padj != "NA",
          log2FoldChange < -1,
          baseMean > 10
        )),
        " Genes <i>Downregulated</i> in <b>",
        condition1, "</b>"
      ),
      hjust = "middle",
      size = 6
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.margin = grid::unit(c(0.5, 0.5, 0.5, 0.5), "in"),
      legend.text = ggplot2::element_text(size = 20),
      legend.position = "top",
      legend.justification = "right",
      axis.title.x = ggplot2::element_text(size = 20),
      axis.title.y = ggplot2::element_text(size = 20),
      axis.text.x = ggplot2::element_text(size = 20),
      axis.text.y = ggplot2::element_text(size = 20),
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = 40
      )
    )
  return(p)
}
