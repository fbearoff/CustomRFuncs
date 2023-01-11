#' Plot MA
#'
#' This function returns a MA plot of the supplied conditions from a DESeq2 object.
#'
#' @param condition1 Reference condition
#' @param condition2 Experimental condition
#' @param top_n_fc Number of DE genes to be highlighted per direction
#' @returns A ggplot of the MA plot
#' @export
plotMA <- function(condition1, condition2, top_n_fc = 12) {
  z <- compDESeq2(condition1, condition2)
  top_fc <- z[padj <= 0.05][baseMean > 10][order(log2FoldChange, decreasing = TRUE)] %>% slice_head(n = top_n_fc)
  top_fc <- rbind(
    top_fc,
    z[padj <= 0.05][baseMean > 10][order(log2FoldChange, decreasing = TRUE)] %>% slice_tail(n = top_n_fc)
  )
  p <- ggplot(
    z,
    aes(
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
  p <- p + geom_point(na.rm = TRUE) +
    geom_text_repel(
      data = top_fc,
      aes(
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
    scale_x_log10(limit = c(0.1, 1e6)) +
    scale_y_continuous(trans = "pseudo_log") +
    scale_color_manual(
      labels = c("padj>0.05", "DE Transcripts", "padj=NA"),
      values = plasma(2, begin = .2, end = .7),
      guide_legend(title = "")
    ) +
    scale_shape_manual(
      guide = "none",
      na.value = 2,
      values = c(0, 1, 2)
    ) +
    scale_size(
      guide = "none"
    ) +
    labs(
      x = "Mean Expression",
      y = "Log2 Fold Change"
    ) +
    guides(color = guide_legend(override.aes = list(
      shape = c(0, 1, 2),
      size = 4
    ))) +
    annotate(
      geom = "richtext",
      y = max(z$log2FoldChange, na.rm = TRUE) / 2,
      x = log(mean(z$baseMean, na.rm = TRUE)) / 10,
      label = paste0(
        nrow(z[padj <= 0.05 & padj != "NA" &
          log2FoldChange > 1 &
          baseMean > 10]),
        " Genes <i>Upregulated</i> in ", "<b>",
        condition1, "</b>"
      ),
      hjust = "middle",
      size = 6
    ) +
    annotate(
      geom = "richtext",
      y = min(z$log2FoldChange, na.rm = TRUE) / 2,
      x = log(mean(z$baseMean, na.rm = TRUE)) / 10,
      label = paste0(
        nrow(z[padj <= 0.05 & padj != "NA" &
          log2FoldChange < -1 &
          baseMean > 10]),
        " Genes <i>Downregulated</i> in <b>",
        condition1, "</b>"
      ),
      hjust = "middle",
      size = 6
    ) +
    theme_minimal() +
    theme(
      plot.margin = grid::unit(c(0.5, 0.5, 0.5, 0.5), "in"),
      legend.text = element_text(size = 20),
      legend.position = "top",
      legend.justification = "right",
      axis.title.x = element_text(size = 20),
      axis.title.y = element_text(size = 20),
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      plot.title = element_text(
        hjust = 0.5,
        size = 40
      )
    )
  return(p)
}
