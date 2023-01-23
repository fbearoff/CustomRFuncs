#' clusteredHeatmap
#'
#' This function returns a clustered heatmap of the top N genes by adjusted p-value
#'
#' @param condition1 Reference condition
#' @param condition2 Experimental condition
#' @param top_count How many genes to plot
#' @returns A ComplexHeatmap object
#' @importFrom dplyr "%>%"
#' @importFrom rlang .data
#' @export
clusteredHeatmap <- function(condition1, condition2, top_count = 50) {
  condition <- NULL
  if (!exists("samples")) {
    stop("`samples` variable not defined in global environment")
  }
  samples <- samples
  if (!exists("txi")) {
    stop("`txi` variable not defined in global environment")
  }
  txi <- txi
  if (!exists("tx2gene")) {
    stop("`tx2gene` variable not defined in global environment")
  }
  tx2gene <- tx2gene

  comparison <- paste0(condition1, " vs ", condition2)
  z <- CustomRFuncs::compDESeq2(condition1, condition2)

  up_genes <- z %>%
    dplyr::arrange(.data$padj) %>%
    dplyr::filter(.data$log2FoldChange > 1) %>%
    dplyr::select(.data$gene_id) %>%
    dplyr::slice_head(n = top_count / 2)

  down_genes <- z %>%
    dplyr::arrange(.data$padj) %>%
    dplyr::filter(.data$log2FoldChange < -1) %>%
    dplyr::select(.data$gene_id) %>%
    dplyr::slice_head(n = top_count / 2)

  top_genes <- dplyr::bind_rows(up_genes, down_genes)
  abund <- txi$abundance[, !(colnames(txi$abundance) %in%
    row.names(subset(
      samples,
      condition != condition1 &
        condition != condition2
    )))]
  abund <- abund[which(rownames(abund) %in% top_genes$gene_id), ]
  rownames(abund) <- tx2gene$gene_symbol[match(
    rownames(abund),
    tx2gene$gene_id
  )]
  abund_scale <- t(scale(t(abund),
    center = TRUE,
    scale = TRUE
  ))
  # anno_col <- samples[which(samples$condition %in% c(condition1, condition2)), ]
  #
  # col_names <- list(condition = stats::setNames(
  #   viridis::viridis(2, begin = 0.2, end = 0.8),
  #   c(condition1, condition2)
  # ))

  ha <- ComplexHeatmap::HeatmapAnnotation(
    group = ComplexHeatmap::anno_block(
      gp = grid::gpar(fill = c(2, "purple")),
      labels = c(condition2, condition1),
      labels_gp = grid::gpar(col = "white")
    )
  )
  # ra <- ComplexHeatmap::rowAnnotation(padj_order = ComplexHeatmap::anno_text(
  #   paste0(
  #     "(",
  #     rank(z[z$gene_symbol %in% rownames(abund), ]$padj,
  #       ties.method = "first"
  #     ),
  #     ") "
  #   ),
  #   gp = grid::gpar(
  #     fontsize = 6,
  #     fontface = "bold"
  #   ),
  #   just = "center",
  #   location = 0.5,
  #   show_name = FALSE
  # ))

  hm <- ComplexHeatmap::Heatmap(abund_scale,
    column_title = paste0("Top ", top_count, " DE Genes in ", comparison),
    column_title_gp = grid::gpar(fontsize = 18),
    top_annotation = ha,
    row_names_gp = grid::gpar(fontsize = 6),
    column_split = 2,
    show_column_names = FALSE,
    show_row_names = TRUE,
    show_parent_dend_line = FALSE,
    col = viridis::plasma(255, direction = -1),
    heatmap_legend_param = list(
      title = "Expression",
      legend_direction = "horizontal",
      title_position = "topcenter",
      labels_gp = grid::gpar(fontsize = 12)
    )
  )

  return(ComplexHeatmap::draw(hm,
    merge_legend = TRUE,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom")
  )
}
