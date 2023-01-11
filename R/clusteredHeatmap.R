#' clusteredHeatmap
#'
#' This function returns a clustered heatmap of the top N genes by adjusted p-value
#'
#' @param condition1 Reference condition
#' @param condition2 Experimental condition
#' @returns A ComplexHeatmap object
#' @export
clusteredHeatmap <- function(condition1, condition2, top_count=50) {
  comparison <- paste0(condition1, " vs ", condition2)
  z <- compDESeq2(condition1, condition2)

  up_genes <- z[order(padj)][log2FoldChange > 1]$gene_id[1:top_count]
  down_genes <- z[order(padj)][log2FoldChange < -1]$gene_id[1:top_count]

  top_genes <- c(up_genes, down_genes)

  abund <- txi$abundance[, !(colnames(txi$abundance) %in%
    row.names(subset(
      samples,
      condition != condition1 &
      condition != condition2
    )))]
  abund <- abund[which(rownames(abund) %in% top_genes), ]
  rownames(abund) <- tx2gene$gene_symbol[match(
    rownames(abund),
    tx2gene$gene_id
  )]
  abund_scale <- t(scale(t(abund),
    center = TRUE,
    scale = TRUE
  ))
  anno_col <- samples[which(samples$condition %in% c(condition1, condition2)), ]

  col_names <- list(condition = setNames(
    viridis(2, begin = 0.2, end = 0.8),
    c(condition1, condition2)
  ))

  ha <- HeatmapAnnotation(
    group = anno_block(
      gp = gpar(fill = c(2, "purple")),
      labels = c(condition2, condition1),
      labels_gp = gpar(col = "white")
    )
  )

  ra <- rowAnnotation(padj_order = anno_text(
    paste0(
      "(",
      rank(z[gene_symbol %in% rownames(abund)]$padj,
        ties.method = "first"
      ),
      ") "
    ),
    gp = gpar(
      fontsize = 6,
      fontface = "bold"
    ),
    just = "center",
    location = 0.5,
    show_name = FALSE
  ))

  hm <- Heatmap(abund_scale,
    column_title = paste0("Top ", top_count, " DE Genes in ", comparison),
    column_title_gp = gpar(fontsize = 18),
    top_annotation = ha,
    row_names_gp = gpar(fontsize = 6),
    column_split = 2,
    show_column_names = FALSE,
    show_row_names = TRUE,
    show_parent_dend_line = FALSE,
    col = plasma(255, direction = -1),
    heatmap_legend_param = list(
      title = "Expression",
      legend_direction = "horizontal",
      title_position = "topcenter",
      labels_gp = gpar(fontsize = 6)
    )
  )

  return(draw(hm,
    merge_legend = TRUE,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom")
  )
}
