#' compDESeq2
#'
#' This function returns a merged data.table object from a DESeq2 results function.
#' `samples`, `txi`, and `tx2gene` objects need to be defined in the global environment
#'
#' @param condition1 Reference condition
#' @param condition2 Experimental condition
#' @returns A data.table object
#' @export
compDESeq2 <- function(condition1, condition2) {
  condition <- NULL
  if (!exists("samples")) {
    stop("`samples` variable not defined in global environment")
  }
  samples <- samples

  if (!exists("tx2gene")) {
    stop("`tx2gene` variable not defined in global environment")
  }
  tx2gene <- tx2gene
  if (!exists("txi")) {
    stop("`txi` variable not defined in global environment")
  }
  txi <- txi
  if (!exists("dds")) {
    stop("`dds` variable not defined in global environment")
  }
  dds <- dds

  res <- DESeq2::results(
    dds,
    contrast = c("condition", condition1, condition2),
    alpha = 0.05
  )

  gene_synonym <- unique(tx2gene[, -1])

  z <- data.frame(res)
  z$gene_symbol <- gene_synonym$gene_symbol[match(rownames(res), gene_synonym$gene_id)]
  z$chr <- gene_synonym$chr[match(rownames(res), gene_synonym$gene_id)]
  z$start <- gene_synonym$start[match(rownames(res), gene_synonym$gene_id)]
  z$end <- gene_synonym$end[match(rownames(res), gene_synonym$gene_id)]
  z$description <- gene_synonym$description[match(rownames(res), gene_synonym$gene_id)]
  z <- merge(z, txi$abundance, by = 0)
  data.table::setnames(z, "Row.names", "gene_id")
  z <-
    z[, !(names(z) %in% row.names(subset(
      samples,
      condition != condition1 &
        condition != condition2
    )))]
  z <- data.table::as.data.table(z)
  z <- z[z[["chr"]] %in% c(1:22, "MT", "X", "Y"), ]
  z <- z[z[["gene_symbol"]] != "", ]
  return(z)
}
