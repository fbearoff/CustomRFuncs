#' Get a tx2gene dataframe for the given species
#'
#' This function returns a tx2gene dataframe object for the queried species
#'
#' @param species A valid ENSEMBL species abbreviation e.g. `mmusculus, rnorvegicus, hsapiesn`
#' @return A dataframe of the tx2gene object
#' @importFrom dplyr "%>%"
#' @export


getTx2gene <- function(species) {
  chromosome_name <- NULL
  mart <- biomaRt::useMart("ensembl")
  ensembl <- biomaRt::useDataset(paste0(species, "_gene_ensembl"), mart)
  tx2gene <- biomaRt::getBM(
    attributes = c(
      "ensembl_transcript_id_version",
      "ensembl_gene_id_version",
      "external_gene_name",
      "description",
      "chromosome_name",
      "start_position",
      "end_position"
    ),
    mart = ensembl
  )
  tx2gene <- tx2gene %>%
    dplyr::filter(chromosome_name %in% c(1:22, "MT", "X", "Y")) %>%
    dplyr::rename(
      "tx_id" = "ensembl_transcript_id_version",
      "gene_id" = "ensembl_gene_id_version",
      "gene_symbol" = "external_gene_name",
      "chr" = "chromosome_name",
      "start" = "start_position",
      "end" = "end_position"
    )
  return(tx2gene)
}
