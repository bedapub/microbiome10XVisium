#' CRC_16 toy data
#'
#' Example output of the decontaminate() function to be used for testing of
#' downstream analyis and visualization, such as spatialPlot(),
#' pseudoBulkProfile() or cooccurrenceNetwork().
#'
#' Original data from Galeano-Niño et al. (Nature 2022) publication.
#' Processed with Roche bioinformatic pipeline and default settings in
#' microbiome10XVisium krakenToMatrix() and decontaminate().
#'
#'
#' @format ## `CRC_16`
#' A list with 3 objects: $seurat_object, $taxid_counts, $matrix
#' \describe{
#'   \item{CRC_16$seurat_object}{An object of class seurat with two assays: GEX (host transcriptomics) and TAXA_G (microbiome)}
#'   \item{CRC_16$taxid_counts}{A dataframe with the column counts, superkingdom, phylum, genus and taxid}
#'   \item{CRC_16$matrix}{A sparce Matrix of class "dgCMatrix", with rows=taxa and columns=spots}
#' }
#' @source <https://doi.org/10.1038/s41586-022-05435-0>
"CRC_16"
