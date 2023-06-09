#' Creates an interactive Seurat SpatialFeaturePlot allowing to probe for different taxa and visualize the UMI counts of taxa on the tissue image.
#'
#' @param object List. Output from decontaminate() function, containing $seurat_object and $taxid_counts.
#' @param taxonomizrDB A string. Path to nameNode.sqlite database required for taxonomic conversions (see README file for how to download).
#'
#' @return a Seurat::SpatialFeaturePlot
#' @export
#'
#' @examples
#' if (FALSE){
#'  interactive_spatialPlot(CRC_16)
#' }
#'
interactive_spatialPlot <- function(object, taxonomizrDB="/projects/site/pred/microbiome/database/taxonomizr_DB/nameNode.sqlite"){

  if(!file.exists(taxonomizrDB)){
    stop("file in taxonomizrDB doesn't exist. Please download SQLite database for taxonomizr by following the instructions in the README.")
  }

  #change dimnames from taxid to taxonomic names
  orig_assay <- object[["seurat_object"]]@assays[["TAXA_G"]]@counts
  taxnames <- taxonomizr::getTaxonomy(rownames(orig_assay), taxonomizrDB, desiredTaxa = "genus")
  rownames(orig_assay) <- taxnames
  orig_assay <- Seurat::CreateAssayObject(counts = orig_assay)
  object$seurat_object[["TAXA_G"]] <- orig_assay

  #create additional superkingdom assays
  bacterial_taxnames <- object$taxid_counts[object$taxid_counts$superkingdom=="Bacteria",]$genus
  if (length(bacterial_taxnames>0)){
    bacterial_assay <- Seurat::CreateAssayObject(counts = orig_assay[bacterial_taxnames,])
    object$seurat_object[["bacterial"]] <- bacterial_assay
  }

  viral_taxnames <- object$taxid_counts[object$taxid_counts$superkingdom=="Viruses",]$genus
  if (length(viral_taxnames>0)){
    viral_assay <- Seurat::CreateAssayObject(counts = orig_assay[viral_taxnames,])
    object$seurat_object[["viral"]] <- viral_assay
  }

  eukaryotic_taxnames <- object$taxid_counts[object$taxid_counts$superkingdom=="Eukaryota",]$genus
  if (length(eukaryotic_taxnames>0)){
    eukaryotic_assay <- Seurat::CreateAssayObject(counts = orig_assay[eukaryotic_taxnames,])
    object$seurat_object[["eukaryotic"]] <- eukaryotic_assay
  }

  #create interactive Shiny plot
  Seurat::SpatialFeaturePlot(object=object$seurat_object, features="nCount_TAXA_G", interactive=TRUE)
}
