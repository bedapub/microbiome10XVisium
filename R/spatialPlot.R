#' Creates a Seurat SpatialFeaturePlot visualizing the number of UMIs or reads of taxa on the tissue image.
#'
#' @param object A list. Output from decontaminate() function.
#' @param taxa A string. Either one of "all", "bacterial", "eukaryotic", "viral", "nTaxa" (i.e. indicating the number of different taxa detected per spot) or a taxa name at genus level, e.g. "Fusobacterium".
#' @param taxonomizrDB A string. Path to nameNode.sqlite database required for taxonomic conversions (see README file for how to download).
#'
#' @return a Seurat::SpatialFeaturePlot
#' @export
#'
#' @examples
#' \dontrun{
#' spatialPlot(CRC_16, taxa="Fusobacterium")
#' spatialPlot(CRC_16, taxa="all")
#' spatialPlot(CRC_16, taxa="eukaryotic")
#' spatialPlot(CRC_16, taxa="nTaxa")
#' }
spatialPlot <- function(object, taxa, taxonomizrDB="/projects/site/pred/microbiome/database/taxonomizr_DB/nameNode.sqlite"){

  if(!file.exists(taxonomizrDB)){
    stop("file in taxonomizrDB doesn't exist. Please download SQLite database for taxonomizr by following the instructions in the README.")
  }

  if(taxa %in% c("all", "bacterial", "eukaryotic", "viral")){
    if (taxa=="all"){
      features="nCount_TAXA_G"
      title <- "total microbiome"
    } else if (taxa=="bacterial"){
      bacterial_taxid <- object$taxid_counts[object$taxid_counts$superkingdom=="Bacteria",]$taxid
      object$seurat_object <- subset(object$seurat_object, features = bacterial_taxid)
      title <- "total bacteria"
    } else if(taxa=="eukaryotic"){
      eukaryotic_taxid <- object$taxid_counts[object$taxid_counts$superkingdom=="Eukaryota",]$taxid
      if (length(eukaryotic_taxid)==0){
        stop("No eukaryotic taxa detected in this sample.")
      }
      object$seurat_object <- subset(object$seurat_object, features = eukaryotic_taxid)
      title <- "total eukaryotes"
    } else if(taxa=="viral"){
      viral_taxid <- object$taxid_counts[object$taxid_counts$superkingdom=="Viruses",]$taxid
      if (length(viral_taxid)==0){
        stop("No viral taxa detected in this sample.")
      }
      object$seurat_object <- subset(object$seurat_object, features = viral_taxid)
      title <- "total viruses"
    }
    features <- "nCount_TAXA_G"
  } else if (taxa=="nTaxa") {
    features="nFeature_TAXA_G"
  } else{
    features=taxonomizr::getId(taxa, taxonomizrDB)
    #select first taxid in case of multiple taxids (first is usually bacterial)
    features <- gsub("\\,.*","",features)

    if (is.na(features)){
      stop("taxa not found. please check spelling of requested taxa.")
    } else {
      #check whether taxa is present in sample
      if (!(taxa %in% object$taxid_counts$genus)){
        stop("taxa not present in sample.")
      }
    }
  }

  p <- suppressWarnings(Seurat::SpatialFeaturePlot(object=object$seurat_object,
                            features=features,
                            alpha=c(0.7,1)))
  p <- p + ggplot2::theme(legend.position = "right", legend.title=ggplot2::element_blank(),
                   legend.text = ggplot2::element_text(size=14))

  if (taxa=="nTaxa"){
    p + ggplot2::ggtitle("Number of taxa detected per spot") + ggplot2::theme(plot.title=ggplot2::element_text(size=22))
  } else if(taxa %in% c("all", "bacterial", "eukaryotic", "viral")){
    p + ggplot2::ggtitle(title) + ggplot2::theme(plot.title=ggplot2::element_text(size=22))
  } else {
    p + ggplot2::ggtitle(taxa) + ggplot2::theme(plot.title=ggplot2::element_text(size=22, face="italic"))
  }

}
