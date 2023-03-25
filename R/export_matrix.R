#' Exports the decontaminated taxa-spot matrix to desired file format for external downstream use.
#'
#' @param object List. Output from decontaminate() function.
#' @param fileFormat String. File format of exported matrix. Defaults to "csv".
#' @param outDir String. Path to desired output directory. Matrix will be saved as either "microbiome_genus_matrix" or "microbiome_taxid_matrix"
#' @param rownames String. Whether rownames of matrix should be taxid ("taxid") or the taxa names at genus level ("genus"). Defaults to "genus".
#' @param taxonomizrDB A string. Path to nameNode.sqlite database required for taxonomic conversions (see README file for how to download).
#'
#' @return none
#' @export
#'
#' @examples
#' export_matrix(object=CRC_16, rownames="genus", fileFormat="csv",
#' outDir=system.file("extdata", "CRC_16/", package="microbiome10XVisium"))
#'
export_matrix <- function(object, rownames="genus", fileFormat="csv", outDir, taxonomizrDB="/projects/site/pred/microbiome/database/taxonomizr_DB/nameNode.sqlite"){

  ##check parameters

  if (!is.list(object)){
    stop("object is not a list and thus probably not the output of the decontaminate() function.")
  }

  if (is.null(rownames)){
    rownames <- "genus"
    warning("rownames set to default: genus")
  } else if (!(rownames %in% c("genus", "taxid"))){
    stop("rownames has to be either read_counts or umi_counts")
  }

  if (is.null(fileFormat)){
    fileFormat <- "csv"
    warning("fileFormat set to default: csv")
  }

  #create output directory if it doesn't exist yet
  if(!dir.exists(outDir)){
    dir.create(outDir, showWarnings = FALSE)
  }

  matrix <- as.matrix(object$matrix)

  if (rownames=="genus"){
    rownames(matrix) <- suppressWarnings(taxonomizr::getTaxonomy(ids=rownames(matrix),sqlFile=taxonomizrDB, desiredTaxa=c("genus")))
  }

  if (fileFormat=="csv"){
    if(rownames=="genus"){
      utils::write.csv(matrix, file=paste0(outDir, "/microbiome_genus_matrix.csv"))
    } else if (rownames=="taxid"){
      utils::write.csv(matrix, file=paste0(outDir, "/microbiome_taxid_matrix.csv"))
    }

  }


}
