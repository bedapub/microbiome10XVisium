#' Performs various decontamination steps for scRNA-seq datasets on the taxid-spot matrix produced by krakenToMatrix().
#'
#' @param sampleName A string. Sample name.
#' @param filePath A string. Path to RDS object created with krakenToMatrix(). (note: either filePath or object has to be supplied)
#' @param object A list. Name of the R object from the output of krakenToMatrix(). (note: either filePath or object has to be supplied)
#' @param outDir A string. Path to output directory. RDS object saved as sampleName_decontaminated.RDS.
#' @param removeSingletons A logical. Whether to remove umi_counts==1 in spots and prevalence(taxid)==1. Only applicable if counts=="umi_counts" in krakenToMatrix(). Defaults to TRUE.
#' @param removeLikelyContaminants A logical. Whether to remove taxa defined as likely contaminants by Poore et al. (Nature 2020). Only applicable if tax_level=="genus" in krakenToMatrix(). Defaults to TRUE.
#' @param removeSpecificTaxa A vector of strings. Vector of genus names of taxa to remove. Defaults to NULL.
#' @param selectGastrointestinal A logical. Whether to only keep taxa defined as oral, fecal or oral/fecal transmitter by Schmidt et al. (2019 eLife). Only applicable if tax_level=="genus" in krakenToMatrix. Defaults to FALSE.
#' @param taxonomizrDB A string. Path to nameNode.sqlite database required for taxonomic conversions (see README file for how to download).
#' @param selectSpecificTaxa A vector of strings. Vector of genus names of taxa to keep in addition to the selectGastrointestinatl. Note: only applicable if selectGastrointestinal==TRUE. Defaults to NULL.
#'
#' @return  A list of $matrix (decontaminated taxid-cell matrix) and $taxid_counts (dataframe with taxids, taxa name and counts).
#' @export
#'
#' @examples
#' \dontrun{
#' decontaminateSingleCell(sampleName="P1_D0", filePath="FILEPATH",
#' outDir="OUTDIR", removeSpecificTaxa=c("Mycobacterium"))
#' }
#'
decontaminateSingleCell <- function(sampleName, filePath=NULL, object=NULL, outDir, removeSingletons=TRUE, removeLikelyContaminants=TRUE, removeSpecificTaxa=NULL, selectGastrointestinal=FALSE, selectSpecificTaxa=TRUE, taxonomizrDB="/projects/site/pred/microbiome/database/taxonomizr_DB/nameNode.sqlite"){

  ##check parameters

  if (is.null(sampleName)){
    stop("please provide sample name in sampleName")
  }

  if (is.null(removeSingletons)){
    removeSingletons <- TRUE
    warning("removeSingletons set to default: TRUE")
  } else if (!is.logical(removeSingletons)){
    stop("removeSingletons has to be a logical (TRUE or FALSE)")
  }

  if (is.null(removeLikelyContaminants)){
    removeLikelyContaminants <- TRUE
    warning("removeLikelyContaminants set to default: TRUE")
  } else if (!is.logical(removeLikelyContaminants)){
    stop("removeLikelyContaminants has to be a logical (TRUE or FALSE)")
  }

  if (!is.logical(selectGastrointestinal)){
    stop("selectGastrointestinal has to be a logical (TRUE or FALSE)")
  }

  if (!is.null(removeSpecificTaxa)){
    for (l in length(removeSpecificTaxa)){
      if (is.na(taxonomizr::getId(removeSpecificTaxa[l], taxonomizrDB)))
        stop("Taxonomic name in removeSpecificTaxa could not be found. Please double-check spelling.")
    }
  }

  if (is.null(filePath) && is.null(object)){
    stop("please provide either filePath to saved RDS object or object (already loaded in R)")
  }

  if(!file.exists(filePath)){
    stop("file in filePath doesn't exist")
  }

  if(!file.exists(taxonomizrDB)){
    stop("file in taxonomizrDB doesn't exist. Please download SQLite database for taxonomizr by following the instructions in the README.")
  }

  #create output directory if it doesn't exist yet
  if(!dir.exists(outDir)){
    dir.create(outDir, showWarnings = FALSE)
  }

  ## load input=output from krakenToMatrix() function
  if (!is.null(filePath)){
    input <- readRDS(file=filePath)
  } else if (!is.null(object)){
    input <- object
  }
  taxid_matrix <- input$matrix
  taxid_matrix_raw <- taxid_matrix
  taxid_counts_raw <- input$taxid_counts

  print(paste0(dim(taxid_matrix)[1], " taxa present before contamination"))
  print(paste0(sum(taxid_matrix), " counts present before contamination"))

  ### DECONTAMINATION STEP#1: remove singleton UMI or read counts (i.e. barcodes with only 1 UMI count for a taxid)
  if(removeSingletons==TRUE){
    print("Decontamination step#1: removing singleton counts")
    taxid_matrix[taxid_matrix==1] <- 0
    #remove taxids with only zeros
    taxid_matrix <- taxid_matrix[Matrix::rowSums(taxid_matrix)>0,]
    print(paste0((dim(taxid_matrix_raw)[1]-dim(taxid_matrix)[1]), " taxa eliminated."))
    print(paste0((sum(taxid_matrix_raw)-sum(taxid_matrix)), " counts eliminated."))
    #print(paste0(sum(taxid_matrix), " counts remaining."))
    taxid_matrix_raw <- taxid_matrix
  }

  ### DECONTAMINATION STEP #2: remove likely contaminants defined by Poore et al. or user-defined
  if(removeLikelyContaminants==TRUE){
    print("Decontamination step#2: removing taxa that are likely contaminants")
    #import list with likely contaminants from Poore et al.
    contaminants <- readr::read_delim(file=system.file("extdata", "contaminants.txt", package="microbiome10XVisium", mustWork=TRUE),
                                      delim = "\t", skip = 1, show_col_types = FALSE) %>%
      dplyr::filter(Category %in% c("LIKELY CONTAMINANT")
      )
    contaminants_taxa <- unique(contaminants$Genera)
    #convert taxa to taxids
    suppressWarnings(contaminants_id <- taxonomizr::getId(taxa=contaminants_taxa, taxonomizrDB))
    contaminants_id <- gsub("\\,.*","",contaminants_id)
    contaminants_id <- contaminants_id[!is.na(contaminants_id)]

    if (!is.null(removeSpecificTaxa)){
      suppressWarnings(additional_contaminants_id <- taxonomizr::getId(taxa=removeSpecificTaxa, taxonomizrDB))
      additional_contaminants_id <- gsub("\\,.*","",additional_contaminants_id)
      additional_contaminants_id <- additional_contaminants_id[!is.na(additional_contaminants_id)]
      contaminants_id <- append(contaminants_id, additional_contaminants_id)
    }

    #remove contaminants from matrix
    contaminant_rownames <- taxid_matrix_raw@Dimnames[[1]] %in% contaminants_id

    #in case all taxa are contaminant taxa (or only one taxa is left) return NULL
    if (length(taxid_matrix@Dimnames[[1]])-sum(contaminant_rownames)==0){
      result <- NULL
      saveRDS(result, file=paste0(outDir, "/", sampleName, "_Seurat_object.RDS"))
      print("No taxa left after decontamination. Returning NULL object.")
      return(result)
    }

    eliminated_taxids <- taxid_matrix@Dimnames[[1]][contaminant_rownames]
    suppressWarnings(eliminated_taxnames <- taxonomizr::getTaxonomy(ids=eliminated_taxids, taxonomizrDB, desiredTaxa = c("genus")))
    taxid_matrix <- taxid_matrix[!contaminant_rownames,]
    #remove taxids with only zeros
    taxid_matrix <- taxid_matrix[Matrix::rowSums(taxid_matrix)>0,]
    print(paste0((dim(taxid_matrix_raw)[1]-dim(taxid_matrix)[1]), " taxa eliminated:"))
    print(paste0(eliminated_taxnames, " genus was removed"))
    eliminated_taxids <- data.frame(taxid=eliminated_taxids)
    eliminated_taxa <- dplyr::left_join(eliminated_taxids,taxid_counts_raw, by="taxid")
    print(paste0(sum(eliminated_taxa$counts), " UMI_counts eliminated"))
    #print(paste0(dim(taxid_matrix)[1], " taxa remaining"))
    #print(paste0(sum(taxid_matrix), " counts remaining"))
    taxid_matrix_raw <- taxid_matrix
  }

  ### DECONTAMINATION STEP#3: selecting only taxa that were defined to be gastrointestinal commensals
  if (selectGastrointestinal==TRUE){
    cat("\n")
    print(paste0("Decontamination step: only keeping taxa that are gastrointestinal commensals"))
    #import list with gastrointestinal taxa at genus level
    commensals <- readr::read_delim(file=system.file("extdata", "schmidt_elife_gastrointestinal_tract_genus_level.txt", package="microbiome10XVisium", mustWork=TRUE),
                                    delim = "\t", col_names = "genus",show_col_types = FALSE)
    #convert taxa to taxids
    suppressWarnings(commensals_id <- taxonomizr::getId(taxa=commensals$genus, taxonomizrDB))
    commensals_id <- gsub("\\,.*","",commensals_id)
    commensals_id <- commensals_id[!is.na(commensals_id)]

    ### DECONTAMINATION STEP#4: selecting additional taxa to keep
    if (!is.null(selectSpecificTaxa)){
      suppressWarnings(additional_commensals_id <- taxonomizr::getId(taxa=selectSpecificTaxa, taxonomizrDB))
      additional_commensals_id <- gsub("\\,.*","",additional_commensals_id)
      additional_commensals_id <- additional_commensals_id[!is.na(additional_commensals_id)]
      commensals_id <- append(commensals_id, additional_commensals_id)
    }

    #only keep commensals in matrix
    commensals_rownames <- taxid_matrix_raw@Dimnames[[1]] %in% commensals_id

    #in case no taxa are commensals  return NULL
    if (length(commensals_rownames)==0){
      result <- NULL
      saveRDS(result, file=paste0(outDir, "/", sampleName, "_Seurat_object.RDS"))
      print("No taxa left after decontamination. Returning NULL object.")
      return(result)
    }

    eliminated_taxids <- taxid_matrix@Dimnames[[1]][!commensals_rownames]
    taxid_matrix <- taxid_matrix[commensals_rownames,]
    #remove taxids with only zeros
    taxid_matrix <- taxid_matrix[Matrix::rowSums(taxid_matrix)>0,]
    print(paste0((dim(taxid_matrix_raw)[1]-dim(taxid_matrix)[1]), " taxa eliminated."))
    eliminated_taxids <- data.frame(taxid=eliminated_taxids)
    eliminated_taxa <- dplyr::left_join(eliminated_taxids,taxid_counts_raw, by="taxid")
    print(paste0(sum(eliminated_taxa$counts), " counts eliminated"))
    #print(paste0(dim(taxid_matrix)[1], " taxa remaining"))
    #print(paste0(sum(taxid_matrix), " counts remaining"))
  }

  if (length(taxid_matrix@Dimnames[[1]])!=0){
    ### update taxid_counts dataframe
    suppressWarnings(dataframe <- data.frame(counts=Matrix::rowSums(taxid_matrix), taxonomizr::getTaxonomy(ids=rownames(taxid_matrix), taxonomizrDB, desiredTaxa = c("superkingdom", "phylum", "genus")),
                            taxid=rownames(taxid_matrix)))
    rownames(dataframe) <- NULL
    dataframe <- dataframe %>% dplyr::filter(counts>0) %>% dplyr::arrange(dplyr::desc(counts))
    print(paste0(nrow(dataframe), " taxa remaining after decontamination"))
    print(paste0(sum(dataframe$counts), " counts remaining after decontamination"))

    ### return and save Seurat object & dataframe
    result <- list("matrix"=taxid_matrix, "taxid_counts"=dataframe)
    saveRDS(result, file=paste0(outDir, "/", sampleName, "_decontaminated.RDS"))
    return(result)
  }
}
