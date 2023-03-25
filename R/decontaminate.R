#' Performs various decontamination steps on the taxid-spot matrix produced by krakenToMatrix().
#'
#' @param sampleName A string. Sample name.
#' @param filePath A string. Path to RDS object created with krakenToMatrix(). (note: either filePath or object has to be supplied)
#' @param object A list. Name of the R object from the output of krakenToMatrix(). (note: either filePath or object has to be supplied)
#' @param spacerangerDir A string. Path to the spaceranger outs folder.
#' @param outDir A string. Path to desired output directory. RDS object will be saved as sampleName_Seurat_object.RDS.
#' @param removeSingletons A logical. Whether to remove counts==1 in spots. Defaults to TRUE.
#' @param removeLikelyContaminants A logical. Whether to remove taxa defined as likely contaminants by Poore et al. (2020 Nature). Only applicable if tax_level=="genus" in krakenToMatrix. Defaults to TRUE.
#' @param spots A string. Which spots to consider, either one of c("tissueOnly", "tissuePlusBordering", "all"). Defaults to "tissueOnly".
#' @param distance An integer. Only required when spots=="tissuePlusBordering". Determines the size of the bordering region surrounding the tissue. Recommended values in range of (2,7). Defaults to 7.
#' @param spatiallyVariable A logical. Whether to perform spatial variability analysis and only keep spatially variable taxa. Defaults to FALSE.
#' @param removeSpecificTaxa A vector of strings. Vector of genus names of taxa to remove. Defaults to NULL.
#' @param selectGastrointestinal A logical. Whether to only keep taxa defined as oral, fecal or oral/fecal transmitter by Schmidt et al. (2019 eLife). Only applicable if tax_level=="genus" in krakenToMatrix. Defaults to FALSE.
#' @param taxonomizrDB A string. Path to nameNode.sqlite database required for taxonomic conversions (see README file for how to download).
#'
#' @return  A list of $seurat_object(Seurat object with decontaminated taxa assay),$taxid_counts (dataframe with taxids, taxonomic names and counts) and $matrix (taxid-spot matrix).
#' @export
#'
#' @examples
#' \dontrun{
#'  decontaminate(sampleName = "CRC_16",
#'  filePath=system.file("extdata", "CRC_16", "genus_umi_counts.RDS", package="microbiome10XVisium"),
#'  spacerangerDir = system.file("extdata", "outs/", package="microbiome10XVisium"),
#'  outDir = system.file("extdata", "CRC_16/", package="microbiome10XVisium"))
#' }
decontaminate <- function(sampleName, filePath=NULL, object=NULL, spacerangerDir, outDir, removeSingletons=TRUE, removeLikelyContaminants=TRUE, selectGastrointestinal=FALSE, spots="tissueOnly", distance=7, spatiallyVariable=FALSE, removeSpecificTaxa=NULL, taxonomizrDB="/projects/site/pred/microbiome/database/taxonomizr_DB/nameNode.sqlite"){

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

  if (is.null(spots)){
    spots <- "tissueOnly"
    warning("spots set to default: tissueOnly")
  } else if (!(spots %in% c("tissueOnly", "tissuePlusBordering", "all"))){
    stop("spots has to be either tissueOnly, tissuePlusBordering or all")
  }

  if (spots=="tissuePlusBordering" && is.null(distance)){
    distance<-7
    warning("distance threshold set to default of 7")
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

  if(!is.null(filePath) && !file.exists(filePath)){
    stop("file in filePath doesn't exist")
  }

  if(!file.exists(taxonomizrDB)){
    stop("file in taxonomizrDB doesn't exist. Please download SQLite database for taxonomizr by following the instructions in the README.")
  }

  if(!dir.exists(spacerangerDir)){
    stop("spaceranger outs directory in spacerangerDir doesn't exist")
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
    cat("\nDecontamination step: removing singleton counts")
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
    cat("\nDecontamination step: removing taxa that are likely contaminants")
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
    print(paste0(sum(eliminated_taxa$counts), " counts eliminated"))
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

  ### DECONTAMINATION STEP#4: create Seurat object and remove taxa in non-selected spots
  cat("\n")
  print(paste0("Decontamination step: only keeping taxa that are present in ", spots, " spots"))
  #creating Seurat object
  if (spots=="all"){
    #loading all spots (not only filtered)
    suppressWarnings(spatial <- Seurat::Load10X_Spatial(
      data.dir = spacerangerDir,
      filename = "raw_feature_bc_matrix.h5",
      slice = "all_spots",
      assay = "GEX",
      filter.matrix = FALSE))
    #in case spatial doesn't contain all 4992 spots
    #subsetting matrix with barcodes present in GEX matrix
    GEX_barcodes <- spatial@assays[["GEX"]]@data@Dimnames[[2]]
    valid_barcodes <- taxid_matrix@Dimnames[[2]] %in% GEX_barcodes
    taxid_matrix <- taxid_matrix[,valid_barcodes]
    # create a new assay to store ADT information
    microbial_assay <- Seurat::CreateAssayObject(counts = taxid_matrix)
    # add this assay to the previously created Seurat object
    suppressWarnings(spatial[["TAXA_G"]] <- microbial_assay)
    suppressWarnings(Seurat::DefaultAssay(spatial) <- "TAXA_G")

  } else if (spots=="tissueOnly"){
    #loading only tissue-covered (filtered) spots
    suppressWarnings(spatial <- Seurat::Load10X_Spatial(
      data.dir = spacerangerDir,
      filename = "filtered_feature_bc_matrix.h5",
      slice = "tissue_spots",
      assay = "GEX",
      filter.matrix = TRUE))
    GEX_barcodes <- spatial@assays[["GEX"]]@data@Dimnames[[2]]
    #subsetting matrix with barcodes present in GEX matrix
    valid_barcodes <- taxid_matrix@Dimnames[[2]] %in% GEX_barcodes
    taxid_matrix <- taxid_matrix[,valid_barcodes]
    #remove taxids with only zeros
    taxid_matrix <- taxid_matrix[Matrix::rowSums(taxid_matrix)>0,]
    #check if any taxa are left
    if (length(taxid_matrix@Dimnames[[1]])==0){
      result <- NULL
      saveRDS(result, file=paste0(outDir, "/", sampleName, "_Seurat_object.RDS"))
      print("No taxa left after decontamination. Returning NULL object.")
      return(result)
    }
    #print(paste0(dim(taxid_matrix)[1], " taxa remaining"))
    #print(paste0(sum(taxid_matrix), " counts remaining"))
    # create a new assay to store ADT information
    microbial_assay <- Seurat::CreateAssayObject(counts = taxid_matrix)
    # add this assay to the previously created Seurat object
    suppressWarnings(spatial[["TAXA_G"]] <- microbial_assay)
    suppressWarnings(Seurat::DefaultAssay(spatial) <- "TAXA_G")

  } else if (spots=="tissuePlusBordering"){
    #loading only tissue-covered spots
    suppressWarnings(spatial_filtered <- Seurat::Load10X_Spatial(
      data.dir = spacerangerDir,
      filename = "filtered_feature_bc_matrix.h5",
      slice = "tissue_spots",
      assay = "GEX",
      filter.matrix = TRUE
    ))
    #loading all spots (not only filtered)
    suppressWarnings(spatial <- Seurat::Load10X_Spatial(
      data.dir = spacerangerDir,
      filename = "raw_feature_bc_matrix.h5",
      slice = "all_spots",
      assay = "GEX",
      filter.matrix = FALSE
    ))
    #identify tissue_only_barcodes
    tissue_only_barcodes <- spatial_filtered@assays[["GEX"]]@data@Dimnames[[2]]
    all_barcodes <- spatial@assays[["GEX"]]@data@Dimnames[[2]]
    tissue_selection <- all_barcodes %in% tissue_only_barcodes
    #add metadata information with tissue/no_tissue
    spatial@meta.data$spot <- "no_tissue"
    spatial@meta.data[tissue_selection,]$spot <- "tissue_only"
    #compute distance matrix between spots
    geometry <- Seurat::GetTissueCoordinates(spatial, cols = c("row", "col"), scale = NULL)
    dist_matrix <- stats::dist(geometry, method = "euclidean") %>% as.matrix()
    tissue_spots <- spatial@meta.data %>% dplyr::filter(spot == "tissue_only") %>% rownames()
    no_tissue_spots <- spatial@meta.data %>% dplyr::filter(spot == "no_tissue") %>% rownames()
    #compute distance from the no_tissue_spots to nearest tissue spot
    spots_distances <- sort(apply(dist_matrix[no_tissue_spots, tissue_spots], 1, min))
    #unique(spots_distances)
    tissue_bordering_spots <- names(spots_distances)[which(spots_distances < distance)]
    #add tissue_bordering factor to meta.data
    bordering_selection <- all_barcodes %in% tissue_bordering_spots
    spatial@meta.data[bordering_selection,]$spot <- "tissue_bordering"
    Seurat::Idents(spatial) <- "spot"
    #subset Seurat object based on Idents (to only include tissue and tissue_bordering spots)
    spatial <- subset(spatial, idents=c("tissue_only", "tissue_bordering"))
    #add microbiome assay
    selected_barcodes <- c(tissue_spots, tissue_bordering_spots)
    #subsetting matrix with barcodes present in GEX matrix
    valid_barcodes <- taxid_matrix@Dimnames[[2]] %in% selected_barcodes
    taxid_matrix <- taxid_matrix[,valid_barcodes]
    #remove taxids with only zeros
    taxid_matrix <- taxid_matrix[Matrix::rowSums(taxid_matrix)>0,]
    # check if any taxa are left
    if (length(taxid_matrix@Dimnames[[1]])==0){
      result <- NULL
      saveRDS(result, file=paste0(outDir, "/", sampleName, "_Seurat_object.RDS"))
      print("No taxa left after decontamination. Returning NULL object.")
      return(result)
    }
    #print(paste0(dim(taxid_matrix)[1], " taxa remaining"))
    #print(paste0(sum(taxid_matrix), " counts remaining"))
    # create a new assay to store ADT information
    microbial_assay <- Seurat::CreateAssayObject(counts = taxid_matrix)
    # add this assay to the previously created Seurat object
    suppressWarnings(spatial[["TAXA_G"]] <- microbial_assay)
    suppressWarnings(Seurat::DefaultAssay(spatial) <- "TAXA_G")
  }


  if (spatiallyVariable==TRUE){
    ### DECONTAMINATION STEP#5: Removing taxa that are not spatially variable (based on Moran's I spatial autocorrelation)
    cat("\nDecontamination step: removing taxa that are not spatially variable")
    spatial <- Seurat::FindSpatiallyVariableFeatures(object=spatial, assay="TAXA_G", slot="counts",selection.method="moransi")
    #filter out taxa based on spatial variability analysis
    moransi_analysis <- spatial@assays[["TAXA_G"]]@meta.features %>% dplyr::filter(MoransI_observed>=0.01, MoransI_p.value<0.001) %>% rownames()
    # create a new assay
    spatially_variable_taxa <- taxid_matrix@Dimnames[[1]] %in% moransi_analysis
    taxid_matrix <- taxid_matrix[spatially_variable_taxa,]
    if (sum(spatially_variable_taxa)==0){
      result <- NULL
      saveRDS(result, file=paste0(outDir, "/", sampleName, "_Seurat_object.RDS"))
      print("No taxa left after decontamination. Returning NULL object.")
      return(result)
    }
    microbial_assay <- Seurat::CreateAssayObject(counts = taxid_matrix)
    # add this assay to the previously created Seurat object
    spatial[["TAXA_G"]] <- microbial_assay
    Seurat::DefaultAssay(spatial) <- "TAXA_G"
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
    result <- list("seurat_object"=spatial, "taxid_counts"=dataframe, "matrix"=taxid_matrix)
    saveRDS(result, file=paste0(outDir, "/", sampleName, "_Seurat_object.RDS"))
    return(result)
  }
}
