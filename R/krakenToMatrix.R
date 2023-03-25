#' Creates a taxid-spot matrix from kraken output file.
#'
#' @param filePath A string. Path to input file (kraken output file from pipeline: SAMPLE_kraken_output.txt).
#' @param counts A character string. Either "umi_counts" or "read_counts". Defaults to "umi_counts".
#' @param tax_level A character string. Either "species", "genus" or "family". It is recommended to use "genus". Defaults to "genus".
#' @param whitelist A character string. Either "10XVisium" (when working with spaceranger output data) or "10XChromium_v3.1" (when working with cellranger and Single Cell 3' v3.1 technology data).
#' @param taxonomizrDB A string. Path to nameNode.sqlite database required for taxonomic conversions (see README file for how to download).
#' @param outDir A character string. Path to desired output directory.
#'
#' @return A list containing $taxid_counts (dataframe with counts, taxonomic names and taxids) and $matrix (matrix with rows=taxids, cols=barcodes).
#' @export
#'
#' @examples
#' #for spatial transcriptomics
#' \dontrun{
#'  krakenToMatrix(filePath=
#'  system.file("extdata", "CRC_16", "CRC_16_kraken_output.txt.gz", package="microbiome10XVisium"),
#'  outDir=system.file("extdata", "CRC_16/", package="microbiome10XVisium"))
#' }
#' #for scRNAseq
#' \dontrun{
#'  krakenToMatrix(filePath="FILEPATH", outDir="OUTDIR", whitelist="10XChromium_v3.1")
#' }
krakenToMatrix <- function(filePath, counts="umi_counts", tax_level="genus", whitelist="10XVisium", taxonomizrDB="/projects/site/pred/microbiome/database/taxonomizr_DB/nameNode.sqlite",outDir){

  ##check parameters
  if (is.null(counts)){
    counts <- "umi_counts"
    warning("counts set to default: umi_counts")
  } else if (!(counts %in% c("read_counts", "umi_counts"))){
    stop("counts has to be either read_counts or umi_counts")
  }

  if (is.null(tax_level)){
    tax_level <- "genus"
    warning("tax_level set to default: genus")
  } else if (!(tax_level %in% c("species", "genus", "family"))){
    stop("tax_level has to be either species, genus or family")
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

  #load the corresponding whitelist file
  if (whitelist=="10XVisium"){
    whitelist.filepath <- system.file("extdata", "whitelist.tsv", package="microbiome10XVisium", mustWork=TRUE)
    whitelist <- utils::read.table(file = whitelist.filepath, header=FALSE)
    colnames(whitelist) <- c("barcode")
  } else if (whitelist=="10XChromium_v3.1") {
    whitelist.filepath <- system.file("extdata", "whitelist_cellranger.tsv.gz", package="microbiome10XVisium", mustWork=TRUE)
    whitelist <- utils::read.table(file = whitelist.filepath, header=FALSE)
    colnames(whitelist) <- c("barcode")
  }

  ## import modified kraken2 output file
  data <- utils::read.table(file=filePath, sep="\t", header=FALSE)
  colnames(data) <- c("barcode", "umi", "taxid")
  #trim whitespaces
  data$barcode <- trimws(data$barcode)
  data$umi <- trimws(data$umi)
  data$taxid <- trimws(data$taxid)

  ## extract valid barcodes from whitelist file (containing the 4992 true spatial barcodes)
  data <- subset(data, data$barcode %in% whitelist$barcode)
  data <- tidyr::as_tibble(data)

  ##remove taxid 131567: cellular organisms and 1:root and 2:Bacteria and 2759:Eukaryota
  data <- data %>% dplyr::filter(!taxid %in% c(131567, 1, 2, 2759))

  if (counts=="read_counts"){

    data <- data %>% dplyr::select(barcode, taxid)

    ## reclassify all taxids to tax_level and remove those taxids at taxonomic level > tax_level
    suppressWarnings(at_tax_level <- taxonomizr::getTaxonomy(ids=data$taxid,sqlFile=taxonomizrDB, desiredTaxa=c(tax_level)))
    data <- data[!is.na(at_tax_level),]
    ## assign taxid at tax_level to all reads
    suppressWarnings(taxa <- taxonomizr::getTaxonomy(ids=data$taxid,sqlFile=taxonomizrDB, desiredTaxa=c(tax_level))[,1])
    suppressWarnings(data$taxid <- taxonomizr::getId(taxa=taxa,sqlFile=taxonomizrDB))
    #select first taxid in case of multiple taxids (first is usually bacterial)
    data$taxid <- gsub("\\,.*","",data$taxid)

    ## collapse all reads with same barcode & taxid --> read counts
    data <- plyr::ddply(data,c("barcode", "taxid"),nrow)
    colnames(data) <- c("barcode", "taxid", "counts")

  } else if (counts=="umi_counts"){

    ## first collapse all reads with identical barcode, umi, taxid
    data <- plyr::ddply(data,c("barcode", "umi", "taxid"),nrow)
    data <- data %>% dplyr::select(barcode, umi, taxid)

    ## reclassify all taxids to tax_level and remove those taxids at taxonomic level > tax_level
    suppressWarnings(at_tax_level <- taxonomizr::getTaxonomy(ids=data$taxid,sqlFile=taxonomizrDB, desiredTaxa=c(tax_level)))
    data <- data[!is.na(at_tax_level),]
    ## assign taxid at tax_level to all reads
    suppressWarnings(taxa <- taxonomizr::getTaxonomy(ids=data$taxid,sqlFile=taxonomizrDB, desiredTaxa=c(tax_level))[,1])
    suppressWarnings(data$taxid <- taxonomizr::getId(taxa=taxa,sqlFile=taxonomizrDB))
    #select first taxid in case of multiple taxids (first is usually bacterial)
    data$taxid <- gsub("\\,.*","",data$taxid)

    ## transform tibble into a matrix with feature=taxid as row and barcode as column
    #drop umi and count column
    data <- plyr::ddply(data,c("barcode", "taxid"),nrow)
    colnames(data) <- c("barcode", "taxid", "counts")
  }

  ## create table with counts, taxids and taxonomic names
  taxa <- data %>% dplyr::select(taxid, counts)
  taxa_abundance <- stats::aggregate(counts ~ taxid, data=taxa, FUN=sum) %>% dplyr::arrange(dplyr::desc(counts))
  print(paste0(sum(taxa_abundance$counts), " total microbial read counts at ", tax_level, " level"))
  suppressWarnings(taxa_names <- taxonomizr::getTaxonomy(ids=taxa_abundance$taxid, taxonomizrDB, desiredTaxa = c("superkingdom", "phylum", "genus")))
  taxa_abundance <- cbind(taxa_abundance, taxa_names)
  col_order <- c("counts", "superkingdom", "phylum", "genus", "taxid")
  taxa_abundance <- taxa_abundance[,col_order]
  rownames(taxa_abundance) <- NULL

  ## convert dataframe to sparse matrix
  data$counts <- as.numeric(data$counts)
  data$taxid <- as.character(data$taxid)
  matrix <- tidyr::pivot_wider(data, names_from=taxid, values_from = counts, values_fill=0)
  #complete barcodes from whitelist with 0s
  matrix <- dplyr::left_join(whitelist, matrix)
  matrix[is.na(matrix)] <- 0
  #add -1 to barcode column
  matrix$barcode <- paste0(matrix$barcode, "-1")
  matrix <- as.matrix(matrix)
  rownames(matrix) <- NULL
  rownames(matrix) <- matrix[,1]
  matrix <- matrix[,2:ncol(matrix)]
  #transpose matrix and convert matrix to dgCMatrix class
  matrix <- t(matrix)
  colnames <- colnames(matrix)
  rownames <- rownames(matrix)
  matrix <- Matrix::Matrix(as.numeric(matrix), nrow=length(rownames), ncol=length(colnames), sparse=TRUE, dimnames=list(rownames, colnames))

  ## return
  result <- list("taxid_counts"=taxa_abundance, "matrix"=matrix)
  saveRDS(result, file=paste0(outDir, "/", tax_level, "_", counts, ".RDS"))
  print(paste0("output saved as RDS file at ", outDir, "/", tax_level, "_", counts, ".RDS" ))
  return(result)

}

