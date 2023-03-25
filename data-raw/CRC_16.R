## code to prepare `CRC_16` dataset goes here

CRC_16 <- readRDS(file=system.file("extdata", "CRC_16", "CRC_16_Seurat_object.RDS", package="microbiome10XVisium"))
usethis::use_data(CRC_16, overwrite = TRUE)
