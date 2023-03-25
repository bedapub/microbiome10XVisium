test_that("throws error when no more taxa left in decontamination step 2 (removeLikelyContaminants)", {
  testthat::expect_equal(decontaminate(sampleName="SN124_A938797_Rep2",
                                       filePath=system.file("extdata", "testing", "genus_umi_counts.RDS", package="microbiome10XVisium"),
                                       spacerangerDir = system.file("extdata", "testing", "outs", package="microbiome10XVisium"),
                                       outDir=system.file("extdata", "testing", package="microbiome10XVisium"),
                                       spots="all", removeSingletons=TRUE, removeLikelyContaminants=TRUE,
                                       removeSpecificTaxa = c("Mycobacterium","Malassezia","Cutibacterium","Anaerostipes","Prevotella","Blautia",
                                                              "Staphylococcus","Fusobacterium","Roseburia","Bacteroides","Acinetobacter","Streptococcus",
                                                              "Mediterraneibacter","Finegoldia","Dysosmobacter","Corynebacterium","Faecalibacillus","Coprococcus",
                                                              "Clostridium","Vescimonas","Lactobacillus","Phocaeicola","Cryptococcus","Ruminococcus",
                                                              "Leishmania","Toxoplasma","Kocuria","Rothia","Aspergillus","Peptoniphilus",
                                                              "Micrococcus","Megasphaera","Hymenobacter")), NULL)
})

test_that("returns NULL when no more taxa left in decontamination step 3 (spots)", {
  testthat::skip_if_not_installed("Matrix", minimum_version = "1.5.0")
  testthat::expect_equal(decontaminate(sampleName="SN124_A938797_Rep2",
                                       filePath=system.file("extdata", "testing", "genus_umi_counts.RDS", package="microbiome10XVisium"),
                                       spacerangerDir = system.file("extdata", "testing", "outs", package="microbiome10XVisium"),
                                       outDir=system.file("extdata", "testing", package="microbiome10XVisium"),
                                       spots="tissueOnly", removeSpecificTaxa = c("Mycobacterium","Malassezia", "Blautia")), NULL)
})

test_that("throws error when no more taxa left in decontamination step 4 (spatiallyVariable)", {
  testthat::skip_if_not_installed("Matrix", minimum_version = "1.5.0")
  testthat::expect_equal(decontaminate(sampleName="SN124_A938797_Rep2",
                                       filePath=system.file("extdata", "testing", "genus_umi_counts.RDS", package="microbiome10XVisium"),
                                       spacerangerDir = system.file("extdata", "testing", "outs", package="microbiome10XVisium"),
                                       outDir=system.file("extdata", "testing", package="microbiome10XVisium"),
                                       spots="tissuePlusBordering", removeSpecificTaxa = c("Mycobacterium","Malassezia", "Roseburia", "Anaerostipes", "Blautia"),
                                       spatiallyVariable=TRUE),NULL)
})

