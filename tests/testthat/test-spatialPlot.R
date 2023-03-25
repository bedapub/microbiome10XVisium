test_that("throws error when no viral reads detected", {
  expect_error(spatialPlot(CRC_16, taxa="viral"))
})

test_that("throws error when taxa misspelled", {
  expect_error(spatialPlot(CRC_16, taxa="Fusobacterim"))
})

test_that("throws error when taxa not detected in sample", {
  expect_error(spatialPlot(CRC_16, taxa="Pseudomonas"))
})
