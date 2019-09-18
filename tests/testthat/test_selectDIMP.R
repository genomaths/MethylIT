library(testthat)
library(GenomicRanges)
library(MethylIT)

context("selectDIMP tests")

test_that("selectDIMP test", {
   data(PS, cutpoint)
   DMPs <- selectDIMP(PS, div.col = 9L, cutpoint = cutpoint$cutpoint)
   expect_true(length(DMPs$T1) <= length(PS$T1))
})

