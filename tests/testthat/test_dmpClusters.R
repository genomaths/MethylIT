library(testthat)
library(GenomicRanges)
library(MethylIT)
library(rtracklayer)


context("dmpClusters tests")

test_that("dmpClusters test", {
   data(PS, cutpoint)
   dmps <- selectDIMP(PS, div.col = 9L, cutpoint = cutpoint$cutpoint)
   x1 = dmpClusters(GR = dmps, maxDist = 6, minNumDMPs = 10, num.cores=2L,
                    verbose = FALSE)
   expect_true(length(x1) == 3)
})


