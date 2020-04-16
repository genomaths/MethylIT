library(testthat)
library(MethylIT)

context("dmpClusters tests")

test_that("dmpClusters test", {
   data(dmps)
   x1 = dmpClusters(GR = dmps, maxDist = 6, minNumDMPs = 10,
                    num.cores=2L, verbose = FALSE)
   expect_true(length(x1) == 13)
})


