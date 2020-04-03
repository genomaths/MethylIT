library(testthat)
library(GenomicRanges)
library(MethylIT)

context("countTest2 tests")

test_that("countTest2 test", {
    data(ds, package = "MethylIT")
    y <- countTest2(DS = ds, num.cores = 1L, maxGrpCV = c(0.6, 0.5),
                    verbose = FALSE)
    expect_true(all(y$adj.pval < 0.1))
})
