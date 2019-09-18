library(testthat)
library(MethylIT)
context("MethylIT getGEOSuppFiles tests")

test_that("getGEOSuppFiles function test", {
    filenames = getGEOSuppFiles(GEO = "GSM881757",pattern = "G_cytosine.txt.gz")
    LR <- readCounts2GRangesList(filenames = filenames, remove = TRUE,
                                sample.id = c("drm2_CG", "drm2_CHG"),
                                columns = c(seqnames = 1, start = 2, mC = 4,
                                        uC = 3),
                                pattern = "^Chr1", verbose = TRUE)
    file.remove(filenames)
    expect_true(length(LR[[2]]) > 1000)
})
