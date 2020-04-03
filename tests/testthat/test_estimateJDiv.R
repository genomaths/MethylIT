library(testthat)
library(GenomicRanges)
library(MethylIT)

context("estimateJDiv tests")

test_that("estimateJDiv using n", {
    jd <- estimateJDiv(p = c(0.5, 0.5))
    expect_true(jd == 0)
})
