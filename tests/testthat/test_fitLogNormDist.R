library(testthat)
library(MethylIT)
context("MethylIT fitLogNormDist tests")

test_that("fitLogNormDist function test", {
    set.seed(126)
    x <- rlnorm(1000, meanlog = 1.03, sdlog = 2.1)
    y <- fitLogNormDist(x)
    expect_true(as.numeric(y$Adj.R.Square[1]) > 0.99)
})
