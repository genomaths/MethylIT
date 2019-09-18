library(testthat)
library(MethylIT)
context("MethylIT fitGGammaDist tests")

test_that("fitGGammaDist function test", {
    set.seed(123)
    x <- rggamma(1000, alpha = 1.03, psi = 0.75, scale = 1.1)
    fit <- fitGGammaDist(x)
    expect_is(fit, "data.frame")
    expect_is(fit$Estimate, "numeric")
    TRUE
    })
