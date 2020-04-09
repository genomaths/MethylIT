library(testthat)
library(GenomicRanges)
library(MethylIT)
context("MethylIT AICmodel tests")

test_that("AICmodel function test", {
    set.seed(77)
    x = runif(100, 1, 5)
    y = 2 * exp(-0.5 * x) + runif(100, 0, 0.1)
    nlm <- nls(Y ~ a * exp( b * X), data = data.frame(X = x, Y = y),
                start = list( a = 1.5, b = -0.7),
                control = nls.control(maxiter = 10^4, tol = 1e-05),
                algorithm = "port")

    expect_equal(AICmodel(nlm), AIC(nlm))
    # Now, using residuals from the fitted model:
    res = y - coef(nlm)[1] * exp(coef(nlm)[2] * x)
    expect_equal(AICmodel(residuals = res, np = 2), AIC(nlm))

    answer <- AICmodel(nlm)
    expect_is(answer, "numeric")
})
