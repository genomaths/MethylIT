library(testthat)
library(GenomicRanges)
library(MethylIT)
context("MethylIT tests")

test_that(".estimateGLM test", {
  data(ds, package = "MethylIT")

  X <- ds$counts[69,]
  baseMeanAndVar <- data.frame(baseMean = mean(X),
                               baseVar = var(X))

  y <- MethylIT:::.estimateGLM(x = X, groups = ds$colData$condition,
                              baseMV = baseMeanAndVar,
                              w = c(1,1), MVrate = 0.95,
                              test = "LRT")
  expect_true(y$pvalue < 0.05)
  expect_true(y$model == "Neg.Binomial")
})
