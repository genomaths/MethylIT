library(testthat)
library(GenomicRanges)
library(MethylIT)

context("pcaLDA tests")

test_that("pcaLDA test", {
  data(iris)
  ld = pcaLDA(formula = Species ~., data = iris, n.pc = 1, max.pc = 2,
              scale = TRUE, center = TRUE)
  set.seed(123)
  idx = sample.int(150, 40)
  newdata = iris[idx, 1:4]
  newdata.prediction = predict(ld, newdata = newdata)

  # The accuracy should be > 0.5
  expect_true(all(newdata.prediction$class[1:2] == "setosa"))
})
