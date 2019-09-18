library(testthat)
library(GenomicRanges)
library(MethylIT)

context("pcaLogisticR tests")

test_that("pcaLogisticR dummy test", {
  data(iris)
  data = iris[ iris$Species != "virginica", ]
  data$Species <- droplevels(data$Species)
  formula = Species ~ Petal.Length + Sepal.Length + Petal.Width
  pca.logistic = pcaLogisticR(formula = formula,
                              data = data, n.pc = 2, scale = TRUE,
                              center = TRUE, max.pc = 2)
  set.seed(123)
  newdata = iris[sample.int(150, 40), 1:4]
  newdata.prediction = predict(pca.logistic, newdata, type = "all")


  # The accuracy should be > 0.5
  expect_true(all(table(newdata.prediction$class) == c(12,28)))
})
