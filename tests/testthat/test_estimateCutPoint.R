library(testthat)
library(GenomicRanges)
library(MethylIT)

context("estimateCutPoint tests")

test_that("estimateCutPoint tests", {
   data(PS)
   cutpoint = estimateCutPoint(LR = PS, simple = FALSE,
                               column = c(hdiv = TRUE, TV = TRUE,
                                           wprob = TRUE, pos = TRUE),
                               classifier1 = "qda",
                               control.names = c("C1", "C2", "C3"),
                               treatment.names = c("T1", "T2", "T3"),
                               tv.cut = 0.92, clas.perf = TRUE, prop = 0.6,
                               div.col = 9L)
  expect_true(round(cutpoint$cutpoint) == 179)
})
