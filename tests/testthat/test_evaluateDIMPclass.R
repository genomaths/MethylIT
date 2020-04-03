library(testthat)
library(MethylIT)

context("evaluateDIMPclass tests")

test_that("evaluateDIMPclass tests", {
    data(dmps, package = "MethylIT")
    perf <- evaluateDIMPclass(LR = dmps,
                            column = c(hdiv = TRUE, TV = TRUE,
                                       wprob = TRUE, pos = TRUE),
                            classifier = 'lda', n.pc = 4L,
                            control.names =  c('C1', 'C2', 'C3'),
                            treatment.names = c('T1', 'T2', 'T3'),
                            center = TRUE, scale = TRUE, prop = 0.6)
    expect_equal(perf$FDR, 0)
    expect_true(perf$Performance$overall[1] == 1)
})
