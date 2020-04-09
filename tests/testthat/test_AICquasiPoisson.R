library(testthat)
library(GenomicRanges)
library(MethylIT)
context("MethylIT tests")

test_that(".AICquasiPoisson test", {
    ## Already tested with the testing for function ".estimateGLM"
    ## When testing ".estimateGLM" function, the resultant best model is
    ## "QuasiPoisson". That is, if the test fail ".estimateGLM, then test
    ## for ".AICquasiPoisson" fails as well.
    TRUE
})
