library(testthat)
library(GenomicRanges)
library(MethylIT)

context("uniqueGRfilterByCov tests")

test_that("uniqueGRfilterByCov with overlapped rows", {
   dfChr1 <- data.frame(chr = "chr1", start = 11:15, end = 11:15,
                       strand = c("+","-","+","*","."), mC = 1:5, uC = 1:5)
   dfChr2 <- data.frame(chr = "chr1", start = 12:18, end = 12:18,
                       strand = '*', mC = 1:7, uC = 1:7)
   gr1 <- makeGRangesFromDataFrame(dfChr1, keep.extra.columns = TRUE)
   gr2 <- makeGRangesFromDataFrame(dfChr2, keep.extra.columns = TRUE)

   r1 <- uniqueGRfilterByCov(gr1, gr2, min.coverage = 4,  high.coverage = 6,
                           percentile = 1, ignore.strand = TRUE,
                           verbose = FALSE)
   r2 <- data.frame(chr = "chr1", start = 12:18, end = 12:18, strand = '*',
                   mC = c(2:5, 0, 0, 0), uC = c(2:5, 0, 0, 0),
                   mC.1 = 1:7, uC = 1:7)
   r2 <- makeGRangesFromDataFrame(r2, keep.extra.columns = TRUE)
   expect_equivalent(r1, r2)
})

