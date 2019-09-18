library(testthat)
library(GenomicRanges)
library(MethylIT)
library(rtracklayer)


context("getDIMPatGenes tests")

test_that("getDIMPatGenes test", {
   genes <- GRanges(seqnames = "1",
                   ranges = IRanges(start = c(3631, 6788, 11649),
                                   end = c(5899, 9130, 13714)),
                   strand = c("+", "-", "-"))
   mcols(genes) <- data.frame(gene_id = c("AT1G01010", "AT1G01020",
                                         "AT1G01030"))

   data(PS, cutpoint)
   DIMPs <- selectDIMP(PS, div.col = 9L, cutpoint = cutpoint$cutpoint)

   DIMR <- getDIMPatGenes(GR = DIMPs$T1, GENES = genes)
   expect_true(length(DIMR) > 0)
})


