#' @rdname computeSizeFactors
#'
#' @title  Estimate size factors
#' @description Size factor estimation provided by DEseq
#' @details A comment about size factor estimation by DEseq based on Simon
#'     Anders (Heidelberg, Germany) work to estimate the library size, simply
#'     taking the total number of (mapped or unmapped) reads is, in our
#'     experience, not a good idea. Sometimes, a few very strongly expressed
#'     genes are differentially expressed, and as they make up a good part of
#'     the total counts, they skew this number. After you divide by total
#'     counts, these few strongly expressed genes become equal, and the whole
#'     rest looks differentially expressed.
#'
#'     The following simple alternative works much better: - Construct a
#'     "reference sample" by taking, for each gene, the geometric mean of the
#'     counts in all samples. - To get the sequencing depth of a sample
#'     relative to the reference, calculate for each gene the quotient of the
#'     counts in your sample divided by the counts of the reference sample. Now
#'     you have, for each gene, an estimate of the depth ratio. - Simply take
#'     the median of all the quotients to get the relative depth of the
#'     library. This is what the 'estimateSizeFactors' function of our DESeq
#'     package does.
#'
#' @param dss DESeqDataSet object
#' @return Estimation of the size factors
#'
#' @references http://seqanswers.com/forums/showpost.php?p=16468&postcount=13
#'
#' @examples
#'     dds <- DESeq2::makeExampleDESeqDataSet(n = 1000, m = 4)
#'     MethylIT:::.computeSizeFactors(dds)
#' @keywords internal
.computeSizeFactors <- function(dss) {
   y <- counts(dss)
   ## cn = colnames( y )
   group <- dss@colData$condition
   lev <- levels(dss@colData$condition)

   x1 <- y[ ,group == lev[1]]
   r1 <- apply(x1, 2, function(v) mean(v, na.rm=TRUE))
   r1 <- r1 / apply(x1, 2, .gmMean)

   x2 <- y[ ,group == lev[2]]
   r2 <- apply(x2, 2, function(v) mean(v, na.rm=TRUE))
   r2 <- r2/apply(x2, 2, .gmMean)
   c(r1, r2)
}
