% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/uniqueGRfilterByCov.R
\name{uniqueGRfilterByCov}
\alias{uniqueGRfilterByCov}
\title{Unique \code{\link[GenomicRanges]{GRanges-class}} of methylation
read counts filtered by coverage.}
\usage{
uniqueGRfilterByCov(
  x,
  y = NULL,
  min.coverage = 4,
  and.min.cov = TRUE,
  min.meth = 0,
  min.umeth = 0,
  min.sitecov = 4,
  percentile = 0.9999,
  min.percentile = TRUE,
  high.coverage = NULL,
  columns = c(mC = 1, uC = 2),
  ignore.strand = FALSE,
  y.centroid = NULL,
  num.cores = 1L,
  tasks = 0L,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{x}{An object from the classes 'GRanges', 'InfDiv', or 'pDMP' with
methylated and unmethylated counts in its meta-column. If the argument 'y'
is not given, then it is assumed that the first four columns of the
\code{\link[GenomicRanges]{GRanges-class}} metadata 'x' are counts:
methylated and unmethylated counts for samples '1' and '2'. If \eqn{x} is a
'InfDiv', or 'pDMP', then \eqn{'y'} is not needed, since samples '1' and '2'
are the first four columns of these objects. In the context of MethylIT
methylation analysis, \eqn{x} is considered a control sample or, more
appropriated (if available) the centroide of control group.}

\item{y}{A \code{\link[GenomicRanges]{GRanges-class}} or a
\code{\link[GenomicRanges]{GRangesList-class}} object with methylated and
unmethylated counts in its meta-column. Default is NULL. If 'y' is a
\code{\link[GenomicRanges]{GRangesList-class}} object, then the first two
metacolumns from 'x' are used in the pair-wise filtering toguether with the
first two metacolumns from each \code{\link[GenomicRanges]{GRanges-class}}
element from 'y'.}

\item{min.coverage}{An integer or an integer vector of length 2. If
'min.coverage' is an integer vector, then the corresponding min coverage
is applied to each sample. Default is 4, i.e., \eqn{min.coverage = c(4,4)}.}

\item{and.min.cov}{Logical. Whether to apply the logical AND to select the
cytosine sites based on \eqn{min.coverage}. If FALSE, then a logical OR is
applied, and cytosine sites where at least one sample hold the
\eqn{min.coverage} are preserved. Default is TRUE.}

\item{min.meth}{An integer or an integer vector of length 2. Cytosine sites
where the numbers of read counts of methylated cytosine in both samples, '1'
and '2', are less than 'min.meth' are discarded. If 'min.meth' is an integer
vector, then the corresponding min number of reads is applied to each
sample. That is, \eqn{min.meth} implement a logical OR.}

\item{min.umeth}{An integer or an integer vector of length 2. Minimum number
of reads to consider cytosine position. Specifically cytosine positions
where \eqn{uC \leq min.umeth} are discarded. Default is \eqn{min.umeth = 0}.}

\item{min.sitecov}{An integer. The minimum total coverage. Only sites where
the total coverage \eqn{cov1 + cov2} is greater than 'min.sitecov' are
considered for downstream analysis, where cov1 and cov2 are the coverages
for samples 1 and 2, respectively.}

\item{percentile}{Threshold to remove the outliers (PCR bias) from each file
and all files stacked. If 'high.coverage = NULL', then the threshold \eqn{q}
will be computed as: \deqn{q1 = quantile(cov1, probs=percentile)}
\deqn{q2 = quantile(cov2, probs=percentile)}.

where \eqn{cov1} and \eqn{cov2} are the coverage vectors from samples 1
and 2, respectively.}

\item{min.percentile}{Logical. Each sample yield a percentile value. The
user must decide whether to use the minimum or the maximum value from these
percentile values. Default is TRUE. Hence, \eqn{min.percentile = FALSE}
implies to use the maximum value.}

\item{high.coverage}{An integer for read counts. Cytosine sites having
higher coverage than this are discarded. Default is NULL.}

\item{columns}{Vector of integer numbers of the columns (from each GRanges
meta-column) where the methylated and unmethylated counts are provided. If
not provided, then the methylated and unmethylated counts are assumed to be
at columns 1 and 2, respectively.}

\item{ignore.strand}{When set to TRUE, the strand information is ignored in
the overlap of \code{\link[GenomicRanges]{GRanges-class}} objects. This is a
parameter passed to \code{\link{uniqueGRanges}} function. Default value:
FALSE.}

\item{y.centroid}{Optional. A \code{\link[GenomicRanges]{GRanges-class}}
object corresponding to the treatment/individual centroide. This information
is applied to reduce the bias originated by missing cytosine sites. The
centroide can be compute applying function \code{\link{poolFromGRlist}}.}

\item{num.cores, tasks}{Parameters for parallel computation using package
\code{\link[BiocParallel]{BiocParallel-package}}: the number of cores to
use, i.e. at most how many child processes will be run simultaneously (see
\code{\link[BiocParallel]{bplapply}} and the number of tasks per job (only
for Linux OS).}

\item{verbose}{if TRUE, prints the function log to stdout}

\item{...}{Additional parameters for \code{\link{uniqueGRanges}} function.}
}
\value{
if 'x' and 'y' are \code{\link[GenomicRanges]{GRanges-class}}
object, then a \code{\link[GenomicRanges]{GRanges-class}} with the columns
of methylated and unmethylated counts filtered for each cytosine position.

A \code{\link[GenomicRanges]{GRangesList-class}} object will returned, if
'y' is a \code{\link[GenomicRanges]{GRangesList-class}} object of same
length as length(y) and named as names(y).
}
\description{
Given two \code{\link[GenomicRanges]{GRanges-class}} objects,
samples '1' and '2', carrying the counts of methylated (mC) and unmethylated
(uC) cytosines in their metacolumns, this function will filter by coverage
each cytosine site from each \code{\link[GenomicRanges]{GRanges-class}}
object.
}
\details{
Cytosine sites with 'coverage' > 'min.coverage' in at least one of
the samples are preserved. Positions with 'coverage' < 'min.coverage' in
both samples, 'x' and 'y', are removed. Positions with 'coverage' >
'percentile' (e.g., 99.9 percentile) are removed as well. It is expected
that the columns of methylated and unmethylated counts are given.

This function is addressed to create pair-wise
\code{\link[GenomicRanges]{GRanges-class}} object with four metacolumns of
count: samples '1' and '2', carrying the counts of methylated (mC) and
unmethylated (uC) cytosines in their metacolumns, respectively. Counts from
sample 1 are typically used as reference counts in computing information
divergences in the downtstream analysis.

The cut-off value to remove PCR bias is computed as:

If \strong{high.coverage} is NULL, then the cut-off point is:
\eqn{q = min(q1, q2)} (if \eqn{min.percentile} = TRUE) or
\eqn{q = max(q1, q2)}. If \strong{high.coverage} is not NULL, then
\deqn{q = max(q, high.coverage)}.

Another source of bias is originated by missing cytosine sites. Missing data
are frequently found in experimental data sets and, in particular, in
bisulfite genomic sequencing data. Typically, in statistical analyses, the
bias originated by missing data (for given variable) is mitigated by using
the mean of the known values for the corresponding variable. That is, in
present case, if the reads for some cytosine site are missed in a sample
from a set of, e.g., three individuals, then the means of reads (methylated
and unmethylated) for such site are applied as an estimation of the best
expected ("guessed") value of missed reads. Obviously, if the reads are
missed in all the samples, then the site is discarded (see examples).

The treatment centroide can be compute applying function
\code{\link{poolFromGRlist}}. Also notice that, since the centroide
correction is only available for the treat group, it is assumed that sample
\eqn{x} carries reads for each (or almost all) cytosine sites are provide.
}
\examples{
## Create new data. It is assumed that sample 'x' carries reads
## for each cytosine sites are provide.
strands <- c("+","-","+","-", "+","-","+","+","+","+","+")
pos <- c(10,11,11,12,13,13,14,15,16,17,18)

x <- data.frame(chr = 'chr1', start = pos, end = pos,
                mC = c(2,3,2,5,10,7,9,11,4,10,7),
                uC = c(2,30,20,4,8,0,10,3,0,8,1),
                strand = strands)

x <- makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
x

## sample y
y <- data.frame(chr = 'chr1', start = 11:18, end = 11:18,
                mC2 = c(4,1,2,1,4,5:7), uC2 = c(0,0,2:7),
                strand = c("+","-","-","+","+","+","+","+"))

y <- makeGRangesFromDataFrame(y, keep.extra.columns = TRUE)
y

## The default settings. Sites where one of the samples has zero methylation
## calling or min.coverage is lesser than 4 in at least one sample are
## discarded. This setting implies a drastic decision and the amount of
## cytosine removed can be lead to strong biased conclusions in the
## downstream analysis.
uniqueGRfilterByCov(x = x, y = y,
                    percentile = 1,
                    ignore.strand = FALSE)

## Setting 'and.min.cov = FALSE' undesired cytosine sites are preseved. For
## example, meaningless situations with  methylation levels
## p = 1/(1 + 0) = 1.
uniqueGRfilterByCov(x = x, y = y,
                    and.min.cov = FALSE,
                    ignore.strand = FALSE,
                    percentile = 1,
                    verbose = FALSE)

## Setting 'min.coverage = 8' does not solves the previous issue, but still
## it preserves cytosine sites with one of the samples with zero or small
## coverage:
uniqueGRfilterByCov(x = x, y = y,
                    and.min.cov = FALSE,
                    min.coverage = 8,
                    percentile = 1,
                    ignore.strand = FALSE,
                    verbose = FALSE)

## A centroid, a vector of means of methylation read for each cytosine site
## from the treatment group can be used as the best estimation to replace
## missing data: 'mC=0' and 'uC=0' or low coverage sites 'mC=1' and 'uC=0'.
y_centroid <- data.frame(chr = 'chr1',
                start = pos, end = pos,
                mC2 = c(8,7,6,7,5,8,1:5), uC2 = 0:10,
                strand = c("+","-","+","-", "+","-","+","+","+","+","+"))
y_centroid <- makeGRangesFromDataFrame(y_centroid,
                                      keep.extra.columns = TRUE)
y_centroid

## The cytosine sites with missing data or low coverage will be still,
## included, using the centroid of the centroinde of th sample group to
## which 'y' belong to.
uniqueGRfilterByCov(x = x, y = y,
                    and.min.cov = FALSE,
                    min.coverage = c(1, 8),
                    ignore.strand = FALSE,
                    y.centroid = y_centroid,
                    min.percentile = FALSE,
                    percentile = 1,
                    verbose = FALSE)
}
