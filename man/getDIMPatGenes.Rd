% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getDIMPatGenes.R
\name{getDIMPatGenes}
\alias{getDIMPatGenes}
\alias{getDIMPatGenes.default}
\alias{getDIMPatGenes.GRanges}
\alias{getDIMPatGenes.pDMP}
\alias{getDIMPatGenes.InfDiv}
\alias{getDIMPatGenes.list}
\title{Count DMPs at gene-body}
\usage{
getDIMPatGenes(
  GR,
  GENES,
  type = "within",
  ignore.strand = TRUE,
  only.hypo = FALSE,
  only.hyper = FALSE,
  output = c("list", "GRanges"),
  by.coord = FALSE,
  gene_id_col = NULL,
  gene_name_col = NULL,
  ...
)

\method{getDIMPatGenes}{default}(
  GR,
  GENES,
  type = "within",
  ignore.strand = TRUE,
  only.hypo = FALSE,
  only.hyper = FALSE,
  output = NULL,
  by.coord = FALSE,
  gene_id_col = NULL,
  gene_name_col = NULL,
  ...
)

\method{getDIMPatGenes}{GRanges}(
  GR,
  GENES,
  type = "within",
  ignore.strand = TRUE,
  only.hypo = FALSE,
  only.hyper = FALSE,
  output = NULL,
  by.coord = FALSE,
  gene_id_col = NULL,
  gene_name_col = NULL,
  ...
)

\method{getDIMPatGenes}{pDMP}(
  GR,
  GENES,
  type = "within",
  ignore.strand = TRUE,
  only.hypo = FALSE,
  only.hyper = FALSE,
  output = c("list", "GRanges"),
  by.coord = FALSE,
  gene_id_col = NULL,
  gene_name_col = NULL,
  ...
)

\method{getDIMPatGenes}{InfDiv}(
  GR,
  GENES,
  type = "within",
  ignore.strand = TRUE,
  only.hypo = FALSE,
  only.hyper = FALSE,
  output = c("list", "GRanges"),
  by.coord = FALSE,
  gene_id_col = NULL,
  gene_name_col = NULL,
  ...
)

\method{getDIMPatGenes}{list}(
  GR,
  GENES,
  type = "within",
  ignore.strand = TRUE,
  only.hypo = FALSE,
  only.hyper = FALSE,
  output = c("list", "GRanges"),
  by.coord = FALSE,
  gene_id_col = NULL,
  gene_name_col = NULL,
  ...
)
}
\arguments{
\item{GR}{An objects object from the any of the classes: 'pDMP', 'InfDiv',
GRangesList, GRanges or a list of GRanges.}

\item{GENES}{A GRanges object with gene coordinates and gene IDs. A
column named \strong{'gene_id'} carrying the gene ids should be included in
the metacolumns. If the meta-column named 'gene_id' is not provided, then
gene (region) ids will be created using the gene (region) coordinates.}

\item{ignore.strand, type}{Same as for
\code{\link[GenomicRanges]{findOverlaps-methods}}.}

\item{only.hypo, only.hyper}{logical(1). Whether to select only
hypo-methylated or hyper-methylated cytosine sites.}

\item{output}{Class of the object to be returned, a "list", or a "GRanges"
object.}

\item{by.coord}{logical(1). If TRUE, then the DMP are count per coordinate
and not per gene id.}

\item{gene_id_col}{Optional. An integer denoting the column from the
\emph{GENES} metacolumn where the gene ids are given.}

\item{gene_name_col}{Optional. An integer denoting the column from the
\emph{GENES} metacolumn where the gene 'name' are given.}

\item{...}{optional arguments for
\code{\link[GenomicRanges]{findOverlaps-methods}}. Users must evaluate
whether specific setting makes sense on each particular context.}
}
\value{
A a list GRanges object.
}
\description{
The function counts DMPs overlapping with gene-body. In fact,
this function also can be used to count DMPs overlapping with any set of
regions given in a GRanges object.
}
\details{
If \strong{\emph{by.coord == FALSE}} and "gene_id" is provided
in \strong{\emph{GENES}} argument, then DMP counts are made per gene-id.
Hence, DMPs from different regions with the same gene-id, say e.g. exons,
will be pooled in the count as they bear the same id.
}
\examples{
## Gene annotation
genes <- GRanges(seqnames = '1',
ranges = IRanges(start = c(3631, 6788, 11649), end = c(5899, 9130, 13714)),
strand = c('+', '-', '-'))
mcols(genes) <- data.frame(gene_id = c('AT1G01010', 'AT1G01020',
'AT1G01030'))

## Get a dataset of potential signals and the estimated cutpoint from the
## package
data(PS, cutpoint)

## The estimated cutpoints are used to discriminate signals from the noise.
## That is, DMPs are selected using the cupoints
DIMPs <- selectDIMP(PS, div.col = 9L, cutpoint = cutpoint$cutpoint)

## Finally DMPs found on genes
DIMR <- getDIMPatGenes(GR = DIMPs$T1, GENES = genes)

}
\seealso{
\code{\link{getDMPatRegions}}
}
