% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lapply.R
\name{lapply}
\alias{lapply}
\alias{lapply.default}
\title{Apply a function over a list-like object preserving its attributes}
\usage{
lapply(x, FUN, ...)

\method{lapply}{default}(x, FUN, keep.attr = FALSE, ...)
}
\arguments{
\item{x}{A list-like or vector-like object}

\item{FUN, ...}{See ?base::\code{\link[base]{lapply}} for a description of
these arguments.}

\item{keep.attr}{Logic. If TRUE, then the original attributes from 'x' are
preserved in the returned list. Default is FALSE.}
}
\value{
Same as in ?base::\code{\link[base]{lapply}} if keep.attr = FALSE.
Otherwise same values preserving original attributes from 'x'.
}
\description{
lapply returns a list of the same length as 'x', each element of
which is the result of applying FUN to the corresponding element of 'x'.
}
\examples{
# Create a list
x <- list(a = 1:10, beta = exp(-3:3), logic = c(TRUE,FALSE,FALSE,TRUE))
class(x) <- 'nice'

# compute the list mean for each list element using 'base::lapply'
class(lapply(x, mean))

# To preserve attributes
class(lapply(x, mean, keep.attr = TRUE))
}
\seealso{
base::\code{\link[base]{lapply}}
}
\keyword{internal}
