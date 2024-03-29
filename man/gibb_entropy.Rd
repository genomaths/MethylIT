% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gibb_entropy.R
\name{gibb_entropy}
\alias{gibb_entropy}
\alias{gibb_entropy,missingORNULL-method}
\alias{gibb_entropy,cdfMODEL-method}
\alias{gibb_entropy,cdfMODELlist-method}
\alias{gibb_entropy,ProbDistrList-method}
\title{Gibbs entropy of Generalized Gamma Distribution}
\usage{
gibb_entropy(model, R = 8.31446261815324, ...)

\S4method{gibb_entropy}{missingORNULL}(
  model,
  R = 8.31446261815324,
  pars,
  log.base = exp(1),
  terms = FALSE
)

\S4method{gibb_entropy}{cdfMODEL}(model, R = 8.31446261815324, log.base = exp(1), terms = FALSE)

\S4method{gibb_entropy}{cdfMODELlist}(model, R = 8.31446261815324, log.base = exp(1), terms = FALSE)

\S4method{gibb_entropy}{ProbDistrList}(model, R = 8.31446261815324, log.base = exp(1), terms = FALSE)
}
\arguments{
\item{model}{An object from any of the classes created in MethylIT
pipeline: \emph{\strong{cdfMODEL}}, \emph{\strong{cdfMODELlist}}, or
\emph{\strong{ProbDistrList}}. If given, then the parameter values are taken
from the model.}

\item{R}{A number or NULL. The gas constant
(\eqn{R = 8.31446261815324 J * (K * mol)^-1}) is given as default value,
which is proportionality constant that relates the energy scale in physics
to the temperature scale and the scale used for amount of substance.}

\item{pars}{Optional. A numerical vector containing the model parameter
values in the given in order: \strong{\emph{alpha, scale}}, and
\strong{\emph{delta}}.}

\item{log.base}{A positive number. The same as in function
\code{\link[base]{log}}.}

\item{terms}{logical(1). If term = TRUE, then a numerical vector with terms
contributing to the Gibb entropy are provided.}
}
\value{
Gibb entropy of the model. Shannon entropy is returned by setting
R = 1 and log.base = 2.
}
\description{
This function computes the Gibbs entropy for member of the
Generalized Gamma (GG) Distribution family. GG density is given as:

\eqn{exp(-y^\alpha) * \alpha*y^(\alpha*\delta - 1)/(scale*\gamma(\delta))}

(see \href{https://goo.gl/1n4kpW}{(Wikipedia)})

A list of the member of the GG distribution family with the corresponding
parameter settings is provided in Table 1 from reference (1). For example,
to compute the Gibbs entropy of a Weibull distribution model, we just set:
\eqn{alpha > 0} and \eqn{delta = 1} (notice that parameter are named
different in reference (1)).
}
\details{
The value of the constant R can be simply 1, which returns the
Shannon entropy in bit units, only carrying informational meaning.
}
\examples{
## Loading the probability distribution models
data(gof, "MethylIT")
## Gibb entropy in J * (K * mol)^-1)
gibb_entropy(gof)
}
\references{
\enumerate{
\item Crooks, Gavin E. (2015) The Amoroso Distribution.
\href{https://arxiv.org/pdf/1005.3274.pdf}{arXiv:1005.3274v2}.
}
}
\seealso{
\code{\link{helmholtz_free_energy}}
}
