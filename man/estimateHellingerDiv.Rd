% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimateHellingerDiv.R
\name{estimateHellingerDiv}
\alias{estimateHellingerDiv}
\title{Hellinger divergence of methylation levels}
\usage{
estimateHellingerDiv(p, n = NULL)
}
\arguments{
\item{p}{A numerical vector of the methylation levels p = c(p1, p2) of
individuals 1 and 2.}

\item{n}{if supplied, it is a vector of integers denoting the coverages used
in the estimation of the methylation levels.}
}
\value{
The Hellinger divergence value for the given methylation levels is
returned
}
\description{
Given a the methylation levels of two individual, the function
computes the information divergence between methylation levels.
}
\details{
The methylation level \eqn{p_ij} for an individual \eqn{i} at
cytosine site \eqn{j} corresponds to a probability vector \eqn{p^ij = (p_ij,
1 - p_ij)}. Then, the information divergence between methylation levels
\eqn{p^1j} and \eqn{p^2j} from individuals 1 and 2 at site \eqn{j} is the
divergence between the vectors \eqn{p^1j = (p_1j, 1 - p_1j)} and \eqn{p^2j =
(p_2j, 1 - p_2j)}. If the vector of coverage is supplied, then the
information divergence is estimated according to the formula:

\deqn{hdiv = 2*(n_1 + 1)*(n_2 + 1)*((sqrt(p_1j) - sqrt(p_2j))^2 +
(sqrt(1 - p_1j) - sqrt(1 - p_2j))^2)/(n_1 + n_2 + 2)}

This formula corresponds to Hellinger divergence as given in the first
formula from Theorem 1 from reference 1. Otherwise:

\deqn{hdiv = (sqrt(p_1j) - sqrt(p_2j))^2 + (sqrt(1 - p_1j) -
sqrt(1 - p_2j))^2}

Missing methylation levels, reported as NA or NaN, are replaced with zero.
}
\examples{
    p <- c(0.5, 0.5)
    estimateHellingerDiv(p)

}
\references{
' 1. Basu  A., Mandal  A., Pardo L (2010) Hypothesis testing for
two discrete populations based on the Hellinger distance. Stat Probab
Lett 80: 206-214.
}
