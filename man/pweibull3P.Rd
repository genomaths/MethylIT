% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pweibull3P.R
\name{pweibull3P}
\alias{pweibull3P}
\title{Weibull distribution with three parameters}
\usage{
pweibull3P(q, shape = 1, scale = 1, mu = 0)
}
\arguments{
\item{q}{vector of quantiles}

\item{shape}{shape parameter, or slope, defaulting to 1}

\item{scale}{scale parameter, or characteristic life,  defaulting to 1}

\item{mu}{location parameter, or failure free life,  defaulting to 0}
}
\value{
3 parameters Weibull distribution
}
\description{
Density, distribution function, quantile function and random
generation for the Weibull distribution with three parameters
}
\examples{
num.samples <- 10000
shape <- 0.75
scale <- 1
x <- rweibull(num.samples, shape = shape, scale = scale)
wei.model <- weibull3P(x)
y <- pweibull3P(x,
                shape = as.numeric(wei.model$Estimate[1]),
                scale = as.numeric(wei.model$Estimate[2]),
                mu = as.numeric(wei.model$Estimate[3]) )

}
