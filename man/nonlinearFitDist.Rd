% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nonlinearFitDist.R
\name{nonlinearFitDist}
\alias{nonlinearFitDist}
\title{Nonlinear fit of Information divergences distribution}
\usage{
nonlinearFitDist(
  LR,
  column = 9,
  dist.name = "Weibull",
  sample.size = 20,
  location.par = FALSE,
  absolute = FALSE,
  min.scale = 1e-05,
  min.mu = 0.01,
  mu.rel.err = 0.89,
  npoints = NULL,
  model = "all",
  maxiter = 1024,
  tol = 1e-12,
  ftol = 1e-12,
  ptol = 1e-12,
  minFactor = 10^-6,
  num.cores = multicoreWorkers(),
  tasks = 0L,
  maxfev = 1e+05,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{LR}{A list of \code{\link[GenomicRanges]{GRanges-class}} objects with
information divergence values in their meta-columns.}

\item{column}{An integer number denoting the index of the GRanges column
where the information divergence is given. Default column = 1}

\item{dist.name}{Name(s) of the distribution to fit. A single character
string or character vector naming the distribution(s): 'Weibull' (default),
gamma with three-parameter (Gamma3P), gamma with two-parameter (Gamma2P),
generalized gamma with three-parameter ('GGamma3P') or four-parameter
('GGamma4P'), and Log-Normal (LogNorm).}

\item{sample.size}{size of the sample}

\item{location.par}{whether to consider the fitting to generalized gamma
distribution (GGamma) including the location parameter, i.e., a GGamma
with four parameters (GGamma4P).}

\item{absolute}{Logic (default, FALSE). Total variation (TV, the difference
of methylation levels) is normally an output in the downstream MethylIT
analysis. If 'absolute = TRUE', then TV is transformed into |TV|, which is an
information divergence that can be fitted to Weibull or to Generalized Gamma
distribution.}

\item{min.scale}{A number. The nonlinear fit of GGamma family distributions
sometimes yields weird estimation of the \emph{scale} parameter with
values < 1e-5. Based on our experience, scale parameter values lower than
min.scale = 1e-5 are probably meaningless. The result of numerical errors,
the nonlinear algorithm are not perfect.}

\item{min.mu, mu.rel.err}{Numbers. Models with location parameter \eqn{\mu}
sometimes reports \eqn{\mu} values close to zero, which generally corresponds
to a non statistically significant model coefficient/parameter. This detail
is not detected by the AIC or the cross-validations correlation coefficient.
In this case the best model is the model without location parameter, where
\eqn{\mu = 0}. If \eqn{\mu < min.mu} and
\eqn{|E(\mu) - mean(idiv)|/mean(idiv) > mu.rel.err}, where \eqn{E(\mu)} is
the estimated value and \eqn{mean(idiv)} is the mean observed value of the
information divergence provided, then the issue is reported in the column
named \emph{'model'}.}

\item{npoints}{number of points to be used in the fit. Default is NULL.}

\item{model}{Optional. Only when dist.name = 'Weibull'. A selection of the
distribution model, two-parameters and three-parameters Weibull model ('2P'
and '3P'). Default is 'all' and the model with the best AIC criterion is
reported. Alternatively, just use dist.name = 'Weibull2P' or dist.name =
'Weibull3P'.}

\item{maxiter}{positive integer. Termination occurs when the number of
iterations reaches maxiter. Default value: 1024}

\item{tol}{A positive numeric value specifying the tolerance level for the
relative offset convergence criterion. Default value: 1e-12,}

\item{ftol}{non-negative numeric. Termination occurs when both the actual
and predicted relative reductions in the sum of squares are at most ftol.
Therefore, ftol measures the relative error desired in the sum of squares.
Default value: 1e-12.}

\item{ptol}{non-negative numeric. Termination occurs when the relative error
between two consecutive iterates is at most ptol. Therefore, ptol measures
the relative error desired in the approximate solution. Default value: 1e-12.}

\item{minFactor}{A positive numeric value specifying the minimum step-size
factor allowed on any step in the iteration. The increment is calculated with
a Gauss-Newton algorithm and successively halved until the residual sum of
squares has been decreased or until the step-size factor has been reduced
below this limit. Default value: 10^-6.}

\item{num.cores}{The number of cores to use, i.e. at most how many child
processes will be run simultaneously (see
\code{\link[BiocParallel]{bpmapply}} function from BiocParallel package).}

\item{tasks}{integer. The number of tasks per job. value must be a scalar
integer >= 0L. In this documentation a job is defined as a single call to a
function, such as bplapply, bpmapply etc. A task is the division of the X
argument into chunks. When tasks == 0 (default), X is divided as evenly as
possible over the number of workers (see
\code{\link[BiocParallel]{MulticoreParam-class}} from BiocParallel package).}

\item{maxfev}{integer; termination occurs when the number of calls to fn has
reached maxfev. Note that nls.lm sets the value of maxfev to 100*(length(par)
\itemize{
\item 
\enumerate{
\item if maxfev = integer(), where par is the list or vector of parameters to
be optimized.
}
}}

\item{verbose}{If TRUE, prints the function log to stdout}

\item{...}{Additional fitting parameters.}
}
\value{
Model table with coeficients and goodness-of-fit results:
Adj.R.Square, deviance, AIC, R.Cross.val, and rho, as well as, the
coefficient covariance matrix.
}
\description{
A wrapper to call functions \code{\link{weibull3P}},
\code{\link{fitGammaDist}}, and \code{\link{fitGGammaDist}} to operate on a
list of \code{\link[GenomicRanges]{GRanges-class}} objects.
}
\details{
The algorithm prepares the information divergence variable
to try fitting Weibull or generalized gamma distribution model to the data.
If Weibull distribution is selected (default: 'Weibull'), function
'Weibull2P' first attempts fitting to the two-parameter Weibull CDF
(Weibull2P). If Weibull2P did not fit, then the algorithm will try to fit
Weibull3P. The Levenberg-Marquardt algorithm implemented in R package
'minpack.lm' is used to perform the nonlinear fit. Cross-validations for the
nonlinear regressions (R.Cross.val) are performed in each methylome as
described in reference (1-2). In addition, Stein's formula for adjusted R
squared (rho) is used as an estimator of the average cross-validation
predictive power (2).

If 'GGamma3P' is selected the call to function \code{\link{fitGGammaDist}}
permits the fitting to the three-parameter GGamma CDF ('GGamma3P'). The fit
to the four-parameter GGamma ('GGamma4P') is also available. GGamma
distribution are fitted using a modification of Levenberg-Marquardt algorithm
implemented in function 'nls.lm' from the 'minpack.lm' R package. Notice that
the fit to GGamma distribution is computationally time consuming (see
?fitGGammaDist for additional information).

Notice that all the estimated parameters/coefficients must be positive.
Negative parameter values would happen and in that case the model must be
dicarded.
}
\examples{
## Load a dataset with Hellinger Divergence of methylation levels on it.
data(HD)

## The nonlinear fit based on three-parameter GGamma distribution
nlms2 <- nonlinearFitDist(HD, npoints = 100, dist.name = 'GGamma3P',
                            verbose = FALSE)

## Weilbull distribution is a particular case of GGamma.
nlms <- nonlinearFitDist(HD, npoints = 100, verbose = FALSE)

## The goodness-of-fit indicators AIC suggests that the best fitted model
## is obtained with GGamma distribution (in this example).
res <- mapply(function(m1,m2) as.numeric(c(Weibull = m1$AIC[1],
                                        GGamma = m2$AIC[1])),
                nlms, nlms2)
rownames(res) <-c('Weibull', 'GGamma')
res

## However, the Cross-validations correlation coefficient is saying that
## the Weibull distribution would be a little better probability
## predictor.
res <- mapply(function(m1,m2) as.numeric(c(Weibull = m1$R.Cross.val[1],
                                        GGamma = m2$R.Cross.val[1])),
                nlms, nlms2)
rownames(res) <-c('Weibull', 'GGamma')
res

}
\references{
\enumerate{
\item R. Sanchez and S. A. Mackenzie, “Information Thermodynamics of
Cytosine DNA Methylation,” PLoS One, vol. 11, no. 3, p. e0150427,
Mar. 2016.
\item Stevens JP. Applied Multivariate Statistics for the Social
Sciences. Fifth Edit. Routledge Academic; 2009.
}
}
\seealso{
\code{\link{gofReport}}
}
\author{
Robersy Sanchez 01/31/2018 \url{https://github.com/genomaths}
}
