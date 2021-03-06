% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimateGLM.R
\name{.estimateGLM}
\alias{.estimateGLM}
\title{Poisson and Negative Binomial regression analysis.}
\usage{
.estimateGLM(
  x,
  groups,
  baseMV,
  w,
  MVrate,
  test = c("Wald", "LRT"),
  p.value = NULL
)
}
\arguments{
\item{x}{Matrix of counts.}

\item{groups}{Groups information derived from a
\code{\link[MethylIT]{glmDataSet}} object.}

\item{baseMV}{Mean and variance of group counts. If
baseMean >= baseVar*MVrate, then the nonlinear fit to 'Poisson' and
'QuasiPoisson' models are performed, otherwise only the nonlinear fit to
'Neg.Binomial' and 'Neg.Binomial with weights' models are performed}

\item{w}{group weights used in glm procedure}

\item{MVrate}{Minimum Mean/Variance rate.}

\item{test}{A character string matching one of 'Wald' or 'LRT'. If test =
'Wald', then the p-value of the Wald test for the coefficient of the
independent variable (\emph{treatment group}) will be reported.
If test = 'LRT', then the p-value from a likelihood ratio test given by
\code{\link[stats]{anova}} function from \emph{stats} packages will be
the reported p-value for the group comparison when the best fitted model
is the negative binomial. As suggested for \code{\link[stats]{glm}}, if
best fitted model is Poisson or quasi-Poisson, then the best test is
'Chi-squared' or 'F-test', respectively. So, for the sake of simplicity,
the corresponding suitable test will be applied when test = 'LRT'.}

\item{p.value}{Cut off p-value to reject the null hypothesis}
}
\value{
GLM model of the group comparison for the given genomic region
}
\description{
This function is called internally by countTest function. You
would need to call it directly only in very special cases.

Perform Poisson and Negative Binomial regression analysis to
compare the counts from different groups, treatment and control x:
vector of counts groups: factor labeling the members from each group
Evaluated models are 'Poisson', 'Quasipoisson', 'Neg.Binomial.W', and
'Neg.Binomial'
}
\examples{
## Get "RangedGlmDataSet" object
data(ds, package = "MethylIT")

X <- ds$counts[69,]
baseMeanAndVar <- data.frame(baseMean = mean(X),
                             baseVar = var(X))

MethylIT:::.estimateGLM(x = X, groups = ds$colData$condition,
                        baseMV = baseMeanAndVar,
                        w = c(1,1), MVrate = 0.95,
                        test = "LRT")

}
\keyword{internal}
