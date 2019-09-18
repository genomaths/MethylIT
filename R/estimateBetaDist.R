#' @rdname estimateBetaDist
#'
#' @title Select the beta distribution that fit specified quantiles
#' @description This function perform a rough estimation of the shape
#'     parameters of beta distribution
#' @details In order to obtain the estimates for shape parameters of beta
#'     distribution, the squared of the difference between the empirical
#'     cumulative distribution function (ecdf) & the theoretical cdf is 
#'     minimized using the Non-Linear Minimization function of 'stats' 
#'     package.
#'
#' @param q prior probabilities
#' @param init.pars  initial parameter values. Defaults to alpha = 1 &
#'     beta = 1, which imply the parsimony pseudo-counts greater than zero.
#'
#' @return the estimated values of the shape parameters of the selected beta
#'    distribution.
#'
#' @importFrom stats ecdf optim nlm pbeta
#' @keywords internal
.estimateBetaDist <- function(q, init.pars=c(1, 1)) {
   ## q: prior probabilities
   ## init.pars: initial parameter values. Defaults to alpha = 1
   ## & beta = 1, which the parsimony pseudo-counts greater than zero.
   Q <- ecdf(q) ## Empirical Cumulative Distribution Function
   dat <- data.frame(Q=Q(q), q=q)
   #### Beta fitting
   ## In order to obtain the estimates for shape paramaters
   ## the squared of the difference between the ecdf &
   ## the theoretical cdf is minimized
   min.RSS <- function(data, inits) {
       alpha <- inits[1] ## This corresponds to the initial
       ##    starting parameter for alpha
       beta <- inits[2] ## This corresponds to the initial
       ##    starting parameter for beta
       ## Because optim minimizes a function, the
       with(data, sum(Q - pbeta(q, alpha, beta)) ^ 2)
   }
   pars <- try(suppressWarnings(nlm(min.RSS, p=init.pars, data=dat)),
               silent=TRUE)
   ## print(pars$minimum)
   ## print(pars$estimate)

   if (inherits(pars, "try-error")) {
       pars <- try(suppressWarnings(optim(par=init.pars, min.RSS, data=dat,
                                           method="BFGS",
                                           control=list(maxit=500,
                                                       abstol=(10 ^ -8)))),
                   silent=TRUE)
       par <- pars$par
       ## print(pars$value)
       ## print(pars$par)
   } else par <- pars$estimate

   if (inherits( pars, "try-error")) {
       par <- c(0, 0)
   }
   par
}
