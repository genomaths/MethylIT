#' @rdname estimateDivergence
#'
#' @title Information Divergences of Methylation Levels
#' @description This function prepares the data for the estimation of
#' information divergences and works as a wrapper calling the functions that
#' compute selected information divergences of methylation levels. In the
#' downstream analysis, the probability distribution of a given information
#' divergence is used in Methyl-IT as the null hypothesis of the noise
#' distribution, which permits, in a further signal detection step, the
#' discrimination of the methylation regulatory signal from the background
#' noise.
#'
#' For the current version, two information divergences of methylation levels
#' are computed by default: 1) Hellinger divergence (\emph{H}) and 2) the total
#' variation distance (\emph{TVD}). In the context of methylation analysis
#' \emph{TVD} corresponds to the absolute difference of methylation levels.
#' Here, although the variable reported is the total variation (\emph{TV}), the
#' variable actually used for the downstream analysis is \emph{TVD}. Once a
#' differentially methylated position (DMP) is identified in the downstream
#' analysis, \emph{TV} is the standard indicator of whether the cytosine
#' position is hyper- or hypo-methylated.
#'
#' The option to compute the J-information divergence (JD) is currently
#' provided. The motivation to introduce this divergence is given in the help of
#' function \code{\link{estimateJDiv}}.
#'
#' @details If read counts are provided, then Hellinger divergence is computed
#' as given in the first formula from Theorem 1 from reference (1). In the
#' present case:
#'
#' \deqn{H = 2 (n_1 + 1) (n_2 + 1)*((sqrt(p_1) - sqrt(p_2))^2 +
#' (sqrt(1-p_2) - sqrt(1-p_2))^2)/(n_1 + n_2 + 2)}
#'
#' where \eqn{n_1} and \eqn{n_2} are the coverage for the control and
#' treatment, respectively. Notice that each row from the matrix of counts
#' correspond to a single cytosine position and has four values corresponding to
#' 'mC1' and 'uC1' (control), and mC2' and 'uC2' for treatment.
#'
#' According with the above equation, to estimate Hellinger divergence, not only
#' the methylation levels are considered in the estimation of H, but also the
#' control and treatment coverage at each given cytosine site. At this point, it
#' is worthy to do mention that if the reference sample is derived with function
#' \code{\link{poolFromGRlist}} using the 'sum' of read counts to compute a
#' methylation pool, then 'min.coverage' parameter value must be used to prevent
#' an over estimation of the divergence for low coverage cytosines sites. For
#' example, if a reference sample is derived as the methylation pool of read
#' count sum from 3 individuals and we want to consider only methylation sites
#' with minimum coverage of 4, then we can set min.coverage = c(12, 4), where
#' the number 12 (3 x 4) is the minimum coverage requested for the each cytosine
#' site in the reference sample.
#'
#' If the methylation levels are provided in place of counts, then the Hellinger
#' divergence is computed as:
#' \deqn{H = (sqrt(p_1) - sqrt(p_2))^2 + (sqrt(1 - p_1) - sqrt(1 - p_2))^2}
#'
#' This formula assumes that the probability vectors derived from the
#' methylation levels  \eqn{p_i = c(p_{i1}, 1 - p_{i2}}) (see
#' \code{\link{estimateHellingerDiv}} are an unbiased estimation of the expected
#' one. The function applies a pairwise filtering after building a single
#' GRanges from the two GRanges objects. Experimentally available cytosine sites
#' are paired using the function 'uniqueGRanges'.
#'
#' It is important to observe that several filtering conditions are provided to
#' select biological meaningful cytosine positions, which prevent to carry
#' experimental errors in the downstream analyses. By filtering the read count
#' we try to remove bad quality data, which would be in the edge of the
#' experimental error originated by the BS-seq sequencing. It is user
#' responsibility to check whether cytosine positions used in the analysis are
#' biological meaningful. For example, a cytosine position with counts mC1 = 10
#' and uC1 = 20 in the 'ref' sample and mC2 = 1 & uC2 = 0 in an 'indv' sample
#' will lead to methylation levels p1 = 0.333 and p2 = 1, respectively, and TV =
#' p2 - p1 = 0.667, which apparently indicates a hypermethylated site. However,
#' there are not enough reads supporting p2 = 1. A Bayesian estimation of TV
#' will reveal that this site would be, in fact, hypomethylated. So, the best
#' practice will be the removing of sites like that. This particular case is
#' removed under the default settings: min.coverage = 4, min.meth = 4, and
#' min.umeth = 0 (see example for function \code{\link{uniqueGRfilterByCov}},
#' called by 'estimateDivergence').
#'
#' @param ref The GRanges object of the reference individual that will be used
#'     in the estimation of the information divergence.
#' @param indiv A list of GRanges objects from the individuals that will be
#'     used in the estimation of the information divergence.
#'@param Bayesian logical(1). Whether to perform the estimations based on
#'     posterior estimations of methylation levels.
#' @param init.pars  initial parameter values. Defaults is NULL and an initial
#' guess is estimated using \code{\link[stats]{optim}} function. If the initial
#' guessing fails initial parameter values are to alpha = 1 &
#' beta = 1, which imply the parsimony pseudo-counts greater than zero.
#'
#' @param columns Vector of one or two integer numbers denoting the indexes of
#'     the columns where the methylated and unmethylated read counts are found
#'     or, if meth.level = TRUE, the columns corresponding to the methylation
#'     levels. If columns = NULL and meth.level = FALSE, then columns = c(1,2)
#'     is assumed. If columns = NULL and meth.level = TRUE, then columns = 1 is
#'     assumed.
#' @param min.coverage An integer or an integer vector of length 2. Cytosine
#'     sites where the coverage in both samples, 'x' and 'y', are less than
#'     'min.coverage' are discarded. The cytosine site is preserved, however, if
#'     the coverage is greater than 'min.coverage' in at least one sample. If
#'     'min.coverage' is an integer vector, then the corresponding min coverage
#'     is applied to each sample.
#' @param min.meth An integer or an integer vector of length 2. Cytosine sites
#'     where the numbers of read counts of methylated cytosine in both samples,
#'     '1' and '2', are less than 'min.meth' are discarded. If 'min.meth' is an
#'     integer vector, then the corresponding min number of reads is applied to
#'     each sample. Default is min.meth = 4.
#' @param min.umeth An integer or an integer vector of length 2
#'     (\eqn{min.umeth = c(min.umeth1, min.umeth2)}). Min number of
#'     reads to consider cytosine position. Specifically cytosine positions
#'     where (uC <= min.umeth) & (mC > 0) & (mC <= min.meth) hold will be
#'     removed, where mC and uC stand for the numbers of methylated and
#'     unmethylated reads. Default is min.umeth = 0.
#' @param min.sitecov An integer. The minimum total coverage. Only sites where
#'     the total coverage (cov1 + cov2) is greater than 'min.sitecov' are
#'     considered for downstream analysis, where cov1 and cov2 are the coverages
#'     for samples 1 and 2, respectively.
#' @param high.coverage An integer for read counts. Cytosine sites having
#'     higher coverage than this are discarded.
#'@param percentile Threshold to remove the outliers from each file and all
#'     files stacked.
#' @param JD Logic (Default:FALSE). Option on whether to add a column with
#'     values of J-information divergence (see \code{\link{estimateJDiv}}).
#'     It is only compute if JD = TRUE and meth.level = FALSE.
#' @param jd.stat logical(1). Whether to compute the \eqn{JD} statistic with
#' asymptotic Chi-squared distribution with one degree of freedom (see
#' \code{\link{estimateJDiv}}).
#' @param num.cores The number of cores to use, i.e. at most how many child
#'     processes will be run simultaneously (see 'bplapply' function from
#'     BiocParallel package).
#' @param tasks integer(1). The number of tasks per job. value must be a scalar
#'     integer >= 0L. In this documentation a job is defined as a single call
#'     to a function, such as bplapply, bpmapply etc. A task is the division of
#'     the X argument into chunks. When tasks == 0 (default), X is divided as
#'     evenly as possible over the number of workers (see MulticoreParam from
#'     BiocParallel package).
#'@param meth.level logical(1) Whether methylation levels are given in place of
#'     counts. Default is FALSE.
#' @param logbase Logarithm base used to compute the JD (if JD = TRUE).
#'     Logarithm base 2 is used as default (bit unit). Use
#'     logbase = \eqn{exp(1)} for natural logarithm.
#' @param verbose if TRUE, prints the function log to stdout
#' @param ... Optional parameters for \code{\link{uniqueGRanges}} function.
#'
#' @return An object from 'infDiv' class with the four columns of counts, the
#' information divergence, and additional columns:
#' \describe{
#' \item{1) \strong{A matrix:}}{The original matrix of methylated
#' \eqn{c_{ij}} and unmethylated \eqn{t_{ij}} read counts from control \eqn{j=1}
#' and treatment \eqn{j=2} samples at positions \eqn{i}.}
#'
#' \item{2) \strong{'p1' and 'p2':}}{methylation levels for control and
#' treatment, respectively. If 'meth.level = FALSE' and 'Bayesian = TRUE'
#' (recommended), 'p1' and 'p2' are estimated following the Bayesian approach
#' described in reference (1).}
#'
#' \item{3) \strong{'bay.TV':}}{total variation TV = p2 - p1}
#'
#' \item{4) \strong{'TV':}}{total variation based on simple counts:
#' \eqn{TV=c_1/(c_1+t_1)-c_2/(c_2+t_2)}, where \eqn{c_i} and \eqn{t_i} denote
#' methylated and unmethylated read counts, respectively.}
#'
#' \item{5) \strong{Hellinger divergence, 'hdiv':}}{If Bayesian = TRUE, the
#' results are based on the posterior estimations of methylation levels. if
#' meth.level = FALSE', then 'hdiv' is computed as given in reference (2),
#' otherwise as: \deqn{hdiv = (sqrt(p_1) - sqrt(p_2))^2 + (sqrt(1 -p_1) - sqrt(1
#' - p_2))^2}}
#' }
#'
#' @references
#' \enumerate{
#' \item Sanchez R, Yang X, Maher T, Mackenzie S. Discrimination of DNA
#' Methylation Signal from Background Variation for Clinical Diagnostics. Int.
#' J. Mol Sci, 2019, 20:5343.
#' \item Basu  A., Mandal  A., Pardo L. Hypothesis testing for two
#' discrete populations based on the Hellinger distance. Stat. Probab.
#' Lett. 2010, 80: 206-214.
#' }
#' @author Robersy Sanchez (\url{https://genomaths.com})
#' @examples
#' ## The read count data are created
#' num.samples <- 250
#' x <- data.frame(chr = 'chr1', start = 1:num.samples,
#'                 end = 1:num.samples,strand = '*',
#'                 mC = rnbinom(size = num.samples, mu = 4, n = 500),
#'                 uC = rnbinom(size = num.samples, mu = 4, n = 500))
#'
#' y <- data.frame(chr = 'chr1', start = 1:num.samples, end = 1:num.samples,
#'                 strand = '*', mC = rnbinom(size = num.samples,
#'                 mu = 4, n = 500),
#' uC = rnbinom(size = num.samples, mu = 4, n = 500))
#'
#' x <- makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
#' y <- makeGRangesFromDataFrame(y, keep.extra.columns = TRUE)
#' hd <- estimateDivergence(ref = x, indiv = list(y), JD = TRUE,
#' verbose = FALSE)[[1]]
#'
#' ## Keep in mind that Hellinger and J divergences are, in general,
#' ## correlated
#' cor.test(x = as.numeric(hd$hdiv), y = as.numeric(hd$jdiv),
#'         method = 'kendall')
#'
#' @importFrom BiocParallel MulticoreParam bplapply SnowParam
#' @importFrom GenomicRanges GRanges GRangesList
#' @seealso \code{\link{estimateBayesianDivergence}}
#' @export
estimateDivergence <- function(
                        ref,
                        indiv,
                        Bayesian = FALSE,
                        init.pars = NULL,
                        columns = NULL,
                        min.coverage = 4,
                        min.meth = 4,
                        min.umeth = 0,
                        min.sitecov = 4,
                        high.coverage = NULL,
                        percentile = 0.999,
                        JD = FALSE,
                        jd.stat = FALSE,
                        num.cores = 1L,
                        tasks = 0L,
                        meth.level = FALSE,
                        logbase = 2,
                        verbose = TRUE,
                        ...) {

    if (is.null(columns) && (!meth.level)) columns <- c(1, 2)
    if (meth.level && (is.null(columns))) columns <- 1
    sn <- names(indiv)

    progressbar = FALSE
    if (verbose) progressbar = TRUE
    if (Sys.info()["sysname"] == "Linux")
        bpparam <- MulticoreParam(workers = num.cores, tasks = tasks,
                                progressbar = progressbar)
    else bpparam <- SnowParam(workers = num.cores, type = "SOCK",
                            progressbar = progressbar)
    if (ncol(mcols(ref)) > 2)
        ref <- ref[, columns]
    indiv <- lapply(indiv, function(x) x[, columns])

    if (meth.level) {
        x = bplapply(seq_len(length(indiv)), function(k, ref, indv, sn) {
            x <- indv[[k]]
            x <- x[, columns]
            x <- uniqueGRanges(
                            ListOfGranges = list(ref, x),
                            num.cores = num.cores,
                            tasks = tasks, type = "equal",
                            verbose = verbose, ...)
            x = estimateBayesianDivergence(
                                        x,
                                        Bayesian = FALSE,
                                        JD = JD,
                                        logbase = logbase,
                                        meth.level = meth.level,
                                        num.cores = num.cores,
                                        tasks = tasks,
                                        verbose = verbose)
            return(x)
        }, BPPARAM = bpparam, ref = ref, indv = indiv, sn = sn)
    } else {
        x = bplapply(seq_len(length(indiv)), function(k,
            ref, indv, sn) {
            if (verbose)
                message("*** Processing sample #", k, " ", sn[k])
            x = uniqueGRfilterByCov(
                            x = ref,
                            y = indv[[k]],
                            min.coverage = min.coverage,
                            min.meth = min.meth,
                            min.umeth = min.umeth,
                            min.sitecov = min.sitecov,
                            percentile = percentile,
                            high.coverage = high.coverage,
                            num.cores = 1L,
                            tasks = tasks,
                            verbose = verbose,
                            type = "equal")
            if (length(x) < 2)
                stop("*** At least two cytosine sites must pass the filtering",
                    " conditions to estimate informations divergences. \n",
                    "The issue was found at sample number: ",
                    k, ", id: ", names(indv)[k])
            x = estimateBayesianDivergence(
                                        x,
                                        Bayesian = Bayesian,
                                        init.pars = init.pars,
                                        JD = JD,
                                        jd.stat = jd.stat,
                                        num.cores = num.cores,
                                        tasks = tasks,
                                        meth.level = meth.level,
                                        logbase = logbase,
                                        verbose = verbose,
                                        ...)
            return(x)
        }, BPPARAM = bpparam, ref = ref, indv = indiv, sn = sn)
    }
    names(x) <- sn
    x <- structure(x, class = c("InfDiv", "list"))
    return(x)
}
