% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimateCutPoint.R
\name{estimateCutPoint}
\alias{estimateCutPoint}
\title{Estimate cutpoints to distinguish the treatment methylation signal
from the control}
\usage{
estimateCutPoint(
  LR,
  control.names,
  treatment.names,
  simple = TRUE,
  column = c(hdiv = TRUE, jdiv = TRUE, jdiv.stat = FALSE, TV = TRUE, bay.TV = FALSE,
    wprob = TRUE, pos = TRUE),
  classifier1 = c("logistic", "pca.logistic", "lda", "qda", "pca.lda", "pca.qda",
    "random_forest"),
  classifier2 = NULL,
  tv.cut = 0.25,
  tv.col = NULL,
  div.col = NULL,
  clas.perf = FALSE,
  post.cut = 0.5,
  prop = 0.6,
  center = FALSE,
  scale = FALSE,
  n.pc = 1,
  interactions = NULL,
  cut.values = NULL,
  stat = 0,
  cutp_data = FALSE,
  maxnodes = NULL,
  ntree = 400,
  nsplit = 1L,
  num.cores = 1L,
  tasks = 0L,
  ...
)
}
\arguments{
\item{LR}{An object from 'pDMP' class. This object is previously obtained
with function \code{\link{getPotentialDIMP}}.}

\item{control.names, treatment.names}{Names/IDs of the control and
treatment samples, which must be include in the variable LR.}

\item{simple}{Logic (default, TRUE). If TRUE, then Youden Index is used to
estimate the cutpoint. If FALSE, the minimum information divergence value
with posterior classification probability greater than \emph{post.cut}
(usually \emph{post.cut} = 0.5) as estimated by \emph{classifier1} will be
the reported cutpoint, except if a better cutpoint is found in the set of
values provided by the user in the parameter \emph{cut.values}.}

\item{column}{a logical vector for column names for the predictor variables
to be used: Hellinger divergence 'hdiv', total variation 'TV', probability
of potential DMP 'wprob', and the relative cytosine site position 'pos' in
respect to the chromosome where it is located. The relative position is
estimated as (x - x.min)/(x.max - x), where x.min and x.max are the maximum
and minimum for the corresponding chromosome, respectively. If 'wprob =
TRUE', then Logarithm base-10 of 'wprob' will be used as predictor in place
of 'wprob'.}

\item{classifier1, classifier2}{Classification model to use. Option
'logistic' applies a logistic regression model; option 'lda' applies a
Linear Discriminant Analysis (LDA); 'qda' applies a Quadratic Discriminant
Analysis (QDA), 'pca.logistic' applies logistic regression model using the
Principal Component (PCs) estimated with Principal Component Analysis (PCA)
as predictor variables. pca.lda' applies LDA using PCs as predictor
variables, and the option 'pca.qda' applies a Quadratic Discriminant
Analysis (QDA) using PCs as predictor variables. If classifier2 is not NULL,
then it will be used to evaluate the classification performance, and the
corresponding best fitted model will be returned.}

\item{tv.cut}{A cutoff for the total variation distance to be applied to
each site/range. Only sites/ranges \emph{k} with \eqn{TVD_{k} > tv.cut} are
are used in the analysis. Its value must be a number. \eqn{0 < tv.cut < 1}.
Default is \eqn{tv.cut = 0.25}.}

\item{tv.col}{Column number for the total variation to be used for filtering
cytosine positions (if provided).}

\item{div.col}{Column number for divergence variable for which the
estimation of the cutpoint will be performed.}

\item{clas.perf}{Logic. Whether to evaluate the classification performance
for the estimated cutpoint using a model classifier when 'simple=TRUE'.
Default, FALSE.}

\item{post.cut}{If 'simple=FALSE', this is posterior probability to decide
whether a DMPs belong to treatment group. Default \emph{post.cut} = 0.5.}

\item{prop}{Proportion to split the dataset used in the logistic regression:
group versus divergence (at DMPs) into two subsets, training and testing.}

\item{center}{A logical value indicating whether the variables should be
shifted to be zero centered.}

\item{scale}{A logical value indicating whether the variables should be}

\item{n.pc}{Number of principal components (PCs) to use if the classifier is
not 'logistic'. In the current case, the maximun number of PCs is 4.}

\item{interactions}{If a logistic classifier is used. Variable interactions
to consider in a logistic regression model. Any pairwise combination of the
variable 'hdiv', 'TV', 'wprob', and 'pos' can be provided. For example:
'hdiv:TV', 'wprob:pos', 'wprob:TV', etc.}

\item{cut.values}{Cut values of the information divergence (ID) specified
in \emph{div.col} where to check the classification performance (0 <
\emph{cut.interval} < max ID). If provided, the search for a cutpoint will
include these values.}

\item{stat}{An integer number indicating the statistic to be used in the
testing when \emph{simple} = FALSE The mapping for statistic names are:
\itemize{
\item 0 = 'Accuracy'
\item 1 = 'Sensitivity'
\item 2 = 'Specificity'
\item 3 = 'Pos Pred Value'
\item 4 = 'Neg Pred Value'
\item 5 = 'Precision'
\item 6 = 'Recall'
\item 7 = 'F1'
\item 8 = 'Prevalence'
\item 9 = 'Detection Rate'
\item 10 = 'Detection Prevalence'
\item 11 = 'Balanced Accuracy'
\item 12 = 'FDR'
}}

\item{cutp_data}{logical(1) (optional). If TRUE, and simple = TRUE, then a
data frame for further analysis or estimation of the optimal cutpoint
based  only on the selected divergence is provided.}

\item{maxnodes, ntree}{Only for Random Forest classifier
(\code{\link[randomForest]{randomForest}}). Maximum number maxnodes of
terminal nodes trees in the forest can have. If not given, trees are grown
to the maximum possible. Parameter ntree stands for the number of trees to
grow. This should not be set to too small a number, to ensure that every
input row gets predicted at least a few times.}

\item{nsplit}{Only for Random Forest classifier. The Random Forest
(\code{\link[randomForest]{randomForest}}) package uses a C+Fortran
implementation which only supports integer indexes, so any dataframe/data
table/matrix with >2^31 elements (limit for integers) gives an error. The
option nsplit is applied to train \eqn{ntrees=floor(ntree/nsplit)} models
with rep(ntrees,nsplit) which are finally combined to obtain a forest with
ntree. Each model would contain ntrees.}

\item{num.cores, tasks}{Parameters for parallel computation using package
\code{\link[BiocParallel]{BiocParallel-package}}: the number of cores to
use, i.e. at most how many child processes will be run simultaneously (see
\code{\link[BiocParallel]{bplapply}} and the number of tasks per job (only
for Linux OS).}

\item{...}{arguments passed to or from other methods.}
}
\value{
Depending the parameter setting will return the following list with
elements:
\enumerate{
\item cutpoint: Cutpoint estimated.
\item testSetPerformance: Performance evaluation on the test set.
\item testSetModel.FDR: False discovery rate on the test set.
\item model: Model used in the performance evaluation.
\item modelConfMatrix: Confusion matrix for the whole dataset
derived applying the model classifier used in the performance
evaluation.
\item initModel: Initial classifier model applied to estimate
posterior classifications used in the cutpoint estimation.
\item postProbCut: Posterior probability used to estimate the
cutpoint
\item classifier: Name of the model classifier used in the
performance evaluation.
\item statistic: Name of the performance statistic used to find the
cutpoint when \emph{simple} =  FALSE.
\item optStatVal: Value of the performance statistic at the
cutpoint.
}
}
\description{
Given a list of two GRanges objects, control and treatment,
carrying the potential signals (prior classification) from controls and
treatments in terms of an information divergence (given the metacolumns),
the function estimates the cutpoints of the control group versus treatment
group.
}
\details{
The function performs an estimation of the optimal cutpoint for the
classification of the differentially methylated (cytosines) positions into
two classes: DMPs from control and DMPs from treatment. The simplest approach
to estimate the cutpoint is based on the application of Youden Index. More
complexes approach based in several machine learning model are provided as
well.

Results of the classification performance resulting from the estimated
cutpoint are normally given, with the exception of those extreme situations
where the statistics to evaluate performance cannot be estimated. More than
one classifier model can be applied. For example, one classifier (logistic
model) can be used to estimate the posterior classification probabilities of
DMP into those from control and those from treatment. These probabilities are
then used to estimate the cutpoint in range of values from, say, 0.5 to 0.8.
Next, a different classifier can be used to evaluate the classification
performance. Different classifier models would yield different performances.
Models are returned and can be used in further prediction with new datasets
from the same batch experiment. This is a machine learning approach to
discriminate the biological regulatory signal naturally generated in the
control from that one induced by the treatment.

Notice that the estimation of an optimal cutpoint based on the application
Youden Index (simple = TRUE) only uses the information provided by the
selected information divergence. As a result, classification results based
only in one variable can be poor or can fail. However, option simple =
FALSE, uses the information from several variables following a
machine-learning (ML) approach.

Nevertheless, when simple = TRUE, still a ML model classifier can be built
using the optimal cutpoint estimated and setting clas.perf = TRUE. Such a ML
model can be used for predictions in further analyses with function
\code{\link{predictDIMPclass}}.
}
\examples{
## Get a dataset of potential signals and the estimates cutpoint
## from the package and performs cutpoint estimation
data(PS)

cutp <- mlCutpoint(LR = PS,
                   column = c(hdiv = TRUE, TV = TRUE,
                              wprob = TRUE, pos = TRUE),
                   classifier1 = 'qda', n.pc = 4,
                   control.names = c('C1', 'C2', 'C3'),
                   treatment.names = c('T1', 'T2', 'T3'),
                   tv.cut = 0.68, prop = 0.6,
                   cut.values = seq(114, 118, 1),
                   div.col = 9L)
cutp

}
\seealso{
\code{\link[MethylIT]{evaluateDIMPclass}}
}
