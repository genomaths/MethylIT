% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pcaLDA.R
\name{pcaLDA}
\alias{pcaLDA}
\alias{predict.pcaLDA}
\title{Linear Discriminant Analysis (LDA) using Principal Component Analysis
(PCA)}
\usage{
pcaLDA(
  formula = NULL,
  data = NULL,
  grouping = NULL,
  n.pc = 1,
  scale = FALSE,
  center = FALSE,
  tol = 1e-04,
  method = "moment",
  max.pc = NULL,
  columns = 9L,
  ...
)

predict.pcaLDA(
  object,
  newdata,
  type = c("lda.pred", "class", "posterior", "scores", "pca.ind.coord", "all"),
  ...
)
}
\arguments{
\item{formula}{Same as in \code{\link[MASS]{lda}} from package 'MASS'.}

\item{data}{Same as in \code{\link[MASS]{lda}} from package 'MASS' or an
object from "pDMP" or "InfDiv" class.}

\item{grouping}{Same as in \code{\link[MASS]{lda}} from package 'MASS'.}

\item{n.pc}{Number of principal components to use in the LDA.}

\item{scale, center, tol, max.pc}{Same as in \code{\link[prcomp]{prcomp}} from
package 'prcomp'.}

\item{method}{Same as in \code{\link[MASS]{lda}} from package 'MASS'.}

\item{columns}{Optional. Only used if 'data' belong to the "pDMP" or
"InfDiv" class. Default is 9L.}

\item{...}{Not in use.}

\item{object}{To use with function 'predict'. A 'pcaLDA' object containing a
list of two objects: 1) an object of class inheriting from 'lda' and 2)
an object of class inheriting from 'prcomp'.}

\item{newdata}{To use with function 'predict'. New data for classification
prediction}

\item{type}{To use with function 'predict'. The type of prediction
required.  The default is 'all' basic predictions: classes and posterior
classification probabilities. Option 'lda.pred' returns the object given by
function 'predict.lda' from MASS package: 'class', 'posterior', 'scores'
(cases scores on discriminant variables, see \code{\link[MASS]{lda}}.}
}
\value{
Function 'pcaLDA' returns an object ('pcaLDA' class) consisting of list
with two objects:
\enumerate{
\item 'lda': an object of class \code{\link[MASS]{lda}} from package
'MASS'.
\item 'pca': an object of class \code{\link[stats]{prcomp}} from package
'stats'.
}
For information on how to use these objects see ?lda and ?prcomp.
}
\description{
The principal components (PCs) for predictor variables provided
as input data are estimated and then the individual coordinates in the
selected PCs are used as predictors in the LDA

Predict using a PCA-LDA model built with function 'pcaLDA'
}
\details{
The principal components (PCs) are obtained using the function
'prcomp' from R package 'stats', while the LDA is performed using the 'lda'
function from R package 'MASS'. The current application only uses basic
functionalities of mentioned functions. As shown in the example, pcaLDA'
function can be used in general classification problems.
}
\examples{
data(iris)
ld1 <- pcaLDA(formula = Species ~ Petal.Length + Sepal.Length + Sepal.Width,
data = iris, n.pc = 1, max.pc = 2, scale = TRUE, center = TRUE)

## ===== Prediction ===== ##
ld2 <- pcaLDA(formula = Species ~., data = iris, n.pc = 1, max.pc = 2,
scale = TRUE, center = TRUE)

set.seed(123)
idx <- sample.int(150, 40)
newdata <- iris[idx, 1:4]
newdata.prediction <- predict(ld2, newdata = newdata)

## ==== The confusion matrix
x <- data.frame(TRUE.class = iris$Species[idx],
PRED.class = newdata.prediction$class)
table(x)

}
\seealso{
\code{\link{pcaQDA}}, \code{\link[MASS]{lda}} and
\code{\link[MASS]{predict.lda}}
}
