\name{crossval,scoopfit-methods}
\docType{methods}
\alias{crossval,scoopfit-method}
\alias{crossval.scoop}
\alias{crossval}
\title{Cross-validation of a scoop model}
\description{Function that computes K-fold cross-validated error of a
  \code{scoop} fit.}
\usage{
\S4method{crossval}{scoopfit}(object,
          K         = 10,
          error     = "deviance",
          threshold = 0.5,
          folds     = NULL,
          verbose   = TRUE)
}
\arguments{
  \item{object}{output of a \code{scoop} run (must be an object of
    class \code{scoopfit}).}
  \item{K}{integer that indicates the number of folds.}
  \item{error}{string that indicates the loss to use for
    cross-validation, either \code{"deviance"} or
    \code{"classification"}. The latter is only available for the
    binomial family. Default is \code{"deviance"}.} 
  \item{threshold}{number between 0 and 1 for the decision threshold when  
    \code{"classification"} error is required.}
  \item{folds}{list of \code{K} vectors that describes the folds to use for the
    cross-validation. If \code{NULL}, they are randomly sampled. Default
    is \code{NULL}.} 
  \item{verbose}{logical; indicates if the progression (the current
    fold number) should be displayed. Default is \code{TRUE}.}
}
\value{Returns an object of class \code{cvscoop} for which a \code{plot}
  method is available.
}
\author{Julien Chiquet, inspired by the function \code{cv.glmnet} of the
  \pkg{glmnet} package by   J. Friedman, T. Hastie and R. Tibshirani.
}

\seealso{See also as \code{\link{plot.cvscoop}},
  \code{\link{scoop}}, \code{\linkS4class{cvscoop}},
  \code{\linkS4class{scoopfit}}.}

\examples{
data(basal)
set.seed(47)
coolasso <- coop.lasso(x, y, group, family="binomial", n.lambda=50)
out.cv <- crossval(coolasso, K=5)
plot(out.cv)
}
\keyword{models}
\keyword{regression}
