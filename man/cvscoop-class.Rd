\name{cvscoop-class}
\Rdversion{1.1}
\docType{class}
\alias{cvscoop-class}

\title{Class "cvscoop"}
\description{Class of object returned by the \code{crossval} method.}
\section{Objects from the Class}{
  Objects can be created by calls to \code{crossval}.}
  
\section{Slots}{
  \describe{

    \item{\code{lambda}:}{vector of penalty levels for which each
      cross-validation has been performed.}

    \item{\code{lambda.min}:}{level of penalty that minimizes the mean
      cross-validated error.}

    \item{\code{lambda.1se}:}{largest level of penalty such has the
      cross-validated error is within 1 standard error of the minimum.}
    
    \item{\code{cv.mean}:}{vector with the same length as \code{lambda} containing
      the mean cross-validated error.}

    \item{\code{cv.error}:}{vector with the same length as \code{lambda} containing
      the estimated standard deviation of \code{cv.mean}.}
    
    \item{\code{folds}:}{list of \code{K} vector indicating the folds used for
      cross-validation.}

    \item{\code{beta.min}:}{vector of parameters corresponding to \code{lambda.min}.}
    
    \item{\code{beta.1se}:}{vector of parameters corresponding to \code{lambda.1se}.}
}
}
\section{Methods}{A \code{plot} method is available and documented.}

\author{Julien Chiquet}

\seealso{See also \code{\link{plot.cvscoop}}.}
\examples{
showClass("cvscoop")
}
\keyword{classes}
