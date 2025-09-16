\name{scoopfit-class}
\Rdversion{1.1}
\docType{class}
\alias{scoopfit-class}
\alias{deviance,scoopfit-method}
\alias{fitted,scoopfit-method}
\alias{predict,scoopfit-method}
\alias{print,scoopfit-method}
\alias{residuals,scoopfit-method}

\title{Class "scoopfit"}
\description{Class of object returned by the \code{scoop} fitting function.}
\section{Objects from the Class}{
  Such objects can be created by calls to \code{scoop}, \code{lasso},
  \code{group.lasso} and \code{coop.lasso}.}

\section{Slots}{
  \describe{
    \item{\code{coefficients}:}{matrix of coefficients with respect to
      the original input (even when normalized). The number of rows of
      \code{coefficients} corresponds the length of \code{lambda}.}
    
    \item{\code{lambda}:}{vector of penalty levels for which the model
      has eventually been fitted.}
    
    \item{\code{wk}:}{vector or lsit which defines how each group is
      weighted regarding the penalizer.}

    \item{\code{group}:}{vector or matrix which defines the grouping of the
      variables.}

    \item{\code{family}:}{response type used. For now, either
      \code{"gaussian"} for linear regression or \code{"binomial"} for
      logistic regression. Default is \code{"gaussian"}.}

    \item{\code{penalty}:}{type of grouping effect used for the
      penalizer: either \code{"coop"}, \code{"group"} or \code{lasso}.}

    \item{\code{monitoring}:}{list wich contains three vectors to
      monitor the optimization process: \code{it.active} gives the
      number of iterations in the active set algorithm for each value of
      lambda; \code{it.optim} gives the number of iterations for each
      call to the optimization function; \code{dual.gap} gives the
      duality gap reached for each lambda.}
    
    \item{\code{x}:}{design matrix saved for further use.}

    \item{\code{y}:}{response vector saved for further use.}

    \item{\code{call}:}{original call to the \code{scoop.default}
      function.}
  }
}

\section{Methods}{This class comes with the usual
  \code{fitted(object, ...)}, \code{predict(object, newx=NULL, ...)},
  \code{residuals(object, ...)}, \code{print(object, ...)} and
  \code{deviance(object, ...)} generic methods.

  Three specific and documented methods are also available: 
  
  \describe{
    \item{crossval}{cross-validates a \code{scoop} object}
    \item{plot}{plots a \code{scoop} object}
    \item{selection}{computes approximated BIC/AIC}
  }
}
\author{Julien Chiquet}

\seealso{See Also \code{\link{scoop}}, \code{link{plot.scoop}}, \code{link{crossval}}, \code{link{selection}}.}

\examples{
showClass("scoopfit")
}
\keyword{classes}
