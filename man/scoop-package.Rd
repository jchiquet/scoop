\name{scoop-package}
\alias{scoop-package}
\docType{package}
\title{Sparse cooperative regression}
\description{%
  This package fits coop-Lasso and group-Lasso paths of solutions for
  linear and logistic regressions. The cooperative-Lasso (in
  short \emph{coop-Lasso}, as introduced in Chiquet et al., 2010, 2011)
  may be viewed as a modification of the group-Lasso penalty that
  promotes sign coherence and allows zeros within groups.
  
  The main function to fit a model is \code{scoop}. It produces an
  object of class \code{scoop} for which \code{predict},
  \code{fitted}, \code{deviance}, \code{residuals}, \code{print} and
  \code{plot} methods exist.
  
  Two other methods are included for model selection purpose (that is,
  to choose an appropriate amount of penalization): the \code{selection}
  function estimates AIC and BIC criteria and relies on an approximation
  of the degrees of freedom of a coop-Lasso fit. The \code{crossval}
  function performs cross-validation and produces an object of class
  \code{cvscoop} for which a \code{plot} method exists.
  
  \emph{Note that this is an early release for peer-review purpose.}
}

\section{Demo available}{%
  \describe{
    \item{\code{demo(basal_tumor)}}{%
      Example on the \code{basal} data set that compares logistic
      coop-Lasso fit and logistic group-Lasso fit on a selected number
      of probe sets. Cross-validation is performed and selected genes
      depicted.  
    }  
  }
}

\author{Julien Chiquet, <julien.chiquet@genopole.cnrs.fr> (also maintainer).}

\references{Julien Chiquet, Yves Grandvalet and Camille Charbonnier
  (2011 - preprint). Sparsity with sign-coherent groups of variables via
  the cooperative-Lasso.
  
  Julien Chiquet, Yves Grandvalet and Christophe Ambroise
  (2010). Inferring  multiple graphical structures, Statistics
  and Computing, \url{http://dx.doi.org/10.1007/s11222-010-9191-2}.
}
\keyword{package}
