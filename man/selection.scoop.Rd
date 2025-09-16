\name{selection,scoopfit-methods}
\docType{methods}
\alias{selection,scoopfit-method}
\alias{selection}
\title{Model selection criteria for scoop models}
\description{Compute estimations of AIC and BIC of a
  \code{scoop} fit for model selection purpose. Only available for the
  Gaussian family (regular linear regression).}
\usage{\S4method{selection}{scoopfit}(object, sigma2 = NULL)}
\arguments{
  \item{object}{output of a \code{scoop} run (must be an object of
    class \code{scoopfit}) 
  }  
  \item{sigma2}{user defined estimator for the variance. If not
    specified, use the OLS estimate (only relevant when
    \if{latex}{\eqn{n>p}}\if{html}{n>p}\if{text}{n>p}). If
    \if{latex}{\eqn{n<p}}\if{html}{n<p}\if{text}{n<p}, the user should
    provide \code{sigma2} otherwise \code{BIC} and \code{AIC} will be
    \code{NULL}.}
}

\details{The following expression are used to estimate the
  BIC (\emph{Bayesian   Information Criterion})
 \if{latex}{\deqn{%
     \widetilde{\textrm{BIC}}(\lambda) =
  \frac{\mathrm{RSS}(\lambda)}{\hat{\sigma}^2} +      \log(n)
  \tilde{\mathrm{df}}(\lambda),
 }}

 \if{html}{\out{
<center>
    BIC<sub>tilde</sub>(&lambda;) =
    RSS((&lambda;)/&sigma;<sub>hat</sub><sup>2</sup> - 
    df<sub>tilde</sub>(&lambda;) log(n),
</center>
}}
\if{text}{BIC.tilde(lambda) = RSS/sigma.hat^2 - log(n) df.tilde(lambda),}

and the AIC (\emph{Akaike Information Criterion})
\if{latex}{\deqn{%
    \widetilde{\textrm{AIC}}(\lambda) =
    \frac{\mathrm{RSS}(\lambda)}{\hat{\sigma}^2} +       2
    \tilde{\mathrm{df}}(\lambda),
  }}  
  \if{html}{\out{
<center>
    AIC<sub>tilde</sub>(&lambda;) =
    RSS((&lambda;)/&sigma;<sub>hat</sub><sup>2</sup> -     2
    df<sub>tilde</sub>(&lambda;),
</center>
}
}
\if{text}{AIC.tilde(lambda) = RSS/sigma.hat^2 - log(n) df.tilde(lambda),}

where
\if{latex}{\eqn{\tilde{\mathrm{df}}}}\if{html}{\out{df<sub>tilde</sub>(&lambda;)}}
\if{text}{df.tilde(lambda)} is an approximation of the degrees of
freedom as described in Yuan and Lin (2006) and Chiquet et al
(submitted). 
}

\value{Returns a list comprising
\item{BIC}{vector containing the values of the BIC for the successive
  values of
  \if{latex}{\eqn{\lambda}}\if{html}{\out{&lambda;}}\if{text}{lambda} used
  along the \code{scoop} fit.} 
\item{beta.BIC}{vector of parameters chosen by minimizing  BIC.}
\item{lambda.BIC}{amount of penalization that minimizes BIC.}
\item{AIC}{vector containing the values of the AIC for the succesive
  values of
  \if{latex}{\eqn{\lambda}}\if{html}{\out{&lambda;}}\if{text}{lambda} used
  along the \code{scoop} fit.} 
\item{beta.AIC}{vector of parameters chosen by minimizing  AIC.}
\item{lambda.AIC}{amount of penalization that minimizes BIC.}
\item{sigma2}{estimated variance.}
\item{df}{vector containing the degrees of freedom for the successive
  values of
  \if{latex}{\eqn{\lambda}}\if{html}{\out{&lambda;}}\if{text}{lambda} used
  along the \code{scoop} fit. }
\item{non.zeros}{vector containing the number of non-zero coefficients
  for the successive  values of
  \if{latex}{\eqn{\lambda}}\if{html}{\out{&lambda;}}\if{text}{lambda} used
  along the \code{scoop}   fit.}
}
\author{Julien Chiquet}
\references{Julien Chiquet, Yves Grandvalet and Camille Charbonnier
  (2011 - preprint). Promotting sign coherence in groups of variables with the
  coop-Lasso, submitted.
  
  M. Yuan and Y. Lin (2006). Model selection and estimation in
  regression with grouped variables, JRSS B, 8, Part 1, pp.49-67.
}

\seealso{See also \code{\link{scoop}}, \code{\linkS4class{scoopfit}}.}
\examples{
\dontrun{
set.seed(87)
## true parameters (including intercept)
beta <- c(2, 0.8, 1.1, 0, -1.5, 0.3, -1.0)
## two groups with three components
group <- c(1, 1, 1, 2, 2, 2)
## random design and additive noise with sd = 2
n <- 60
p <- 6 
x <- matrix(rnorm(p * n), nrow = n)
y <- cbind(1, x) \%*\% beta + rnorm(n,0,2)
## fit coop-lasso and look for BIC selection capability
coolasso <- scoop(x, y, group, lambda.min=0.1)
scoo     <- selection(coolasso)
## BIC choice
scoo$beta.BIC}
}
\keyword{models}
\keyword{regression}
