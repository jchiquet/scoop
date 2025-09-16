\name{plot,cvscoop-methods}
\docType{methods}
\alias{plot,cvscoop-method}
\alias{plot.cvscoop}
\title{Plot method for cross validated error of a scoop model}
\description{Produce a plot of the cross validated error of a \code{scoop} model.}
\usage{
\S4method{plot}{cvscoop}(x, y,
        log.scale = TRUE,
        xlab = ifelse(log.scale, "lambda (log scale)", "lambda"),
        ylab = "Mean CV error", ...)}

\arguments{
  \item{x}{output of a \code{crossval} run (must be of class
  \code{cvscoop}).}
  \item{y}{used for S4 compatibility}
  \item{xlab}{title for the x axis. If \code{NULL}, will be set to
  \code{"lambda"} or \code{"lambda (log-scale)"} according to \code{log.scale}.}
  \item{ylab}{title for the y axis}
  \item{log.scale}{logical; indicates if a log-scale should be used.}
  \item{\dots}{additional arguments for generic plot.}
}

\author{Julien Chiquet}
\seealso{See also \code{\linkS4class{cvscoop}}, \code{\link{crossval.scoop}}.}
\examples{
data(basal)
set.seed(47)
cat("\nmight take about half a minute...")
coolasso <- coop.lasso(x, y, group, family="binomial", n.lambda=50)
out.cv <- crossval(coolasso, K=5)
plot(out.cv)
}
\keyword{hplot}
