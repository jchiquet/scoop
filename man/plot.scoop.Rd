\name{plot,scoopfit-methods}
\docType{methods}
\alias{plot,scoopfit-method}
\alias{plot.scoop}
\title{Plot method for scoop object}
\description{Produce a plot of the solution path of a \code{scoop} fit.}
\usage{%
\S4method{plot}{scoopfit}(x, y,
            xvar = "lambda", yvar="coefficients",
            main = paste(x@penalty," path (", yvar, ")", sep=""),
            crit = NULL, log.scale = TRUE, labels = NULL,
            col = c(1+switch(yvar,
              group = c(1:nlevels(as.factor(x@group))),
              c(x@group))),
            lty = switch(yvar,
              group = c(1:nlevels(as.factor(x@group))),
              c(x@group)), ...)
}
\arguments{
  \item{x}{output of a \code{scoop} fit (must be of class \code{scoopfit}).}
  \item{y}{used for S4 compatibility.}
  \item{xvar}{variable to plot on the X-axis: either \code{"lambda"},
    \code{"fraction"} or \code{"df"}. Default is set to \code{"lambda"}.}
  \item{yvar}{variable to plot on the Y-axis: either \code{"coefficients"}
    or \code{"group"} for plotting the group-norm. Default is set to \code{"coefficients"}.}
  \item{main}{main title. Default is set to the model name followed by what
    is on the Y-axis.}
  \item{crit}{vector depicting the chosen criterion with the same number of
    entries as the number of \code{lambda} values in \code{x} (typically,
    BIC, AIC, CV-error or test-error). If specified, a vertical line is
    plotted to represent the value of \code{lambda} which minimizes this
    criterion.}
  \item{log.scale}{logical; indicates if a log-scale should be used
    when \code{xvar="lambda"}.}
  \item{labels}{vector indicating the names associated to the plotted
    variables (either the names of the coefficients or the names of the
    groups, depending on \code{yvar}). When specified, the names appear
    on the Y-axis in order to identify each variable in the
    solution path. Remind that the intercept does not count. Default is
    \code{NULL}.}
  \item{col}{specification for the plotting color. Default is to use
    the same color within a group.}
  \item{lty}{specification for the line type. Default is to use
    the same type within a group.}
  \item{\dots}{additional arguments for generic \code{plot}.}
}
\author{Julien Chiquet}
\seealso{See also \code{\link{scoop}}, \code{\linkS4class{scoopfit}}.}
\examples{
set.seed(87)
## true parameters (including intercept)
beta <- c(2, 0.8, 1.1, 0, -1.5, 0.3, -1.0)
## two groups with three components
group <- c(1, 1, 1, 2, 2, 2)
## Labels for coefficients and groups
labels.coef <- c("v11","v12","v13","v21","v22", "v23")
labels.grp  <- c("g1","g2")
## random design and additive noise with sd = 2
n <- 60
p <- 6 
x <- matrix(rnorm(p * n), nrow = n)
y <- cbind(1, x) \%*\% beta + rnorm(n,0,2)
## fit coop-lasso and look for BIC selection capability
coolasso <- coop.lasso(x, y, group, lambda.min=0.1)
scoo     <- selection(coolasso)
## plot the coefficient path with chosen model
par(mfrow=c(1,2))
plot(coolasso, crit=scoo$BIC, labels=labels.coef)
plot(coolasso, crit=scoo$BIC, yvar="group", labels=labels.grp)
}
\keyword{hplot}
