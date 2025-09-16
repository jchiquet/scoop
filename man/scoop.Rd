\name{scoop}
\alias{scoop}
\alias{lasso}
\alias{group.lasso}
\alias{coop.lasso}
\alias{sparse.group.lasso}
\alias{sparse.coop.lasso}
\alias{tree.group.lasso}
\alias{tree.coop.lasso}
\title{Function to fit a penalized regression problem with various
  grouping effects.}
\description{Fit a penalized regression problem with possible various
  grouping effects. Implements variants such as Lasso, group-Lasso, coop-Lasso, 
  sparse group-Lasso, sparse coop-Lasso, tree group-Lasso and tree
  coop-Lasso penalties.
  
  \code{lasso}, \code{group.lasso}, \code{coop.lasso} and
  \code{sparse.group.lasso}  are aliases of   \code{scoop} for which the
  type of penalty does not need to be specified.
}
\usage{
lasso(x, y, family = "gaussian", intercept = TRUE,
       normalize    = ifelse(family=="gaussian",TRUE,FALSE),
       wk           = rep(1,ncol(x)),
       lambda       = NULL, n.lambda = 100, lambda.min = 1e-2,
       verbose      = FALSE, eps = 1e-5, max.iter = 2*ncol(x),
       optim.method = ifelse(family=="gaussian","fista","bfgs"))

group.lasso(x, y, group, family = "gaussian", intercept = TRUE,
       normalize    = ifelse(family=="gaussian",TRUE,FALSE),
       wk           = sqrt(tabulate(group)),
       lambda       = NULL, n.lambda = 100, lambda.min = 1e-2,
       verbose      = FALSE, eps = 1e-5, max.iter = 2*ncol(x),
       optim.method = ifelse(family=="gaussian","fista","bfgs"))

coop.lasso(x, y, group, family = "gaussian", intercept = TRUE,
       normalize    = ifelse(family=="gaussian",TRUE,FALSE),
       wk           = sqrt(tabulate(group)),
       lambda       = NULL, n.lambda = 100, lambda.min = 1e-2,
       verbose      = FALSE, eps = 1e-5, max.iter = 2*ncol(x),
       optim.method = ifelse(family=="gaussian","fista","bfgs"))

sparse.group.lasso(x, y, group, family = "gaussian", intercept = TRUE,
       normalize    = ifelse(family=="gaussian",TRUE,FALSE),
       wk           = sqrt(tabulate(group)),
       lambda       = NULL, n.lambda = 100, lambda.min = 1e-2,
       verbose      = FALSE, eps = 1e-5, max.iter = 2*ncol(x))

sparse.coop.lasso(x, y, group, family = "gaussian", intercept = TRUE,
       normalize    = ifelse(family=="gaussian",TRUE,FALSE),
       wk           = sqrt(tabulate(group)),
       lambda       = NULL, n.lambda = 100, lambda.min = 1e-2,
       verbose      = FALSE, eps = 1e-5, max.iter = 2*ncol(x))

tree.group.lasso(x, y, group, family = "gaussian", intercept = TRUE,
       normalize    = ifelse(family=="gaussian",TRUE,FALSE),
       wk           = lapply(apply(group,1,tabulate),sqrt),
       lambda       = NULL, n.lambda = 100, lambda.min = 1e-2,
       verbose      = FALSE, eps = 1e-5, max.iter = 2*ncol(x))

tree.coop.lasso(x, y, group, family = "gaussian", intercept = TRUE,
       normalize    = ifelse(family=="gaussian",TRUE,FALSE),
       wk           = lapply(apply(group,1,tabulate),sqrt),
       lambda       = NULL, n.lambda = 100, lambda.min = 1e-2,
       verbose      = FALSE, eps = 1e-5, max.iter = 2*ncol(x))

}

\arguments{
  \item{x}{design matrix (do NOT include intercept).}
  \item{y}{response vector.}
  \item{group}{defines how variables should be grouped. Components
  sharing the same number define a group. Should be a vector for
  funciont \code{lasso}, \code{group.lasso}, \code{coop.lasso},
  \code{sparse.group.lasso} and \code{sparse.coop.lasso}. Should be a
  matrix for \code{tree.group.lasso}, \code{tree.coop.lasso} with as
  much line as there are levels in the tree. See details.
  }
  \item{family}{response type. For now, either \code{"gaussian"} for linear
    regression or \code{"binomial"} for logistic regression. Default is
    \code{"gaussian"}.}
  \item{intercept}{logical; indicates if an intercept should be included in
    the model. Default is \code{TRUE}.}
  \item{normalize}{logical; indicates if variables should be
    normalized to have unit L2 norm before fitting. Default is
    \code{TRUE} for linear regression, \code{FALSE} for logistic
    regression.}
  \item{lambda}{vector of decreasing penalty levels. If \code{NULL}
    (the default), an appropriate vector will be generated with
    \code{n.lambda}  entries, starting from a level of penalty where only
    the intercept is included, and then shrinked to \code{lambda.min}.
  }
  \item{n.lambda}{integer that indicates the number of values to put
    in the \code{lambda} vector.  Ignored if \code{lambda} is provided.
  }

  \item{lambda.min}{minimal value of penalty that will be tried. Default
   is \code{1e-3}. Ignored if \code{lambda} is provided.}
   \item{wk}{an object with real postitive values that weight the
   penalty of each group of features. The default weights each feature
   according to the square root of the size of the group it belongs
   to. Should be a vector with \code{nlevels(group)} elements for
   \code{lasso}, \code{group.lasso}, \code{coop.lasso},
   \code{sparse.group.lasso} and \code{sparse.coop.lasso}. Should be a
   list with as much entries as there are levels in the tree for
   \code{tree.group.lasso}, \code{tree.coop.lasso}. Each element in the
   list contains as many elements as groups at the present level.
 }
 \item{verbose}{logical; indicates verbose mode to display
  progression. Default is \code{FALSE}.} 
  \item{eps}{tolerance used to stop the optimization algorithm. Default
  is \code{1e-4}.}
  \item{max.iter}{maximum number of iterations before the whole optimization
    algorithm stops (for a fixed value of lambda).
  }
  \item{optim.method}{The optimization method used to solve the underlying
    optimization problem on the currently activated set of variables. Can
    be \code{"bfgs"} (a BFGS quasi-Newton method with box constraints),
    \code{"ista"} or   \code{"fista"} (proximal methods). Default is
    \code{"fista"} for linear regression and \code{"bfgs"} for logistic
    regression (the fastest in each category of problem). Only
    \code{"fista"} is available for \code{sparse.group.lasso},
   \code{sparse.coop.lasso}, \code{tree.group.lasso} and
   \code{tree.coop.lasso}.}
}

\details{The general strategy of the algorithm relies on maintaining a
  working set of variable, starting from a vector of zero and entering
  the features by groups of variables. The underlying optimization
  problem is solved only on the activated groups, thus solving small
  smooth problems with increasing size. The working set algorithm stops
  when the optimality conditions are met, up to a tolerance level
  defined by the first element of the \code{eps} argument. This also
  controls the convergence of the underlying, smooth optimization
  problem on the currently activated variable. Various methods are
  implemented to solve this latter problem, which is controlled by
  \code{optim.method}. The Lasso is implemented here and enjoys the
  tools of the \pkg{scoop} package, yet note that it is not as
  performant as in \pkg{lars} or \pkg{glmnet}.
}

\value{Returns an object of class \code{scoopfit} (see the documentation).}
\author{Julien Chiquet}

\seealso{See also \code{\linkS4class{scoopfit}}, \code{\link{plot.scoop}},
  \code{\link{selection}} and \code{\link{crossval}}.}
\examples{
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
coolasso <- coop.lasso(x, y, group, lambda.min=0.1)
scoo     <- selection(coolasso)
## plot the coefficient path with the chosen model
plot(coolasso, crit=scoo$BIC)
}
\keyword{models}
\keyword{regression}
