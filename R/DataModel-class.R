DataModel <- R6::R6Class(
  classname = "DataModel",
  private = list(
    names     = NA
  ),
  public = list(
    ## model-related fields
    data = NA,
    initialize = function(covariates, outcome, group, intercept=TRUE, standardize=TRUE) {

      if (is.null(colnames(covariates))) colnames(covariates) <- 1:ncol(covariates)
      self$data  <- list(X = covariates, y = outcome, g = group)

      ## group labels MUST start from 1
      ## sorting the groups and the columns of the design matrix
      self$data$o <- order(self$data$g,decreasing=FALSE)
      self$data$X <- self$data$X[, self$data$o]
      self$data$g <- self$data$g[self$data$o]
      self$data$s <- tabulate(self$data$g)
      
      if (intercept) {
        private$names <- c("intercept", colnames(self$data$X))
        x_bar <- colMeans(self$data$X)
        self$data$X <- scale(self$data$X, x_bar, FALSE)
        attr(self$data$X, "x_bar") <- x_bar
      } else {
        private$names <- colnames(self$data$X)
      }
      
      ## normalizing the data
      if (standardize) {
        norm   <- sqrt(drop(colSums(self$data$X^2)))
        self$data$X <- scale(self$data$X, FALSE, norm)
        attr(self$data$X, "norm") <- norm
      }

    }
  ), 
  active = list(
    d = function() ncol(self$data$X),
    n = function() nrow(self$data$X),
    k = function() length(self$data$s),
    has_intercept = function() !is.null(attr(self$data$X, "x_bar")),
    is_standardized = function() !is.null(attr(self$data$X, "norm")),
    varnames = function() private$names
  )
)

#' @export
GaussianModel <- R6::R6Class(
  classname = "GaussianModel",
  inherit = DataModel,
  public = list(
    initialize = function(covariates, outcome, group, 
                          intercept = TRUE, standardize = TRUE) {
      super$initialize(covariates, outcome, group, intercept, standardize)
      if (intercept) {
        y_bar <- mean(self$data$y)
        self$data$y <- self$data$y - y_bar
        attr(self$data$y, "y_bar") <- y_bar
      }
    },
    loss = function(theta) {
      y_hat <- self$data$X %*% theta
      res <- .5 * mean( (self$data$y - y_hat)^2 )
      attr(res, "grad") <- crossprod(self$data$X, y_hat - self$data$y)
      res
    }
  ),
  active = list(
    name = function() "Gaussian response (Linear Regression)"
  )
)

#' @export
BinaryModel <- R6::R6Class(
  classname = "BinaryModel",
  inherit = DataModel,
  public = list(
    loss = function(theta) {
      eta <- self$data$X %*% theta
      res <- sum(log(1 + exp(eta)) - self$data$y * eta)
      attr(res, "grad") <- crossprod(self$data$X, .sigmoid(eta) - self$data$y)
      res
    }
  ),
  active = list(
    name = function() "Binary response (Logistic Regression)"
  )
)
