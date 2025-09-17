DataModel <- R6::R6Class(
  classname = "DataModel",
  private = list(
    transform = NA,
    names     = NA
  ),
  public = list(
    ## model-related fields
    data = NA,  
    loss = NA,
    initialize = function(covariates, response, group, intercept, standardize) {

      if (is.null(colnames(x))) colnames(x) <- 1:ncol(x)
      self$transform <- list(intercept = intercept, standardize = standardize)
      self$data      <- list(x = covariates, y = responses, g = group)
      self$pretreatment()    

    },
    pretreatment = function() {
      if (private$transform$intercept) {
        private$names <- c("intercept", colnames(self$data$x))
        x_bar <- colMeans(self$data$x)
        self$data$x <- scale(self$data$x, x_bar, FALSE) 
        attribute(self$data$x, "x_bar") <- x_bar
      } else {
        private$names <- colnames(self$data$x)
      }
      
      ## normalizing the data
      if (private$transform$standardize) {
        norm   <- sqrt(drop(colSums(data$x^2)))
        self$data$x <- scale(self$data$x, FALSE, norm)
        attribute(data$x, "norm") <- norm
      }
      
      ## group labels MUST start from 1
      ## sorting the groups and the columns of the design matrix
      self$data$o <- order(self$data$g,decreasing=FALSE)
      self$data$x <- self$data$x[, self$data$o]
      self$data$g <- self$data$group[self$data$o]
      self$data$s <- tabulate(self$data$g)
    }
  ), 
  active = list(
    d = function() ncol(self$data$x),
    n = function() nrow(self$data$x),
    varnames = function() private$names
  )
)

#' @export
GaussianModel <- R6::R6Class(
  classname = "GaussianModel",
  inherit = DataModel,
  public = list(
    pretreatment = function() {
      super$pretreatment()
      self$data$Xty <- crossprod(X,y)
      if (private$transform$intercept) {
        y_bar <- mean(self$data$y)
        self$data$y <- self$data$y - y_bar
        attribute(self$data$y, "y_bar") <- y_bar
      }
    },
    loss = function(theta) {
      res <- .5 * mean( (selfdata$y - crossprod(self$data$x, theta))^2 )
      attributes(res, "grad") <- drop(crossprod(self$data$XtX, theta)) - self$data$Xty
      res
    }
  )
)

#' @export
BinaryModel <- R6::R6Class(
  classname = "BinaryModel",
  inherit = DataModel,
  public = list(
    loss = function(theta) {
      eta <- crossprod(t(data$X), theta)
      res <- -sum( data$y * eta - log(1 + exp(eta)))
      attributes(res, "grad") <- -crossprod(data$X, data$y - .sigmoid(eta) )
      res
    }
  )
)

