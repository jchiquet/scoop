Loss <- R6::R6Class(
  classname = "Loss",
  public = list(
     ## model-related fields
     value = function() {}
  )
)

MSELoss <- R6::R6Class(
  classname = "MSELoss",
  inherit = Loss,
  public = list(
    initialize = function(data) {
      self$value <- function(param) {
        res <- .5 * sum( (data$y - crossprod(data$X, param))^2 )
        attributes(res, "grad") <- drop(crossprod(data$XtX, beta)) - data$Xty
      }
    }
  )
)

LogisticLoss <- R6::R6Class(
  classname = "LogisticLoss",
  inherit = Loss,
  public = list(
    initialize = function(data) {
      self$func <- function(param) {
        eta <- crossprod(t(data$X), param)
        res <- - sum( data$y * eta - log(1 + exp(eta)))
        attributes(res, "grad") <- - crossprod(data$X, data$y - .sigmoid(eta) )
        res
      }
    }
  )
)

