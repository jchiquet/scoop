Penalty <- R6::R6Class(
  classname = "Penalty",
  public = list(
    initialize = function(group, weights) {
      private$g <- tabulate(group)
      private$w <- weights
    },
    ## function to compute penalty term (sum of the weighted group norm)
    norm = function() {},
    ## function to compute the norm at the group/element level
    elt_norm = function() {},
    ## function to compute the norm in the subdifferential associated to
    ## the  coefficient ( linked to the dual norm - componentwisely)
    dual_norm = function() {},
    ## function to compute the proximal operator of the current penalty
    proximal = function() {}
  ),
  private = list(
    w = NA,
    g = NA
  ), 
  active = list(
    name = function() {private$id},
    weights = function() {private$w},
    group = function() {private$g}
  )
)

#' @export
Penalty_Lasso <- R6::R6Class(
  classname = "Penalty_Lasso",
  inherit = Penalty,
  public = list(
    initialize = function(weights) {
      super$initialize(group = 1:length(weights), weights)
      private$id <- "L1 (Lasso)"
    },
    norm      = function(theta) pen_norm_L1(theta, private$w),
    elt_norm  = function(theta) elt_norm_L1(theta),
    dual_norm = function(theta) dual_norm_L1(theta),
    proximal  = function(theta, lambda) proximal_L1(theta, private$w, lambda)
  )
)

#' @export
Penalty_Bounded <- R6::R6Class(
  classname = "Penalty_Bounded",
  inherit = Penalty,
  public = list(
    initialize = function(weights) {
      super$initialize(group = 1:length(weights), weights)
      private$id <- "L-infty (Bounded)"
    },
    norm      = function(theta) pen_norm_LINF(theta, private$w),
    elt_norm  = function(theta) elt_norm_LINF(theta),
    dual_norm = function(theta) dual_norm_LINF(theta),
    proximal  = function(theta, lambda) proximal_LINF(theta, private$w, lambda)
  )
)

#' @export    
Penalty_GroupLasso <- R6::R6Class(
  classname = "Penalty_GroupLasso",
  inherit = Penalty,
  public = list(
    initialize = function(group, weights) {
      super$initialize(group = group, weights)
      private$id <- "L1/L2 (Standard Group-Lasso)"
    },
    norm      = function(theta) pen_norm_L1L2(theta, private$g, private$w),
    elt_norm  = function(theta) elt_norm_L1L2(theta, private$g),
    dual_norm = function(theta) dual_norm_L1L2(theta, private$g),
    proximal  = function(theta, lambda) proximal_L1L2(theta, private$g, private$w, lambda)
  )
)

#' @export
Penalty_GroupLassoInf <- R6::R6Class(
  classname = "Penalty_GroupLassoInf",
  inherit = Penalty,
  public = list(
    initialize = function(group, weights) {
      super$initialize(group = group, weights)
      private$id <- "L1/LInfty (Variant of Group-Lasso)"
    },
    norm      = function(theta) pen_norm_L1LINF(theta, private$g, private$w),
    elt_norm  = function(theta) elt_norm_L1LINF(theta, private$g),
    dual_norm = function(theta) dual_norm_L1LINF(theta, private$g),
    proximal  = function(theta, lambda) proximal_L1LINF(theta, private$g, private$w, lambda)
  )
)

#' @export
Penalty_CoopLasso <- R6::R6Class(
  classname = "Penalty_CoopLasso",
  inherit = Penalty,
  public = list(
    initialize = function(group, weights) {
      super$initialize(group = group, weights)
      private$id <- "L1/Signed-L2 (Cooperative-Lasso)"
    },
    norm      = function(theta) pen_norm_COOP(theta, private$g, private$w),
    elt_norm  = function(theta) elt_norm_COOP(theta, private$g),
    dual_norm = function(theta) dual_norm_COOP(theta, private$g),
    proximal  = function(theta, lambda) proximal_COOP(theta, private$g, private$w, lambda)    
  )
)
