Penalty <- R6::R6Class(
  classname = "Penalty",
  public = list(
    initialize = function(group, weights) {
      private$g <- group
      private$w <- weights
    },
    ## the function to compute penalty term (sum of the weighted group norm)
    value = function() {},
    ## the function to compute the norm at the group level according to the penalty
    group_norm = function() {},
    ## the function to compute the norm in the subdifferential associated to
    ## the  coefficient ( linked to the dual norm - componentwisely)
    dual_norm = function() {},
    ## the function to compute how much a non zero coefficient
    ## is shrunhen in the subgradient term of the penalty (componentwisely)
    subgrad_norm = function() {}
  ),
  private = list(
    w = NA,
    g = NA
  ), 
  active = list(
    weights = function() {private$w},
    group = function() {private$g}
  )
)

Penalty_Lasso <- R6::R6Class(
  classname = "Penalty_Lasso",
  inherit = Penalty,
  public = list(
    value        = function(theta) sum(private$w * abs(theta)),
    group_norm   = function(theta) abs(theta),
    dual_norm    = function(theta) abs(theta),
    subgrad_norm = function(theta) abs(theta)  
  )
)
    
Penalty_GroupLasso <- R6::R6Class(
  classname = "Penalty_GroupLasso",
  inherit = Penalty,
  public = list(
    value        = function(theta) sum(private$w * group.norm(theta, private$g)),
    group_norm   = group.norm.rep,
    dual_norm    = group.norm.rep,
    subgrad_norm = group.norm.rep
  )
)

Penalty_CoopLasso <- R6::R6Class(
  classname = "Penalty_CoopLasso",
  inherit = Penalty,
  public = list(
    value        = function(theta) sum(private$w * coop.norm(theta, private$g)),
    group_norm   = coop.norm.rep,
    dual_norm    = function(theta) {
          pmax.int(group.norm.rep(pmax.int(0, theta),private$g), 
                   group.norm.rep(pmax.int(0,-theta),private$g))
      },
    subgrad_norm = function(theta) {
          a.pos <- group.norm.rep(pmax.int(0, theta), private$g)
          a.neg <- group.norm.rep(pmax.int(0,-theta), private$g)
          a.pos[x < 0] <- 0
          a.neg[x > 0] <- 0
          a.pos + a.neg
      }
  )
)
