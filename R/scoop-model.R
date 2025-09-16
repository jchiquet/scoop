setClass("scoop-model", representation = representation(
   family        = "character",
   penalty       = "character",
   nloglik       = "function",
   gradient      = "function",
   optim         = "function",
   name          = "character",
   grp.norm      = "ANY",
   pen.norm      = "ANY",
   subgrad.norm  = "ANY",
   subgrad.dual  = "ANY",
   dual.norm     = "ANY")
)

new.scoop.model <- function(family, penalty, optim.method) {
  
  new("scoop-model",
      family    = family,
      penalty   = penalty,

      ## the negative log-likelihood of the data-fitting term
      nloglik   = switch(family,
        "gaussian" = nloglik.gaussian,
        "binomial" = nloglik.binomial),
      ## the gradient of the negative log-likelihood of the data-fitting term
      gradient  = switch(family,
        "gaussian" = gradient.gaussian,
        "binomial" = gradient.binomial),

      ## the function to compute the norm at the group level according to the penalty
      grp.norm = switch(penalty,
        "lasso" = function(x,pk=NULL) {abs(x)},
        "group" = group.norm.rep,
        "coop"  = coop.norm.rep ,
        NULL),

      ## the function to compute penalty term (sum of the weighted group norm)
      pen.norm  = switch(penalty,
        "lasso" = function(x, g, wk) {sum(wk * abs(x))},
        "group" = function(x, g, wk) {sum(wk * group.norm(x,g))},
        "coop"  = function(x, g, wk) {sum(wk * coop.norm(x,g))} , NULL),
      
      ## the function to compute how much a non zero coefficient
      ## is shrunhen in the subgradient term of the penalty (componentwisely)
      subgrad.norm  = switch(penalty,
        "lasso" = function(x,pk=NULL) {abs(x)},
        "group" = group.norm.rep,
        "coop"  = function(x, g) {
          a.pos <- group.norm.rep(pmax.int(0, x),g)
          a.neg <- group.norm.rep(pmax.int(0,-x),g)
          a.pos[x < 0] <- 0
          a.neg[x > 0] <- 0
          return(a.pos + a.neg)
        },
        NULL),

      ## the function to compute the norm in the subdifferential associated to
      ## the  coefficient ( linked to the dual norm - componentwisely)
      subgrad.dual  = switch(penalty,
        "lasso" = function(x,pk=NULL) {abs(x)},
        "group" = group.norm.rep,
        "coop"  = function(x, g) {
          pmax.int(group.norm.rep(pmax.int(0, x),g), 
                   group.norm.rep(pmax.int(0,-x),g))
        },
        NULL),
      
      ## the optimiszation procedure used
      optim     = switch(optim.method,
        "bfgs" = bfgs, "ista" = ista, "fista" = fista),

      ## the model name 
      name      = paste(switch(family,
        "gaussian" = "Linear Regression Model",
        "binomial" = "Logistic Regression Model"),
        switch(penalty,
               "lasso" = "Lasso penalty",
               "group" = "group-Lasso penalty",
               "coop"  = "coop-Lasso penalty",
               "treegroup"  = "tree group-Lasso penalty",
               "treecoop"  = "tree coop-Lasso penalty"), sep = " with "))
}

setGeneric("lambda.grid",
   function(object,x,y,pk,wk,n,lmin){standardGeneric ("lambda.grid")}
)

setMethod("lambda.grid", "scoop-model",
   function(object,x,y,pk,wk,n,lmin) {
     if (object@penalty %in% c("treegroup","treecoop")) {
       lmax <- findLambdaMaxTree(x, y, pk, wk, object@penalty)
     } else {
       lmax <- max(object@grp.norm(crossprod(y,x),pk)/rep.int(wk,pk))
     }
     lambda <- 10^(seq(from=log10(lmax),to=log10(lmin), length=n))
   }
)

setGeneric("optimization",
   function(model, x0, fx, fy, lambda, pk, wk, lower=-Inf, upper=+Inf, L0, eps) {
     standardGeneric ("optimization")
   }
)

setMethod("optimization", "scoop-model",
   function(model, x0, fx, fy, lambda, pk, wk, lower=-Inf, upper=+Inf, L0, eps) {
      return(model@optim(model, x0, fx, fy, lambda, pk, wk, lower, upper, L0, eps))
   }
)

setGeneric("solver",
    function(object, fx, fy, group, lambda, wk, beta0, eps, iter) {
      standardGeneric ("solver")
   }
)

setMethod("solver", "scoop-model",
   function(object, fx, fy, group, lambda, wk, beta0, eps, iter) {
     return(switch(object@penalty,
    "lasso" = working.set(object, fx, fy, group, lambda, wk, beta0, eps, iter),
    "group" = working.set(object, fx, fy, group, lambda, wk, beta0, eps, iter),
    "coop"  = solver.coop (object, fx, fy, group, lambda, wk, beta0, eps, iter),
    "treegroup" = solver.tree (object, fx, fy, group, lambda, wk, beta0, eps, iter),
    "treecoop"  = solver.tree (object, fx, fy, group, lambda, wk, beta0, eps, iter)))
   }
)

