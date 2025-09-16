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

nloglik.gaussian <- function(XtX, Xty, beta) {
  return(.5*crossprod(beta,crossprod(XtX,beta))-crossprod(beta,Xty))
}

nloglik.binomial <-  function(X, y, beta) {
  eta <- crossprod(t(X), beta)
  return(- sum( y * eta - log(1 + exp(eta))))
}

gradient.gaussian <- function(XtX, Xty, beta) {
  return(drop(crossprod(XtX,beta)) - Xty)
}

gradient.binomial <- function(X, y, beta) {
  return(- crossprod(X,y-1/(1+exp(-crossprod(t(X), beta)))))
}

group.norm <- function(x,g) {
  p <- length(x)
  K <- length(g)
  return(c(.C("groupnorm",
              as.integer(p),
              as.integer(K),
              as.integer(g),
              as.double(x),
              out=as.double(rep(0,K))
              ##,PACKAGE="scoop"
              ))$out)
}

group.norm.rep <- function(x,g) {
  p <- length(x)
  K <- length(g)
  return(c(.C("groupnormrep",
              as.integer(p),
              as.integer(K),
              as.integer(g),
              as.double(x),
              out=as.double(rep(0,p))
              #,PACKAGE="scoop"
              ))$out)
}

coop.norm <- function(x,g) {
  p <- length(x)
  K <- length(g)
  return(c(.C("coopnorm",
              as.integer(p),
              as.integer(K),
              as.integer(g),
              as.double(x),
              out=as.double(rep(0,K))
              #,PACKAGE="scoop"
              ))$out)
}

coop.norm.rep <- function(x,g) {
  return(rep(coop.norm(x,g),g))
}

proximal.tree.group <- function(u, lambda, tpk, Ks) {
  return(c(.C("proximal_tree_grp_standalone",
              as.integer(length(u)),
              as.integer(length(Ks)),
              as.integer(length(tpk)),
              as.integer(Ks),
              as.integer(tpk),
              as.double(lambda),
              as.double(1),              
              out=as.double(u)))$out)
}

proximal.tree.coop <- function(u, lambda, tpk, Ks) {
  return(c(.C("proximal_tree_coo_standalone",
              as.integer(length(u)),
              as.integer(length(Ks)),
              as.integer(length(tpk)),
              as.integer(Ks),
              as.integer(tpk),
              as.double(lambda),
              as.double(1),              
              out=as.double(u)))$out)
}

## TREE-STRUCTURED GROUP LASSO
## Finding the minimal value of lambda which zeroes everything
findLambdaMaxTree <- function (x, y, pk, wk, pen, lambda0 = 100, eps=10e-8) {

  proximal.tree <- switch(pen,
                          "treegroup" = proximal.tree.group,
                          "treecoop"  = proximal.tree.coop)
  
  lambda <- lambda0

  pk <- rev(pk)
  Ks  <- sapply(pk,length)
  pk <- unlist(pk)
  wk  <- unlist(rev(wk))
  
  ## l = - sous-gradient en 0, ie quand lambda est à l'infini = -Xty
  l0 <- crossprod(y,x)
  
  ## si lambda est tel que les betas sont tous nuls
  if (all(proximal.tree(l0,lambda*wk,pk,Ks)==0)) {
    ## alors on l'affecte à lambda2 
    lambda2 <- lambda
    
    ## et on cherche lambda1 = le plus grand lambda tel que betas non tous nuls
    lambda1 <- lambda/2
    l <- proximal.tree(l0,lambda1*wk,pk,Ks)
      
    while (all(l == 0)) {
      lambda1 <- lambda1/2
      l <- proximal.tree(l0,lambda1*wk,pk,Ks)
    }
  } else { ## sinon, on cherche le premier lambda qui annule tous les betas
    ## on sait que lambda n'annule pas tous les betas => on le met dans lambda1
    lambda1 <- lambda
    
    ## on augmente lambda jusqu'à ce qu'il annule tous les betas
    lambda2 <- lambda*2
    l <- proximal.tree(l0,lambda2*wk,pk,Ks)
    
    while (any(l != 0)) {
      lambda2 <- lambda2*2
      l <- proximal.tree(l0,lambda2*wk,pk,Ks)
    }
  }
  
  while (lambda2 - lambda1 > eps) {
    lambda <- (lambda2 + lambda1) / 2
    if (all(proximal.tree(l0,lambda*wk,pk,Ks) == 0)) {
      lambda2 <- lambda
    } else {
      lambda1 <- lambda
    }
  }
  return(lambda)
}
