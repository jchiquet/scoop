.sigmoid <- function(eta) {
  1 / (1 + exp(-eta))
}

.reg_fitting <- function(x, y, intercept, family) {
  
  fit <- NULL  
  if (intercept) {
    my_formula <- y ~ x + 1
  } else {
    my_formula <- y ~ x + 0
  }

  ## OLS in the gaussian case
  if (family == "gaussian") {        
    if (nrow(x) <= ncol(x)) {
      try(fit <- lm.ridge(my_formula, lambda=0), silent=TRUE)
    } else {
      try(fit <- lm(my_formula), silent=TRUE)
    }
  }
  
  ## Logistic in the binomial case
  if (family == "binomial") {    
    if (nrow(x) > ncol(x)) {
      try(fit <- glm(my_formula, family="binomial"), silent=TRUE)
    }
  }

  fit
}

nloglik_gaussian <- function(XtX, Xty, beta) {
  return(.5 * crossprod(beta, crossprod(XtX, beta)) - crossprod(beta, Xty))
}

nloglik_binomial <-  function(X, y, beta) {
  eta <- crossprod(t(X), beta)
  return(- sum( y * eta - log(1 + exp(eta))))
}

gradient_gaussian <- function(XtX, Xty, beta) {
  return(drop(crossprod(XtX,beta)) - Xty)
}

gradient_binomial <- function(X, y, beta) {
  return(- crossprod(X,y - 1 /(1 + exp(-crossprod(t(X), beta)))))
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
