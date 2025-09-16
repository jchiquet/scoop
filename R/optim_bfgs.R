bfgs <- function(model, x0, fx, fy, lambda, pk, wk, lower=-Inf, upper=+Inf, L0, eps) {

  L  <- function(beta, fx, fy, lambda, pk, wk) {
    model@nloglik(fx,fy,beta) + lambda * model@pen.norm(beta,pk,wk)
  }
  
  dL <- function(beta, fx, fy, lambda, pk, wk) {
    
    dL      <- rep.int(0,length(beta))
    nabla.f <- model@gradient(fx,fy,beta)
    lambda  <- lambda * rep.int(wk,pk)

    n.beta <- model@grp.norm(beta,pk)
    zero <- n.beta == 0
    
    ## terms whose norm1 is different from zero
    dL[!zero] <- nabla.f[!zero] + lambda[!zero] * beta[!zero]/n.beta[!zero]
    
    ## terms whose norm2 is zero
    n.nabla.f <- model@grp.norm(nabla.f,pk)
    shrunken  <- zero & n.nabla.f > lambda

    ## shrink toward zero
    shrink <- lambda[shrunken]/n.nabla.f[shrunken]
    dL[shrunken] <- nabla.f[shrunken] * (1 - shrink)

    return(dL) 
  }

  dL.coop <- function(beta,fx,fy,lambda,pk,wk) {
    
    null.vector <- rep.int(0,length(beta))
    dL <- null.vector
    lambda <- lambda * rep.int(wk,pk)
    nabla.f <- model@gradient(fx,fy,beta)

    a.pos <- group.norm.rep(pmax.int(0, beta),pk)
    a.neg <- group.norm.rep(pmax.int(0,-beta),pk)
    
    ## manage cases where alpha+ > 0 or alpha- < 0
    neg <- beta < 0
    pos <- beta > 0
    dL[pos] <- nabla.f[pos] + lambda[pos]*pmax.int(null.vector[pos], beta[pos]) / a.pos[pos]
    dL[neg] <- nabla.f[neg] - lambda[neg]*pmax.int(null.vector[neg],-beta[neg]) / a.neg[neg]
    
    ## manage cases where beta = 0 and (alpha- = 0 and/or alpha+ = 0)
    n2.theta.pos <- group.norm.rep(pmax.int(0, nabla.f),pk)
    n2.theta.neg <- group.norm.rep(pmax.int(0,-nabla.f),pk)
    
    pos.shrunken  <- abs(beta) == 0 & nabla.f >= 0 & n2.theta.pos > lambda
    neg.shrunken  <- abs(beta) == 0 & nabla.f <= 0 & n2.theta.neg > lambda
    
    shrink.pos <- lambda[pos.shrunken]/n2.theta.pos[pos.shrunken]
    shrink.neg <- lambda[neg.shrunken]/n2.theta.neg[neg.shrunken]

    dL[pos.shrunken] <- nabla.f[pos.shrunken] * (1-shrink.pos)
    dL[neg.shrunken] <- nabla.f[neg.shrunken] * (1-shrink.neg)
    
    ##    dL[zero.pos][pos.to.zero] <- 0  ## redondant
    ##    dL[zero.neg][neg.to.zero] <- 0  ## redondant
    
    return(dL)
  }
  
  res <- optim(x0, method="L-BFGS-B",
               fn     = L,
               gr     = switch(model@penalty, "lasso"=dL, "group"=dL, "coop"=dL.coop, NULL),
               fx     = fx,
               fy     = fy,
               lambda = lambda,
               pk     = pk,
               wk     = wk,
               lower  = lower,
               upper  = upper,
               control= list(pgtol=eps/length(x0), factr=eps*length(x0), maxit=50))
  i <- res$counts[1]
  names(i) <- NULL
  return(list(i=i,xk=res$par))
}


## bfgs <- function(model, x0, fx, fy, lambda, pk, wk, lower=-Inf, upper=+Inf, L0, eps) {

##   L  <- function(beta, fx, fy, lambda, pk, wk) {
##     model@nloglik(fx,fy,beta) + lambda * model@pen.norm(beta,pk,wk)
##   }
  
##   dL <- function(beta, fx, fy, lambda, pk, wk) {
    
##     ## This manages case where the term goes to zero
##     dL      <- rep.int(0,length(beta))
##     nabla.f <- model@gradient(fx,fy,beta)
##     lambda  <- lambda * rep.int(wk,pk)

##     ## Which are the non null guys
##     phi.beta <- model@subgrad.norm(beta,pk)
##     zero <- phi.beta == 0
    
##     ## terms whose norm is different from zero
##     dL[!zero] <- nabla.f[!zero] + lambda[!zero] * beta[!zero]/phi.beta[!zero]
    
##     ## terms shrunken toward zero
##     phi.nabla.f <- model@subgrad.norm(nabla.f,pk)
##     ##    shrunken    <- zero & n.nabla.f > lambda
##     shrunken    <-  zero & phi.nabla.f > lambda

##     dL[shrunken] <- nabla.f[shrunken] * (1 - lambda[shrunken]/phi.nabla.f[shrunken])
    
##     return(dL) 
##   }
  
##   res <- optim(x0, method="L-BFGS-B",
##                fn     = L,
##                gr     = dL,
##                fx     = fx,
##                fy     = fy,
##                lambda = lambda,
##                pk     = pk,
##                wk     = wk,
##                lower  = lower,
##                upper  = upper,
##                control= list(pgtol=eps/10,factr=eps/10))
##   i <- res$counts[1]
##   names(i) <- NULL
##   return(list(i=i,xk=res$par))
## }
