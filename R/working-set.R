working.set <- function(model, fx, fy, group, lambda, wk, beta, eps, max.it) {
  ## _____________________________________________________________
  ##
  ## INITIALIZATION
  ## _____________________________________________________________
  ##
  iter <- c()
  pk  <- tabulate(group)
  rwk <- rep.int(wk,pk)
  active <- model@grp.norm(beta, pk) > 0
  ## _____________________________________________________________
  ## 
  ## CHECK IF OPTIMALITY IS REACHED AT THE STARTING POINT
  ## _____________________________________________________________
  ##
  nabla.f   <- model@gradient(fx,fy,beta)
  n.nabla.f <- model@subgrad.dual(nabla.f,pk)
  dual.all  <- n.nabla.f - lambda * rwk
  gap.old   <- max(0,dual.all)
  if (gap.old < eps) {
    return(list(beta=beta, iter=iter, gap=gap.old, status="converged"))
  } else {
    ## UPDATE THE ACTIVE SET
    new <- which(dual.all == gap.old)
    active[new] <- TRUE
  }
  ##
  while (1) {
    ## _____________________________________________________________
    ##
    ## (1) OPTIMIZATION OVER THE CURRENTLY ACTIVATED VARIABLES
    ## _____________________________________________________________
    ##
    fx.active <- switch (model@family,
                         "gaussian" = fx[active,active],
                         "binomial" = fx[, active])
    fy.active <- switch (model@family,
                         "gaussian" = fy[active],
                         "binomial" = fy)
    L0 <- switch (model@family,
                  "gaussian" = max(eigen(fx.active)$values,TRUE,TRUE),
                  "binomial" = ifelse(is.null(dim(fx.active)),
                    sum(fx.active^2),sum(rowSums(fx.active^2))))                  
    ## Get the currently activated group of variables
    pk.active <- tabulate(group[active])
    wk.active <- wk[which(pk.active != 0)]
    pk.active <- pk.active[pk.active != 0]
    ##
    res <- optimization(model, beta[active], fx.active, fy.active, lambda,
                        pk.active, wk.active, L0=L0, eps=eps)
    ##L0 <- res$L
    iter <- c(iter,res$i)
    ## update the appropriate variables on the active set
    beta[active]  <- res$xk
    nabla.f[active]   <- model@gradient(fx.active, fy.active, beta[active])
    n.nabla.f[active] <- model@grp.norm(nabla.f[active], pk.active)
    ## _____________________________________________________________
    ##
    ## (2) GROUP DELETION IF APPLICABLE
    ## _____________________________________________________________
    ##
    zeroed <- which(active)[model@grp.norm(beta[active],pk.active) < eps &
                            n.nabla.f[active] < lambda * rwk[active] + eps]
    if (length(zeroed) > 0) {
      active[zeroed]  <- FALSE
      beta[zeroed]      <- 0
      fx.zeroed <- switch (model@family,
                           "gaussian" = fx[zeroed,zeroed],
                           "binomial" = fx[, zeroed])
      fy.zeroed <- switch (model@family,
                           "gaussian" = fy[zeroed],
                           "binomial" = fy)
      nabla.f[zeroed]   <- model@gradient(fx.zeroed, fy.zeroed, beta[zeroed])
      n.nabla.f[zeroed] <- lambda * rwk[zeroed]
    }
    ## _____________________________________________________________
    ##
    ## (3) OPTIMALITY TESTING AND GROUP ACTIVATION IF APPLICABLE
    ## _____________________________________________________________
    ##
    dual.all <- n.nabla.f - lambda * rwk
    gap      <- max(0,dual.all)
    if (length(iter) == max.it) {
      return(list(beta=beta, iter=iter, gap=gap, status="max # of iterates reached"))
    }
    if (abs(gap - gap.old) < eps) {
      return(list(beta=beta, iter=iter, gap=gap, status="no further improvement..."))
    }
    if (gap < eps) {
      return(list(beta=beta, iter=iter, gap=gap, status="converged"))
    } else {
      ## UPDATE THE ACTIVE SET / ACTIVE GROUP
      gap.old <- gap
      new <- which(dual.all == gap.old)
      active[new] <- TRUE
    }
  }
}

solver.coop <- function(model, fx, fy, group, lambda, wk, beta, eps, max.it){
  ## _____________________________________________________________
  ##
  ## INITIALIZATION
  ## _____________________________________________________________
  ##
  iter <- c()
  pk  <- tabulate(group)
  rwk <- rep.int(wk,pk)
  p <- length(beta)
  active.pos <- group.norm.rep(pmax.int(0, beta),pk) > 0
  active.neg <- group.norm.rep(pmax.int(0,-beta),pk) > 0
  ## _____________________________________________________________
  ## 
  ## CHECK IF OPTIMALITY IS REACHED AT THE STARTING POINT
  ## _____________________________________________________________
  ##
  nabla.f <- model@gradient(fx, fy, beta)
  n2.nabla.f.pos  <- group.norm.rep(pmax.int(0, nabla.f),pk)
  n2.nabla.f.neg  <- group.norm.rep(pmax.int(0,-nabla.f),pk)
  n.nabla.f <- model@subgrad.dual(nabla.f,pk)
  dual.all <- n.nabla.f - lambda * rwk
  gap.old  <- max(0,dual.all)
  if (gap.old < eps) {
    return(list(beta = beta, iter=iter, gap=gap.old, status="converged"))
  } else {
    ## UPDATE THE ACTIVE SET / ACTIVE GROUP
    new <- which(dual.all == gap.old)
    if (length(new) > 0) {
      if (n2.nabla.f.pos[new][1] > n2.nabla.f.neg[new][1]) {
        active.neg[new] <- TRUE
      } else {
        active.pos[new] <- TRUE
      }
    }
  }
  
  while (1) {
    
    ## (1) OPTIMIZATION OVER A
    ## _____________________________________________________________
    lower.bound <- rep(-Inf,p)
    upper.bound <- rep(+Inf,p)
    lower.bound[active.pos & !active.neg] <- 0
    upper.bound[active.neg & !active.pos] <- 0
    active <- active.pos | active.neg

    fx.active <- switch (model@family,
                         "gaussian" = fx[active,active],
                         "binomial" = fx[, active])
    fy.active <- switch (model@family,
                         "gaussian" = fy[active],
                         "binomial" = fy)
    L0 <- switch (model@family,
                  "gaussian" = max(eigen(fx.active)$values,TRUE,TRUE),
                  "binomial" = ifelse(is.null(dim(fx.active)),
                    sum(fx.active^2),sum(rowSums(fx.active^2))))
    pk.active <- tabulate(group[active])
    wk.active <- wk[which(pk.active != 0)]
    pk.active <- pk.active[pk.active != 0]

    res <- optimization(model, beta[active], fx.active, fy.active, lambda,
                        pk.active, wk.active,
                        lower = lower.bound[active],upper=upper.bound[active],
                        L0=L0,eps=eps)
    L0 <- res$L
    iter <- c(iter,res$i)
    beta[active] <- res$xk
    nabla.f[active] <- model@gradient(fx.active, fy.active, beta[active])
    n.nabla.f[active] <- model@subgrad.dual(nabla.f[active], pk.active)
    
    ## (2) TEST GROUP DELETION IF APPLICABLE
    ## _____________________________________________________________
    go.zero.pos <- intersect(which(group.norm(pmax.int(0,beta),pk)<eps),group[active.pos])
    if (length(go.zero.pos) > 0) {
      zeroed <- group %in% go.zero.pos & group.norm.rep(pmax.int(0,-nabla.f),pk) < lambda * rwk + eps & beta > 0
      active.pos[zeroed] <- FALSE
      beta[zeroed]       <- 0
      nabla.f[zeroed] <- lambda + rwk[zeroed]
    }

    go.zero.neg <- intersect(which(group.norm(pmax.int(0,-beta),pk) < eps),group[active.neg])
    if (length(go.zero.neg) > 0) {
      zeroed <- group %in% go.zero.neg & group.norm.rep(pmax.int(0,nabla.f),pk) < lambda * rwk + eps & beta<0
      active.neg[zeroed] <- FALSE
      beta[zeroed]       <- 0
      nabla.f[zeroed] <- lambda + rwk[zeroed]
    }
    
    ## (3) OPTIMALITY TESTING AND GROUP ACTIVATION IF APPLICABLE
    ## _____________________________________________________________
    ## Get the entry for which the constraint is the most violated
    n2.nabla.f.pos  <- group.norm.rep(pmax.int(0, nabla.f),pk)
    n2.nabla.f.neg  <- group.norm.rep(pmax.int(0,-nabla.f),pk)
    dual.all <- n.nabla.f - lambda * rwk
    gap <- max(0,dual.all)
    if (length(iter) == max.it) {
      return(list(beta=beta, iter=iter, gap=gap, status="max # of iterates reached"))
    }
    if (abs(gap - gap.old) < eps) {
      return(list(beta=beta, iter=iter, gap=gap, status="no further improvement..."))
    }
    if (gap < eps) {
      return(list(beta=beta, iter=iter, gap=gap, status="converged"))
    } else {
      gap.old <- gap
      ## UPDATE THE ACTIVE SET / ACTIVE GROUP
      new <- which(dual.all == gap.old)
      if (length(new) > 0) {
        if (n2.nabla.f.pos[new][1] > n2.nabla.f.neg[new][1]) {
          active.neg[new] <- TRUE
        } else {
          active.pos[new] <- TRUE
        }
      }
    }
  }
}
