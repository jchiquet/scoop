lasso <- function(x,
                  y,
                  family       = "gaussian",
                  intercept    = TRUE,
                  normalize    = ifelse(family=="gaussian",TRUE,FALSE),
                  wk           = rep(1,ncol(x)),
                  lambda       = NULL,
                  n.lambda     = 100,
                  lambda.min   = 1e-2,
                  verbose      = FALSE,
                  eps          = 1e-5,
                  max.iter     = 2*ncol(x),
                  optim.method = ifelse(family=="gaussian","fista","bfgs")) {

  out <- scoop(x = x,
               y = y,
               group        = 1:ncol(x),
               family       = family,
               penalty      = "lasso", 
               intercept    = intercept,
               normalize    = normalize,
               wk           = wk,
               lambda       = lambda,
               n.lambda     = n.lambda,
               lambda.min   = lambda.min,
               verbose      = verbose,
               eps          = eps,
               max.iter     = max.iter,
               optim.method = optim.method,
               call         = match.call())
}

group.lasso <- function(x,
                        y,
                        group,
                        family       = "gaussian",
                        intercept    = TRUE,
                        normalize    = ifelse(family=="gaussian",TRUE,FALSE),
                        wk           = sqrt(tabulate(group)),
                        lambda       = NULL,
                        n.lambda     = 100,
                        lambda.min   = 1e-2,
                        verbose      = FALSE,
                        eps          = 1e-5,
                        max.iter     = 2*ncol(x),
                        optim.method = ifelse(family=="gaussian","fista","bfgs")) {
  
  out <- scoop(x = x,
               y = y,
               group        = group,
               family       = family,
               penalty      = "group", 
               intercept    = intercept,
               normalize    = normalize,
               wk           = wk,
               lambda       = lambda,
               n.lambda     = n.lambda,
               lambda.min   = lambda.min,
               verbose      = verbose,
               eps          = eps,
               max.iter     = max.iter,
               optim.method = optim.method,
               call         = match.call())
}

coop.lasso <- function(x,y,
                       group,
                       family       = "gaussian",
                       intercept    = TRUE,
                       normalize    = ifelse(family=="gaussian",TRUE,FALSE),
                       wk           = sqrt(tabulate(group)),
                       lambda       = NULL,
                       n.lambda     = 100,
                       lambda.min   = 1e-2,
                       verbose      = FALSE,
                       eps          = 1e-5,
                       max.iter     = 2*ncol(x),
                       optim.method = ifelse(family=="gaussian","fista","bfgs")) {
  
  out <- scoop(x = x,
               y = y,
               group        = group,
               family       = family,
               penalty      = "coop", 
               intercept    = intercept,
               normalize    = normalize,
               wk           = wk,
               lambda       = lambda,
               n.lambda     = n.lambda,
               lambda.min   = lambda.min,
               verbose      = verbose,
               eps          = eps,
               max.iter     = max.iter,
               optim.method = optim.method,
               call         = match.call())

}

sparse.group.lasso <- function(x,y,
                               group,
                               family       = "gaussian",
                               intercept    = TRUE,
                               normalize    = ifelse(family=="gaussian",TRUE,FALSE),
                               wk           = sqrt(tabulate(group)),
                               lambda       = NULL,
                               n.lambda     = 100,
                               lambda.min   = 1e-2,
                               verbose      = FALSE,
                               eps          = 1e-5,
                               max.iter     = 2*ncol(x)) {

  ## sparse group lasso is a special case of tree structured group lasso
  out <- scoop(x = x,
               y = y,
               group        = rbind(group,1:ncol(x)),
               family       = family,
               penalty      = "treegroup", 
               intercept    = intercept,
               normalize    = normalize,
               wk           = list(wk,rep(1,ncol(x))),
               lambda       = lambda,
               n.lambda     = n.lambda,
               lambda.min   = lambda.min,
               verbose      = verbose,
               eps          = eps,
               max.iter     = max.iter,
               optim.method = "fista",
               call         = match.call())

}

sparse.coop.lasso <- function(x,y,
                              group,
                              family       = "gaussian",
                              intercept    = TRUE,
                              normalize    = ifelse(family=="gaussian",TRUE,FALSE),
                              wk           = sqrt(tabulate(group)),
                              lambda       = NULL,
                              n.lambda     = 100,
                              lambda.min   = 1e-2,
                              verbose      = FALSE,
                              eps          = 1e-5,
                              max.iter     = 2*ncol(x)) {
  
  ## sparse coop lasso is a special case of tree structured coop lasso
  out <- scoop(x = x,
               y = y,
               group        = rbind(group,1:ncol(x)),
               family       = family,
               penalty      = "treecoop", 
               intercept    = intercept,
               normalize    = normalize,
               wk           = list(wk,rep(1,ncol(x))),
               lambda       = lambda,
               n.lambda     = n.lambda,
               lambda.min   = lambda.min,
               verbose      = verbose,
               eps          = eps,
               max.iter     = max.iter,
               optim.method = "fista",
               call         = match.call())

}

tree.group.lasso <- function(x,y,
                             group,
                             family       = "gaussian",
                             intercept    = TRUE,
                             normalize    = ifelse(family=="gaussian",TRUE,FALSE),
                             wk           = lapply(apply(group,1,tabulate),sqrt),
                             lambda       = NULL,
                             n.lambda     = 100,
                             lambda.min   = 1e-2,
                             verbose      = FALSE,
                             eps          = 1e-5,
                             max.iter     = 2*ncol(x)) {
  
  out <- scoop(x = x,
               y = y,
               group        = group,
               family       = family,
               penalty      = "treegroup", 
               intercept    = intercept,
               normalize    = normalize,
               wk           = wk,
               lambda       = lambda,
               n.lambda     = n.lambda,
               lambda.min   = lambda.min,
               verbose      = verbose,
               eps          = eps,
               max.iter     = max.iter,
               optim.method = "fista",
               call         = match.call())

}

tree.coop.lasso <- function(x,y,
                             group,
                             family       = "gaussian",
                             intercept    = TRUE,
                             normalize    = ifelse(family=="gaussian",TRUE,FALSE),
                             wk           = lapply(apply(group,1,tabulate),sqrt),
                             lambda       = NULL,
                             n.lambda     = 100,
                             lambda.min   = 1e-2,
                             verbose      = FALSE,
                             eps          = 1e-5,
                             max.iter     = 2*ncol(x)) {
  
  out <- scoop(x = x,
               y = y,
               group        = group,
               family       = family,
               penalty      = "treecoop", 
               intercept    = intercept,
               normalize    = normalize,
               wk           = wk,
               lambda       = lambda,
               n.lambda     = n.lambda,
               lambda.min   = lambda.min,
               verbose      = verbose,
               eps          = eps,
               max.iter     = max.iter,
               optim.method = "fista",
               call         = match.call())
 }

scoop <- function(x,
                  y,
                  group       ,
                  family      , 
                  penalty     , 
                  intercept   , 
                  normalize   , 
                  wk          , 
                  lambda      , 
                  n.lambda    , 
                  lambda.min  , 
                  verbose     , 
                  eps         , 
                  max.iter    , 
                  optim.method,
                  call) {

  ## Keep the colnames for later
  if (is.null(colnames(x))) {colnames(x) <- 1:ncol(x)}
  names <- colnames(x)
  
  ## INITIALIZATION OF MAIN OBJECTS 
  ## fit's construction
  scoop.fit <- new("scoopfit",
                   coefficients = matrix(),
                   lambda       = numeric(),
                   group        = group,
                   wk           = wk,
                   family       = family,
                   penalty      = penalty,
                   x            = x,
                   y            = y,
                   monitoring   = list(),
                   call         = call)

  ## model's construction
  scoop.mdl <- new.scoop.model(family, penalty, optim.method)
  
  ## ======================================================
  ## INTERCEPT TREATMENT
  mu0 <- NULL
  if (intercept) {
    names <- c("intercept",names)
    x.bar <- colMeans(x)
    x     <- scale(x,x.bar,FALSE) 
    ## Gaussian family
    if (family == "gaussian") {
      y.bar <- mean(y)
      y <- y-y.bar
    }
    ## Binomial family
    if (family == "binomial") {
      mu0 <- log(sum(y==1)/sum(y==0))
    }   
  }
  
  ## ======================================================
  ## NORMALIZATION, REORDERING

  ## normalizing the data
  if (normalize) {
    norm <- sqrt(drop(colSums(x^2)))
    x    <- scale(x, FALSE, norm)
  }
  ## group numerotation MUST start from 1
  ## sorting the groups and the columns of the design matrix
  if (!(penalty %in% c("treegroup","treecoop"))) {
    o       <- order(group,decreasing=FALSE)
    x       <- x[,o]
    group   <- group[o]
    tpk     <- tabulate(group)
  } else {
    tpk     <- apply(group,1,tabulate)
  }
  
  ## ======================================================
  ## GENERATE A GRID OF LAMBDA IF NONE HAS BEEN PROVIDED
  if (is.null(lambda)) {
    y.tilde <- scale(y, (family == "binomial" & intercept), FALSE)
    lambda <- lambda.grid(scoop.mdl, x, y.tilde, tpk, wk, n.lambda, lambda.min)
  }
  
  ## ========================================================
  ## LOOP OVER LAMBDA 

  ## quantities which won't move along the regularization path
  ## differs weither the family is gaussian or binomial
  fx <- switch(family, "gaussian" = crossprod(x,x)   , "binomial" = x)
  fy <- switch(family, "gaussian" = c(crossprod(y,x)), "binomial" = y)

  ## Some initializations
  beta  <- c()
  beta0 <- c(mu0,rep(0,ncol(x))) ## mu0 is NULL if no intercept required
  last.beta   <- beta0-2*eps
  real.lambda <- c()
  monitoring  <- list(it.active=c(), it.optim=c(), dual.gap=c())

  ## Managing intercept for the binomial case
  if (intercept & family == "binomial") {
    fx <- cbind(1,fx)
    if (penalty %in% c("treegroup","treecoop")) {
      wk <- lapply(wk, function(x) c(0,x))
      group  <- cbind(0,group) + 1
    } else {
      wk <- c(0,wk)
      group  <- c(0,group) + 1
    }
  }

  for (l in lambda) {
    if (verbose) cat("\nlambda =",l)
    
    ## CALL TO THE FUNCTION THAT DOES THE WORK
    out <- solver(scoop.mdl, fx, fy, group, l, wk, beta0, eps, max.iter)

    monitoring$it.active <- c(monitoring$it.active, length(out$iter))
    monitoring$it.optim  <- c(monitoring$it.optim , out$iter)
    monitoring$dual.gap  <- c(monitoring$dual.gap , out$gap)
    
    ## When the path of solutions has stabilized, try to earn some time
    ## by avoiding irreleavant computations
##    if (sqrt(sum((last.beta-out$beta)^2)) < eps & !all(last.beta == 0)) {
    ## if (sqrt(sum((last.beta-out$beta)^2)) < eps & length(out$iter) > 0 & (!(penalty %in% c("treegroup","treecoop")))) {
    ##   missing <- n.lambda - length(real.lambda)
    ##   beta    <- rbind(beta,t(matrix(rep(out$beta,missing),ncol=missing)))
    ##   real.lambda <- lambda
    ##   break
    ## }
    
    ## Move to the next lambda
    beta  <- rbind(beta,out$beta)
    real.lambda <- c(real.lambda,l)
    last.beta <- out$beta
    
    ## When the algorithm does not converge, send a message and start over from
    ## zero for the next iterate
    if (verbose & out$status != "converged") {
      cat("\nscoop:",out$status, "with max.min||subgradient|| =", out$gap)
      beta0 <-  c(mu0,rep(0,ncol(x)))
    } else {
      beta0 <- out$beta
    }

  } ## END OF THE LOOP OVER LAMBDA

  ## Managing intercept for the binomial case
  if (intercept & family == "binomial") {
    mu   <- beta[, 1]
    beta <- beta[,-1]
  }
  
  scoop.fit@lambda     <- real.lambda
  scoop.fit@monitoring <- monitoring
  
  ## ======================================================
  ## RESCALING, REORDERING
  ## unsorting beta
  if (!(penalty %in% c("treegroup","treecoop"))) {
    if (length(lambda) > 1) {
      beta  <- beta[, order(o)]
    } else {
      beta  <- beta[order(o)]
    }
  }
  ## re-scaling beta
  if (normalize) {
    if (length(lambda) > 1) {
      beta <- scale(beta, FALSE, norm)
    } else {
      beta <- beta/norm
    }
  }
  
  ## ======================================================
  ## INTERCEPT TREATMENT
  beta <- matrix(beta,ncol=ncol(x), nrow=length(lambda))
  if (intercept) {
    mu <- switch(family,
                 "gaussian" = y.bar -  beta %*% x.bar,
                 "binomial" = mu - apply(t(beta) * x.bar, 2, sum))
    beta <- cbind(mu,beta)    
  }
  colnames(beta) <- names
  rownames(beta) <- round(lambda,2)
  
  scoop.fit@coefficients <- beta

  return(scoop.fit)

}
