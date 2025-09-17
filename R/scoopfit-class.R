#' @export
setClass("scoopfit",
  representation = representation(
     coefficients  = "matrix",
     lambda        = "numeric",
     group         = "ANY",
     wk            = "ANY",
     family        = "character",
     penalty       = "character",
     x             = "matrix",
     y             = "ANY",
     monitoring    = "list",
     call          = "call")
)

# exportClasses(scoopfit, cvscoop)
# exportMethods(print       ,
#               plot        ,
#               fitted      ,
#               predict     ,
#               residuals   ,
#               deviance    ,
#               selection   ,
#               crossval)


#' @export
setMethod("print", "scoopfit", definition =
   function(x, ...) {
     ncoef <- ncol(x@coefficients)
     if (is.null(ncoef)) {ncoef <- length(x@coefficients)}
     cat("- number of coefficients:", ncoef, "\n")
     cat("- number of groups:", length(unique(x@group)), "\n")
     cat("- penalty parameter lambda:", length(x@lambda), "points from",
         format(max(x@lambda), digits = 3),"to",
         format(min(x@lambda), digits = 3))
     cat("\n")
     invisible(x)
   }
)

group.norm.unordered <- function(beta,group) {
  n2 <- c(sqrt(rowsum(beta^2,group,reorder=FALSE)))
  return(n2)
}

#' @export
setMethod("plot", "scoopfit", definition =
   function(x, y,
            xvar = "lambda", yvar="coefficients",
            main = paste(x@penalty," path (", yvar, ")", sep=""),
            crit = NULL, log.scale = TRUE, labels = NULL,
            col = c(1+switch(yvar,
              group = c(1:nlevels(as.factor(x@group))),
              c(x@group))),
            lty = switch(yvar,
              group = c(1:nlevels(as.factor(x@group))),
              c(x@group)), ...) {
    
     if (length(x@lambda) == 1) {
       stop("Not available when length(lambda) == 1")
     }

     ## removing intercept for plotting
     beta   <- x@coefficients[,!is.intercept(x)]
     if (yvar == "group") {
       beta <- t(apply(beta,1,group.norm.unordered,group=as.factor(x@group)))
     }
     
     xv <- switch(xvar,
                  "df"       = df.scoop(x),
                  "fraction" = apply(abs(beta),1,sum)/max(apply(abs(beta),1,sum)),
                  ifelse(rep(log.scale,length(x@lambda)),log10(x@lambda),x@lambda))
     xlab <- switch(xvar,
                    "df"       = "degrees of freedom",
                    "fraction" = "fraction",
                    ifelse(log.scale, "lambda (log scale)", "lambda"))
     
     matplot(xv, beta, xlab = xlab, ylab = yvar, col=col, lty=lty, main=main,
             type="l", yaxt = switch(!is.null(labels), T="n", "s"), ...)
     abline(h = 0, lty = 3)
     
     if (!is.null(labels)) { 
       axis(switch(xvar,"lambda"=2,4), at = beta[which.max(rowSums(abs(beta))),],
            labels = labels, cex.axis=.65, las=2)
     }
     
     if (!is.null(crit)) {
       x.star <- xv[max(which(crit == min(crit, na.rm=TRUE)), na.rm=TRUE)]
       abline(v=x.star)
     }
   }
)


#' @export
setMethod("fitted", "scoopfit", definition =
   function(object, ...) {
     return(predict(object))
   }
)

#' @export
setMethod("predict", "scoopfit", definition =
   function (object, newx=NULL, ...)  {

     sigmoid <- function(eta) {
       return( 1/(1+exp(-eta)) )
     }
     
     if (is.null(newx)) {newx <- object@x}
  
     if (object@family == "gaussian") {
       if (sum(is.intercept(object))) {
         intercept <- matrix(rep(t(object@coefficients[,1]), nrow(newx)),
                             nrow=nrow(newx), byrow=TRUE)
         return(intercept + newx %*% t(object@coefficients[,-1]))
       } else {
         return(newx %*% t(object@coefficients))
       }
     }

     if (object@family == "binomial") {
       if (sum(is.intercept(object))) {
         intercept <- matrix(rep(t(object@coefficients[,1]), nrow(newx)),
                             nrow=nrow(newx), byrow=TRUE)
         eta <- as.matrix(newx %*% t(object@coefficients[,-1]))
         return(sigmoid(intercept + eta))
       } else {
         return(sigmoid(newx %*% t(object@coefficients)))
       }    
     }  
   }
)

#' @export
setMethod("residuals", "scoopfit", definition =
   function(object, ...) {
     n <- length(object@lambda)
       return(matrix(rep(object@y, n), ncol=n) - fitted(object))
   }
)

#' @export
setMethod("deviance", "scoopfit", definition =
   function(object, ...) {

     if (object@family == "gaussian") {
       dev <- apply(residuals(object)^2,2,sum)
     }
     if (object@family == "binomial") {
       y.hat <- predict(object)
       y     <- matrix(rep(as.matrix(object@y),ncol(y.hat)),ncol=ncol(y.hat))
       dev <- 2*colSums((log(1-y))^(1-y) + (log(y))^y - (log(1-y.hat))^(1-y)
                        - (log(y.hat))^y)
     }
     return(dev)
     
   }
)

setGeneric("selection", function(object, sigma2=NULL)
           {standardGeneric("selection")})

#' @export
setMethod("selection", "scoopfit", definition =
   function(object, sigma2=NULL) {
     
     stopifnot(object@family == "gaussian" |
               (object@family == "binomial" & object@penalty == "lasso"))
     
     Beta   <- object@coefficients
     lambda <- object@lambda
     n  <- nrow(object@x)
     df <- df.scoop(object)
     L  <- length(lambda)
     p  <- ncol(Beta)
     if (sum(is.intercept(object))) {
       p <- p-1
     }
     
     BIC <- beta.BIC <- lambda.BIC <- df.BIC <- NA
     AIC <- beta.AIC <- lambda.AIC <- df.AIC <- NA
     
     RSS <- deviance(object) 
     
     if (length(RSS) < p & (is.null(sigma2))) {
       sigma2 <- NA
       stop("Cannot compute OLS nor estimate sigma2 in the high dimensional context.")
     }
    
     if (is.null(sigma2)) {
       reg.fit <- suppressWarnings(reg.fitting(object@x, object@y, sum(is.intercept(object)), object@family))
       sigma2 <- deviance(reg.fit)/df.residual(reg.fit)
     }
     
     ## BIC choice
     BIC    <- RSS / sigma2  + log(n) * df
     choice <- which.min(BIC)
     beta.BIC   <- Beta[choice,]
     lambda.BIC <- lambda[choice]
     
     ## AIC choice
     AIC    <- RSS / sigma2 + 2*df
     choice <- which.min(AIC)
     beta.AIC   <- Beta[choice,]
     lambda.AIC <- lambda[choice]
     
     return(list(BIC=BIC, beta.BIC=beta.BIC, lambda.BIC=lambda.BIC, 
                 AIC=AIC, beta.AIC=beta.AIC, lambda.AIC=lambda.AIC, 
                 sigma2 = sigma2, df = df,
                 non.zeros = rowSums(object@coefficients != 0)))
     
   }
)

setGeneric("crossval", function(object,
                                K          = 10,
                                error      = "deviance",
                                threshold  = 0.5,
                                folds      = NULL,
                                verbose    = TRUE)
           {standardGeneric("crossval")})

#' @export
setMethod("crossval", "scoopfit", definition =
   function (object,
             K          = 10,
             error      = "deviance",
             threshold  = 0.5,
             folds      = NULL,
             verbose    = TRUE) {

     ## =====================================================================
     ## INITIALIZATION
     ## extract all variables needed from object scoop
     x         <- object@x
     y         <- object@y
     group     <- object@group
     wk        <- object@wk
     family    <- object@family
     penalty   <- object@penalty
     intercept <- sum(is.intercept(object)) > 0
     lambda    <- object@lambda
     call      <- object@call
     ## some more work to extract eps and normalize from the call to scoop
     call.eps <- match(c("eps"), names(call), 0L)
     if (call.eps == 0) {
       eps <- 1e-4
     } else {
       if(is.numeric(call[call.eps][[1]])) {
         eps <- call[call.eps][[1]]
       } else {
         eps <- 1e-4
       }
     }
     call.normalize <- match(c("normalize"), names(call), 0L)
     if (call.normalize == 0) {
       normalize <- ifelse(family=="gaussian",TRUE,FALSE)
     } else {
       if (is.logical(call[call.normalize][[1]])) {
         normalize <- call[call.normalize][[1]]
       } else {
         normalize <- ifelse(family=="gaussian",TRUE,FALSE)
       }
     }
     call.optim <- match(c("optim.method"), names(call), 0L)
     if (call.optim == 0) {
       optim.method <- "fista"
     } else {
       if (is.character(call[call.optim][[1]])) {
         optim.method <- call[call.optim][[1]]
       } else {
         optim.method <- "fista"
       }
     }
     call.max.iter <- match(c("max.iter"), names(call), 0L)
     if (call.max.iter == 0) {
       max.iter <- 2*ncol(x)
     } else {
       if(is.numeric(call[call.max.iter][[1]])) {
         max.iter <- call[call.max.iter][[1]]
       } else {
         max.iter <- 2*ncol(x)
       }
     }
     
     ## =====================================================================
     ## CROSS-VALIDATION WORK

     ## if folds are specified by the user
     if (is.null(folds)) {
       folds <- split(sample(1:length(y)), rep(1:K, length = length(y)))
     } else {
       K <- length(folds)
     }
     
     ## now loop on the fold and cross-validate
     if(verbose) {cat("omitting bloc")}
     err <- matrix(0, length(lambda), K)

     func <- function(i){
       if(verbose) {cat("",i)}
       omit <- folds[[i]]
       fit <- scoop(x[-omit,], y[-omit], group, penalty = penalty, family=family,
                    intercept = intercept, normalize = normalize, lambda = lambda,
                    eps = eps, wk = wk, optim.method = optim.method,
                    verbose = FALSE, max.iter = max.iter, call = call)
       if (error == "class" & family == "binomial") {
         y.hat <- predict(fit,x[omit,])     
         y.hat[y.hat  > threshold] <- 1
         y.hat[y.hat <= threshold] <- 0
         return(apply( (y[omit] - y.hat) ^2,2,mean))
       } else {
         return(err.deviance(fit, x[omit,], y[omit]) / length(omit))
       }
     }
     err <- sapply(1:K, func)
     
     ## average over folds
     cv.mean     <- apply(err, 1, mean)
     cv.error    <- sqrt(apply(err, 1, var)/K)
     
     lambda.min <- max(lambda[cv.mean <= min(cv.mean)], na.rm=TRUE)
     lambda.1se <- max(lambda[cv.mean<=(cv.mean+cv.error)[match(lambda.min,lambda)]],
                       na.rm=TRUE)
     
     beta.min <- object@coefficients[match(lambda.min, lambda),]
     beta.1se <- object@coefficients[match(lambda.1se, lambda),]

     return(new("cvscoop",
                lambda     = lambda,
                lambda.min = lambda.min,
                lambda.1se = lambda.1se,
                cv.mean    = cv.mean,
                cv.error   = cv.error,
                folds      = folds,
                beta.min   = beta.min,
                beta.1se   = beta.1se))
   }
)


setGeneric("err.deviance", function(o, x, y) {standardGeneric("err.deviance")})

setMethod("err.deviance", "scoopfit", definition =
   function(o, x, y) {

     if (o@family == "gaussian") {
       dev <- apply((y-predict(o,x))^2,2,sum)
     }
     if (o@family == "binomial") {
       y.hat <- predict(o, x)
       y     <- matrix(rep(as.matrix(y),ncol(y.hat)),ncol=ncol(y.hat))    
       dev <- 2*colSums((log(1-y))^(1-y) + (log(y))^y - (log(1-y.hat))^(1-y)
                        - (log(y.hat))^y)
     }
     return(dev)
   }
)

setGeneric("is.intercept", function(o){standardGeneric("is.intercept")})

setMethod("is.intercept", "scoopfit", definition =
   function(o) {
     ncoef <- ncol(o@coefficients)
     if (is.null(ncoef)) {ncoef <- length(o@coefficients)}
     if (is.null(dim(o@group))) {
       inter <- ifelse (length(o@group) == ncoef, FALSE, TRUE) 
     } else {
       inter <- ifelse (ncol(o@group) == ncoef, FALSE, TRUE) 
     }
     if (!inter) {
       return(rep(FALSE,ncoef))
     } else {
       return(c(TRUE,rep(FALSE,ncoef-1)))
     }
   }
)

setGeneric("df.scoop", function(o){standardGeneric("df.scoop")})

setMethod("df.scoop", "scoopfit", definition =
   function(o) {

     if (o@family == "binomial") {
       warning("df for binomial family not implemented - returning NULL.")
       return(NULL)
     }

     reg.fit <- suppressWarnings(reg.fitting(o@x, o@y, sum(is.intercept(o)), o@family))
     
     if (sum(is.intercept(o))) {
       beta   <- o@coefficients[,-1]
       beta.reg <- coefficients(reg.fit)[-1]
     } else {
       beta   <- o@coefficients
       beta.reg <- coefficients(reg.fit)
     }
     group <- o@group
     
     if (is.null(reg.fit$lambda)) {
       gamma <- 0
     } else {
       gamma <- reg.fit$lambda
     }
     
     L <- nrow(beta)
     g <- length(table(group))
     
     if (o@penalty == "lasso") {
       df <- rowSums(beta != 0)  
     } else
     if (o@penalty == "group") {
       
       n2.reg <- group.norm.unordered(beta.reg,group)
       pg <- tapply(beta.reg != 0,group,sum)
       
       n2       <- apply(as.matrix(beta),1,group.norm.unordered,group)
       n2.scale <- t(scale(t(n2),FALSE,n2.reg/((pg-1) / (1+gamma) )))
       
       I <- matrix(0,ncol=L,nrow=g)
       I[n2>0] <- 1
       
       df <- colSums( I * (1 + n2.scale) )
       
     } else
     if (o@penalty == "coop") {
       
       n2.reg.p <- group.norm.unordered(pmax(0,beta.reg),group)
       n2.reg.m <- group.norm.unordered(pmin(0,beta.reg),group)
       pg.p     <- tapply(beta.reg > 0, group,sum)
       pg.m     <- tapply(beta.reg < 0, group,sum)
       
       ## Manage 0/0 = 0
       n2.reg.p[n2.reg.p == 0] <- Inf
       n2.reg.m[n2.reg.m == 0] <- Inf
       
       n2.p <- apply(matrix(pmax(0,beta),nrow=L), 1, group.norm.unordered, group)
       n2.m <- apply(matrix(pmin(0,beta),nrow=L), 1, group.norm.unordered, group)
       
       n2.scale.p  <- t(scale(t(n2.p), FALSE, n2.reg.p / ( (pg.p-1)/(1+gamma) ) ))
       n2.scale.m  <- t(scale(t(n2.m), FALSE, n2.reg.m / ( (pg.m-1)/(1+gamma) ) ))
       
       Ip  <- matrix(0,ncol=L,g)
       Im <- matrix(0,ncol=L,g)
       Ip[n2.p>0] <- 1
       Im[n2.m>0] <- 1
       df <- colSums( Ip*(1 + n2.scale.p) + Im*(1 + n2.scale.m ) )
     }
     
     if (sum(is.intercept(o))) {
       df <- df + 1
     }
     names(df) <- round(o@lambda,2)
     
     return(df)
   }
)

reg.fitting <- function(x, y, intercept, family) {
  
  reg.fit <- NULL  
  if (intercept) {
    my.formula <- y ~ x + 1
  } else {
    my.formula <- y ~ x + 0
  }
  ## - OLS in the gaussian case
  if (family == "gaussian") {        
    if (nrow(x) <= ncol(x)) {
      try(reg.fit <- lm.ridge(my.formula, lambda=0), silent=TRUE)
    } else {
      try(reg.fit <- lm(my.formula), silent=TRUE)
    }
  }
  ## - Logistic in the binomial case
  if (family == "binomial") {    
    if (nrow(x) > ncol(x)) {
      try(reg.fit <- glm(my.formula, family="binomial"), silent=TRUE)
    }
  }

  return(reg.fit)
}

