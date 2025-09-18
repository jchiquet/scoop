#' R6 abstract class for a generic sparse modelmodel

#' @export
SCOOP <- R6::R6Class(

  classname = "SCOOP",

  ## ______________________________________________________
  ##
  ## PRIVATE MEMBERS
  ##

  private = list(
    ## fit-related fields
    theta         = NA,
    lambda        = NA,
    ## optim-related fields
    optimizer     = NA,
    monitoring    = NA
  ),

  ## ______________________________________________________
  ##
  ## PUBLIC MEMBERS
  ##
  public = list(
    model   = NA,
    penalty = NA,

    initialize = function(data_model, penalty_model) {
       self$model   <- data_model
       self$penalty <- penalty_model
     },

    plot = function(x, y,
            xvar = "lambda", yvar="coefficients",
            main = paste(self$penalty_type," path (", yvar, ")", sep=""),
            crit = NULL, log_scale = TRUE, labels = NULL,
            col = c(1 + switch(yvar,group = 1:self$k, self$grouping)),
            lty = switch(yvar, group = self$k, self$grouping), ...) {

     if (length(private$lambda) == 1) {
       stop("Not available when length(lambda) == 1")
     }

     ## removing intercept for plotting
     beta   <- private$beta[,!is.intercept(x)]
     if (yvar == "group") {
       beta <- t(apply(beta,1,group.norm.unordered,group=as.factor(private$group)))
     }

     xv <- switch(xvar,
                  "df"       = df.scoop(x),
                  "fraction" = apply(abs(beta),1,sum)/max(apply(abs(beta),1,sum)),
                  ifelse(rep(log.scale,length(x@lambda)),log10(private$lambda),private$lambda))
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
    },

    predict = function (newx = self$data$x, ...)  {

     if (self$family == "gaussian") {
       if (sum(is.intercept(object))) {
         intercept <- matrix(rep(t(private$beta[,1]), nrow(newx)),
                             nrow=nrow(newx), byrow=TRUE)
         return(intercept + newx %*% t(private$beta[,-1]))
       } else {
         return(newx %*% t(object@coefficients))
       }
     }

     if (self$family == "binomial") {
       if (sum(is.intercept(object))) {
         intercept <- matrix(rep(t(private$beta[,1]), nrow(newx)),
                             nrow=nrow(newx), byrow=TRUE)
         eta <- as.matrix(newx %*% t(private$beta[,-1]))
         return(.sigmoid(intercept + eta))
       } else {
         return(.sigmoid(newx %*% t(private$beta)))
       }
     }
    },

    show = function() {
     cat("- number of coefficients:", self$d, "\n")
     cat("- number of groups:", self$k, "\n")
     cat("- penalty parameter lambda:", length(private$lambda), "points from",
         format(max(private$lambda), digits = 3),"to",
         format(min(private$lambda), digits = 3))
     cat("\n")
     invisible(self)
     }

  ),

  ## ______________________________________________________
  ##
  ## ACTIVE BINDINGS
  ##

  active = list(
    #' @field d number of variables
    d = function() self$model$d,
    #' @field k number of groups
    k = function() self$model$k,
    #' @field k number of groups
    n = function() self$model$n,
    #' @field group vector of penalty grouping
    grouping = function() self$penalty$group,
    #' @field data type of data model
    model_type = function() self$model$name,
    #' @field penalty type of penalty function
    penalty_type = function() self$penalty$name,
    #' @field penalty_level vector of amounts of penalty applied (aka lambda)
    penalty_levels = function() private$lambda,
    #' @field penalty_weights vector of weights applied at the group level (or coefficients level if no grouping)
    penalty_weights = function() self$penalty$weights
  #'   #' @field fitted a matrix: fitted values
  #'   fitted = function() {self$predict()},
  #'   #' @filed residuals matrix of residuals
  #'   residuals = function() {
  #'     matrix(data$y - self$fitted)
  #'   },
  #'   #' @field deviance
  #'   deviance = function() {
  #'     if (self$family == "gaussian") {
  #'       dev <- apply(residuals(object)^2,2,sum)
  #'     }
  #'     if (self$family == "binomial") {
  #'       y_hat <- self$fitted
  #'       y <- matrix(rep(data$y,ncol(y_hat)),ncol=ncol(y_hat))
  #'       dev <- 2*colSums((log(1-y))^(1-y) + (log(y))^y - (log(1-y.hat))^(1-y) - (log(y.hat))^y)
  #'     }
  #'     dev
  #'   }
  )
)

LASSO <- R6::R6Class(
  
  classname = "LASSO",
  inherit = SCOOP,
  
  ## ______________________________________________________
  ##
  ## PRIVATE MEMBERS
  ##
  
  ## ______________________________________________________
  ##
  ## PUBLIC MEMBERS
  ##
  public = list(
    solver = function(theta, lambda, eps, max.it) {
      ## _____________________________________________________________
      ##
      ## INITIALIZATION
      ## _____________________________________________________________
      ##
      iter <- c()
      pk  <- self$penalty$group_size
      rwk <- rep.int(self$penalty$weights,pk)
      active <- self$penalty$elt_norm(theta) > 0
      ## _____________________________________________________________
      ## 
      ## CHECK IF OPTIMALITY IS REACHED AT THE STARTING POINT
      ## _____________________________________________________________
      ##
      nabla_f <- self$model$gradient(fx,fy,beta)
      n.nabla_f <- model@subgrad.dual(nabla.f,pk)
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
  )
)

GROUP_LASSO <- R6::R6Class(
  
  classname = "GROUP_LASSO",
  inherit = SCOOP,
  
  ## ______________________________________________________
  ##
  ## PRIVATE MEMBERS
  ##
  
  ## ______________________________________________________
  ##
  ## PUBLIC MEMBERS
  ##
  public = list(
    initialize = function() {
      optimizer <- list(main = working.set)
      
    }
  )
)

COOPERATIVE_LASSO <- R6::R6Class(
  
  classname = "COOPERATIVE_LASSO",
  inherit = SCOOP,
  
  ## ______________________________________________________
  ##
  ## PRIVATE MEMBERS
  ##
  
  ## ______________________________________________________
  ##
  ## PUBLIC MEMBERS
  ##
  public = list(
  )
)

# "lasso" = working.set(object, fx, fy, group, lambda, wk, beta0, eps, iter),
# "group" = working.set(object, fx, fy, group, lambda, wk, beta0, eps, iter),
# "coop"  = solver.coop (object, fx, fy, group, lambda, wk, beta0, eps, iter),
