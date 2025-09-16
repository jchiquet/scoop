#' R6 abstract class for a generic sparse modelmodel

SCOOP <- R6::R6Class(
  
  classname = "SCOOP",

  ## ______________________________________________________
  ##
  ## PRIVATE MEMBERS
  ##  
 
  private = list(
     name          = NA,
     formula       = NA,
     beta          = NA,
     ## model-related fields
     family        = NA,
     nloglik       = NA,
     gradient      = NA,
     ## model-related fields
     penalty       = NA,
     group         = NA,
     lambda        = NA,
     wk            = NA,
     grp_norm      = NA,
     pen_norm      = NA,
     subgrad_norm  = NA,
     subgrad_dual  = NA,
     dual_norm     = NA,
    ## optim-related fields
     solver        = NA,
     
     monitoring    = NA,     
  ),
  
  ## ______________________________________________________
  ##
  ## PUBLIC MEMBERS
  ##  
  public = list(
    data = NA,
     
    initialize = function(data, family) {
       
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
    d = function() ncol(private$beta),
    #' @field k number of groups
    k = function() length(private$wk),
    #' @field group vector of penalty grouping, when applicable
    grouping = function() private$group,
    #' @field family type of loss (aka family)
    family = function() private$family,
    #' @field penalty type of penalty function
    penalty_type = function() private$penalty,
    #' @field penalty_level vector of amounts of penalty applied (aka lambda)
    penalty_levels = function() private$lambda,
    #' @field penalty_weights vector of weights applied at the group level (or coefficients level if no grouping)
    penalty_weights = function() private$wk,
    #' @field fitted a matrix: fitted values
    fitted = function() {self$predict()},
    #' @filed residuals matrix of residuals
    residuals = function() {
      matrix(data$y - self$fitted)
    },
    #' @field deviance
    deviance = function() {
      if (self$family == "gaussian") {
        dev <- apply(residuals(object)^2,2,sum)
      }
      if (self$family == "binomial") {
        y_hat <- self$fitted
        y <- matrix(rep(data$y,ncol(y_hat)),ncol=ncol(y_hat))
        dev <- 2*colSums((log(1-y))^(1-y) + (log(y))^y - (log(1-y.hat))^(1-y) - (log(y.hat))^y)
      }
      dev
    }
  )
)

