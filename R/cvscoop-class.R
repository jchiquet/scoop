setClass("cvscoop",
   representation = representation(
     lambda     = "numeric",
     lambda.min = "numeric",
     lambda.1se = "numeric",
     cv.mean    = "numeric",
     cv.error   = "numeric",
     folds      = "list",
     beta.min   = "numeric",
     beta.1se   = "numeric")
)

setMethod("plot", "cvscoop", definition =
   function(x, y,
            log.scale = TRUE,
            xlab = ifelse(log.scale, "lambda (log scale)", "lambda"),
            ylab = "Mean CV error", ...) {

     if(log.scale){
       xv <- log10(x@lambda)
       l.min <- log10(x@lambda.min)
       l.1se <- log10(x@lambda.1se)
     } else {
       xv <- x@lambda
       l.min <- x@lambda.min
       l.1se <- x@lambda.min
     }
     
     upper <- x@cv.mean + x@cv.error
     lower <- x@cv.mean - x@cv.error
     
     plot(xv, x@cv.mean, type="b", xlab=xlab, ylab=ylab,
          ylim=range(upper,lower, na.rm=TRUE), col="red", ...)
     barw <- diff(range(xv)) / length(xv)
     segments(xv, upper, xv, lower, col="gray")
     segments(xv - barw, upper, xv + barw, upper, col="gray")
     segments(xv - barw, lower, xv + barw, lower, col="gray")
     range(upper, lower)
     abline(v=l.min)
     axis(3, at = l.min, labels = "cross val min", cex.axis=0.75)
     abline(v=l.1se)
     axis(3, at = l.1se, labels = "cross val 1se", cex.axis=0.75)
   }
)
