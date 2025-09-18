set.seed(1234)

d <- 30
k <- 3
n <- 100

theta <- rnorm(d)
group <- rep(1:k, each = d/k)
weights_elt <- rep(1,d)
weights_grp <- rep(sqrt(d/k), k)

X <- matrix(rnorm(n * d), n, d)
epsilon <- rnorm(n)
mu <- 5


test_that("DataModel for Gaussian response (Linear Regression) works",  {
  
  y <- mu + X %*% theta + epsilon
  
  myModel   <- GaussianModel$new(X, y, group)
  myPenalty <- Penalty_GroupLasso$new(group, weights_grp)  

  Scoop <- SCOOP$new(myModel, myPenalty)
})

