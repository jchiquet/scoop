set.seed(1234)

d <- 30
k <- 3
n <- 100

theta <- rnorm(d)
group <- rep(1:k, each=d/k)
weights_elt <- rep(1,d)
weights_grp <- rep(sqrt(d/k), k)

X <- matrix(rnorm(n * d), n, d)
epsilon <- rnorm(n)
mu <- 5


test_that("DataModel for Gaussian response (Linear Regression) works",  {

  y <- mu + X %*% theta + epsilon
  
  Data <- GaussianModel$new(X, y, group)
  expect_true(Data$has_intercept)
  expect_true(Data$is_standardized)
  expect_equal(Data$d, d)
  expect_equal(Data$n, n)
  expect_equal(Data$k, k)
  expect_true(Data$loss(theta) > 0)
  
})

test_that("DataModel for Binary response (Logistic Regression) works",  {
  
  y <- round(scoop:::.sigmoid(mu + X %*% theta))
  
  Data <- BinaryModel$new(X, y, group)
  expect_true(Data$has_intercept)
  expect_true(Data$is_standardized)
  expect_equal(Data$d, d)
  expect_equal(Data$n, n)
  expect_equal(Data$k, k)
  expect_true(Data$loss(theta) > 0)
  
})