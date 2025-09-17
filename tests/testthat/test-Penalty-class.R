
set.seed(1234)
theta <- rnorm(30)
group <- rep(1:3, each=10)
weights_elt <- rep(1,30)
weights_grp <- c(1,1,1)

test_that("Penalty: Lasso",  {

  lasso <- Penalty_Lasso$new(weights_elt)  
  expect_equal(lasso$norm(theta), sum(abs(theta)))
  expect_equal(sum(lasso$elt_norm(theta)), sum(abs(theta)))
  expect_equal(lasso$dual_norm(theta), max(abs(theta)))
  expect_equal(
    as.numeric(lasso$proximal(theta, 0.5)),
    pmax(0, abs(theta) - 0.5) * sign(theta)
  )
})

test_that("Penalty: Bounded",  {
  
  bounded <- Penalty_Bounded$new(weights_elt)  
  expect_equal(bounded$norm(theta), max(abs(theta)))
  expect_equal(as.numeric(bounded$elt_norm(theta)), abs(theta))
  expect_equal(bounded$dual_norm(theta), lasso$norm(theta))
  
})

test_that("Penalty: Group-Lasso L1/L2",  {
  group_lasso <- Penalty_GroupLasso$new(group, weights_grp)  
  group_norm <- as.numeric(sqrt(rowsum(theta^2, group)))
  
  expect_equal(group_lasso$norm(theta), sum(group_norm))
  expect_equal(as.numeric(group_lasso$elt_norm(theta)),  group_norm)
  expect_equal(group_lasso$dual_norm(theta), max(group_norm))
  expect_equal(
    as.numeric(group_lasso$proximal(theta, 3)),
    pmax(0, 1 - 3/rep(sqrt(rowsum(theta^2, group)), each=10)) * theta
  )
})

test_that("Penalty: Group-Lasso L1/LINF",  {

  group_norm <- as.numeric(tapply(abs(theta), group, max))
  group_lasso <- Penalty_GroupLassoInf$new(group, weights_grp)  
  group_norm_dual <- as.numeric(tapply(abs(theta), group, sum))
  
  expect_equal(group_lasso$norm(theta), sum(group_norm))
  expect_equal(as.numeric(group_lasso$elt_norm(theta)),  group_norm)
  expect_equal(group_lasso$dual_norm(theta), max(group_norm_dual))
  # group_lasso$proximal(theta, 3)

})

test_that("Penalty: Coop-Lasso",  {
  
  coop_lasso <- Penalty_CoopLasso$new(group, weights_grp)  
  coop_lasso$norm(theta)
  coop_lasso$elt_norm(theta)
  coop_lasso$dual_norm(theta)
  coop_lasso$proximal(theta, 3)

})
