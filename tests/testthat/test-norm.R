set.seed(1234)

d <- 30
k <- 3
n <- 100

theta <- rnorm(d)
group <- rep(1:k, each=d/k)
grp_s <- tabulate(group)
weights_elt <- rep(1,d)
weights_grp <- rep(sqrt(d/k), k)

theta_prox1 <- scoop:::proximal_L1(theta, weights_elt, 0.5)

scoop:::proximal_las(rep(0.5,d) , 1, theta)

print(sum( (theta_prox1 - theta)^2 ))

ggplot2::autoplot(microbenchmark::microbenchmark(
  scoop:::proximal_las(weights_elt * 0.5, 1, theta),
  scoop:::proximal_L1(theta, weights_elt, 0.5)
))
