solver.tree <- function (model, fx, fy, group, lambda, wk, beta, eps, iter) {

  ## INITIALIZATION
  ## forme compact de liste des groupes
  pk <- rev(apply(group,1,tabulate))

  wk <- rev(wk)

  L0 <- switch (model@family,
                "gaussian" = max(eigen(fx)$values,TRUE,TRUE),
                "binomial" = sum(rowSums(fx^2)))
  
  res <- optimization(model, beta, fx, fy, lambda, pk, wk, L0=L0, eps=eps)
  
  return(list(beta=res$xk, iter=res$i, gap=eps, status="converged"))
  
}

