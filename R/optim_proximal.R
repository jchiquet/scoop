ista <- function(model, x0, fx, fy, lambda, pk, wk, lower=NULL, upper=NULL, L0, eps) {
  if (model@penalty %in% c("treegroup","treecoop")) {
    Ks <- sapply(pk,length)
    pk <- unlist(pk)
    wk <- unlist(wk)
  } else {
    Ks <- 0
  }

  it <- 0
  switch(model@family,
  "gaussian" =
    ISTA_LM(Ks, pk, fx, fy, lambda*wk, L0, 10000L, eps/10, x0, x0, 
            switch(model@penalty, "lasso"=1, "group"=2, "coop"=3, "treegroup"=4,"treecoop"=5), it),
  "binomial" = 
    ISTA_LRM(Ks, pk, fx, fy, lambda*wk, L0, 10000L, eps/length(x0), x0, x0,
            switch(model@penalty, "lasso"=1, "group"=2, "coop"=3, "treegroup"=4,"treecoop"=5), it)
  )
  list(xk = x0, i = it)
}

fista <- function(model, x0, fx, fy, lambda, pk, wk, lower=NULL, upper=NULL, L0, eps) {
  if (model@penalty %in% c("treegroup","treecoop")) {
    Ks <- sapply(pk,length)
    pk <- unlist(pk)
    wk <- unlist(wk)
  } else {
    Ks <- 0
  }
  it <- 0
  switch(model@family,
    "gaussian" =
    FISTA_LM(Ks, pk, fx, fy, lambda*wk, L0, 10000L, eps/10, x0, x0, 
            switch(model@penalty, "lasso"=1, "group"=2, "coop"=3, "treegroup"=4,"treecoop"=5), 0),
    "binomial" = 
    FISTA_LRM(Ks, pk, fx, fy, lambda*wk, L0, 10000L, eps/length(x0), x0, x0,
            switch(model@penalty, "lasso"=1, "group"=2, "coop"=3, "treegroup"=4,"treecoop"=5), 0)
  )
  list(xk = x0, i = it)
}
