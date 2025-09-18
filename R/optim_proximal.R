# ista <- function(model, x0, fx, fy, lambda, pk, wk, lower=NULL, upper=NULL, L0, eps) {
#   if (model@penalty %in% c("treegroup","treecoop")) {
#     Ks <- sapply(pk,length)
#     pk <- unlist(pk)
#     wk <- unlist(wk)
#   } else {
#     Ks <- 0
#   }
#   return(switch(model@family,
#                 "gaussian" = .C("ISTA_LM",
#                   as.integer(length(x0)),
#                   as.integer(length(Ks)),
#                   as.integer(length(pk)),
#                   as.integer(Ks),
#                   as.integer(pk),
#                   as.double(fx),
#                   as.double(fy),
#                   as.double(lambda * wk),
#                   as.double(L0),
#                   as.integer(10000), # max # of iteration
#                   as.double(eps/10),
#                   x0 = as.double(x0),
#                   xk = as.double(x0),
#                   as.integer(switch(model@penalty, "lasso"=1, "group"=2, "coop"=3,
#                                     "treegroup"=4,"treecoop"=5)),
#                   i  = as.integer(0)),
#                 "binomial" = .C("ISTA_LRM",
#                   as.integer(length(x0)),
#                   as.integer(length(Ks)),
#                   as.integer(length(pk)),
#                   as.integer(Ks),
#                   as.integer(pk),
#                   as.integer(length(fy)),
#                   as.double(fx),
#                   as.double(fy),
#                   as.double(lambda * wk),
#                   L = as.double(L0),
#                   as.integer(10000),
#                   as.double(eps/length(x0)),
#                   x0 = as.double(x0),
#                   xk = as.double(x0),
#                   as.integer(switch(model@penalty, "lasso"=1, "group"=2, "coop"=3,
#                                     "treegroup"=4, "treecoop"=5)),
#                   i  = as.integer(0))))
# }
# 
# fista <- function(model, x0, fx, fy, lambda, pk, wk, lower=NULL, upper=NULL, L0, eps) {
#   if (model@penalty %in% c("treegroup","treecoop")) {
#     Ks <- sapply(pk,length)
#     pk <- unlist(pk)
#     wk <- unlist(wk)
#   } else {
#     Ks <- 0
#   }
#   return(switch(model@family,
#                 "gaussian" = .C("FISTA_LM",
#                   as.integer(length(x0)),
#                   as.integer(length(Ks)),
#                   as.integer(length(pk)),
#                   as.integer(Ks),
#                   as.integer(pk),
#                   as.double(fx),
#                   as.double(fy),
#                   as.double(lambda * wk),
#                   as.double(L0),
#                   as.integer(10000), # max # of iteration
#                   as.double(eps/10),
#                   x0 = as.double(x0),
#                   xk = as.double(x0),
#                   as.integer(switch(model@penalty, "lasso"=1, "group"=2, "coop"=3,
#                                     "treegroup"=4, "treecoop"=5)),
#                   i  = as.integer(0)),
#                 "binomial" = .C("FISTA_LRM",
#                   as.integer(length(x0)),
#                   as.integer(length(Ks)),
#                   as.integer(length(pk)),
#                   as.integer(Ks),
#                   as.integer(pk),
#                   as.integer(length(fy)),
#                   as.double(fx),
#                   as.double(fy),
#                   as.double(lambda * wk),
#                   L = as.double(L0),
#                   as.integer(10000),
#                   as.double(eps/length(x0)),
#                   x0 = as.double(x0),
#                   xk = as.double(x0),
#                   as.integer(switch(model@penalty, "lasso"=1, "group"=2, "coop"=3,
#                                     "treegroup"=4, "treecoop"=5)),
#                   i  = as.integer(0))))
# }
