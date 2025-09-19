# load gglasso library
library(gglasso)

# load bardet data set
data(bardet)

# define group index
group_bardet <- rep(1:20,each=5)

# fit group lasso penalized least squares
m_gglasso <- gglasso(x=bardet$x, y=bardet$y, group=group_bardet, loss="ls")

m_scoop <- scoop::group.lasso(x=bardet$x, y=bardet$y, group=group_bardet, 
                              family = "gaussian", optim.method = "bfgs")
