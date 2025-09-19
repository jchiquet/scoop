data(basal)

cat("\nExecuting the whole script is a matter of minutes... \n")
## =================================================================
## FITTING COOP, GROUP AND LASSO MODELS
cat("\nMODEL FITTING \n")

cat("\nComputing sparse group-Lasso path...")
lasso <- lasso(x, y, family="binomial")

cat("\nComputing sparse group-Lasso path...")
spglasso <- sparse.group.lasso(x, y, group, family="binomial")

cat("\nComputing group-Lasso path...")
grplasso <- group.lasso(x, y, group, family="binomial", optim.method = "fista")

cat("\nComputing cooperative-Lasso path...")
coolasso <- coop.lasso(x, y, group, family="binomial")

K    <- 5
n    <- length(y)
seed <- 5555 #save seed for reproductibility
set.seed(seed)
folds <- split(sample(1:n), rep(1:K, length = n))

cat("\n5-fold CV for sparse group-Lasso - ")
cv.spg <- crossval(spglasso, 5, error="class", folds=folds)
cat("\n5-fold CV for group-Lasso - ")
cv.grp <- crossval(grplasso, 5, error="class", folds=folds)
cat("\n5-fold CV for coop-Lasso - ")
cv.coo <- crossval(coolasso, 5, error="class", folds=folds)

id.genes.spg <- sqrt(rowsum(cv.spg@beta.min[-1]^2, group)) > 0
id.genes.grp <- sqrt(rowsum(cv.spg@beta.min[-1]^2, group)) > 0
id.genes.coo <- sqrt(rowsum(cv.coo@beta.min[-1]^2, group)) > 0

genes.spg <- genes[id.genes.spg]
genes.grp <- genes[id.genes.grp]
genes.coo <- genes[id.genes.coo]

par(mfrow=c(2,3))
plot(grplasso, col=1+group, lty=group)
plot(spglasso, col=1+group, lty=group)
plot(coolasso, col=1+group, lty=group)
plot(cv.grp)
plot(cv.spg)
plot(cv.coo)

readline("Press for next plot")

id.genes <- id.genes.coo | id.genes.grp | id.genes.spg

all.genes <- genes[id.genes]
grp.genes <- tabulate(group)[id.genes]
names(grp.genes) <- as.character(all.genes)

beta.coo <- cv.coo@beta.min[-1][rep(id.genes, tabulate(group))]
names(beta.coo) <- as.character(rep(all.genes, tabulate(group)[id.genes]))

beta.grp <- cv.grp@beta.min[-1][rep(id.genes, tabulate(group))]
names(beta.grp) <- as.character(rep(all.genes, tabulate(group)[id.genes]))

beta.spg <- cv.spg@beta.min[-1][rep(id.genes, tabulate(group))]
names(beta.spg) <- as.character(rep(all.genes, tabulate(group)[id.genes]))

couleur <- rep(c(1,2,3,4), 3)[-1]
symbol  <- c(0,2,16,5,6,15,3,17,18,4,8)
cols <- rep(couleur, grp.genes)
pch  <- rep(symbol , grp.genes)

layout(mat=matrix(c(1:3,4,4,4),ncol=3,byrow=TRUE),
       width=c(1,1,1,3), height=c(1,1,1,1/3))
plot(beta.grp, col=cols, pch=pch, xaxt="n",
     ylim=range(beta.coo,beta.grp,beta.spg),lwd=2,
     xlab='probe sets',ylab="coefficients", main="group-Lasso")
abline(h=0, lty="dotted")
plot(beta.coo, col=cols ,pch=pch, xaxt="n",
     ylim=range(beta.coo,beta.grp,beta.spg),lwd=2,
     xlab='probe sets',ylab="coefficients", main="coop-Lasso")
abline(h=0, lty="dotted")
plot(beta.spg, col=cols, pch=pch, xaxt="n",
     ylim=range(beta.coo,beta.grp,beta.spg),lwd=2,
     xlab='probe sets',ylab="coefficients", main="sparse group-Lasso")
abline(h=0, lty="dotted")
plot(1:11, rep(0,11), xlim=c(0.95, 11.05), pch=symbol, lwd=2, col=couleur,
     axes=FALSE,  ylab="", xlab="gene names (group)")
text(1:11, rep(.5,11),  labels=all.genes, cex=.75)
