setwd("~/treevalues/")
devtools::load_all()
n <- 200
p <- 10
sigma_y <- 5
nTrials <- 1000
pvals <- rep(0, nTrials)
pvals2 <- rep(0, nTrials)

for (i in 1:nTrials) {
  set.seed(i)
  print(i)

  X <- MASS::mvrnorm(n, mu=rep(0,p), Sigma=diag(1, nrow=p, ncol=p))
  beta = 0
  mu_y_1 <- beta*I(X[,1] > 0)
  mu_y_2 <- mu_y_1 + 2*beta*(I(X[,1] > 0 & X[,2] > 0)) - 2*beta*(I(X[,1] < 0 & X[,2] > 0))
  mu_y <- mu_y_2 + beta*I(X[,3] > 0 & X[,2] > 0 & X[,1] > 0) - beta*I(X[,3] > 0 & X[,2] < 0 & X[,1] < 0)
  
  
  
  y <- rnorm(n, mu_y, sigma_y)

  ### This turns out to be necessary to reading the split rules.
  dat <- data.frame(y=y,X=X)
  nameX <- sapply(1:p, function(u) paste0("X",u))
  names(dat) = c("y", nameX)

  ### Build an rpart of depth d
  base_tree <- rpart::rpart(y~., data=dat, control=rpart.control(maxdepth = 3,
                                                                 minsplit=1, minbucket=1,
                                                                 cp=-1, maxcompete=0,
                                                                 maxsurrogate=0), model=TRUE)
  

  terminalNodes <- sort(unique(base_tree$where))
  locTest <- terminalNodes[1:2]
  if (locTest[2] == locTest[1]+1) {
    pvals[i] <- getSplitPval(base_tree, locTest, sigma_y)
    pvals2[i] <- getSplitPval_CAT(base_tree, locTest, sigma_y)
  } else {
  pvals[i] <- NA
  pvals2[i] <- NA
  }
}


mean(pvals < 0.05, na.rm=TRUE)

qqsample <- sort(pvals[!is.na(pvals)])
qqtheory <- qunif(seq(0,1,length.out=length(pvals[!is.na(pvals)])))
plot(qqsample, qqtheory)
abline(0,1, col="red")

qqsample2 <- sort(pvals2[!is.na(pvals2)])
qqtheory2 <- qunif(seq(0,1,length.out=length(pvals2[!is.na(pvals2)])))
plot(qqsample2, qqtheory2)
abline(0,1, col="red")

