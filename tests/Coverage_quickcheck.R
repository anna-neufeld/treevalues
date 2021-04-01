library(treevalues)
library(rpart)
n <- 100
p <- 10
sigma_y <- 5
nTrials <- 3
pvals <- rep(0, nTrials)

CIs <- matrix(NA, nrow=nTrials, ncol=2)
truths <- rep(NA, nTrials)

for (i in 1:nTrials) {
  set.seed(i)
  print(i)

  X <- MASS::mvrnorm(n, rep(0,p), diag(rep(1,p)))
  beta = 10
  mu_y <- beta*I(X[,1] > 0)
  #mu_y_2 <- mu_y_1 + 2*beta*(I(X[,1] > 0 & X[,2] > 0)) - 2*beta*(I(X[,1] < 0 & X[,2] > 0))
  #mu_y <- mu_y_2 + beta*I(X[,3] > 0 & X[,2] > 0 & X[,1] > 0) - beta*I(X[,3] > 0 & X[,2] < 0 & X[,1] < 0)

  y <- rnorm(n, mu_y, sigma_y)

  ### This turns out to be necessary to reading the split rules.
  dat <- data.frame(y=y,X=X)
  nameX <- sapply(1:p, function(u) paste0("X",u))
  names(dat) = c("y", nameX)

  ### Build an rpart of depth d
  base_tree <- rpart::rpart(y~., data=dat, control=rpart.control(maxdepth = 1,
                                                                 minsplit=1, minbucket=1,
                                                                 cp=0.15, maxcompete=0,
                                                                 maxsurrogate=0), model=TRUE)


  terminalNodes <- sort(unique(base_tree$where))
  if (length(terminalNodes) > 1) {
  node <- sample(terminalNodes, 1)
  #splits <- getAncestors(base_tree, node)
  CIs[i,] <- nodeInference(base_tree, node, sigma_y)$confint
  truths[i] <- mean(mu_y[base_tree$where==node])
  }

}

mean(CIs[,1] < truths & CIs[,2] > truths)
#qqsample <- sort(pvals[pvals!=0])
#qqtheory <- qunif(seq(0,1,length.out=sum(pvals!=0)))
#plot(qqsample, qqtheory)
#abline(0,1, col="red")
