setwd("~/treevalues")
devtools::load_all()
library(rpart)


nTrials <- 10
pvals1 <- rep(NA, nTrials)

for (i in 1:nTrials) {
  print(i)
  set.seed(i)

  n <- 150
  p <- 8
  sigma_y <- 5
  X <- cbind(MASS::mvrnorm(n, rep(0,p/2), diag(rep(1,p/2))), matrix(rbinom(n*p/2,size=1,p=0.3), ncol=p/2))
  mu_y <- 0
  y <- rnorm(n, mu_y, sigma_y)

  dat <- data.frame(y=y,X=X)
  nameX <- sapply(1:p, function(u) paste0("X",u))
  names(dat) = c("y", nameX)

  ### Build an rpart of depth d
  base_tree <- rpart::rpart(y~., data=dat, model=TRUE, control=rpart.control(maxdepth = 4, minbucket=10,cp=0.01, maxcompete=0,maxsurrogate=0))
  if (NROW(base_tree$frame) > 1) {

  region <- sample(row.names(base_tree$frame)[-1],size=1)
  splits <- getBranch(base_tree, region)
  membership <-  getRegion(base_tree,region)
  nu <- (as.numeric(membership))/sum(membership)

  fullInt1 <- getInterval(base_tree,nu, splits,grow=FALSE)
  #fullInt2 <- getInterval2(base_tree,nu, splits,grow=FALSE)

  pvals1[i] <- correctPVal(fullInt1, nu, y, sigma_y)
  #pvals2[i] <- correctPVal(fullInt2, nu, y, sigma_y)
  }
}

#all.equal(pvals1,pvals2)

 #par(mfrow=c(1,2))
 qqsample <- sort(pvals1[!is.na(pvals1)])
 qqtheory <- qunif(seq(0,1,length.out=length(pvals1[!is.na(pvals1)])))
 plot(qqsample, qqtheory)
 abline(0,1, col="red")
#
 #qqsample2 <- sort(pvals2[!is.na(pvals2)])
 #qqtheory2 <- qunif(seq(0,1,length.out=length(pvals2[!is.na(pvals2)])))
 #plot(qqsample2, qqtheory2)
#abline(0,1, col="red")
#
