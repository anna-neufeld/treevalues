library(rpart)
library(treevalues)
library(intervals)

n <- 110
p <- 7
sigma_y <- 5
nTrials <- 1000
pvals <- rep(0, nTrials)
pvals2 <- rep(0, nTrials)
lengthS <- rep(0, nTrials)
lengthS2 <- rep(0, nTrials)
CIs <- rep(0, nTrials)
CIs2 <- rep(0, nTrials)
cov <- rep(0, nTrials)
cov2 <- rep(0, nTrials)

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
  base_tree <- rpart::rpart(y~., data=dat, control=rpart.control(maxdepth = 4,
                                                                 minsplit=1, minbucket=1,
                                                                 cp=0.02, maxcompete=0,
                                                                 maxsurrogate=0), model=TRUE)
  if (length(unique(base_tree$where))>1){
    region <- sample(row.names(base_tree$frame)[-1], size=1)
    res1 <- regInference(base_tree,region,permutation=TRUE,sigma_y=sigma_y,computeCI=TRUE)
    res2 <- regInference(base_tree,region,permutation=FALSE,sigma_y=sigma_y,computeCI=TRUE)
    pvals[i] <- res1$pval
    pvals2[i] <- res2$pval
    lengthS[i] <- sum(size(res1$condset))
    lengthS2[i] <- sum(size(res2$condset))
    CIs[i] <- res1$confint[2]-res1$confint[1]
    CIs2[i] <- res2$confint[2]-res2$confint[1]
    cov[i] <- res1$confint[2] > 0 & res1$confint[1] < 0
    cov2[i] <-  res2$confint[2] > 0 & res2$confint[1] < 0
  }
}

all.equal(pvals,pvals2)

par(mfrow=c(1,2))
qqsample <- sort(pvals[!is.na(pvals)])
qqtheory <- qunif(seq(0,1,length.out=length(pvals[!is.na(pvals)])))
plot(qqsample, qqtheory)
abline(0,1, col="red")

qqsample2 <- sort(pvals2[!is.na(pvals2)])
qqtheory2 <- qunif(seq(0,1,length.out=length(pvals2[!is.na(pvals2)])))
plot(qqsample2, qqtheory2)
abline(0,1, col="red")

