library(rpart)
library(treevalues)
library(intervals)

n <- 150
p <- 8
sigma_y <- 5
nTrials <- 1000
pvals <- rep(0, nTrials)
pvals2 <- rep(0, nTrials)

for (i in 1:nTrials) {
  set.seed(i+nTrials)
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
  base_tree <- rpart::rpart(y~., data=dat, control=rpart.control(maxdepth = 5,
                                                                 minsplit=1, minbucket=1,
                                                                 cp=0.03, maxcompete=0,
                                                                 maxsurrogate=0), model=TRUE)
  base_tree2 <- rpart::rpart(y~., data=dat, control=rpart.control(maxdepth = 5,
                                                                 minsplit=1, minbucket=1,
                                                                 cp=0.03), model=TRUE)
  if (length(unique(base_tree$where))>1){
    region <- sample(row.names(base_tree$frame)[-1], size=1)
    branch <-getBranch(base_tree, region)
    pvals[i] <- branchInference(base_tree,branch,type="sib",sigma_y=sigma_y,computeCI=FALSE)$pval

    temp_tree <- base_tree
    temp_tree$cp <- 0
    temp_tree$control$cp <- 0
    pvals2[i] <- branchInference(temp_tree,branch,type="sib",sigma_y=sigma_y,computeCI=FALSE)$pval

  }
}


#par(mfrow=c(1,2))
#qqsample <- sort(pvals[!is.na(pvals)])
#qqtheory <- qunif(seq(0,1,length.out=length(pvals[!is.na(pvals)])))
#plot(qqsample, qqtheory)
#abline(0,1, col="red")

#qqsample2 <- sort(pvals2[!is.na(pvals2)])
#qqtheory2 <- qunif(seq(0,1,length.out=length(pvals2[!is.na(pvals2)])))
#plot(qqsample2, qqtheory2)
#abline(0,1, col="red")

