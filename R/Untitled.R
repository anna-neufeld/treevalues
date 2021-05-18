# library(rpart)
# library(treevalues)
# library(intervals)
#
# n <- 120
# p <- 7
# sigma_y <- 5
# nTrials <- 10
# pvals <- rep(0, nTrials)
# pvals2 <- rep(0, nTrials)
# lengthS <- rep(0, nTrials)
# lengthS2 <- rep(0, nTrials)
#
# for (i in 1:nTrials) {
#   set.seed(i)
#   print(i)
#
#   X <- MASS::mvrnorm(n, mu=rep(0,p), Sigma=diag(1, nrow=p, ncol=p))
#   y <- rnorm(n, 0, sigma_y)
#
#   ### This turns out to be necessary to reading the split rules.
#   dat <- data.frame(y=y,X=X)
#   nameX <- sapply(1:p, function(u) paste0("X",u))
#   names(dat) = c("y", nameX)
#
#   ### Build an rpart of depth d
#   base_tree <- rpart::rpart(y~., data=dat, control=rpart.control(maxdepth = 2,
#                                                                  minsplit=1, minbucket=1,
#                                                                  cp=0.0, maxcompete=0,
#                                                                  maxsurrogate=0), model=TRUE)
#
#   terms <- row.names(base_tree$frame)[base_tree$frame$var == "<leaf>"]
#
#   region <- sample(terms, size=1)
#   branch <- getBranch(base_tree, region)
#
#   #print(branchInference(base_tree,branch)$condset)
#   #print(branchInference(base_tree,branch[2:1])$condset)
#   print(regInference(base_tree,region,FALSE,computeCI=FALSE)$condset)
#   print(regInference(base_tree,region,TRUE,computeCI=FALSE)$condset)
# }
#
# all.equal(pvals,pvals2)
#
# par(mfrow=c(1,2))
# qqsample <- sort(pvals[!is.na(pvals)])
# qqtheory <- qunif(seq(0,1,length.out=length(pvals[!is.na(pvals)])))
# plot(qqsample, qqtheory)
# abline(0,1, col="red")
#
# qqsample2 <- sort(pvals2[!is.na(pvals2)])
# qqtheory2 <- qunif(seq(0,1,length.out=length(pvals2[!is.na(pvals2)])))
# plot(qqsample2, qqtheory2)
# abline(0,1, col="red")
#
