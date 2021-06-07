context("basic testing")
library(treevalues)
library(rpart)

test_that("Region: do different methods match??", {

  n <- 150
  p <- 8
  sigma_y <- 5
  X <- MASS::mvrnorm(n, rep(0,p), diag(rep(1,p)))
  mu_y <- 0
  y <- rnorm(n, mu_y, sigma_y)

  dat <- data.frame(y=y,X=X)
  nameX <- sapply(1:p, function(u) paste0("X",u))
  names(dat) = c("y", nameX)

  ### Build an rpart of depth d
  base_tree <- rpart(y~., data=dat, model=TRUE, control=rpart.control(maxdepth = 4, minsplit=2, minbucket=1,cp=0.01))
  unpruned_tree <- rpart(y~., data=dat, model=TRUE, control=rpart.control(maxdepth = 4, minsplit=2, minbucket=1,cp=0.0))

  region <- sample(row.names(base_tree$frame)[-1],size=1)
  splits <- getBranch(base_tree, region)
  membership <-  getRegion(base_tree,region)
  nu <- (as.numeric(membership))/sum(membership)

  fullInt <- getInterval(base_tree,nu, splits,grow=FALSE)
  growInt <- getInterval(base_tree,nu, splits,grow=TRUE)

  growInt2 <- getInterval(unpruned_tree,nu, splits,grow=FALSE)

  pruneInt <- suppressWarnings(interval_difference(growInt,fullInt))

  expect_true(length(suppressWarnings(interval_difference(fullInt,growInt)))==0)
  expect_true(length(suppressWarnings(interval_difference(growInt2,growInt)))==0)
  expect_true(length(suppressWarnings(interval_difference(growInt,growInt2)))==0)

  }
)

test_that("Basic Hypothesis Tests; Siblings under Null!!!!", {
  set.seed(153)
  n <- 150
  p <- 8
  sigma_y <- 5
  X <- MASS::mvrnorm(n, rep(0,p), diag(rep(1,p)))
  mu_y <- 0*X[,1]
  y <- rnorm(n, mu_y, sigma_y)

  dat <- data.frame(y=y,X=X)
  nameX <- sapply(1:p, function(u) paste0("X",u))
  names(dat) = c("y", nameX)

  ### Build an rpart of depth d
  base_tree <- rpart(y~X1+X3+X5+X7, data=dat, control=rpart.control(cp=0.03, maxcompete=0, maxsurrogate=0), model=TRUE)

  ### For each pair of terminal nodes, compute results.
  region <- sample(row.names(base_tree$frame)[-1],size=1)
  if (as.numeric(region)%%2==0) {
    sib = as.numeric(region)+1
  } else { sib = as.numeric(region)-1}


  splits <- getBranch(base_tree, region)
  membership <-  getRegion(base_tree,region)
  membershipsib <-  getRegion(base_tree,sib)


  nu <- (as.numeric(membership))/sum(membership) - (as.numeric(membershipsib))/sum(membershipsib)

  phi_bounds1 <- getInterval(base_tree,nu, splits)
  pTree <- branchInference(base_tree, splits, "sib", sigma_y=sigma_y)$pval
  CI1 <- branchInference(base_tree, splits, "sib", sigma_y=sigma_y)$confint
  expect_true(pTree==correctPVal(phi_bounds1, nu, y, sigma_y))
  expect_true(all.equal(CI1, computeCI(nu, y, sigma_y,phi_bounds1, alpha=0.05)))
})


test_that("Permutation!!!", {
  n <- 150
  p <- 8
  sigma_y <- 5
  X <- MASS::mvrnorm(n, rep(0,p), diag(rep(1,p)))
  mu_y <- 0
  y <- rnorm(n, mu_y, sigma_y)

  dat <- data.frame(y=y,X=X)
  nameX <- sapply(1:p, function(u) paste0("X",u))
  names(dat) = c("y", nameX)

  ### Build an rpart of depth d
  base_tree <- rpart(y~., data=dat, model=TRUE, control=rpart.control(maxdepth = 4, minsplit=2, minbucket=1,cp=0.01))

  region <- sample(row.names(base_tree$frame)[-1],size=1)
  splits <- getBranch(base_tree, region)

  res <- branchInference(base_tree, splits, type="reg", permute=TRUE)
  res2 <- branchInference(base_tree, splits, type="reg", permute=FALSE)

  expect_true(length(suppressWarnings(interval_difference(res2$condset,res$condset)))==0)
  expect_true((res$confint[2]-res$confint[1]) <= (res2$confint[2]-res2$confint[1]))
})
