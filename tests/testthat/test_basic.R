context("basic testing")
library(treevalues)

test_that("Region: do different methods match??", {
  #library(rpart)
  #library(intervals)
  #library(rpart.plot)
  #library(partykit)

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
  base_tree <- rpart::rpart(y~., data=dat, model=TRUE, control=rpart.control(maxdepth = 4, minsplit=2, minbucket=1,cp=0.05, maxcompete=0,maxsurrogate=0))
  unpruned_tree <- rpart::rpart(y~., data=dat, model=TRUE, control=rpart.control(maxdepth = 4, minsplit=2, minbucket=1,cp=0.0, maxcompete=0,maxsurrogate=0))

  region <- sample(row.names(base_tree$frame)[-1],size=1)
  splits <- getBranch(base_tree, region)
  membership <-  getRegion(base_tree,region)
  nu <- (as.numeric(membership))/sum(membership)

  fullInt <- getInterval_full(base_tree,nu, splits,grow=FALSE)
  growInt <- getInterval_full(base_tree,nu, splits,grow=TRUE)

  growInt2 <- getInterval_full(unpruned_tree,nu, splits,grow=FALSE)

  pruneInt <- suppressWarnings(interval_difference(growInt,fullInt))

  expect_true(length(suppressWarnings(interval_difference(fullInt,growInt)))==0)
  expect_true(length(suppressWarnings(interval_difference(growInt2,growInt)))==0)
  expect_true(length(suppressWarnings(interval_difference(growInt,growInt2)))==0)

  }
)

test_that("Basic Hypothesis Tests; Null Model", {
  set.seed(153)
  n <- 150
  p <- 8
  sigma_y <- 5
  X <- MASS::mvrnorm(n, rep(0,p), diag(rep(1,p)))
  mu_y <- 0*X[,1]
  y <- rnorm(n, mu_y, sigma_y)

  ### This turns out to be necessary to reading the split rules.
  dat <- data.frame(y=y,X=X)
  nameX <- sapply(1:p, function(u) paste0("X",u))
  names(dat) = c("y", nameX)

  ### Build an rpart of depth d
  base_tree <- rpart::rpart(y~., data=dat, control=rpart.control(maxdepth = 2,
                                                                 minsplit=2, minbucket=1,
                                                                 cp=-1, maxcompete=0,
                                                                 maxsurrogate=0), model=TRUE)

  ### For each pair of terminal nodes, compute results.
  terminalNodes <- sort(unique(base_tree$where))



  locTest = c(terminalNodes[1], terminalNodes[2])
  splits <- getBranch(base_tree, locTest[1])
  y1 <- y[base_tree$where==locTest[1]]
  y2 <- y[base_tree$where==locTest[2]]
  nu <- (base_tree$where==locTest[1])/sum((base_tree$where==locTest[1])) - (base_tree$where==locTest[2])/sum(base_tree$where==locTest[2])
  true_signal <- abs(t(nu)%*%mu_y)
  phi_bounds1 <- getInterval_full(base_tree,nu, splits)
  #pTree <- splitInference(base_tree, locTest, sigma_y)$pval
  #expect_true((pTree-0.8460218)<1e-6)
  #expect_true(pTree==correctPVal(phi_bounds1, nu, y, sigma_y))
})


test_that("Full Inference", {
  set.seed(1)
  n <- 200
  p <- 10
  sigma_y <- 5
  X <- MASS::mvrnorm(n, rep(0,p), diag(rep(1,p)))
  mu_y <- 0*X[,1]
  y <- rnorm(n, mu_y, sigma_y)

  ### This turns out to be necessary to reading the split rules.
  dat <- data.frame(y=y,X=X)
  nameX <- sapply(1:p, function(u) paste0("X",u))
  names(dat) = c("y", nameX)

  ### Build an rpart of depth d
  base_tree <- rpart::rpart(y~., data=dat, model=TRUE, control=rpart.control(maxdepth = 3,
                                                                             minsplit=2, minbucket=1,
                                                                             cp=-1, maxcompete=0,
                                                                             maxsurrogate=0))
  mat <- fullTreeInference(base_tree, sigma_y)
  expect_true(NROW(mat)==length(unique(base_tree$where)))

  treeval.plot(base_tree, mat, sigma_y)
}
)
