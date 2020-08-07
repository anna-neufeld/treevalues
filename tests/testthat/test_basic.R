context("basic testing")
library(treevalues)

test_that("Basic Conf Int Terminal Node", {
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
  base_tree <- rpart::rpart(y~., data=dat, model=TRUE, control=rpart.control(maxdepth = 2,
                                                                 minsplit=2, minbucket=1,
                                                                 cp=-1, maxcompete=0,
                                                                 maxsurrogate=0))
  node <- base_tree$where[1]
  splits <- getAncestors(base_tree, node)
  y1 <- y[base_tree$where==node]
  nu <- (base_tree$where==node)/sum((base_tree$where==node))
  true_signal <- abs(t(nu)%*%mu_y)
  phi_bounds <- getInterval(base_tree,nu, splits)
  bounds2 <- Intervals(as.matrix(phi_bounds)/(sigma_y*sqrt(sum(nu^2))))
  CI <- getNodeInterval(base_tree, node, sigma_y)
  expect_true(mean(y1) > CI[1] & mean(y1) < CI[2])
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
  splits <- getAncestors(base_tree, locTest[1])
  y1 <- y[base_tree$where==locTest[1]]
  y2 <- y[base_tree$where==locTest[2]]
  nu <- (base_tree$where==locTest[1])/sum((base_tree$where==locTest[1])) - (base_tree$where==locTest[2])/sum(base_tree$where==locTest[2])
  true_signal <- abs(t(nu)%*%mu_y)
  phi_bounds <- getInterval(base_tree,nu, splits)
  pTree <- getSplitPval(base_tree, locTest, sigma_y)
  expect_true(pTree-0.8460218<1e-6)
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

  #expect_equal(sum(mat$lower < mat$effectSize &  mat$effectSize < mat$upper), NROW(mat))
  expect_equal(sum(mat$branch1lower < mat$branch1mean &  mat$branch1mean < mat$branch1upper), NROW(mat))
  expect_equal(sum(mat$branch2lower < mat$branch2mean &  mat$branch2mean < mat$branch2upper), NROW(mat))

  treeval.plot(base_tree, mat, sigma_y)
}
)


test_that("Two Versions", {
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
  base_tree <- rpart::rpart(y~., data=dat, model=TRUE, control=rpart.control(maxdepth = 4,
                                                                             minsplit=2, minbucket=1,
                                                                             cp=-1, maxcompete=0,
                                                                             maxsurrogate=0))

  terminalNodes <- sort(unique(base_tree$where))
  locTest = c(terminalNodes[1], terminalNodes[2])
  splits <- getAncestors(base_tree, locTest[1])
  nu <- (base_tree$where==locTest[1])/sum((base_tree$where==locTest[1])) - (base_tree$where==locTest[2])/sum(base_tree$where==locTest[2])
  int1 <- getInterval(base_tree, nu, splits)
  int2 <- getInterval_OLD(base_tree, nu, splits)
  expect_true(sum(as.matrix(int1)-as.matrix(int2)) < 1e-8)

}
)
