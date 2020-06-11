context("basic testing")
library(treevalues)

test_that("Basic Hypothesis Tests; Null Model", {
  set.seed(2)
  reject <- 0
  for (j in 1:10) {
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
                                                                 maxsurrogate=0))

  ### For each pair of terminal nodes, compute results.
  terminalNodes <- sort(unique(base_tree$where))
  if (length(terminalNodes)==4) {
    locs = c(1,3)
  } else{
    locs = c(1)
  }
  for (i in locs) {
    locTest = c(terminalNodes[i], terminalNodes[i+1])
    splits <- getSplits(base_tree, locTest)
    y1 <- y[base_tree$where==locTest[1]]
    y2 <- y[base_tree$where==locTest[2]]
    nu <- (base_tree$where==locTest[1])/sum((base_tree$where==locTest[1])) - (base_tree$where==locTest[2])/sum(base_tree$where==locTest[2])
    true_signal <- abs(t(nu)%*%mu_y)
    phi_bounds <- getInterval(base_tree, X,y, locTest, splits, dat)
    pTree <- correctPVal(phi_bounds, nu, y, sigma_y)
    if (pTree < 0.1) {
      reject = reject + 1
    }
  }
  }
  expect_equal(reject,2)
})
