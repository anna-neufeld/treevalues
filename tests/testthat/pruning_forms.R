setwd("~/treevalues")
devtools::load_all()
library(rpart)


n <- 200
p <- 10
sigma_y <- 5
X <- MASS::mvrnorm(n, rep(0,p), diag(rep(1,p)))
mu_y <- 5*(X[,1] <= 0) + 2*(X[,2]>0)
y <- rnorm(n, mu_y, sigma_y)
  
dat <- data.frame(y=y,X=X)
nameX <- sapply(1:p, function(u) paste0("X",u))
names(dat) = c("y", nameX)
  
  ### Build an rpart of depth d
base_tree <- rpart::rpart(y~., data=dat, model=TRUE, control=rpart.control(maxdepth = 4,
 minsplit=2, minbucket=1,cp=0.03, maxcompete=0, maxsurrogate=0))
node <- sort(base_tree$where)[1]
if ((node+1) %in% unique(base_tree$where)) {
nu <- (base_tree$where==node)/sum(base_tree$where==node)
nu2 <- (base_tree$where==node)/sum(base_tree$where==node) - (base_tree$where==(node+1))/sum(base_tree$where==(node+1))
splits <- getAncestors(base_tree, node)
print(getInterval_full(base_tree,nu,splits))
print(getInterval_full(base_tree,nu2,splits))
mean(y[base_tree$where==(node+1)])
}


  CI1 <- getNodeInterval(base_tree,node,sigma_y,permutations=FALSE)
  CI2 <- getNodeInterval(base_tree,node,sigma_y,permutations=TRUE)
  length2 <- CI2[2]-CI2[1]
  length1 <- CI1[2]-CI1[1]
  expect_true(length1 >= length2)
  
  splits <- getAncestors(base_tree, node)
  y1 <- y[base_tree$where==node]
  nu <- (base_tree$where==node)/sum((base_tree$where==node))
  true_signal <- abs(t(nu)%*%mu_y)
  phi_bounds_direct <- getInterval_full(base_tree,nu, splits)
  phi_bounds_indirect <- suppressWarnings(interval_intersection(
    getInterval_build(base_tree,nu, splits),
    phi_bounds_prune <- getInterval_prune(base_tree,nu, splits)
  ))
  expect_true(all.equal(phi_bounds_direct,phi_bounds_indirect))




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
  phi_bounds1 <- getInterval_full(base_tree,nu, splits)
  phi_bounds2 <- suppressWarnings(interval_intersection(
    getInterval_build(base_tree,nu, splits),
    getInterval_prune(base_tree,nu, splits)
  ))
  pTree <- getSplitPval(base_tree, locTest, sigma_y)
  expect_true((pTree-0.8460218)<1e-6)
  expect_true(pTree==correctPVal(phi_bounds1, nu, y, sigma_y))




  