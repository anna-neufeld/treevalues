get_tree_design_matrix <- function(tree, node) {
  splits <- getAncestors(tree, node)
  dat <- tree$model
  design.mat <- matrix(1, ncol=1, nrow=nrow(dat))
  for (l in 1:length(splits)) {
    cy<- eval(parse(text = splits[l]))*design.mat[,l]
    design.mat <- cbind(design.mat, cy)
  }
  design.mat
}


one_rep_T <- function(seed, beta=0, depth=3) {
  set.seed(seed)
  n <- 200
  p <- 10
  sigma_y= 5
  X <- MASS::mvrnorm(n, rep(0,p), diag(rep(1,p)))
  mu_y <- beta*I(X[,1] > 0)
  mu_y <- mu_y + 2*beta*(I(X[,1] > 0 & X[,2] > 0)) - 2*beta*(I(X[,1] < 0 & X[,2] > 0))
  mu_y <- mu_y + beta*I(X[,3] > 0 & X[,2] > 0 & X[,1] > 0) - beta*I(X[,3] > 0 & X[,2] < 0 & X[,1] < 0)
  y <- rnorm(n, mu_y, sigma_y)
  dat <- data.frame(y=y,X=X)
  nameX <- sapply(1:p, function(u) paste0("X",u))
  names(dat) = c("y", nameX)

  tree1 <- rpart::rpart(y~., data=dat, model=TRUE, control=rpart.control(maxdepth = depth,
                                                                       minsplit=2, minbucket=1,
                                                                       cp=-1, maxcompete=0,maxsurrogate=0))

  terminalNodes <- sort(unique(tree1$where))
  locTest <- terminalNodes[1:2]
  nu <- (tree1$where ==locTest[1])/sum(tree1$where==locTest[1]) -  (tree1$where==locTest[2])/sum(tree1$where==locTest[2])
  splits <- getAncestors(tree1, locTest[1])


  Pi_perp <- diag(rep(1,n)) - nu%*%t(nu)/sum(nu^2)

  C <- get_tree_design_matrix(tree1, locTest[1])
  Clast <- C[,1:(NCOL(C)-1)]
  Hl <- diag(1,n) - C%*%solve(t(C)%*%C)%*%t(C)
  Hlast <- diag(1,n) - Clast%*%solve(t(Clast)%*%Clast)%*%t(Clast)
  c <- C[,NCOL(C)]

  alpha <- sqrt(as.numeric(t(y)%*%Hl%*%y)/((n-p-1)*as.numeric(t(c)%*%Hlast%*%c)))

  phi_bounds <- getInterval_T(tree1, nu, splits, alpha)

  ### Once we have the bounds: TRUNCATED T DISTRIBUTION.
  ### probably need to deal with some numerical stability things??

  sampS <- as.numeric(sqrt(1/(t(c)%*%Hlast%*%c)*t(y)%*%Hl%*%y/(n-p-1)))
  samp_sig <- as.numeric(t(nu)%*%y)

  true_sig <- as.numeric(t(nu)%*%mu_y)

  samp_T <- samp_sig/sampS
  test_interval <- interval_complement(Intervals(c(-abs(samp_T), abs(samp_T))))
  numerator <- interval_intersection(test_interval, phi_bounds)

  areaNum <- sum(apply(numerator, 1, function(u) pt(u[2], df=n-p-1)-pt(u[1], df=n-p-1)))
  areaDenom <- sum(apply(phi_bounds, 1, function(u) pt(u[2], df=n-p-1)-pt(u[1], df=n-p-1)))
  pval = areaNum/areaDenom
  return(c(pval, true_sig))


  #t <- samp_T
  #yT <- Pi_perp%*%y + nu/sum(nu^2)*alpha*t

  #treeT <- rpart::rpart(yT~X, model=TRUE, control=rpart.control(maxdepth = depth,
  #                                                                       minsplit=2, minbucket=1,
  #                                                                       cp=-1, maxcompete=0,maxsurrogate=0))

  #all.equal(tree1$where, treeT$where)

}







