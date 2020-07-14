for (i in 2:500) {
  set.seed(i)
  n <- sample(100:300, size=1)
  p <- sample(3:20, size=1)
  depth <- sample(2:4, size=1)
  sigma_y=5
  X <- MASS::mvrnorm(n, rep(0,p), diag(rep(1,p)))
  beta <- 5

  # A mechanism that includes mostly real splits and some non-real splits
  mu_y <- beta*I(X[,1] > 0)
  mu_y <- mu_y + 2*beta*(I(X[,1] > 0 & X[,2] > 0)) - 2*beta*(I(X[,1] < 0 & X[,2] > 0))
  mu_y <- mu_y + beta*I(X[,3] > 0 & X[,2] > 0 & X[,1] > 0) - beta*I(X[,3] > 0 & X[,2] < 0 & X[,1] < 0)


y <- rnorm(n, mu_y, sigma_y)

# This helps read the split rules
dat <- data.frame(y=y,X=X)
nameX <- sapply(1:p, function(u) paste0("X",u))
names(dat) = c("y", nameX)

base_tree <- rpart::rpart(y~., data=dat, model=TRUE,
                          control=rpart.control(maxdepth = depth,
                                                minsplit=2, minbucket=1,
                                                cp=-1, maxcompete=0,maxsurrogate=0))

terminalNodes <- sort(unique(base_tree$where))
locTest <- terminalNodes[1:2]
if(locTest[2] != locTest[1]+1) {
  next
}
splits <- getAncestors(base_tree, locTest[1])
nu <- (base_tree$where ==locTest[1])/sum(base_tree$where==locTest[1]) -  (base_tree$where==locTest[2])/sum(base_tree$where==locTest[2])

if (all.equal(getInterval(base_tree, nu, splits),getInterval_EXPERIMENT(base_tree, nu, splits))==TRUE) {
print(i)
  } else {
    print(getInterval(base_tree, nu, splits))
    print(getInterval_EXPERIMENT(base_tree, nu, splits))
  }
}



for (n in seq(100,2000, by=100)) {
  for (p in seq(5,30, by=5)) {
    #set.seed(i)
    #n <- sample(50:1000, size=1)
    #p <- sample(3:20, size=1)
    depth <- 3
    sigma_y=5
    X <- MASS::mvrnorm(n, rep(0,p), diag(rep(1,p)))
    beta <- 5

    # A mechanism that includes mostly real splits and some non-real splits
    mu_y <- beta*I(X[,1] > 0)
    mu_y <- mu_y + 2*beta*(I(X[,1] > 0 & X[,2] > 0)) - 2*beta*(I(X[,1] < 0 & X[,2] > 0))
    mu_y <- mu_y + beta*I(X[,3] > 0 & X[,2] > 0 & X[,1] > 0) - beta*I(X[,3] > 0 & X[,2] < 0 & X[,1] < 0)


    y <- rnorm(n, mu_y, sigma_y)

    # This helps read the split rules
    dat <- data.frame(y=y,X=X)
    nameX <- sapply(1:p, function(u) paste0("X",u))
    names(dat) = c("y", nameX)

    base_tree <- rpart::rpart(y~., data=dat, model=TRUE,
                              control=rpart.control(maxdepth = depth,
                                                    minsplit=2, minbucket=1,
                                                    cp=-1, maxcompete=0,maxsurrogate=0))

    terminalNodes <- sort(unique(base_tree$where))
    locTest <- terminalNodes[1:2]
    if(locTest[2] != locTest[1]+1) {
      next
    }
    splits <- getAncestors(base_tree, locTest[1])
    nu <- (base_tree$where ==locTest[1])/sum(base_tree$where==locTest[1]) -  (base_tree$where==locTest[2])/sum(base_tree$where==locTest[2])

    print(c(n,p))
    print(system.time(getInterval_EXPERIMENT(base_tree, nu, splits))[1]-system.time(getInterval(base_tree, nu, splits))[1] > 0)
  }
}
