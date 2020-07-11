n=200
p=10
sigma_y=5
X <- MASS::mvrnorm(n, rep(0,p), diag(rep(1,p)))

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
                          control=rpart.control(maxdepth = 3,
                                                minsplit=2, minbucket=1,
                                                cp=-1, maxcompete=0,maxsurrogate=0))

terminalNodes <- sort(unique(base_tree$where))
locTest <- terminalNodes[1:2]
splits <- getAncestors(base_tree, locTest[1])

nu <- (base_tree$where ==locTest[1])/sum(base_tree$where==locTest[1]) -  (base_tree$where==locTest[2])/sum(base_tree$where==locTest[2])


print(all.equal(getInterval(base_tree, nu, splits),getInterval_OLD(base_tree, nu, splits)))
print(system.time(getInterval_OLD(base_tree, nu, splits)))
print(system.time(getInterval(base_tree, nu, splits)))

