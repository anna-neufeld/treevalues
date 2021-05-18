n <- 150
p <- 8
sigma_y <- 5

X <- MASS::mvrnorm(n, mu=rep(0,p), Sigma=diag(1, nrow=p, ncol=p))
beta = 5
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
    