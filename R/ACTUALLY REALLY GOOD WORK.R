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
                           control=rpart.control(maxdepth = depth,
                                                 minsplit=2, minbucket=1,
                                                 cp=-1, maxcompete=0,maxsurrogate=0))

locTest = c(11,12)
terminalNodes <- sort(unique(base_tree$where))
splits <- getAncestors(base_tree, 4)


  
  ### rpart orderes things nicely YAY
  dat <- base_tree$model
  y <- dat[,1]
  X <- dat[,-1]
  n <- nrow(X)
  p <- NCOL(X)
  
  nu <- (base_tree$where ==4)/sum(base_tree$where==5) -  (base_tree$where==4)/sum(base_tree$where==5)
  
  normnusq <- sum(nu^2)
  
  Pi_perp <- diag(rep(1,n)) - nu%*%t(nu)/sum(nu^2)
  C <- (Pi_perp%*%y)%*%t((Pi_perp%*%y))
  B <- (Pi_perp%*%y%*%t(nu) + nu%*%t(y)%*%Pi_perp)/sum(nu^2)
  A <- (nu%*%t(nu))/sum(nu^2)^2
  
  #### Overall Indicator Vars
  cs <- matrix(NA, nrow=n*p, ncol=n)
  for (j in 1:p) {
    cs[((j-1)*n+1):(j*n),] <- sapply(X[,j], function(u) (u <= X[,j]))
  }
  ### Remove places that will lead us to divide by 0. They are trivial splits.
  cs <- cs[rowSums(cs)!=n,]
  
  
  ### In the first layer, nothing gets zeroed out
  nulled <- rep(1, n)
  C_prev <- rep(1,n)
  
  ### This has such a simple form. I should write it in this way. 
  H_prev <- matrix(-1/n, nrow=n, ncol=n)
  diag(H_prev) <- 1-1/n
  #H_prev <- diag(1,n) - C_prev%*%solve(t(C_prev)%*%C_prev)%*%t(C_prev)
  
  full_interval = Intervals(c(-Inf, Inf))
  
  for (i in 1:length(splits)) {
    cy<- eval(parse(text = splits[i]))*nulled
    
    
    vec <- runif(200)
    if (NCOL(C_prev) > 1) {
    regions <- rowSums(C_prev)
    } else{regions <- C_prev}
    means <- sapply(1:max(regions), function(u) mean(vec[regions==u]))
    resids <- vec - means[regions]
    print(all.equal(as.numeric(H_prev%*%vec), as.numeric(resids)))
    
    
    c_now <- t(apply(cs, 1, function(u) u*nulled))
    c_now <- c_now[rowSums(c_now)!=0,]
    c_now <- c_now[rowSums(c_now)!=sum(nulled),]
    
    
    ### Faster version
    denoms2 <- apply(c_now, 1, function(u) colSums(u*(H_prev%*%u)))
    
    leafSize <- sum(cy)
    num1s <- rowSums(c_now) ## O(np*n)
    #num2s <- rowSums(nulled-c_now) ## O(np*n)
    
    denoms <- num1s*(leafSize-num1s)/(leafSize) ## O(np*1)
    
    all.equal(denoms, denoms2)
    ### denomoniator is same for everyone. It's leaf size. 
    
    
    avec <- apply(c_now, 1, function(u) colSums(u*(hAh%*%u)))/denoms - partialA
    bvec <- apply(c_now, 1, function(u) colSums(u*(hBh%*%u)))/denoms - partialB
    cvec <- apply(c_now, 1, function(u) colSums(u*(hCh%*%u)))/denoms - partialC

    #### SOMEHOW COMBINE THIS RES WITH PREVIOUS RES.
    

    #### If this inversion is slow could use the rank 1 update!!!
    #### Update these things for next level.
    C_prev <- cbind(C_prev, cy)
    
    ### OH FIX THIS. DUH. LOW HANGING FRUIT. 
    Hc <- H_prev%*%cy
    #H_prev <- H_prev - ()
      
      
    H_prev <- diag(1,n) - C_prev%*%solve(t(C_prev)%*%C_prev)%*%t(C_prev)
    
    ### Shoot. If we want to set "nulled" in this way- we really need to be only working with a single brach of
    ### tree. Cannot condition on extra branches.
    nulled <- cy
  }
 
