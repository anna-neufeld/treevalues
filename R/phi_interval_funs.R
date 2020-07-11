#' getPhiInterval
#' Should work for CIs or Hyp Tests.
#'
#' Also make this play nice with other people's real data
#'
#' Add warnings to ensure that their base_tree was build properly!!!
#' Would be a lot nicer if people only ha to pass dat. Or, better yet, only base_tree,
#' because they will use model=TRUE
getInterval_OLD <- function(base_tree, nu, splits) {

  ### rpart orderes things nicely YAY
  dat <- base_tree$model
  y <- dat[,1]
  X <- dat[,-1]
  n <- nrow(X)
  p <- NCOL(X)

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
  H_prev <- diag(1,n) - C_prev%*%solve(t(C_prev)%*%C_prev)%*%t(C_prev)
  full_interval = Intervals(c(-Inf, Inf))

  for (i in 1:length(splits)) {
    cy<- eval(parse(text = splits[i]))*nulled

    c_now <- t(apply(cs, 1, function(u) u*nulled))
    c_now <- c_now[rowSums(c_now)!=0,]
    c_now <- c_now[rowSums(c_now)!=sum(nulled),]


    hAh <- H_prev%*%A%*%H_prev
    hBh <- H_prev%*%B%*%H_prev
    hCh <- H_prev%*%C%*%H_prev

    denom <- colSums(cy*(H_prev%*%cy))
    partialA <- colSums(cy*(hAh%*%cy))/denom
    partialB <- colSums(cy*(hBh%*%cy))/denom
    partialC <- colSums(cy*(hCh%*%cy))/denom

    ### Faster version
    denoms <- apply(c_now, 1, function(u) colSums(u*(H_prev%*%u)))
    avec <- apply(c_now, 1, function(u) colSums(u*(hAh%*%u)))/denoms - partialA
    bvec <- apply(c_now, 1, function(u) colSums(u*(hBh%*%u)))/denoms - partialB
    cvec <- apply(c_now, 1, function(u) colSums(u*(hCh%*%u)))/denoms - partialC

    coeffs <- cbind(avec,bvec,cvec)
    phi_bounds <-  t(apply(coeffs, 1, getBounds))

    ### INSIDE INTERVALS
    numIn <- nrow(phi_bounds[phi_bounds[,3]==1,])
    if (!is.null(numIn)) {
      inside_comp_mat <- matrix(c(-Inf, Inf), nrow=2*numIn, ncol=2, byrow=TRUE)
      inside_comp_mat[1:numIn,2] <- phi_bounds[phi_bounds[,3]==1,1]
      inside_comp_mat[(numIn+1):(2*numIn), 1] <- phi_bounds[phi_bounds[,3]==1,2]
      ints_comp_inside <- Intervals(inside_comp_mat, closed=c(TRUE, TRUE))
      intersection1 <- interval_complement(interval_union(ints_comp_inside))
    } else {
      intersection1 <- Intervals(c(-Inf, Inf))
    }

    ### OUTSIDE INTERVALS
    ints_outside <- Intervals(phi_bounds[phi_bounds[,3]==0,1:2], closed=c(TRUE, TRUE))
    intersection2 <- interval_complement(interval_union(ints_outside))
    if(length(intersection1)==0) {
      res1 <- intersection2
    } else {
      res1 <- suppressWarnings(interval_intersection(intersection1, intersection2))
    }

    #### SOMEHOW COMBINE THIS RES WITH PREVIOUS RES.

    full_interval <- suppressWarnings(interval_intersection(full_interval, res1))

    #### If this inversion is slow could use the rank 1 update!!!
    #### Update these things for next level.
    C_prev <- cbind(C_prev, cy)
    H_prev <- diag(1,n) - C_prev%*%solve(t(C_prev)%*%C_prev)%*%t(C_prev)

    ### Shoot. If we want to set "nulled" in this way- we really need to be only working with a single brach of
    ### tree. Cannot condition on extra branches.
    nulled <- cy
  }

  return(full_interval)
}


#' getPhiInterval
#' Should work for CIs or Hyp Tests.
#'
#' Also make this play nice with other people's real data
#'
#' Add warnings to ensure that their base_tree was build properly!!!
#' Would be a lot nicer if people only ha to pass dat. Or, better yet, only base_tree,
#' because they will use model=TRUE
getInterval <- function(base_tree, nu, splits) {

  ### rpart orderes things nicely YAY
  dat <- base_tree$model
  y <- dat[,1]
  X <- dat[,-1]
  n <- nrow(X)
  p <- NCOL(X)

  Pi_perp <- diag(rep(1,n)) - nu%*%t(nu)/sum(nu^2)

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
  H_prev <- diag(1,n) - C_prev%*%solve(t(C_prev)%*%C_prev)%*%t(C_prev)
  full_interval = Intervals(c(-Inf, Inf))
  leafSize=n

  norm_nu <- sum(nu^2)
  y_proj <- (nu%*%t(nu)/norm_nu)%*%y
  nut_y <- t(nu)%*%y

  nunu <- nu%*%t(nu)/sum(nu^2)

  for (i in 1:length(splits)) {

    if (NCOL(C_prev) > 1) {
      regions <- rowSums(C_prev)
    } else{regions <- C_prev}

    Hnu <- nu - (sapply(1:max(regions), function(u) mean(nu[regions==u])))[regions]
    Hy <- y - (sapply(1:max(regions), function(u) mean(y[regions==u])))[regions]
    Hyproj <- y_proj - (sapply(1:max(regions), function(u) mean(y_proj[regions==u])))[regions]
    Hminus <- Hy-Hyproj

    cy<- eval(parse(text = splits[i]))*nulled

    Hcy <- cy - (sapply(1:max(regions), function(u) mean(cy[regions==u])))[regions]


    c_now <- t(apply(cs, 1, function(u) u*nulled))
    c_now <- c_now[rowSums(c_now)!=0,]
    c_now <- c_now[rowSums(c_now)!=sum(nulled),]

    denom <- colSums(cy*(H_prev%*%cy))

    indices <- which(cy==1)
    partialA <- (sum(Hnu[indices])^2/norm_nu^2)/denom
    partialB <- as.numeric((2*sum(Hy[indices])*(t(nu)%*%H_prev%*%cy)/
      norm_nu - 2*(t(cy)%*%Hnu)^2%*%nut_y/
      norm_nu^2)/denom)
    partialC <- (sum(Hminus[indices]))^2/denom

    num1s <- rowSums(c_now) ## O(np*n)
    denoms <- num1s*(leafSize-num1s)/(leafSize)

    ### Faster version
    other_As <- apply(c_now, 1, function(u) sum(Hnu[u==1])^2/norm_nu^2)
    other_Bs <- apply(c_now, 1, function(u) 2*(sum(Hy[u==1])%*%t(u)%*%Hnu)/
                        norm_nu
         - 2*((t(u)%*%Hnu)^2%*%nut_y)/
                        norm_nu^2)
    other_Cs <- apply(c_now, 1, function(u) sum(Hminus[u==1])^2)

    #bvec <- apply(c_now, 1, function(u) colSums(u*(hBh%*%u)))/denoms - partialB

    avec <- other_As/denoms-partialA
    bvec <- other_Bs/denoms-partialB
    cvec <- other_Cs/denoms-partialC

    coeffs <- cbind(avec,bvec,cvec)
    phi_bounds <-  t(apply(coeffs, 1, getBounds))

    ### INSIDE INTERVALS
    numIn <- nrow(phi_bounds[phi_bounds[,3]==1,])
    if (!is.null(numIn)) {
      inside_comp_mat <- matrix(c(-Inf, Inf), nrow=2*numIn, ncol=2, byrow=TRUE)
      inside_comp_mat[1:numIn,2] <- phi_bounds[phi_bounds[,3]==1,1]
      inside_comp_mat[(numIn+1):(2*numIn), 1] <- phi_bounds[phi_bounds[,3]==1,2]
      ints_comp_inside <- Intervals(inside_comp_mat, closed=c(TRUE, TRUE))
      intersection1 <- interval_complement(interval_union(ints_comp_inside))
    } else {
      intersection1 <- Intervals(c(-Inf, Inf))
    }

    ### OUTSIDE INTERVALS
    ints_outside <- Intervals(phi_bounds[phi_bounds[,3]==0,1:2], closed=c(TRUE, TRUE))
    intersection2 <- interval_complement(interval_union(ints_outside))
    if(length(intersection1)==0) {
      res1 <- intersection2
    } else {
      res1 <- suppressWarnings(interval_intersection(intersection1, intersection2))
    }

    #### SOMEHOW COMBINE THIS RES WITH PREVIOUS RES.

    full_interval <- suppressWarnings(interval_intersection(full_interval, res1))

    #### If this inversion is slow could use the rank 1 update!!!
    #### Update these things for next level.
    C_prev <- cbind(C_prev, cy)


    H_prev <- H_prev - Hcy%*%t(Hcy)/sum(Hcy[cy==1])

    leafSize <- sum(cy)
    ### Shoot. If we want to set "nulled" in this way- we really need to be only working with a single brach of
    ### tree. Cannot condition on extra branches.
    nulled <- cy
  }

  return(full_interval)
}



### Get all the ancestor splits of a certain node.
### You need this to be able to call phibounds
### TREE MUST HAVE BEEN BUILT WITH MODEL=TRUE
### So we can access the data.
getAncestors <- function(tree, node)
{
  try1 <- try({
    object <- partykit::as.party(tree)
    rules <- partykit:::.list.rules.party(object)
    relevantRules <- rules[[as.character(node)]]
    relevantRules <- strsplit(relevantRules, '&')[[1]]
    relevantRules <- sapply(relevantRules, trimws)
    relevantRules <- sapply(relevantRules, function(u) paste0("dat$", u))
  })
  if (class(try1)=="try-error") {
    return(errorCondition("Could not Parse Rule Set"))
  } else {
    return(relevantRules)
  }
}


### This is when we want a specific name attached
getAncestors_ALT <- function(tree, node, name)
{
  try1 <- try({
    object <- partykit::as.party(tree)
    rules <- partykit:::.list.rules.party(object)
    relevantRules <- rules[[as.character(node)]]
    relevantRules <- strsplit(relevantRules, '&')[[1]]
    relevantRules <- sapply(relevantRules, trimws)
    relevantRules <- sapply(relevantRules, function(u) paste0(name, "$", u))
  })
  if (class(try1)=="try-error") {
    return(errorCondition("Could not Parse Rule Set"))
  } else {
    return(relevantRules)
  }
}



#' Turn quadratic coefficients into an interval
#'
#' @param vec stores a,b,c coeffs of a quadratic equation
getBounds <- function(vec) {
  a <- vec[1]
  b <- vec[2]
  c <- vec[3]
  tol <- 5*10^(-10)
  if (abs(a) < tol & abs(b) < tol & abs(c) < tol) {
    return(c(-Inf,Inf, 1))
  }

  ## If we just have a straight line
  if (abs(a) < tol) {
    root = -c/b
    if (b>0) {
      return(c(-Inf,root,1))
    } else{
      return(c(root, Inf, 1))
    }
  }

  det <- b^2-4*a*c
  if (det < 0) {
    if (a > 0) {
      #return(simpleError("STOP! EMPTY INTERVAL"))
      return(c(0,0,0)) ## TO AVOID ERRORS FOR NOW PLZ CHANGE LATER
    }
    return(c(-Inf, Inf, 1))
  }
  if (det==0) {
    if (a > 0) {c(-b/(2*a), -b/(2*a), 1)}
    return(c(-Inf, Inf, 1))
  }
  ## Down here determinant is positive: 2 roots
  root1 <- (-b + sqrt(det))/(2*a)
  root2 <- (-b - sqrt(det))/(2*a)
  if (a > 0) {
    return(c(min(root1, root2), max(root1,root2), 1))
  } else {
    return(c(min(root1, root2), max(root2,root1), 0))
  }
}
