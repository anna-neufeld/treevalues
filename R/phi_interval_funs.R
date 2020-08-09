# getInterval <- function(base_tree, nu, splits) {
#   ### rpart orderes things nicely YAY
#   dat <- base_tree$model
#   ### Add an error condition reminding people to call rpart with model=TRUE
#   y <- dat[,1]
#   X <- dat[,-1]
#   n <- nrow(X)
#   p <- NCOL(X)
#
#
#   ### In the first layer, nothing gets zeroed out
#   nulled <- rep(TRUE, n)
#   #H_prev <- diag(1,n,n) - matrix(1/n,n,n)
#   full_interval = Intervals(c(-Inf, Inf))
#   leafSize=n
#
#   ### O(n^2*p) work up front to get the inequalities??
#   ### Because for each of n*p possible splits, it is O(n) work to fill up the vector.
#   ### Consider optimizing this
#   ### If you made the C vecs as you went using an ordered version of the Xs, this might be faster
#   #OH yeah DUH. If you sort the Xs first, the cs have a predictable form
#   # This work is likely avoidable
#   cs<- array(NA, c(p,n,n))
#   for (j in 1:p) {
#     cs[j,,]  <-  sapply(X[,j], function(u) (u <= X[,j]))
#   }
#
#   # O(n)
#   norm_nu <- sum(nu^2)
#
#   # O(n)
#   nut_y <- sum(nu*y)
#
#   #nunu <- nu%*%t(nu)/norm_nu
#   y_proj <-nu*(nut_y/norm_nu)
#
#   #yproj2 <-
#
#
#   #sum(nunu==nunu2)/n
#
#   regions <- rep(1,n)
#
#   L <- length(splits)
#   cur <- 1
#
#   coeffs_full <- matrix(0, nrow=n*p*L, ncol=3)
#
#   ### MULTIPLY EVERYTHING IN LOOP BY O(L)
#   for (l in 1:L) {
#
#     #if (NCOL(C_prev) > 1) {
#     #  regions <- rowSums(C_prev) ## O(n)
#     #} else{regions <- C_prev}
#
#     ### Each of these should be O(n). Techically number of regions can scale with L?? So maybe O(nL)??
#     Hnu <- nu - (sapply(1:max(regions), function(u) mean(nu[regions==u])))[regions]
#     Hy <- y - (sapply(1:max(regions), function(u) mean(y[regions==u])))[regions]
#     Hyproj <- y_proj - (sapply(1:max(regions), function(u) mean(y_proj[regions==u])))[regions]
#     Hminus <- Hy-Hyproj
#
#     ### All of these are O(n) too
#     Hnu2 <- Hnu[nulled]
#     Hy2 <- Hy[nulled]
#     Hminus2 <- Hminus[nulled]
#
#     cy<- eval(parse(text = splits[l]))*nulled
#
#     ### This is an interesting choice.
#     ### If the region could have been formed with only one split instead of two,
#     ### say that no interval is possible??
#     if (sum(cy)==sum(leafSize)) {
#       coeffs_full[cur:(cur+1),] <- c(1,1,1)
#       break
#     }
#
#     Hcy <- cy - (sapply(1:max(regions), function(u) mean(cy[regions==u])))[regions]
#
#     ### I'm doing O(n)
#     denom <- sum(cy)*(leafSize - sum(cy))/leafSize
#
#
#
#     indices <- which(cy==1)
#     partialA <- (sum(Hnu[indices])^2/norm_nu^2)/denom
#     partialB <- as.numeric((2*sum(Hy[indices])*(t(nu)%*%Hcy)/
#                               norm_nu - 2*(t(cy)%*%Hnu)^2%*%nut_y/
#                               norm_nu^2)/denom)
#     partialC <- (sum(Hminus[indices]))^2/denom
#
#
#
#
#     for (j in 1:p) {
#       c_now <- cs[j,nulled, nulled]
#       hardwork <- sort(rowSums(c_now),index.return=TRUE) ## I hope this is O(n)!!!!
#       num1sSORT <- hardwork$x
#       sorted <- hardwork$ix
#       len <- length(sorted)
#
#       denomsSORT <- (num1sSORT*(leafSize-num1sSORT)/leafSize)
#
#
#
#       ### These should be O(n)
#       unus <- cumsum(Hnu2[sorted])[1:(sum(nulled)-1)]
#       uys <- cumsum(Hy2[sorted])[1:(sum(nulled)-1)]
#       umin <- cumsum(Hminus2[sorted])[1:(sum(nulled)-1)]
#       denomsSORT <- denomsSORT[1:(sum(nulled)-1)]
#
#       ### This is the clever thing where we hope we are only paying O(n) because
#       ### we are not multiplying matricies
#       other_As <- unus^2/norm_nu^2
#       other_Bs <- 2*uys*unus/norm_nu- 2*(unus^2*(nut_y))/norm_nu^2
#       other_Cs <- umin^2
#
#
#       avec <- other_As/denomsSORT-partialA
#       bvec <- other_Bs/denomsSORT-partialB
#       cvec <- other_Cs/denomsSORT-partialC
#
#       num_coeffs <- length(avec[ denomsSORT != 0])
#       if (num_coeffs != 0) {
#       coeffs_full[cur:(cur+num_coeffs-1),] <- cbind(avec,bvec,cvec)[ denomsSORT != 0,]
#       }
#       cur <- cur+num_coeffs
#
#     }
#
#     #### If this inversion is slow could use the rank 1 update!!!
#     #### Update these things for next level.
#     regions <- regions+cy
#     #H_prev <- H_prev - Hcy%*%t(Hcy)/sum(Hcy[cy==1])
#     leafSize <- sum(cy)
#     nulled <- cy==1 ## requires that we are testing SIBLING nodes.
#
#   }
#
#   phi_bounds <-  t(apply(coeffs_full, 1, getBounds))
#
#   ### INSIDE INTERVALS
#   numIn <- nrow(phi_bounds[phi_bounds[,3]==1,])
#   if (!is.null(numIn)) {
#     inside_comp_mat <- matrix(c(-Inf, Inf), nrow=2*numIn, ncol=2, byrow=TRUE)
#     inside_comp_mat[1:numIn,2] <- phi_bounds[phi_bounds[,3]==1,1]
#     inside_comp_mat[(numIn+1):(2*numIn), 1] <- phi_bounds[phi_bounds[,3]==1,2]
#     ints_comp_inside <- Intervals(inside_comp_mat, closed=c(TRUE, TRUE))
#     intersection1 <- interval_complement(interval_union(ints_comp_inside))
#   } else {
#     intersection1 <- Intervals(c(-Inf, Inf))
#   }
#
#   ### OUTSIDE INTERVALS
#   ints_outside <- Intervals(phi_bounds[phi_bounds[,3]==0,1:2], closed=c(TRUE, TRUE))
#   intersection2 <- interval_complement(interval_union(ints_outside))
#   #if(length(intersection1)==0) {
#     #res1 <- intersection2
#   #} else {
#   res1 <- suppressWarnings(interval_intersection(intersection1, intersection2))
#   #}
#
#   return(res1)
# }
#
#

#' getPhiInterval
#' Should work for CIs or Hyp Tests.
#'
#' Also make this play nice with other people's real data
#'
#' Add warnings to ensure that their base_tree was build properly!!!
#' Would be a lot nicer if people only ha to pass dat. Or, better yet, only base_tree,
#' because they will use model=TRUE
getInterval_MEDIUM_OLD <- function(base_tree, nu, splits) {

  if (NROW(base_tree$model) > 450) {
    return(getInterval_EXPERIMENT(base_tree, nu, splits))
  } else {
  ### rpart orderes things nicely YAY
  dat <- base_tree$model
  y <- dat[,1]
  X <- dat[,-1]
  n <- nrow(X)
  p <- NCOL(X)

  #### Overall Indicator Vars
  ### this is o(n^2*p) work. Is there an way to do better??
  cs <- matrix(NA, nrow=n*p, ncol=n)
  for (j in 1:p) {
    cs[((j-1)*n+1):(j*n),] <- sapply(X[,j], function(u) (u <= X[,j]))
  }
  ### Remove places that will lead us to divide by 0. They are trivial splits.
  ### O(n^2)
  cs <- cs[rowSums(cs)!=n,]


  ### In the first layer, nothing gets zeroed out
  nulled <- rep(1, n)
  C_prev <- rep(1,n)


  ### I guess O(n^2) just to make the matrix.
  H_prev <- diag(1,n,n) - matrix(1/n,n,n)
  full_interval = Intervals(c(-Inf, Inf))

  ### This stuff is tecnhically O(n^2) too, but could be optimized??
  norm_nu <- sum(nu^2)
  nunu <- nu%*%t(nu)/sum(nu^2)
  y_proj <- nunu%*%y
  nut_y <- t(nu)%*%y



  for (i in 1:length(splits)) {
    leafSize <- sum(nulled)

    if (NCOL(C_prev) > 1) {
      regions <- rowSums(C_prev)
    } else{regions <- C_prev}

    Hnu <- nu - (sapply(1:max(regions), function(u) mean(nu[regions==u])))[regions]
    Hy <- y - (sapply(1:max(regions), function(u) mean(y[regions==u])))[regions]
    Hyproj <- y_proj - (sapply(1:max(regions), function(u) mean(y_proj[regions==u])))[regions]
    Hminus <- Hy-Hyproj

    cy<- eval(parse(text = splits[i]))*nulled

    Hcy <- cy - (sapply(1:max(regions), function(u) mean(cy[regions==u])))[regions]


    ### THIS IS STILL SLOW. What a dumb thing to spend time on, right??
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

    other_Bs <- apply(c_now, 1, function(u) 2*sum(Hy[u==1])*sum(Hnu[u==1])/
                       norm_nu
         - 2*(sum(Hnu[u==1])^2%*%nut_y)/
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

    ### Shoot. If we want to set "nulled" in this way- we really need to be only working with a single brach of
    ### tree. Cannot condition on extra branches.
    nulled <- cy
  }

  return(full_interval)
  }
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

### Get all the ancestor splits of a certain node.
### You need this to be able to call phibounds
### TREE MUST HAVE BEEN BUILT WITH MODEL=TRUE
### So we can access the data.
getAncestors_ROUND <- function(tree, node, altsplits)
{
  tree$splits[,4] <- altsplits
  try1 <- try({
    object <- partykit::as.party(tree)
    rules <- partykit:::.list.rules.party(object)
    relevantRules <- rules[[as.character(node)]]
    relevantRules <- strsplit(relevantRules, '&')[[1]]
    relevantRules <- sapply(relevantRules, trimws)
    relevantRules <- sapply(relevantRules, function(u) substr(u,1,11))
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

  #vec[abs(vec) < 1e-15] <- 0
  a <- vec[1]
  b <- vec[2]
  c <- vec[3]
  tol <- 1*10^(-10)
  if (abs(a) < tol & abs(b) < tol & abs(c) < tol) {
    return(c(-Inf,Inf, 1))
  }

  ## If we just have a straight line
  if (abs(a) < tol) {
    if (abs(b) < tol) {
      if (c < 0) {return(c(-Inf,Inf, 1))
        } else {return(c(-Inf, Inf, 0))}
    }
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
      return(c(-Inf,Inf,0)) ## TO AVOID ERRORS FOR NOW PLZ CHANGE LATER
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


#### A certain phi gives u the same region.
#### But if u had started with that phi, the first split to win would not have matched
### the first split that we wanted to win. Need to draw a pic of this.
