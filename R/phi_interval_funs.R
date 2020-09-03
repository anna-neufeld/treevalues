#' Workhorse function for computing conditioning sets.
#'
#' Returns the interval for \phi such that Tree(y'(phi)) (where the definition of y'(phi) depends on nu)
#' contains the series of splits "splits".Either base_tree with model=TRUE or X and y must be provided, not both.
#'
#' @param base_tree An rpart object. Must have been built with model=TRUE arguement unless X and y are provided separately
#' @param nu The contrast vector that defines the parameter.
#' @param splits A vector of strings specifying the splits.
#' @param X the covariate matrix for the tree.
#' @param y the response vector for the tree.
#'
#' @return an object of class Interval that defines the set S.
#' @export
getInterval_PROB <- function(base_tree, nu, splits) {
  minbucket <- base_tree$control$minbucket

  dat <- base_tree$model
  ### Add an error condition reminding people to call rpart with model=TRUE
  y <- dat[,1]
  X <- dat[,-1]
  n <- nrow(X)
  p <- NCOL(X)


  ### In the first layer, nothing gets zeroed out
  nulled <- rep(TRUE, n)
  #H_prev <- diag(1,n,n) - matrix(1/n,n,n)
  full_interval = Intervals(c(-Inf, Inf))
  leafSize=n

  ### O(n^2*p) work up front to get the inequalities??
  ### Because for each of n*p possible splits, it is O(n) work to fill up the vector.
  ### Consider optimizing this
  ### If you made the C vecs as you went using an ordered version of the Xs, this might be faster
  #OH yeah DUH. If you sort the Xs first, the cs have a predictable form
  # This work is likely avoidable
  cs<- array(NA, c(p,n,n))
  for (j in 1:p) {
    cs[j,,]  <-  sapply(X[,j], function(u) (u <= X[,j]))
  }

  # O(n)
  norm_nu <- sum(nu^2)

  # O(n)
  nut_y <- sum(nu*y)

  #nunu <- nu%*%t(nu)/norm_nu
  y_proj <-nu*(nut_y/norm_nu)



  regions <- rep(1,n)

  L <- length(splits)
  cur <- 1

  coeffs_full <- matrix(0, nrow=n*p*L, ncol=3)

  ### MULTIPLY EVERYTHING IN LOOP BY O(L)
  for (l in 1:L) {

    #if (NCOL(C_prev) > 1) {
    #  regions <- rowSums(C_prev) ## O(n)
    #} else{regions <- C_prev}

    ### Each of these should be O(n). Techically number of regions can scale with L?? So maybe O(nL)??
    Hnu <- nu - (sapply(1:max(regions), function(u) mean(nu[regions==u])))[regions]
    Hy <- y - (sapply(1:max(regions), function(u) mean(y[regions==u])))[regions]
    Hyproj <- y_proj - (sapply(1:max(regions), function(u) mean(y_proj[regions==u])))[regions]
    Hminus <- Hy-Hyproj

    ### All of these are O(n) too
    Hnu2 <- Hnu[nulled]
    Hy2 <- Hy[nulled]
    Hminus2 <- Hminus[nulled]

    cy<- eval(parse(text = splits[l]))*nulled

    ### This is an interesting choice.
    ### If the region could have been formed with only one split instead of two,
    ### say that no interval is possible??
    if (sum(cy)==sum(leafSize)) {
      coeffs_full[cur:(cur+1),] <- c(1,1,1)
      break
    }

    Hcy <- cy - (sapply(1:max(regions), function(u) mean(cy[regions==u])))[regions]

    ### I'm doing O(n)
    denom <- sum(cy)*(leafSize - sum(cy))/leafSize



    indices <- which(cy==1)
    partialA <- (sum(Hnu[indices])^2/norm_nu^2)/denom
    partialB <- as.numeric((2*sum(Hy[indices])*(t(nu)%*%Hcy)/
                              norm_nu - 2*(t(cy)%*%Hnu)^2%*%nut_y/
                              norm_nu^2)/denom)
    partialC <- (sum(Hminus[indices]))^2/denom




    for (j in 1:p) {
      c_now <- cs[j,nulled, nulled]
      hardwork <- sort(rowSums(c_now),index.return=TRUE) ## I hope this is O(n)!!!!
      num1sSORT <- hardwork$x
      sorted <- hardwork$ix
      len <- length(sorted)

      denomsSORT <- (num1sSORT*(leafSize-num1sSORT)/leafSize)
      minDenomSORT <- ((minbucket-1)*(leafSize-(minbucket-1))/leafSize)

      ### These should be O(n)
      unus <- cumsum(Hnu2[sorted])[1:(sum(nulled)-1)]
      uys <- cumsum(Hy2[sorted])[1:(sum(nulled)-1)]
      umin <- cumsum(Hminus2[sorted])[1:(sum(nulled)-1)]
      denomsSORT <- denomsSORT[1:(sum(nulled)-1)]

      ### This is the clever thing where we hope we are only paying O(n) because
      ### we are not multiplying matricies
      other_As <- unus^2/norm_nu^2
      other_Bs <- 2*uys*unus/norm_nu- 2*(unus^2*(nut_y))/norm_nu^2
      other_Cs <- umin^2


      avec <- other_As/denomsSORT-partialA
      bvec <- other_Bs/denomsSORT-partialB
      cvec <- other_Cs/denomsSORT-partialC


      ### SHOOT THIS IS GOOD FOR CATEGORICAL. BUT SOMEHOW NOT FOR REGULAR???
      #### PUT BACK SOON
      uniqueindices <- c(cumsum(table(hardwork$x)))[1:(length(unique(hardwork$x))-1)]
      uniqueindices <- uniqueindices[denomsSORT[uniqueindices] > minDenomSORT]

      ###### NOTE TO SELD: PUT BACK IN!!!
      #uniqueindices <- uniqueindices[num1sSORT[uniqueindices] >= minbucket &
      #                                (leafSize - num1sSORT[uniqueindices]) >= minbucket]

      #uniqueindices <- (1:length(denomsSORT))[denomsSORT !=0]

      #print(uniqueindices)
      #print(which(denomsSORT==0))
      #print("---------")

      num_coeffs <- length(avec[uniqueindices])
      #print(num_coeffs)
      if (num_coeffs !=0 & !is.na(uniqueindices[1])) {
        coeffs_full[cur:(cur+num_coeffs-1),] <- cbind(avec,bvec,cvec)[uniqueindices,]
      }
      cur <- cur+num_coeffs
      if (cur > 2835) {
        break
      }
    }

    #### If this inversion is slow could use the rank 1 update!!!
    #### Update these things for next level.
    regions <- regions+cy
    #H_prev <- H_prev - Hcy%*%t(Hcy)/sum(Hcy[cy==1])
    leafSize <- sum(cy)
    nulled <- cy==1 ## requires that we are testing SIBLING nodes.

  }

  coeffs_full <- coeffs_full[1:cur,]
  phi_bounds <-  t(apply(coeffs_full, 1, getBounds))

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
  #if(length(intersection1)==0) {
  #res1 <- intersection2
  #} else {
  res1 <- suppressWarnings(interval_intersection(intersection1, intersection2))
  #}

  return(res1)
}

#' Obtain the list of splits that refine a certain region.
#' Necessary precursor to calling "getInterval"
#'
#' @param tree An rpart object. Must have been built with model=TRUE
#' @param node A number identifying the region that you want the ancestors of. Numbers correspond to elements in tree$where.
#'
#' @return A vector of strings describing the splits that define the node in the tree.
#' @export
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
    return(errorCondition("Could not obtain splits from tree. Make sure that the rpart tree
                          was built with model=TRUE parameter."))
  } else {
    return(relevantRules)
  }
}


#' An internal function that was just just used for some simulations.
#'
#' It is the same as getAncestors but you can paste a different name onto the inequalities
#' for parsing purposes.
#' @keywords internal
#' @noRd
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
#' @return A vector describing the regions where the quadratic is negative.
#' @keywords internal
#' @noRd
getBounds <- function(vec) {

  #vec[abs(vec) < 1e-15] <- 0
  a <- vec[1]
  b <- vec[2]
  c <- vec[3]
  tol <- 10*10^(-10)
  if (abs(a) < tol & abs(b) < tol & abs(c) < tol) {
    return(c(-Inf,Inf, 1))
  }

  ## If we just have a straight line
  if (abs(a) < tol) {
    if (abs(b) < tol) {
      if (c < 1e-6) {return(c(-Inf,Inf, 1))
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


#' IF ANYTHING SEEMS BROKEN, TURN THIS BACK ON!!!!!!!!!!!!!!!!!!!
getInterval <- function(base_tree = NULL, nu, splits, X=NULL, y=NULL,
                                  minbucket = base_tree$control$minbucket) {

  #if (is.null(minbucket)) {
  #  if (is.null(base_tree$call$control$minsplit)) {
  #  minbucket = round(20/3)
  #  } else {
  #    minbucket = round(base_tree$call$control$minsplit/3)
  #  }
  #}
  ### rpart orderes things nicely YAY
  if (!is.null(base_tree$model)) {
    dat <- base_tree$model
    ### Add an error condition reminding people to call rpart with model=TRUE
    y <- dat[,1]
    X <- dat[,-1]
  } else {
    if (is.null(X) | is.null(y)) {
      stop('Must provide data X and y, or must ensure that base_tree is build with model=TRUE')
    }
  }

  n <- nrow(X)
  p <- NCOL(X)

  ### In the first layer, nothing gets zeroed out
  nulled <- rep(TRUE, n)
  full_interval = Intervals(c(-Inf, Inf))
  leafSize=n

  ### O(n^2*p) work up front to get the inequalities??
  ### Because for each of n*p possible splits, it is O(n) work to fill up the vector.
  ### If you sort the Xs first, the cs have a predictable form
  ### Can likely avoid/ optimize this step.
  cs<- array(NA, c(p,n,n))
  for (j in 1:p) {
    cs[j,,]  <-  sapply(X[,j], function(u) (u <= X[,j]))
  }

  # O(n)
  norm_nu <- sum(nu^2)

  # O(n)
  nut_y <- sum(nu*y)

  y_proj <-nu*(nut_y/norm_nu)

  regions <- rep(1,n)

  L <- length(splits)
  cur <- 1

  coeffs_full <- matrix(0, nrow=n*p*L, ncol=3)

  ### MULTIPLY EVERYTHING IN LOOP BY O(L)
  for (l in 1:L) {

    ### Each of these should be O(n). Techically number of regions can scale with L?? So maybe O(nL)??
    Hnu <- nu - (sapply(1:max(regions), function(u) mean(nu[regions==u])))[regions]
    Hy <- y - (sapply(1:max(regions), function(u) mean(y[regions==u])))[regions]
    Hyproj <- y_proj - (sapply(1:max(regions), function(u) mean(y_proj[regions==u])))[regions]
    Hminus <- Hy-Hyproj

    ### All of these are O(n) too
    Hnu2 <- Hnu[nulled]
    Hy2 <- Hy[nulled]
    Hminus2 <- Hminus[nulled]

    cy<- eval(parse(text = splits[l]))*nulled

    ### This is an interesting choice.
    ### If the region could have been formed with only one split instead of two,
    ### say that no interval is possible??
    if (sum(cy)==sum(leafSize)) {
      coeffs_full[cur:(cur+1),] <- c(1,1,1)
      break
    }

    Hcy <- cy - (sapply(1:max(regions), function(u) mean(cy[regions==u])))[regions]

    ### I'm doing O(n)
    denom <- sum(cy)*(leafSize - sum(cy))/leafSize



    indices <- which(cy==1)
    partialA <- (sum(Hnu[indices])^2/norm_nu^2)/denom
    partialB <- as.numeric((2*sum(Hy[indices])*(t(nu)%*%Hcy)/
                              norm_nu - 2*(t(cy)%*%Hnu)^2%*%nut_y/
                              norm_nu^2)/denom)
    partialC <- (sum(Hminus[indices]))^2/denom




    for (j in 1:p) {
      c_now <- cs[j,nulled, nulled]
      hardwork <- sort(rowSums(c_now),index.return=TRUE) ## I hope this is O(n)!!!!
      num1sSORT <- hardwork$x
      sorted <- hardwork$ix
      len <- length(sorted)

      denomsSORT <- (num1sSORT*(leafSize-num1sSORT)/leafSize)

      minDenomSORT <- ((minbucket-1)*(leafSize-(minbucket-1))/leafSize)



      ### These should be O(n)
      unus <- cumsum(Hnu2[sorted])[1:(sum(nulled)-1)]
      uys <- cumsum(Hy2[sorted])[1:(sum(nulled)-1)]
      umin <- cumsum(Hminus2[sorted])[1:(sum(nulled)-1)]
      denomsSORT <- denomsSORT[1:(sum(nulled)-1)]

      ### This is the clever thing where we hope we are only paying O(n) because
      ### we are not multiplying matricies
      other_As <- unus^2/norm_nu^2
      other_Bs <- 2*uys*unus/norm_nu- 2*(unus^2*(nut_y))/norm_nu^2
      other_Cs <- umin^2


      avec <- other_As/denomsSORT-partialA
      bvec <- other_Bs/denomsSORT-partialB
      cvec <- other_Cs/denomsSORT-partialC

      indices <- which(denomsSORT > minDenomSORT)
      num_coeffs <- length(indices)


      if (num_coeffs != 0) {
        coeffs_full[cur:(cur+num_coeffs-1),] <- cbind(avec,bvec,cvec)[ indices,]
      }
      cur <- cur+num_coeffs

    }

    #### Update these things for next level.
    regions <- regions+cy
    leafSize <- sum(cy)
    nulled <- cy==1
  }

  phi_bounds <-  t(apply(coeffs_full, 1, getBounds))

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
  #if(length(intersection1)==0) {
  #res1 <- intersection2
  #} else {
  res1 <- suppressWarnings(interval_intersection(intersection1, intersection2))
  #}

  return(res1)
}
