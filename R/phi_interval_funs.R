#' Workhorse function for computing conditioning sets.
#'
#' Returns the interval for phi such that Tree(y'(phi)) (where the definition of y'(phi) depends on nu)
#' contains the series of branch "branch". An rpart object build with model=TRUE must be provided.
#'
#' @param tree An rpart object. Must have been built with model=TRUE arguement unless X and y are provided separately
#' @param nu The contrast vector that defines the parameter.
#' @param branch A vector of strings specifying the branch.
#' @param sib If you are doing inference and nu=nu_sib, this takes advantage of some computational speedups!
#' @param grow Set this to true if you only want to compute S_grow
#' @param prune Set this to true if you only want to see what gets REMOVED from S_grow during pruning
#' @return an object of class Interval that defines the set S.
#' @importFrom stats predict
#' @export
getInterval_full <- function(tree, nu, branch,sib=FALSE,grow=FALSE,prune=FALSE) {

  branch <- paste("dat$", branch)
  minbucket <- tree$control$minbucket
  cp <- tree$control$cp
  if (is.null(tree$model)) {
    stop("rpart tree must be built with model=TRUE")
  }
  dat <- tree$model

  y <- dat[,1]
  X <- dat[,-1]
  n <- nrow(X)
  p <- NCOL(X)


  # Go through and compute the h() and g() function for the subtrees that don't depend on phi!!!
  # This work isonly needed for s_prune.
  if (length(branch) > 1) {

    vecChildren <- eval(parse(text = paste(branch, collapse=" & ")))
    childrenGains <-   sum((y[as.logical(vecChildren)]-mean(y[as.logical(vecChildren)]))^2)  - sum((y[as.logical(vecChildren)]-predict(tree)[as.logical(vecChildren)])^2)
    childrenSize <- length(unique(tree$where[as.logical(vecChildren)]))-1

    inequalities <- sapply(branch, function(u) eval(parse(text = u)))
    counts <- apply(inequalities, 1, cumsum)
    membership <- t(sapply(1:NROW(counts), function(u) counts[u,]==u))

    off_subgroup_membership <- matrix(0, nrow=NROW(membership), ncol=NCOL(membership))
    off_subgroup_membership[1,] <- !membership[1,]
    for (c in 2:NROW(off_subgroup_membership)) {
      off_subgroup_membership[c,] <- membership[c,]==0 & membership[(c-1),]==1
    }

    SSE_full_tree <- sum((y-predict(tree))^2)

    subgroup_GAINZ <- function(vec) {
      sum((y[as.logical(vec)]-mean(y[as.logical(vec)]))^2)  - sum((y[as.logical(vec)]-predict(tree)[as.logical(vec)])^2)
    }

    numbranch <- function(vec) {
      length(unique(tree$where[as.logical(vec)])) - 1
    }

    gainZ <- apply(off_subgroup_membership, 1, subgroup_GAINZ)
    gainZfull <- cumsum(gainZ[length(gainZ):1])

    ### There is a chance that these are all 0 always because of the way I generated
    ### my branch ... more later.
    sizeZ <- apply(off_subgroup_membership, 1, numbranch)
    sizeZfull <- cumsum(sizeZ[length(sizeZ):1])
  } else {

    off_subgroup_membership <- !eval(parse(text = branch))
    vec <- off_subgroup_membership
    gainZfull <-   sum((y[as.logical(vec)]-mean(y[as.logical(vec)]))^2)  - sum((y[as.logical(vec)]-predict(tree)[as.logical(vec)])^2)
    sizeZfull <- length(unique(tree$where[as.logical(vec)])) - 1

    vecChildren <- eval(parse(text = branch))
    childrenGains <-   sum((y[as.logical(vecChildren)]-mean(y[as.logical(vecChildren)]))^2)  - sum((y[as.logical(vecChildren)]-predict(tree)[as.logical(vecChildren)])^2)
    childrenSize <- length(unique(tree$where[as.logical(vecChildren)])) - 1 ### DO I MINUS 1???? no, right?? NOT AGAIN???

  }

  cp_cutoff <- cp*(sum((y-mean(y))^2))

  ### In the first layer, nothing gets zeroed out
  nulled <- rep(TRUE, n)
  #H_prev <- diag(1,n,n) - matrix(1/n,n,n)
  full_interval = Intervals(c(-Inf, Inf))
  leafSize=n

  ### This work is actually the slow part!! Essentially sorting the obsertions for each covariate.
  cs<- array(NA, c(p,n,n))
  for (j in 1:p) {
    cs[j,,]  <-  sapply(X[,j], function(u) (u <= X[,j]))
  }

  # O(n)
  norm_nu <- sum(nu^2)
  nut_y <- sum(nu*y)
  y_proj <-nu*(nut_y/norm_nu)


  regions <- rep(1,n)

  L <- length(branch)
  cur <- 1

  coeffs_full <- matrix(0, nrow=n*p*L, ncol=3)

  gainsA <- rep(0,L)
  gainsB <- rep(0,L)
  gainsC <- rep(0,L)

  ### Someday: should I redo with Daniela's notation?? To match the proofs??
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

    cy<- eval(parse(text = branch[l]))*nulled

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

    gainsA[L-l+1] <- partialA
    gainsB[L-l+1] <- partialB
    gainsC[L-l+1] <- partialC

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



      uniqueindices <- c(cumsum(table(hardwork$x)))[1:(length(unique(hardwork$x))-1)]
      uniqueindices <- uniqueindices[denomsSORT[uniqueindices] > minDenomSORT]

      num_coeffs <- length(avec[uniqueindices])
      if (num_coeffs !=0 & !is.na(uniqueindices[1])) {
        coeffs_full[cur:(cur+num_coeffs-1),] <- cbind(avec,bvec,cvec)[uniqueindices,]
      }
      cur <- cur+num_coeffs
    }


    regions <- regions+cy
    leafSize <- sum(cy)
    nulled <- cy==1

  }

  coeffs_full <- coeffs_full[1:cur,]
  phi_bounds <-  t(apply(coeffs_full, 1, getBounds))


  ### INSIDE INTERVALS
  if (sib) {
    #### FIX
    inside <- phi_bounds[phi_bounds[,3]==1,]
    inside <- inside[inside[,1] != -Inf,]
    if (length(inside)>=1) {
      intersection1 <- Intervals(c(max(inside[,1]), min(inside[,2])))
    } else {
      intersection1 <- Intervals(c(-Inf, Inf))
    }
    ### OUTSIDE INTERVALS
    outside <- phi_bounds[phi_bounds[,3]==0,]
    intersection2 <- interval_complement(Intervals(c(min(outside[,1]), max(outside[,2]))))
  } else {
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
  }


  res1 <- suppressWarnings(interval_intersection(intersection1, intersection2))

  #### Can ignore all the rest of the steps in this case!!!!!
  if(grow) {return(res1)}

  if (length(res1)==0) {
    return(res1)
  } else{


    sizes <- (1:L +sizeZfull + childrenSize) + 1 - 1

    cpMat <- cbind(cumsum(gainsA)/(sizes), cumsum(gainsB)/(sizes), cumsum(gainsC)/(sizes) - cp_cutoff + gainZfull/sizes + childrenGains/sizes)
    cpMat <- cpMat[sizes != 0, ]

    if (NCOL(cpMat) == 1) {
      cpBounds <- getBounds(cpMat)
      if (cpBounds[3]==0) {return(res1)}
    } else {
      cpBounds <- t(apply(cpMat, 1, getBounds))
      cpBounds <- cpBounds[cpBounds[,3]==1,]
    }

    cpInterval <- Intervals(c(-Inf, Inf))
    if (cp_cutoff==0) {
      cpInterval = Intervals(c(-Inf,Inf))
    } else {
      if (NCOL(cpBounds) == 1) {
        cpInterval <-  suppressWarnings(interval_intersection(cpInterval, interval_complement(Intervals(cpBounds[1:2]))))
      } else {
        if (NROW(cpBounds) != 0) {
          for (i in 1:NROW(cpBounds)) {
            cpInterval <- suppressWarnings(interval_intersection(cpInterval, interval_complement(Intervals(cpBounds[i,1:2]))))
          }
        }
      }
    }

    #print(res1)
    #print(cpInterval)

    return(suppressWarnings(interval_intersection(res1, cpInterval)))
  }
}




#' #' Workhorse function for computing conditioning sets; prune only!!
#' #'
#' #' Returns the interval for phi such that Tree(y'(phi)) (where the definition of y'(phi) depends on nu)
#' #' contains the series of branch "branch".Either tree with model=TRUE or X and y must be provided, not both.
#' #'
#' #' @param tree An rpart object. Must have been built with model=TRUE arguement unless X and y are provided separately
#' #' @param nu The contrast vector that defines the parameter.
#' #' @param branch A vector of strings specifying the branch.
#' #'
#' #' @return an object of class Interval that defines the set S.
#' #' @export
#' #' @importFrom stats predict
#' getInterval_prune <- function(tree, nu, branch,sibs=FALSE) {
#'
#'   branch <- paste("dat$", branch)
#'   minbucket <- tree$control$minbucket
#'   cp <- tree$control$cp
#'   if (is.null(tree$model)) {
#'     stop("rpart tree must be built with model=TRUE")
#'   }
#'   dat <- tree$model
#'   ### Add an error condition reminding people to call rpart with model=TRUE
#'   y <- dat[,1]
#'   X <- dat[,-1]
#'   n <- nrow(X)
#'   p <- NCOL(X)
#'
#'
#'   ###### improvement factors for the OTHER parts of the subtree created by each split.
#'   if (length(branch) > 1) {
#'
#'     vecChildren <- eval(parse(text = paste(branch, collapse=" & ")))
#'     childrenGains <-   sum((y[as.logical(vecChildren)]-mean(y[as.logical(vecChildren)]))^2)  - sum((y[as.logical(vecChildren)]-predict(tree)[as.logical(vecChildren)])^2)
#'     childrenSize <- length(unique(tree$where[as.logical(vecChildren)]))-1  #DO WE MINUS 1???
#'
#'     inequalities <- sapply(branch, function(u) eval(parse(text = u)))
#'     counts <- apply(inequalities, 1, cumsum)
#'     membership <- t(sapply(1:NROW(counts), function(u) counts[u,]==u))
#'
#'     off_subgroup_membership <- matrix(0, nrow=NROW(membership), ncol=NCOL(membership))
#'     off_subgroup_membership[1,] <- !membership[1,]
#'     for (c in 2:NROW(off_subgroup_membership)) {
#'       off_subgroup_membership[c,] <- membership[c,]==0 & membership[(c-1),]==1
#'     }
#'
#'     SSE_full_tree <- sum((y-predict(tree))^2)
#'
#'     subgroup_GAINZ <- function(vec) {
#'       sum((y[as.logical(vec)]-mean(y[as.logical(vec)]))^2)  - sum((y[as.logical(vec)]-predict(tree)[as.logical(vec)])^2)
#'     }
#'
#'     numbranch <- function(vec) {
#'       length(unique(tree$where[as.logical(vec)])) - 1
#'     }
#'
#'     gainZ <- apply(off_subgroup_membership, 1, subgroup_GAINZ)
#'     gainZfull <- cumsum(gainZ[length(gainZ):1])
#'
#'     ### There is a chance that these are all 0 always because of the way I generated
#'     ### my branch ... more later.
#'     sizeZ <- apply(off_subgroup_membership, 1, numbranch)
#'     sizeZfull <- cumsum(sizeZ[length(sizeZ):1])
#'   } else {
#'
#'     ### Wait if there is
#'     off_subgroup_membership <- !eval(parse(text = branch))
#'     vec <- off_subgroup_membership
#'     gainZfull <-   sum((y[as.logical(vec)]-mean(y[as.logical(vec)]))^2)  - sum((y[as.logical(vec)]-predict(tree)[as.logical(vec)])^2)
#'     sizeZfull <- length(unique(tree$where[as.logical(vec)])) - 1
#'
#'     vecChildren <- eval(parse(text = branch))
#'     childrenGains <-   sum((y[as.logical(vecChildren)]-mean(y[as.logical(vecChildren)]))^2)  - sum((y[as.logical(vecChildren)]-predict(tree)[as.logical(vecChildren)])^2)
#'     childrenSize <- length(unique(tree$where[as.logical(vecChildren)])) - 1 ### DO I MINUS 1???? no, right?? NOT AGAIN???
#'
#'   }
#'
#'   cp_cutoff <- cp*(sum((y-mean(y))^2))
#'
#'   ### In the first layer, nothing gets zeroed out
#'   nulled <- rep(TRUE, n)
#'   H_prev <- diag(1,n,n) - matrix(1/n,n,n)
#'   full_interval = Intervals(c(-Inf, Inf))
#'   leafSize=n
#'
#'   ### O(n^2*p) work up front to get the inequalities??
#'   ### Because for each of n*p possible branch, it is O(n) work to fill up the vector.
#'   ### Consider optimizing this
#'   ### If you made the C vecs as you went using an ordered version of the Xs, this might be faster
#'   #OH yeah DUH. If you sort the Xs first, the cs have a predictable form
#'   # This work is likely avoidable
#'
#'
#'   #cs<- array(NA, c(p,n,n))
#'   #for (j in 1:p) {
#'   #  cs[j,,]  <-  sapply(X[,j], function(u) (u <= X[,j]))
#'   #}
#'
#'   # O(n)
#'   norm_nu <- sum(nu^2)
#'
#'   # O(n)
#'   nut_y <- sum(nu*y)
#'
#'   #nunu <- nu%*%t(nu)/norm_nu
#'   y_proj <-nu*(nut_y/norm_nu)
#'
#'
#'
#'   regions <- rep(1,n)
#'
#'   L <- length(branch)
#'   cur <- 1
#'
#'   #coeffs_full <- matrix(0, nrow=n*p*L, ncol=3)
#'
#'   gainsA <- rep(0,L)
#'   gainsB <- rep(0,L)
#'   gainsC <- rep(0,L)
#'
#'   ### MULTIPLY EVERYTHING IN LOOP BY O(L)
#'   for (l in 1:L) {
#'
#'     #if (NCOL(C_prev) > 1) {
#'     #  regions <- rowSums(C_prev) ## O(n)
#'     # } else{regions <- C_prev}
#'
#'     ### Each of these should be O(n). Techically number of regions can scale with L?? So maybe O(nL)??
#'     Hnu <- nu - (sapply(1:max(regions), function(u) mean(nu[regions==u])))[regions]
#'     Hy <- y - (sapply(1:max(regions), function(u) mean(y[regions==u])))[regions]
#'     Hyproj <- y_proj - (sapply(1:max(regions), function(u) mean(y_proj[regions==u])))[regions]
#'     Hminus <- Hy-Hyproj
#'
#'     ### All of these are O(n) too
#'     Hnu2 <- Hnu[nulled]
#'     Hy2 <- Hy[nulled]
#'     Hminus2 <- Hminus[nulled]
#'
#'     cy<- eval(parse(text = branch[l]))*nulled
#'
#'     ### This is an interesting choice.
#'     ### If the region could have been formed with only one split instead of two,
#'     ### say that no interval is possible??
#'     #if (sum(cy)==sum(leafSize)) {
#'     # coeffs_full[cur:(cur+1),] <- c(1,1,1)
#'     #  break
#'     #}
#'
#'     Hcy <- cy - (sapply(1:max(regions), function(u) mean(cy[regions==u])))[regions]
#'
#'     ### I'm doing O(n)
#'     denom <- sum(cy)*(leafSize - sum(cy))/leafSize
#'
#'
#'
#'     indices <- which(cy==1)
#'     partialA <- (sum(Hnu[indices])^2/norm_nu^2)/denom
#'     partialB <- as.numeric((2*sum(Hy[indices])*(t(nu)%*%Hcy)/
#'                               norm_nu - 2*(t(cy)%*%Hnu)^2%*%nut_y/
#'                               norm_nu^2)/denom)
#'     partialC <- (sum(Hminus[indices]))^2/denom
#'     #print("---------")
#'     #print(partialA)
#'     #print(partialB)
#'     #print(partialC)
#'
#'
#'     gainsA[L-l+1] <- partialA
#'     gainsB[L-l+1] <- partialB
#'     gainsC[L-l+1] <- partialC
#'
#'     #for (j in 1:p) {
#'     # c_now <- cs[j,nulled, nulled]
#'     #hardwork <- sort(rowSums(c_now),index.return=TRUE) ## I hope this is O(n)!!!!
#'     #num1sSORT <- hardwork$x
#'     #sorted <- hardwork$ix
#'     #len <- length(sorted)
#'
#'     #denomsSORT <- (num1sSORT*(leafSize-num1sSORT)/leafSize)
#'     #minDenomSORT <- ((minbucket-1)*(leafSize-(minbucket-1))/leafSize)
#'
#'     ### These should be O(n)
#'     #unus <- cumsum(Hnu2[sorted])[1:(sum(nulled)-1)]
#'     #uys <- cumsum(Hy2[sorted])[1:(sum(nulled)-1)]
#'     #umin <- cumsum(Hminus2[sorted])[1:(sum(nulled)-1)]
#'     #denomsSORT <- denomsSORT[1:(sum(nulled)-1)]
#'
#'     ### This is the clever thing where we hope we are only paying O(n) because
#'     ### we are not multiplying matricies
#'     #other_As <- unus^2/norm_nu^2
#'     #other_Bs <- 2*uys*unus/norm_nu- 2*(unus^2*(nut_y))/norm_nu^2
#'     #other_Cs <- umin^2
#'
#'
#'
#'
#'     #avec <- other_As/denomsSORT-partialA
#'     # bvec <- other_Bs/denomsSORT-partialB
#'     #cvec <- other_Cs/denomsSORT-partialC
#'
#'
#'     ### SHOOT THIS IS GOOD FOR CATEGORICAL. BUT SOMEHOW NOT FOR REGULAR???
#'     #### PUT BACK SOON
#'     #uniqueindices <- c(cumsum(table(hardwork$x)))[1:(length(unique(hardwork$x))-1)]
#'     #uniqueindices <- uniqueindices[denomsSORT[uniqueindices] > minDenomSORT]
#'
#'     ###### NOTE TO SELD: PUT BACK IN!!!
#'     #uniqueindices <- uniqueindices[num1sSORT[uniqueindices] >= minbucket &
#'     #                                (leafSize - num1sSORT[uniqueindices]) >= minbucket]
#'
#'     #uniqueindices <- (1:length(denomsSORT))[denomsSORT !=0]
#'
#'     #print(uniqueindices)
#'     #print(which(denomsSORT==0))
#'     #print("---------")
#'
#'     #num_coeffs <- length(avec[uniqueindices])
#'     #if (num_coeffs !=0 & !is.na(uniqueindices[1])) {
#'     #  coeffs_full[cur:(cur+num_coeffs-1),] <- cbind(avec,bvec,cvec)[uniqueindices,]
#'     #}
#'     #cur <- cur+num_coeffs
#'     #}
#'
#'     #### If this inversion is slow could use the rank 1 update!!!
#'     #### Update these things for next level.
#'     regions <- regions+cy
#'     #H_prev <- H_prev - Hcy%*%t(Hcy)/sum(Hcy[cy==1])
#'     leafSize <- sum(cy)
#'     nulled <- cy==1 ## requires that we are testing SIBLING nodes.
#'
#'   }

  #coeffs_full <- coeffs_full[1:cur,]
  #phi_bounds <-  t(apply(coeffs_full, 1, getBounds))





  ### INSIDE INTERVALS
  #numIn <- nrow(phi_bounds[phi_bounds[,3]==1,])
  #if (!is.null(numIn)) {
  #  inside_comp_mat <- matrix(c(-Inf, Inf), nrow=2*numIn, ncol=2, byrow=TRUE)
  #  inside_comp_mat[1:numIn,2] <- phi_bounds[phi_bounds[,3]==1,1]
  #  inside_comp_mat[(numIn+1):(2*numIn), 1] <- phi_bounds[phi_bounds[,3]==1,2]
  #  ints_comp_inside <- Intervals(inside_comp_mat, closed=c(TRUE, TRUE))
  #  intersection1 <- interval_complement(interval_union(ints_comp_inside))
  #} else {
  #  intersection1 <- Intervals(c(-Inf, Inf))
  #}

  ### OUTSIDE INTERVALS
  #ints_outside <- Intervals(phi_bounds[phi_bounds[,3]==0,1:2], closed=c(TRUE, TRUE))
  #intersection2 <- interval_complement(interval_union(ints_outside))
  #if(length(intersection1)==0) {
  #res1 <- intersection2
  #} else {
  #res1 <- suppressWarnings(interval_intersection(intersection1, intersection2))

  #if (length(res1)==0) {
  # return(res1)
  #} else{
  #}

  ### Just TRY the minus 1 from the stack overflow post
#'   sizes <- (1:L +sizeZfull + childrenSize) + 1 - 1
#'
#'
#'   #### NEW AND EXPERIMENTAL
#'   ### CHILD_GAINS gets added to EVERY SPLIT. NO one can separate from CHILDGAINS
#'   cpMat <- cbind(cumsum(gainsA)/(sizes), cumsum(gainsB)/(sizes), cumsum(gainsC)/(sizes) - cp_cutoff + gainZfull/sizes + childrenGains/sizes)
#'   cpMat <- cpMat[sizes != 0, ]
#'
#'   if (NCOL(cpMat) == 1) {
#'     cpBounds <- getBounds(cpMat)
#'     if (cpBounds[3]==0) {return(interval_complement(Intervals(c(-Inf, Inf))))}
#'   } else {
#'     cpBounds <- t(apply(cpMat, 1, getBounds))
#'     #cpBounds <- cpBounds[cpBounds[,3]==1,]
#'   }
#'
#'   cpInterval <- Intervals(c(-Inf, Inf))
#'   if (cp_cutoff==0) {
#'     cpInterval = Intervals(c(-Inf,Inf))
#'   } else {
#'     if (NCOL(cpBounds) == 1) {
#'       cpInterval <-  suppressWarnings(interval_intersection(cpInterval, interval_complement(Intervals(cpBounds[1:2]))))
#'
#'     } else {
#'       if (NROW(cpBounds) != 0) {
#'         for (i in 1:NROW(cpBounds)) {
#'
#'           cpInterval <- suppressWarnings(interval_intersection(cpInterval, interval_complement(Intervals(cpBounds[i,1:2]))))
#'         }
#'       }
#'     }
#'   }
#'
#'   return(cpInterval)
#' }
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#' #' Workhorse function for computing conditioning sets.
#' #'
#' #' Returns the interval for phi such that Tree(y'(phi)) (where the definition of y'(phi) depends on nu)
#' #' contains the series of branch "branch".Either tree with model=TRUE or X and y must be provided, not both.
#' #'
#' #' @param tree An rpart object. Must have been built with model=TRUE arguement unless X and y are provided separately
#' #' @param nu The contrast vector that defines the parameter.
#' #' @param branch A vector of strings specifying the branch.
#' #' @param sibs does nu have the form nu_sib? If so, we can do some computational speedups!!!
#' #'
#' #' @return an object of class Interval that defines the set S.
#' @export
#' getInterval_grow <- function(tree, nu, branch, sibs=FALSE) {
#'
#'   branch <- paste("dat$", branch)
#'   minbucket <- tree$control$minbucket
#'
#'   dat <- tree$model
#'   if (is.null(dat)) {stop("Must build rpart tree with parameter model=TRUE")}
#'   y <- dat[,1]
#'   X <- dat[,-1]
#'   n <- nrow(X)
#'   p <- NCOL(X)
#'
#'
#'   ### In the first layer, nothing gets zeroed out
#'   nulled <- rep(TRUE, n)
#'   #H_prev <- diag(1,n,n) - matrix(1/n,n,n)
#'   full_interval = Intervals(c(-Inf, Inf))
#'   leafSize=n
#'
#'   ### O(n^2*p) work up front to get the inequalities??
#'   ### Because for each of n*p possible branch, it is O(n) work to fill up the vector.
#'   ### Consider optimizing this
#'   ### If you made the C vecs as you went using an ordered version of the Xs, this might be faster
#'   #OH yeah DUH. If you sort the Xs first, the cs have a predictable form
#'   # This work is likely avoidable
#'   cs<- array(NA, c(p,n,n))
#'   for (j in 1:p) {
#'     cs[j,,]  <-  sapply(X[,j], function(u) (u <= X[,j]))
#'   }
#'
#'   # O(n)
#'   norm_nu <- sum(nu^2)
#'
#'   # O(n)
#'   nut_y <- sum(nu*y)
#'
#'   #nunu <- nu%*%t(nu)/norm_nu
#'   y_proj <-nu*(nut_y/norm_nu)
#'
#'
#'
#'   regions <- rep(1,n)
#'
#'   L <- length(branch)
#'   cur <- 1
#'
#'   coeffs_full <- matrix(0, nrow=n*p*L, ncol=3)
#'
#'   ### MULTIPLY EVERYTHING IN LOOP BY O(L)
#'   for (l in 1:L) {
#'
#'     #if (NCOL(C_prev) > 1) {
#'     #  regions <- rowSums(C_prev) ## O(n)
#'     #} else{regions <- C_prev}
#'
#'     ### Each of these should be O(n). Techically number of regions can scale with L?? So maybe O(nL)??
#'     Hnu <- nu - (sapply(1:max(regions), function(u) mean(nu[regions==u])))[regions]
#'     Hy <- y - (sapply(1:max(regions), function(u) mean(y[regions==u])))[regions]
#'     Hyproj <- y_proj - (sapply(1:max(regions), function(u) mean(y_proj[regions==u])))[regions]
#'     Hminus <- Hy-Hyproj
#'
#'     ### All of these are O(n) too
#'     Hnu2 <- Hnu[nulled]
#'     Hy2 <- Hy[nulled]
#'     Hminus2 <- Hminus[nulled]
#'
#'     cy<- eval(parse(text = branch[l]))*nulled
#'     if (sum(cy) < minbucket | sum(cy) > (leafSize-minbucket)) {
#'       coeffs_full[cur:(cur+1),] <- c(1,1,1)
#'       break
#'     }
#'
#'     ### This is an interesting choice.
#'     ### If the region could have been formed with only one split instead of two,
#'     ### say that no interval is possible??
#'     if (sum(cy)==sum(leafSize)) {
#'       coeffs_full[cur:(cur+1),] <- c(1,1,1)
#'       break
#'     }
#'
#'     Hcy <- cy - (sapply(1:max(regions), function(u) mean(cy[regions==u])))[regions]
#'
#'     ### I'm doing O(n)
#'     denom <- sum(cy)*(leafSize - sum(cy))/leafSize
#'
#'
#'
#'     indices <- which(cy==1)
#'     partialA <- (sum(Hnu[indices])^2/norm_nu^2)/denom
#'     partialB <- as.numeric((2*sum(Hy[indices])*(t(nu)%*%Hcy)/
#'                               norm_nu - 2*(t(cy)%*%Hnu)^2%*%nut_y/
#'                               norm_nu^2)/denom)
#'     partialC <- (sum(Hminus[indices]))^2/denom
#'
#'
#'
#'
#'     for (j in 1:p) {
#'       c_now <- cs[j,nulled, nulled]
#'       hardwork <- sort(rowSums(c_now),index.return=TRUE) ## I hope this is O(n)!!!!
#'       num1sSORT <- hardwork$x
#'       sorted <- hardwork$ix
#'       len <- length(sorted)
#'
#'       denomsSORT <- (num1sSORT*(leafSize-num1sSORT)/leafSize)
#'       minDenomSORT <- ((minbucket-1)*(leafSize-(minbucket-1))/leafSize)
#'
#'       ### These should be O(n)
#'       unus <- cumsum(Hnu2[sorted])[1:(sum(nulled)-1)]
#'       uys <- cumsum(Hy2[sorted])[1:(sum(nulled)-1)]
#'       umin <- cumsum(Hminus2[sorted])[1:(sum(nulled)-1)]
#'       denomsSORT <- denomsSORT[1:(sum(nulled)-1)]
#'
#'       ### This is the clever thing where we hope we are only paying O(n) because
#'       ### we are not multiplying matricies
#'       other_As <- unus^2/norm_nu^2
#'       other_Bs <- 2*uys*unus/norm_nu- 2*(unus^2*(nut_y))/norm_nu^2
#'       other_Cs <- umin^2
#'
#'
#'       avec <- other_As/denomsSORT-partialA
#'       bvec <- other_Bs/denomsSORT-partialB
#'       cvec <- other_Cs/denomsSORT-partialC
#'
#'
#'       ### SHOOT THIS IS GOOD FOR CATEGORICAL. BUT SOMEHOW NOT FOR REGULAR???
#'       #### PUT BACK SOON
#'       uniqueindices <- c(cumsum(table(hardwork$x)))[1:(length(unique(hardwork$x))-1)]
#'       uniqueindices <- uniqueindices[denomsSORT[uniqueindices] > minDenomSORT]
#'
#'       ###### NOTE TO SELD: PUT BACK IN!!!
#'       num_coeffs <- length(avec[uniqueindices])
#'       if (num_coeffs !=0 & !is.na(uniqueindices[1])) {
#'         coeffs_full[cur:(cur+num_coeffs-1),] <- cbind(avec,bvec,cvec)[uniqueindices,]
#'       }
#'       cur <- cur+num_coeffs
#'     }
#'
#'     #### If this inversion is slow could use the rank 1 update!!!
#'     #### Update these things for next level.
#'     regions <- regions+cy
#'     #H_prev <- H_prev - Hcy%*%t(Hcy)/sum(Hcy[cy==1])
#'     leafSize <- sum(cy)
#'     nulled <- cy==1 ## requires that we are testing SIBLING nodes.
#'
#'   }
#'
#'   coeffs_full <- coeffs_full[1:cur,]
#'   phi_bounds <-  t(apply(coeffs_full, 1, getBounds))
#'
#'   ### INSIDE INTERVALS
#'   numIn <- nrow(phi_bounds[phi_bounds[,3]==1,])
#'   if (!is.null(numIn)) {
#'     inside_comp_mat <- matrix(c(-Inf, Inf), nrow=2*numIn, ncol=2, byrow=TRUE)
#'     inside_comp_mat[1:numIn,2] <- phi_bounds[phi_bounds[,3]==1,1]
#'     inside_comp_mat[(numIn+1):(2*numIn), 1] <- phi_bounds[phi_bounds[,3]==1,2]
#'     ints_comp_inside <- Intervals(inside_comp_mat, closed=c(TRUE, TRUE))
#'     intersection1 <- interval_complement(interval_union(ints_comp_inside))
#'   } else {
#'     intersection1 <- Intervals(c(-Inf, Inf))
#'   }
#'
#'   ### OUTSIDE INTERVALS
#'   ints_outside <- Intervals(phi_bounds[phi_bounds[,3]==0,1:2], closed=c(TRUE, TRUE))
#'   intersection2 <- interval_complement(interval_union(ints_outside))
#'   #if(length(intersection1)==0) {
#'   #res1 <- intersection2
#'   #} else {
#'   res1 <- suppressWarnings(interval_intersection(intersection1, intersection2))
#'   #}
#'
#'   return(res1)
#' }
#'
#'
#'
#'
#'
#' #' Turn quadratic coefficients into an interval
#' #'
#' #' @param vec stores a,b,c coeffs of a quadratic equation
#' #' @return A vector describing the regions where the quadratic is negative.
#' #' @keywords internal
#' #' @noRd
#' getBounds <- function(vec) {
#'
#'   #vec[abs(vec) < 1e-15] <- 0
#'   a <- vec[1]
#'   b <- vec[2]
#'   c <- vec[3]
#'   tol <- 10*10^(-10)
#'   if (abs(a) < tol & abs(b) < tol & abs(c) < tol) {
#'     return(c(-Inf,Inf, 1))
#'   }
#'
#'   ## If we just have a straight line
#'   if (abs(a) < tol) {
#'     if (abs(b) < tol) {
#'       if (c < 1e-6) {return(c(-Inf,Inf, 1))
#'         } else {return(c(-Inf, Inf, 0))}
#'     }
#'     root = -c/b
#'     if (b>0) {
#'       return(c(-Inf,root,1))
#'     } else{
#'       return(c(root, Inf, 1))
#'     }
#'   }
#'
#'   det <- b^2-4*a*c
#'   if (det < 0) {
#'     if (a > 0) {
#'       #return(simpleError("STOP! EMPTY INTERVAL"))
#'       return(c(-Inf,Inf,0)) ## TO AVOID ERRORS FOR NOW PLZ CHANGE LATER
#'     }
#'     return(c(-Inf, Inf, 1))
#'   }
#'   if (det==0) {
#'     if (a > 0) {c(-b/(2*a), -b/(2*a), 1)}
#'     return(c(-Inf, Inf, 1))
#'   }
#'   ## Down here determinant is positive: 2 roots
#'   root1 <- (-b + sqrt(det))/(2*a)
#'   root2 <- (-b - sqrt(det))/(2*a)
#'   if (a > 0) {
#'     return(c(min(root1, root2), max(root1,root2), 1))
#'   } else {
#'     return(c(min(root1, root2), max(root2,root1), 0))
#'   }
#' }
