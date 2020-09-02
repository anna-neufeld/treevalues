getInterval_TEMP <- function(base_tree, nu, splits) {
  ### rpart orderes things nicely YAY
  minbucket <- base_tree$control$minbucket - 1

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

      #### FIND ONLY THE UNIQUE ONES. SHIT DIDN"T WE SORT?? SO THIS WON"T WORK??
      ### OH wai no the sorting is PERF.
      ### Need to make sure that this didn't break my usual simulations.

      ### SHOOT THIS IS GOOD FOR CATEGORICAL. BUT SOMEHOW NOT FOR REGULAR???
      #### PUT BACK SOON
      #uniqueindices <- c(cumsum(table(hardwork$x)))[1:(length(unique(hardwork$x))-1)]
      #uniqueindices <- uniqueindices[denomsSORT[uniqueindices]!=0]

      ###### NOTE TO SELD: PUT BACK IN!!!
      #uniqueindices <- uniqueindices[num1sSORT[uniqueindices] >= minbucket &
       #                                (leafSize - num1sSORT[uniqueindices]) >= minbucket]

      uniqueindices <- (1:length(denomsSORT))[denomsSORT !=0]

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

