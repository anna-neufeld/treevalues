getInterval_EXPERIMENT <- function(base_tree, nu, splits) {
  
  ### rpart orderes things nicely YAY
  dat <- base_tree$model
  y <- dat[,1]
  X <- dat[,-1]
  n <- nrow(X)
  p <- NCOL(X)
  
  #Pi_perp <- diag(rep(1,n)) - nu%*%t(nu)/sum(nu^2)
  
  
  cs <- array(NA, c(p,n,n))
  for (j in 1:p) {
    mat <-  sapply(X[,j], function(u) (u <= X[,j]))
    cs[j,,] <- mat
  }
  
  ### In the first layer, nothing gets zeroed out
  nulled <- rep(1, n)
  C_prev <- rep(1,n)
  H_prev <- diag(1,n) - C_prev%*%solve(t(C_prev)%*%C_prev)%*%t(C_prev)
  full_interval = Intervals(c(-Inf, Inf))
  leafSize=n
  
  norm_nu <- sum(nu^2)
  y_proj <- (nu%*%t(nu)/norm_nu)%*%y
  nut_y <- as.numeric(t(nu)%*%y)
  
  nunu <- nu%*%t(nu)/sum(nu^2)
  
  for (l in 1:length(splits)) {
    X2 <- X[as.logical(nulled),]
    y2 <- y[as.logical(nulled)]
    n2 <- sum(nulled)
    
    cs2<- array(NA, c(p,n2,n2))
    for (j in 1:p) {
      mat <-  sapply(X2[,j], function(u) (u <= X2[,j]))
      cs2[j,,] <- mat
    }
    
    if (NCOL(C_prev) > 1) {
      regions <- rowSums(C_prev)
    } else{regions <- C_prev}
    
    Hnu <- nu - (sapply(1:max(regions), function(u) mean(nu[regions==u])))[regions]
    Hy <- y - (sapply(1:max(regions), function(u) mean(y[regions==u])))[regions]
    Hyproj <- y_proj - (sapply(1:max(regions), function(u) mean(y_proj[regions==u])))[regions]
    Hminus <- Hy-Hyproj
    
    Hnu2 <- Hnu[as.logical(nulled)]
    Hy2 <- Hy[as.logical(nulled)]
    Hminus2 <- Hminus[as.logical(nulled)]
    
    cy<- eval(parse(text = splits[l]))*nulled
    
    Hcy <- cy - (sapply(1:max(regions), function(u) mean(cy[regions==u])))[regions]
    
    denom <- colSums(cy*(H_prev%*%cy))
    
    indices <- which(cy==1)
    partialA <- (sum(Hnu[indices])^2/norm_nu^2)/denom
    partialB <- as.numeric((2*sum(Hy[indices])*(t(nu)%*%H_prev%*%cy)/
                              norm_nu - 2*(t(cy)%*%Hnu)^2%*%nut_y/
                              norm_nu^2)/denom)
    partialC <- (sum(Hminus[indices]))^2/denom
    
  
    
    
    for (j in 1:p) {
      c_now <- t(apply(cs2[j,,], 1, function(u) u))
      c_now3 <- c_now[rowSums(c_now)!=sum(nulled),]
      num1s <- rowSums(c_now)
      sorted <- sort(num1s, index.return=TRUE)$ix ## I hope this is O(n)!!!!
      len <- length(sorted)
      num1sSORT <- sort(rowSums(c_now))
      denomsSORT <- (num1sSORT*(leafSize-num1sSORT)/leafSize)
      
      
      denoms <- (num1s2*(leafSize-num1s2)/leafSize)[(rowSums(c_now2)!=sum(nulled))]
      num1s3 <- rowSums(c_now3)
      denoms3 <- (num1s3*(leafSize-num1s3)/leafSize)[(rowSums(c_now3)!=sum(nulled))]
      
      
      
      #other_As2 <- apply(c_now3, 1, function(u) sum(Hnu2[u==1])^2/norm_nu^2)
      #other_Bs2 <- apply(c_now3, 1, function(u) 2*sum(Hy2[u==1])*sum(Hnu2[u==1])/
       #                    norm_nu
        #                 - 2*(sum(Hnu2[u==1])^2%*%nut_y)/
         #                  norm_nu^2)
      #other_Cs2 <- apply(c_now3, 1, function(u) sum(Hminus2[u==1])^2)
      
      
      
      unus <- cumsum(Hnu2[sorted])[1:(sum(nulled)-1)]
      uys <- cumsum(Hy2[sorted])[1:(sum(nulled)-1)]
      umin <- cumsum(Hminus2[sorted])[1:(sum(nulled)-1)]
      denomsSORT <- denomsSORT[1:(sum(nulled)-1)]
      
      other_As <- unus^2/norm_nu^2
      other_Bs <- 2*uys*unus/norm_nu- 2*(unus^2*(nut_y))/norm_nu^2
      other_Cs <- umin^2
      
      #print(l)
      #print(all.equal(as.numeric(sort(other_As)), as.numeric(sort(other_As2))))
      #print(all.equal(as.numeric(sort(other_Bs)), as.numeric(sort(other_Bs2))))
      #print(all.equal(as.numeric(sort(other_Cs)), as.numeric(sort(other_Cs2))))
      
      indices2 <- denomsSORT != 0
      avec <- other_As/denomsSORT-partialA
      bvec <- other_Bs/denomsSORT-partialB
      cvec <- other_Cs/denomsSORT-partialC
      
      #indices3 <- denoms3 != 0
      #avec2 <- other_As2/denoms3-partialA
      #bvec2 <- other_Bs2/denoms3-partialB
      #cvec2 <- other_Cs2/denoms3-partialC
      
      #print(all.equal(as.numeric(sort(avec2[indices3])), as.numeric(sort(avec[indices2]))))
      #print(all.equal(as.numeric(sort(bvec2[indices3])), as.numeric(sort(bvec[indices2]))))
      #print(all.equal(as.numeric(sort(cvec2[indices3])), as.numeric(sort(cvec[indices2]))))
      
      coeffs <- cbind(avec,bvec,cvec)
      phi_bounds <-  t(apply(coeffs[indices2,], 1, getBounds))
      
      ### INSIDE INTERVALS
      numIn <- nrow(phi_bounds[phi_bounds[,3]==1,])
      if (!is.null(numIn)) {
        if (numIn > 0) {
          inside_comp_mat <- matrix(c(-Inf, Inf), nrow=2*numIn, ncol=2, byrow=TRUE)
          inside_comp_mat[1:numIn,2] <- phi_bounds[phi_bounds[,3]==1,1]
          inside_comp_mat[(numIn+1):(2*numIn), 1] <- phi_bounds[phi_bounds[,3]==1,2]
          ints_comp_inside <- Intervals(inside_comp_mat, closed=c(TRUE, TRUE))
          intersection1 <- interval_complement(interval_union(ints_comp_inside))
        }
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
    }
    
    #### If this inversion is slow could use the rank 1 update!!!
    #### Update these things for next level.
    C_prev <- cbind(C_prev, cy)
    H_prev <- H_prev - Hcy%*%t(Hcy)/sum(Hcy[cy==1])
    leafSize <- sum(cy)
    nulled <- cy ## requires that we are testing SIBLING nodes.
  
  }
  
  return(full_interval)
}
