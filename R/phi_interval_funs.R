#' Workhorse function for computing conditioning sets. Not typically needed by users.
#'
#' Returns the interval for phi such that Tree(y'(phi,nu))
#' contains the set of regions induced by branch \code{branch}. An rpart object built with
#' \code{model=TRUE} must be provided. This function shouldn't be
#' needed by most users (it is called internally by \code{branchInference}), but is needed to reproduce our paper simulations.
#'
#' @param tree An rpart object. Must have been built with \code{model=TRUE} argument.
#' @param nu The vector that defines the parameter nu^T mu.
#' @param branch A vector of strings specifying the branch.
#' @param sib If you are doing inference and nu=nu_sib, this takes advantage of some computational speedups!
#' @param grow Set this to true if you only want to compute S_grow (ignore cost complexity pruning).
#' @param prune Set this to true if you only want to see what gets REMOVED from S_grow during pruning
#' @return an object of class Interval that defines the set S.
#' @importFrom stats predict
#' @export
#' @examples
#' data(blsdata, package="treevalues")
#' bls.tree <-rpart::rpart(kcal24h0~hunger+disinhibition+resteating+rrvfood+liking+wanting,
#'       model = TRUE, data = blsdata, cp=0.02)
#' branch <- getBranch(bls.tree, 2)
#' left_child <- getRegion(bls.tree,2)
#' right_child <- getRegion(bls.tree,3)
#' nu_sib <- left_child/sum(left_child) -  right_child/sum(right_child)
#' S_sib <- getInterval(bls.tree, nu_sib,branch)
#' S_sib
#' branchInference(bls.tree, branch, type="sib")$condset
getInterval <- function(tree, nu, branch,sib=FALSE,grow=FALSE,prune=FALSE) {

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
  # This work is only needed for s_prune, and so if for some reason someone has asked for only s_{grow}
  # you can skip it!!
  if (length(branch) > 1 & !grow) {

    ### Total SSE for all descendents of branch??
    vecChildren <- eval(parse(text = paste(branch, collapse=" & ")))
    childrenGains <-   sum((y[as.logical(vecChildren)]-mean(y[as.logical(vecChildren)]))^2)  - sum((y[as.logical(vecChildren)]-predict(tree)[as.logical(vecChildren)])^2)
    childrenSize <- length(unique(tree$where[as.logical(vecChildren)]))-1

    inequalities <- sapply(branch, function(u) eval(parse(text = u)))
    counts <- apply(inequalities, 1, cumsum)
    membership <- t(sapply(1:NROW(counts), function(u) counts[u,]==u))

    ### Deals with the half of the descendents who do not depend on PHI.
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
    ### Deals with necessary pruning calculations in case where branch has length 1.
    if (!grow) {
      off_subgroup_membership <- !eval(parse(text = branch))
      vec <- off_subgroup_membership
      gainZfull <-   sum((y[as.logical(vec)]-mean(y[as.logical(vec)]))^2)  - sum((y[as.logical(vec)]-predict(tree)[as.logical(vec)])^2)
      sizeZfull <- length(unique(tree$where[as.logical(vec)])) - 1

      vecChildren <- eval(parse(text = branch))
      childrenGains <-   sum((y[as.logical(vecChildren)]-mean(y[as.logical(vecChildren)]))^2)  - sum((y[as.logical(vecChildren)]-predict(tree)[as.logical(vecChildren)])^2)
      childrenSize <- length(unique(tree$where[as.logical(vecChildren)])) - 1
    }
  }

  cp_cutoff <- cp*(sum((y-mean(y))^2))

  ### In the first layer, nothing gets zeroed out
  nulled <- rep(TRUE, n)
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
  y_perp <- y - nu*(nut_y/norm_nu)



  regions <- rep(1,n)

  L <- length(branch)
  cur <- 1

  coeffs_full <- matrix(0, nrow=n*p*L, ncol=3)

  gainsA <- rep(0,L)
  gainsB <- rep(0,L)
  gainsC <- rep(0,L)

  for (l in 1:L) {

    nu2 <- nu[nulled]
    y_perp2 <- y_perp[nulled]

    ### Fix the comments Anna.
    cy<- eval(parse(text = branch[l]))*nulled
    indices_left <- as.logical(cy==1)
    indices_right <- as.logical(cy==0 & nulled==1)
    n_left <- sum(indices_left)
    n_right <- sum(indices_right)



    ### This is a trivial case where the split did like an (n-0) split. In this case, say that
    ### no interval is possible??
    if (sum(cy)==sum(leafSize)) {
      coeffs_full[cur:(cur+1),] <- c(1,1,1)
      break
    }


    ##### NOTE TO SELF--- the last terms should always cancel and so we should be able to ignore??
    partialC <- (sum(y_perp[indices_left])^2/n_left + sum(y_perp[indices_right])^2/n_right)
    partialA <- 1/(norm_nu^2)*(sum(nu[indices_left])^2/n_left + sum(nu[indices_right])^2/n_right)
    partialB <- 2*1/(norm_nu)*(
      (sum(nu[indices_left])*sum(y_perp[indices_left]))/n_left +
        (sum(nu[indices_right])*sum(y_perp[indices_right]))/n_right)

    #### Save these because they get used during pruning.
    gainsA[L-l+1] <- partialA - 1/(norm_nu^2)*sum(nu[as.logical(nulled)])^2/leafSize
    gainsB[L-l+1] <- partialB - 2/(norm_nu)*(sum(nu[as.logical(nulled)])*sum(y_perp[as.logical(nulled)]))/leafSize
    gainsC[L-l+1] <- partialC - sum(y_perp[as.logical(nulled)])^2/leafSize

    for (j in 1:p) {
      c_now <- cs[j,nulled, nulled]
      hardwork <- sort(rowSums(c_now),index.return=TRUE) ## I hope this is O(n)!!!!
      num1sSORT <- hardwork$x
      sorted <- hardwork$ix
      reverse_sorted <- hardwork$ix[length(hardwork$ix):1]
      len <- length(sorted)

      ### I definitely mostly forget why this minbucket hack worked at all :(.
      denomsSORT <- (num1sSORT*(leafSize-num1sSORT)/leafSize)
      minDenomSORT <- ((minbucket-1)*(leafSize-(minbucket-1))/leafSize)

      ### These should be O(n)
      Inus_left <- cumsum(nu2[sorted])[1:(sum(nulled)-1)]
      Iys_left <- cumsum(y_perp2[sorted])[1:(sum(nulled)-1)]
      Inus_right <- cumsum(nu2[reverse_sorted])[((sum(nulled)-1)):1]
      Iys_right <- cumsum(y_perp2[reverse_sorted])[((sum(nulled)-1)):1]

      denoms_left <- num1sSORT[1:(sum(nulled)-1)]
      denoms_right <- leafSize-num1sSORT[1:(sum(nulled)-1)]


      other_As <-1/(norm_nu^2)*(Inus_left^2/denoms_left + Inus_right^2/denoms_right)
      other_Bs <- 2/(norm_nu)*(Inus_left*Iys_left/denoms_left + Inus_right*Iys_right/denoms_right)
      other_Cs <- (Iys_left^2/denoms_left + Iys_right^2/denoms_right)


      avec <- other_As-partialA
      bvec <- other_Bs-partialB
      cvec <- other_Cs-partialC



      ### ANNA BE SURE TO EVENTUALLY CHECK ON DISCRETE OR MINBUCK DATA
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

    return(suppressWarnings(interval_intersection(res1, cpInterval)))
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
