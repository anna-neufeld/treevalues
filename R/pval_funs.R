#' Get a valid pvalue for a test
#'
#' Based on the similar function from the Outference package,
#' which accomplishes a similar task.
#'
#' @param phiInterval the selective inference interval
#' @param nu the contrast vector for the hypothesis test
#' @param y the response data y
#' @param sigma the (assumed known) noise standard deviation
correctPVal <- function(phiInterval, nu, y, sigma) {
  delta1 <- phiInterval
  delta2 <- interval_complement(Intervals(c(-abs(t(nu)%*%y), abs(t(nu)%*%y))))
  numerator_region <- suppressWarnings(interval_intersection(delta1, delta2))

  ## This really should never happen, because delta1 should always contain
  ## nu^T y
  if(nrow(numerator_region)==0) {return(0)}

  # try exact calculation
  denom <- 0
  for (i in 1:nrow(delta1)) {
    denom = denom + myGetNormProb(delta1[i,1], delta1[i,2], sqrt(sum(nu^2))*sigma)
  }
  numer <- 0
  for (i in 1:nrow(numerator_region)) {
    numer <- numer + myGetNormProb(numerator_region[i,1], numerator_region[i,2], sqrt(sum(nu^2))*sigma)
  }

  ### If the p-value is approximately 0/0, have numeric instability
  ### So use the approximation suggested in the Outference package.
  if (denom < 1e-100 || numer < 1e-100) {
    ### The magic function assumes we are working with a standard normal
    ### So divide by the (known) standard deviation everywhere.
    delta1 <- delta1/(sqrt(sum(nu^2))*sigma)
    numerator_region <- numerator_region/(sqrt(sum(nu^2))*sigma)

    ### standardized pivot
    q <- ((t(nu)%*%y)^2)/(sqrt(sum(nu^2))*sigma)
    res <- myTNRatioApprox(numerator_region, delta1, scale =q*q)
    if (is.nan(res)) { # try approximation one more time
      res <- myTNRatioApprox(numerator_region, delta1, scale = NULL)
    }
    return(max(0, min(1, res)))
  }
  # we know denom and num are both reasonably >  0
  res <- numer / denom
  #print(numer/denom)
  # force the result to lie in [0, 1]
  return(max(0, min(1, res)))
}

#' Approximate the survival function of a truncated normal
#'
#' Slight adaptation of the code from Outference
#'
#' @export
#'
#' @param E1 numerator interval
#' @param E2 denominator interval
myTNRatioApprox <- function(E1, E2, scale = NULL) {

  #E1 <- mySortE(E1)
  #E2 <- mySortE(E2)
  E1 <- SortE(E1)
  E2 <- SortE(E2)

  #### If scale not provided, try a bunch until you find something OK.
  if (is.null(scale)) {
    temp <- (c(E1, E2))^2
    scale.grid <- stats::quantile(temp, probs = seq(0, 1, 0.2))

    for(scale in scale.grid) {
      if(scale==Inf) {break}
      temp <- myTNRatioApprox(E1, E2, scale = scale)
      if (!is.na(temp)) {
        return(temp)
      }
      # if temp is NaN, proceed to the next loop
    }
    # if all scale.grid does not work, then return NaN
    return(NaN)
  }

  scale <- as.vector(scale) ## Avoids warnings

  num1 <- magicfun((E1[, 1])) * exp(-(E1[, 1]^2 - scale)/2)

  num2 <- magicfun((E1[, 2])) * exp(-(E1[, 2]^2 - scale)/2)
  num2[is.nan(num2)] <- 0

  denom1 <- magicfun((E2[, 1])) * exp(-(E2[, 1]^2 - scale)/2)

  denom2 <- magicfun((E2[, 2])) * exp(-(E2[, 2]^2 - scale)/2)
  denom2[is.nan(denom2)] <- 0

  res <- sum(num1-num2)/sum(denom1-denom2)
  return(res)
}

#' Approximate a normal density for the purpose of taking ratios
#'
#' Add a nice citation like they have in Outference.

magicfun = function(z){
  z2 <- z*z
  z3 <- z*z*z
  temp <- (z2 + 5.575192695 * z + 12.77436324) /
    (sqrt(2*pi) * z3 + 14.38718147*z2 + 31.53531977*z + 2*12.77436324)
  return(temp)
}


#' normalCDF
#'
#' I think this is also from outference and I'm not convinced that I need it
#'
#' Hi!!!!
myGetNormProb <- function(lo, up,sd) {
  if (up == Inf) {
    return(stats::pnorm(lo, 0, sd, lower.tail = FALSE))
  }
  if (lo==-Inf) {
    return(stats::pnorm(up, 0, sd, lower.tail = TRUE))
  }
  # we know up < Inf, want P(lo <= X <= up).
  # we want small - small (big - big will mask small numbers),
  try1 <- stats::pnorm(lo, 0, sd, lower.tail = FALSE) - stats::pnorm(up, 0, sd, lower.tail = FALSE)
  if (try1 != 0) return(try1)
  try2 <- stats::pnorm(up, 0, sd, lower.tail = TRUE) - stats::pnorm(lo, 0, sd, lower.tail = TRUE)
  return(try2)
}



#mySortE <- function(E) {
#  E.mySortEd <- lapply(1:nrow(E), function(i){
#    temp <- as.numeric(E[i, ])
#    if (temp[1] <= 0 & temp[2] <= 0) {
#      return(sort(-temp))
#    }
#    if (temp[1] >= 0 & temp[2] >= 0) {
#      return(sort(temp))
#    }
    # we know temp[1] < 0, temp[2] > 0 OR temp[1] > 0, temp[2] < 0
#    temp <- abs(temp)
#    return(rbind(c(0, temp[1]), c(0, temp[2])))
#  })
#  E.mySortEd <- do.call(rbind, E.mySortEd)
  # in order to use the approximation, we translate Inf to a large number
 # return((E.mySortEd))
#}

myFiniteE <- function(E) {
  ind.inf <- which(E == Inf)
  if (length(ind.inf) == 0) return(E)
  # we know there are some infinite entries
  E.max <- max(E[-ind.inf])
  E[which(E == Inf)] <- max(10000, E.max * 2)
  return(E)
}



