### Almost all of the code (and documentation) in this file comes directly from the ``Outference" package,
#### written by Shuxiao Chen.
### The code and the documentation can be found at https://github.com/shuxiaoc/outference.

#' Get a pvalue for a selective Z-test using a truncated normal distribution. Not typically needed by users.
#'
#' Based on the similar function from the Outference package,
#' which accomplishes a similar task. Code modified from code found at
#' https://github.com/shuxiaoc/outference.
#' This function shouldn't be
#' needed by most users (it is called internally by \code{branchInference()}), but is needed to reproduce our paper simulations.
#'
#' @param phiInterval the conditioning set (truncation interval). An object of class \code{Interval}, where the rows represent the union of
#' disjoint intervals on the real line.
#' @param nu the vector that defines the parameter of interest. We are testing the hypothesis that nu^T mu = 0.
#' @param y the response data y.
#' @param sigma the (assumed known) noise standard deviation. The assumption is that y_i ~ N(mu_i, sigma^2).
#' @return a p-value.
#' @importFrom intervals interval_complement
#' @importFrom intervals interval_intersection
#' @export
#' @examples
#' data(blsdata, package="treevalues")
#' bls.tree <-rpart::rpart(kcal24h0~hunger+disinhibition+resteating+rrvfood+liking+wanting,
#'     model = TRUE,  data = blsdata, cp=0.02)
#' branch <- getBranch(bls.tree, 2)
#' left_child <- getRegion(bls.tree,2)
#' right_child <- getRegion(bls.tree,3)
#' nu_sib <- left_child/sum(left_child) -  right_child/sum(right_child)
#' S_sib <- getInterval(bls.tree, nu_sib,branch)
#' correctPVal(S_sib, nu_sib, blsdata$kcal24h0, sd(blsdata$kcal24h0))
#' # Same answer as using branchInference()
#' branchInference(bls.tree, branch, type="sib")$pval
correctPVal <- function(phiInterval, nu, y, sigma) {
  delta1 <- phiInterval
  delta2 <- interval_complement(Intervals(c(-abs(t(nu)%*%y), abs(t(nu)%*%y))))
  numerator_region <- suppressWarnings(interval_intersection(delta1, delta2))


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
#' Slight adaptation of the code from the Outference package. Code modified from code found at
#' https://github.com/shuxiaoc/outference.
#'
#' @keywords internal
#' @noRd
#'
#' @param E1 numerator interval
#' @param E2 denominator interval
#' @param scale determines the grid search
#' @keywords internal
#' @noRd
myTNRatioApprox <- function(E1, E2, scale = NULL) {
  E1 <- sortE(E1)
  E2 <- sortE(E2)

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

#' normalCDF helper function. Computes P(low <= X <= up), where X ~ N(0, sd).
#' Code take from https://github.com/shuxiaoc/outference
#' @param lo
#' @keywords internal
#' @noRd
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
