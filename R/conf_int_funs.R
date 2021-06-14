### Almost all of the code (and documentation) in this file comes directly from the ``Outference" package,
#### written by Shuxiao Chen.
### The code and the documentation can be found at https://github.com/shuxiaoc/outference.


# ----- Truncated Normal Distribution -----
#' Survival function of truncated normal distribution. This comes from the Outference Package https://github.com/shuxiaoc/outference.
#'
#' This function returns the upper tail probability of a truncated normal distribution
#'     at quantile \code{q}.
#'
#' Let \eqn{X} be a normal random variable with \code{mean} and \code{sd}. Truncating
#'     \eqn{X} to the set \eqn{E} is equivalent to conditioning on \eqn{{X \in E}}. So this function
#'     returns \eqn{P(X \ge q | X \in E)}.
#'
#' @keywords internal
#' @noRd
#'
#' @param q the quantile.
#' @param mean the mean parameter
#' @param sd the standard deviation
#' @param E the truncation set, an "Intervals" object or a matrix where rows represents
#'     a union of disjoint intervals.
#' @param approx should the approximation algorithm be used? Default is \code{FALSE},
#'     where the approximation is not used in the first place. But when the result is wacky,
#'     the approximation will be used.
#'
#' @return This function returns the value of the survival function evaluated at quantile \code{q}.
#'
#' @references Bryc, Wlodzimierz. "A uniform approximation to the right normal tail integral."
#'     Applied mathematics and computation 127.2 (2002): 365-374.
TNSurv <- function(q, mean, sd, E, approx = FALSE) {
  # check if truncation is empty (i.e. truncated to the empty set)
  if (nrow(E) == 0) {
    stop("truncation set is empty")
  }

  # check if truncation is the whole real line
  if (isSameIntervals(E, intervals::Intervals(c(-Inf, Inf)))) {
    return(stats::pnorm(q, mean, sd, lower.tail = FALSE))
  }

  # E is not empty and is not the whole real line,
  # i.e. 0 < P(X in E) < 1

  # we want P(X > q | X in E) = P(X >= q AND X in E) / P(X in E)
  # {X >= q} = {Z >= (q-mean)/sd}
  # {X in E} = {Z in (E-mean)/sd}
  # Z ~ N(0, 1)
  q <- (q-mean)/sd
  E <- (E-mean)/sd
  mean <- 0
  sd <- 1
  q2 <- q*q
  region <- suppressWarnings(intervals::interval_intersection(E, intervals::Intervals(c(q, Inf))))
  # check if the result is 0 or 1
  if(nrow(region) == 0) return(0)
  if (isSameIntervals(E, region)) return(1)

  # transform region and E so that intervals have positive endpoints
  region <- sortE(region)
  E <- sortE(E)

  # we want P(Z in region) / P(Z in E)
  # try approximate calculation
  if (approx) {
    res <- TNRatioApprox(region, E, scale = q2)
    if (is.nan(res)) { # try approximation one more time
      res <- TNRatioApprox(region, E, scale = NULL)
    }
    return(max(0, min(1, res)))
  }

  # try exact calculation
  denom <- TNProb(E)
  num <- TNProb(region)

  if (denom < 1e-100 || num < 1e-100) {
    res <- TNRatioApprox(region, E, scale = q2)
    if (is.nan(res)) { # try approximation one more time
      res <- TNRatioApprox(region, E, scale = NULL)
    }
    return(max(0, min(1, res)))
  }

  # we know denom and num are both reasonably > 0

  res <- num / denom
  # force the result to lie in [0, 1]
  return(max(0, min(1, res)))
}


#' A helper function for approximating normal tail probabilities
#' This comes directly from the Outference Package.
#' The code and the documentation can be found at https://github.com/shuxiaoc/outference.
#'
#' For \eqn{Z ~ N(0, 1)}, we have the approximation
#'     \eqn{P(Z \ge z) \approx }\code{magicfun(z)*exp(-z^2/2)}.
#'
#'
#'
#'
#' @param z, the number where the function is evaluated.
#'
#' @return This function returns the value of the function evaluated at \code{z}.
#'
#' @references Bryc, Wlodzimierz. "A uniform approximation to the right normal tail integral."
#'     Applied mathematics and computation 127.2 (2002): 365-374.
#' @keywords internal
#' @noRd
magicfun = function(z){
  z2 <- z*z
  z3 <- z*z*z
  temp <- (z2 + 5.575192695 * z + 12.77436324) /
    (sqrt(2*pi) * z3 + 14.38718147*z2 + 31.53531977*z + 2*12.77436324)
  return(temp)
}


#' Make endpoints of intervals finite
#'
#' This function modifies a union of intervals with positive but possibly infinite endpoints
#'    into a union of intervals with positive and \emph{finite} endpoints, while ensuring
#'    the probability of a \eqn{N(0, 1)} falling into it numerically the same.
#'
#'
#' @param E an "Intervals" object or a matrix where rows represents
#'     a union of intervals with \emph{positive} but possibly infinite endpoints.
#'
#' @return This function returns an "Intervals" object or a matrix depending on the input.
#' @keywords internal
#' @noRd
finiteE <- function(E) {
  ind.inf <- which(E == Inf)
  if (length(ind.inf) == 0) return(E)
  # we know there are some infinite entries
  E.max <- max(E[-ind.inf])
  E[which(E == Inf)] <- max(10000, E.max * 2)
  return(E)
}


#' Make endpoints of intervals positive
#'
#' This function modifies a union of intervals with possibly negative enpoints
#'     into a union of intervals with \emph{positive} endpoints, while ensuring
#'    the probability of a \eqn{N(0, 1)} falling into it numerically the same.
#'
#' @keywords internal
#' @noRd
#' @param E an "Intervals" object or a matrix where rows represents
#'     a union of intervals with \emph{positive} but possibly infinite endpoints.
#'
#' @return This function returns an "Intervals" object or a matrix depending on the input.
#' @keywords internal
#' @noRd
sortE <- function(E) {
  E.sorted <- lapply(1:nrow(E), function(i){
    temp <- as.numeric(E[i, ])
    if (temp[1] <= 0 & temp[2] <= 0) {
      return(sort(-temp))
    }
    if (temp[1] >= 0 & temp[2] >= 0) {
      return(sort(temp))
    }
    # we know temp[1] < 0, temp[2] > 0 OR temp[1] > 0, temp[2] < 0
    temp <- abs(temp)
    return(rbind(c(0, temp[1]), c(0, temp[2])))
  })
  E.sorted <- do.call(rbind, E.sorted)
  # in order to use the approximation, we translate Inf to a large number
  return(finiteE(E.sorted))
}


#' Approximation of the ratio of two normal probabilities
#'
#' This function returns an approximation of \eqn{P(Z \in E1)/P(Z \in E2)}, where \eqn{Z ~ N(0, 1)}.
#'
#' @keywords internal
#'
#' @param E1,E2 "Intervals" objects or matrices where rows represents
#'     a union of intervals with \emph{positive and finite} endpoints.
#' @param scale scaling parameter.
#'
#' @return This function returns the value of the approximation.
#'
#' @references Bryc, Wlodzimierz. "A uniform approximation to the right normal tail integral."
#'     Applied mathematics and computation 127.2 (2002): 365-374.
#' @keywords internal
#' @noRd
TNRatioApprox <- function(E1, E2, scale = NULL) {

  if (is.null(scale)) {
    temp <- (c(E1, E2))^2
    scale.grid <- stats::quantile(temp, probs = seq(0, 1, 0.2))

    for(scale in scale.grid) {
      temp <- TNRatioApprox(E1, E2, scale = scale)
      if (!is.na(temp)) {
        return(temp)
      }
      # if temp is NaN, proceed to the next loop
    }

    # if all scale.grid does not work, then return NaN
    return(NaN)
  }
  num1 <- magicfun(E1[, 1]) * exp(-(E1[, 1]^2 - scale)/2)
  num2 <- magicfun(E1[, 2]) * exp(-(E1[, 2]^2 - scale)/2)
  denom1 <- magicfun(E2[, 1]) * exp(-(E2[, 1]^2 - scale)/2)
  denom2 <- magicfun(E2[, 2]) * exp(-(E2[, 2]^2 - scale)/2)
  res <- sum(num1-num2)/sum(denom1-denom2)
  return(res)
}



#' Probability of a standard normal in a single interval
#'
#' This function returns \eqn{P(lo \le Z \le up)}, where \eqn{Z ~ N(0, 1)}.
#'
#' @keywords internal
#' @noRd
#'
#' @param lo,up quantiles.
#'
#' @return This function returns the desired probability.
TNProbEachInt <- function(lo, up) {
  if (up == Inf) {
    return(stats::pnorm(lo, 0, 1, lower.tail = FALSE))
  }
  # we know up < Inf, want P(lo <= X <= up).
  # we want small - small (big - big will mask small numbers),

  try1 <- stats::pnorm(lo, 0, 1, lower.tail = FALSE) - stats::pnorm(up, 0, 1, lower.tail = FALSE)
  if (try1 != 0) return(try1)

  try2 <- stats::pnorm(up, 0, 1, lower.tail = TRUE) - stats::pnorm(lo, 0, 1, lower.tail = TRUE)
  return(try2)

}

#' Probability of a standard normal in a union of intervals
#'
#' This function returns \eqn{P(Z \in E)}, where \eqn{Z ~ N(0, 1)}.
#'
#' @keywords internal
#'
#' @param E an "Intervals" object or a matrix where rows represents
#'     a union of disjoint intervals.
#'
#' @return This function returns the desired probability.
TNProb <- function(E) {
  # sum cdf over each disjoint interval of E
  res <- sum(sapply(1:nrow(E), function(v) {
    return(TNProbEachInt(E[v, 1], E[v, 2]))
  }))
  return(res)
}

# ----- Helper Functions ------

#' Comparison between two intervals
#'
#' This functions returns \code{TRUE} if and only if two intervals are the same.
#'
#' @keywords internal
#'
#' @param int1,int2 "Intervals" objects.
#'
#' @return This function returns the desired logical result.
#' @keywords internal
#' @noRd
isSameIntervals <- function(int1, int2) {

  # first make int1, int2 to the default order
  int1 <- intervals::reduce(int1)
  int2 <- intervals::reduce(int2)

  if (nrow(int1) != nrow(int2)) return(FALSE)

  # int1 and int2 has the same number of intervals

  if (sum(int1 != int2) > 0) return(FALSE)

  # int1 and int2 has the same elements
  return(TRUE)
}


# ----- computing the confidence intervals -----

#' Return selective confidence interval for v^T mu (a single parameter).
#'
#' Compute selective confidence interval for parameter v^T mu based on a truncated normal distribution. A slight modification of code found in the
#' Outference package, available at https://github.com/shuxiaoc/outference.
#' This function shouldn't be
#' needed by most users (it is called internally by \code{branchInference}), but is needed to reproduce our paper simulations.
#' @export
#'
#' @param truncation, the truncation set for the statistic v'y.
#' Computes a confidence interval for the mean of a truncated normal distribution.
#' @param v the vector that defines the parameter of interest; v^T mu
#' @param y the observed response vector
#' @param sigma The known noise standard deviation. If unknown, we recommend a conservative estimate. If it
#' is left blank, we automatically use a conservative estimate.
#' @param alpha, the significance level.
#'
#' @return This function returns a vector of lower and upper confidence limits.
#' @examples
#' data(blsdata, package="treevalues")
#' bls.tree <- rpart::rpart(kcal24h0~hunger+disinhibition+resteating+rrvfood+liking+wanting,
#'     model = TRUE, data = blsdata, cp=0.02)
#' branch <- getBranch(bls.tree, 2)
#' full_result <- branchInference(bls.tree, branch, type="sib")
#' left_child <- getRegion(bls.tree,2)
#' right_child <- getRegion(bls.tree,3)
#' nu_sib <- left_child/sum(left_child) -  right_child/sum(right_child)
#' S_sib <- getInterval(bls.tree, nu_sib,branch)
#'
#' computeCI(nu_sib, blsdata$kcal24h0, sd(blsdata$kcal24h0),S_sib)
#' full_result$confint
computeCI <- function(v, y, sigma=NULL, truncation, alpha=0.05) {
  ### Conservative guess
  if (is.null(sigma)) {
    sigma <- sd(y)
  }

  ### Standardization expected.
  truncation <- Intervals(as.matrix(truncation)/(sigma*sqrt(sum(v^2))))



  #browser()
  vTv <- sum(v*v)
  scale <- sigma * sqrt(vTv)
  q <- sum(v*y) / scale

  fun <- function(x) {
    return(TNSurv(q, x/scale, 1, truncation))
  }

  # L: fun.L(L) = 0
  fun.L <- function(x) {
    return(fun(x) - alpha/2)
  }
  # U: fun.U(U) = 0
  fun.U <- function(x) {
    return(fun(x) - (1-alpha/2))
  }

  # find the starting point (x1, x2) such that
  # fun.L(x1), fun.U(x1) <= 0 AND fun.L(x2), fun.U(x2) >= 0.
  # i.e. fun(x1) <= alpha/2 AND fun(x2) >= 1-alpha/2.

  # find x1 s.t. fun(x1) <= alpha/2
  # what we know:
  # fun is monotone incresing;
  # fun(x) = NaN if x too small;
  # fun(x) > alpha/2 if x too big.
  # so we can do a modified bisection search to find x1.
  step <- 0
  x1.up <- q * scale + scale
  x1 <- q * scale - 10 * scale
  f1 <- fun(x1)
  while(step <= 20) {
    if (is.na(f1)) { # x1 is too small
      x1 <- (x1 + x1.up) / 2
      f1 <- fun(x1)
    }
    else if (f1 > alpha/2) { # x1 is too big
      x1.up <- x1
      x1 <- x1 - 10 * 1.4^step
      f1 <- fun(x1)
      step <- step + 1
    }
    else { # fun(x1) <= alpha/2, excited!
      break
    }
  }

  # find x2 s.t. fun(x2) <= 1 - alpha/2
  # what we know:
  # fun is monotone incresing;
  # fun(x) = NaN if x too big;
  # fun(x) < 1 - alpha/2 if x too small.
  # again can do a modified bisection search to find x2.
  step <- 0
  x2 = q * scale + 10 * scale
  x2.lo = q * scale - scale
  f2 = fun(x2)
  while(step <= 20) {
    if (is.na(f2)) { # x2 is too big
      x2 <- (x2 + x2.lo) / 2
      f2 <- fun(x2)
    }
    else if (f2 < 1 - alpha/2) { # x2 is too small
      x2.lo <- x2
      x2 <- x2 + 10 * 1.4^step
      f2 <- fun(x2)
      step <- step + 1
    }
    else { # fun(x2) >= 1 - alpha/2, excited!
      break
    }
  }

  # if the above search does not work, set up a grid search
  # for starting points
  if (is.na(f1)||(f1 > alpha/2)||is.na(f2)||(f2 < 1-alpha/2)) {
    grid <- seq(from = q * scale - 1000*scale, to = q*scale + 1000*scale)
    value <- sapply(grid, fun)
    # want max x1: fun(x1) <= alpha/2
    ind1 <- rev(which(value <= alpha/2))[1]
    x1 <- grid[ind1]
    #f1 = value[ind1]
    # want min x2: fun(x2) >= 1-alpha/2
    ind2 <- which(value >= 1 - alpha/2)[1]
    x2 <- grid[ind2]
    #f2 = value[ind2]
  }

  # if the above fails, then either x1, x2 = NA, so uniroot() will throw error,
  # in which case we set (-Inf, Inf) as the CI

  # we know the functions are increasing

  L <- tryCatch({
    stats::uniroot(fun.L, c(x1, x2), extendInt = "upX", tol = 1e-5)$root
  }, error = function(e) {
    -Inf
  })


  U <- tryCatch({
    stats::uniroot(fun.U, c(x1, x2), extendInt = "upX", tol = 1e-5)$root
  }, error = function(e) {
    Inf
  })

  return(c(L, U))
}
