### NOTE THAT EVERYTHING HERE CAME FROM OUTFERENCE

# ----- Truncated Normal Distribution -----
#' Survival function of truncated normal distribution
#'
#' This function returns the upper tail probability of a truncated normal distribution
#'     at quantile \code{q}.
#'
#' Let \eqn{X} be a normal random variable with \code{mean} and \code{sd}. Truncating
#'     \eqn{X} to the set \eqn{E} is equivalent to conditioning on \eqn{{X \in E}}. So this function
#'     returns \eqn{P(X \ge q | X \in E)}.
#'
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
#'
#' For \eqn{Z ~ N(0, 1)}, we have the approximation
#'     \eqn{P(Z \ge z) \approx }\code{magicfun(z)*exp(-z^2/2)}.
#'
#'
#' @param z, the number where the function is evaluated.
#'
#' @return This function returns the value of the function evaluated at \code{z}.
#'
#' @references Bryc, Wlodzimierz. "A uniform approximation to the right normal tail integral."
#'     Applied mathematics and computation 127.2 (2002): 365-374.
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
#'
#' @param E an "Intervals" object or a matrix where rows represents
#'     a union of intervals with \emph{positive} but possibly infinite endpoints.
#'
#' @return This function returns an "Intervals" object or a matrix depending on the input.
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





# ----- Truncated Chi-square Distribution -----

#' Survival function of truncated central chi-squared distribution
#'
#' This function returns the upper tail probability of a truncated central chi-squared distribution
#'     at quantile \code{q}.
#'
#' Let \eqn{X} be a central chi-squared random variable with \code{df} degrees of freedom. Truncating
#'     \eqn{X} to the set \eqn{E} is equivalent to conditioning on \eqn{{X \in E}}. So this function
#'     returns \eqn{P(X \ge q | X \in E)}.
#'
#' @keywords internal
#'
#' @param q the quantile.
#' @param df the degree of freedom.
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
#' @references Canal, Luisa. "A normal approximation for the chi-square distribution."
#'     Computational statistics & data analysis 48.4 (2005): 803-808.
TChisqSurv <- function(q, df, E, approx = FALSE) {
  # check if truncation is empty (i.e. truncated to the empty set)
  if (nrow(E) == 0) {
    stop("truncation set is empty")
  }

  # check if truncation is the whole support: [0, Inf]
  if (isSameIntervals(E, intervals::Intervals(c(0, Inf)))) {
    return(stats::pchisq(q = q, df = df, lower.tail = FALSE))
  }

  # E is not empty and is not the whole support
  # i.e. 0 < P(X in E) < 1
  region <- suppressWarnings(intervals::interval_intersection(E, intervals::Intervals(c(q, Inf))))

  # check if the result is 0 or 1
  if(nrow(region) == 0) return(0)
  if (isSameIntervals(E, region)) return(1)

  # now we are in a nontrivial situation

  if (approx) {
    res = TChisqRatioApprox(df = df, E1 = region, E2 = E)
    return(max(0, min(1, res)))

  }


  # first try exact calculation
  denom <- TChisqProb(df = df, E = E)
  num <- TChisqProb(df = df, E = region)

  if (denom < 1e-100 || num < 1e-100) {
    res = TChisqRatioApprox(df = df, E1 = region, E2 = E)
    return(max(0, min(1, res)))
  }

  # we know denom and num are reasonably > 0
  return(max(0, min(1, num/denom)))
}



#' Approximation of the ratio of two chi-squared probabilities
#'
#' This function returns an approximation of \eqn{P(X \in E1)/P(X \in E2)}, where
#'     \eqn{X} is a central chi-squared random variable with \code{df} degrees of freedom.
#'
#' @keywords internal
#'
#' @param df degree of freedom of the chi-squared random variable.
#' @param E1,E2 "Intervals" objects or matrices where rows represents
#'     a union of intervals with \emph{positive and finite} endpoints.
#'
#' @return This function returns the value of the approximation.
#'
#' @references Bryc, Wlodzimierz. "A uniform approximation to the right normal tail integral."
#'     Applied mathematics and computation 127.2 (2002): 365-374.
#' @references Canal, Luisa. "A normal approximation for the chi-square distribution."
#'     Computational statistics & data analysis 48.4 (2005): 803-808.
TChisqRatioApprox <- function(df, E1, E2) {

  # the transform that makes x into a N(0, 1) r.v. such that
  # P(X >= x) = P(Z >= Chisq2N(x)), X ~ chisq(df), Z ~ N(0, 1)
  # this function can take either scaler, vector or matrix
  Chisq2N <- function(x, df, tol = 1e-6) {

    if (is.numeric(x) && length(x) == 1) {
      if (x <= tol) { # x <= 0
        return(-Inf)
      }
      if (x == Inf) {
        return(Inf)
      }
      # we know x > 0 and x is finite
      x <- (x/df)^(1/6) - (1/2) * (x/df)^(1/3) + (1/3) * (x/df)^(1/2)
      mu <- 5/6 - 1/(9*df) - 7/(648*df^2) + 25/(2187*df^3)
      sig <- sqrt(1/(18*df) + 1/(162*df^2) - 37/(11664*df^3))
      return((x-mu)/sig)
    }

    if (is.vector(x)) {
      return(vapply(X = x, FUN = Chisq2N, FUN.VALUE = 0.1, df = df))
    }

    if (is.matrix(x)) {
      return(structure(vapply(X = x, FUN = Chisq2N, FUN.VALUE = 0.1, df = df), dim = dim(x)))
    }

    return(intervals::Intervals())

  }


  E1 <- Chisq2N(E1, df)
  E1 <- sortE(E1) # notice that Chisq2N can be negative
  E2 <- Chisq2N(E2, df)
  E2 <- sortE(E2)

  # now we want P(Z in E1) / P(Z in E2), Z ~ N(0, 1)
  return(TNRatioApprox(E1, E2))
}



#' Probability of a chi-squared random variable in a single interval
#'
#' This function returns \eqn{P(lo \le X \le up)}, where \eqn{X} is a
#'     central chi-squared random variable with \code{df} degrees of freedom.
#'
#' @keywords internal
#'
#' @param df degree of freedom of the chi-squared random variable.
#' @param lo,up quantiles.
#'
#' @return This function returns the desired probability.
TChisqProbEachInt <- function(df, lo, up) {

  if (up == Inf) { # P(X > lo)
    return(stats::pchisq(lo, df, lower.tail = FALSE))
  }

  # we know up < Inf, want P(lo <= X <= up).
  # we want small - small (big - big will mask small numbers),
  # since chisq is not centered at zero, we try two kinds of calculations
  try1 <- stats::pchisq(lo, df, lower.tail = FALSE) - stats::pchisq(up, df, lower.tail = FALSE)
  if (try1 != 0) return(try1)

  # we know try1 gives 0
  try2 <- stats::pchisq(up, df, lower.tail = TRUE) - stats::pchisq(lo, df, lower.tail = TRUE)
  return(try2)

}


#' Probability of a chi-squared random variable in a union of intervals
#'
#' This function returns \eqn{P(X \in E)}, where \eqn{X} is a
#'     central chi-squared random variable with \code{df} degrees of freedom.
#'
#' @keywords internal
#'
#' @param df degree of freedom of the chi-squared random variable.
#' @param E an "Intervals" object or a matrix where rows represents
#'     a union of disjoint intervals.
#'
#' @return This function returns the desired probability.
TChisqProb <- function(df, E) {

  # sum cdf over each disjoint interval of E
  res <- sum(sapply(1:nrow(E), function(v) {
    return(TChisqProbEachInt(df, E[v, 1], E[v, 2]))
  }))
  return(res)
}













# ----- Truncated t and F distribution -----

#' Approximation of the ratio of two F probabilities
#'
#' This function returns an approximation of \eqn{P(X \in E1)/P(X \in E2)}, where
#'     \eqn{X} is a central F random variable with \code{df1, df2} degrees of freedom.
#'
#' @keywords internal
#'
#' @param df1,df2 degrees of freedom of the F random variable.
#' @param E1,E2 "Intervals" objects or matrices where rows represents
#'     a union of intervals with \emph{positive and finite} endpoints.
#'
#' @return This function returns the value of the approximation.
#'
#' @references Bryc, Wlodzimierz. "A uniform approximation to the right normal tail integral."
#'     Applied mathematics and computation 127.2 (2002): 365-374.
#' @references Peizer, David B., and John W. Pratt. "A normal approximation for binomial,
#'     F, beta, and other common, related tail probabilities, I." Journal of the American
#'     Statistical Association 63.324 (1968): 1416-1456.
TFRatioApprox <- function(df1, df2, E1, E2) {
  if (df1 == 1) {

    # reduces to P(T^2 in E1) / P(T^2 in E2), where T ~ t(df2),
    # which further reduces to
    # (2*P(T in sqrt(E1))) / (2*P(T in sqrt(E2))) = P(T in sqrt(E1)) / P(T in sqrt(E2))
    E1 <- sqrt(E1)
    E2 <- sqrt(E2)

    # By "A Study of the Accuracy of Some Approximations for t, chi-squre, and F Tail Probabilities",
    # Equation (2.6),
    # P(T >= t) = P( Z >= (df2 - 2/3 + 1/(10*df2)) * sqrt(log(1 + t^2/df2) / (df2 - 5/6)) )

    # Since E is a union of intervals, then we can simply use
    # P(T in E) = P( Z in (df2 - 2/3 + 1/(10*df2)) * sqrt(log(1 + E^2/df2) / (df2 - 5/6)) )

    # recall after taking sqrt(), we want P(T in E1) / P(T in E2)

    # For x >= 0, P(T >= x) = P(Z >= T2N(x)), T ~ t(df), Z ~ N(0, 1)
    # this function can take vector or matrix
    T2N <- function(x, df, tol = 1e-6) {

      # when x is a scaler
      if (is.numeric(x) && length(x) == 1) {
        if (x < -tol) stop("x must be >= 0")
        if (abs(x) <= tol) return(0)
        const <- (df - 2/3 + 1/(10*df)) / sqrt(df - 5/6)
        x <- const * sqrt(log(1 + x^2/df))
        return(x)
      }

      if (is.vector(x)) {
        return(vapply(X = x, FUN = T2N, FUN.VALUE = 0.1, df = df))
      }

      if (is.matrix(x)) {
        return(structure(vapply(X = x, FUN = T2N, FUN.VALUE = 0.1, df = df), dim = dim(x)))
      }

      stop("x must be a scaler, vector or matrix")


    }

    # notice that T2N(Inf) = Inf, and T2N >= 0 always holds
    E1 <- T2N(E1, df = df2)
    E1 <- finiteE(E1)
    E2 <- T2N(E2, df = df2)
    E2 <- finiteE(E2)

    # now we want P(Z in region) / P(Z in E), Z ~ N(0, 1)
    return(TNRatioApprox(E1, E2))
  }

  # we know df1 > 1, so it does not reduce to the square of t-distribution
  # recall we want P(X in E1) / P(X in E2), X ~ F(df1, df2)

  # By "A Study of the Accuracy of Some Approximations for t, chi-square, and F Tail Probabilities",
  # Equation 4.3,
  # P(X >= x) = P( Z >= d * sqrt( (1 + q * g(S / (n*p)) + p * g(R / (n*q))) / ((n + 1/6) * p * q) ) ),
  # where S = (df2 - 1)/2, R = (df1 - 1)/2, n = (df1 + df2 - 2)/2,
  # p = df2/(df1 * x + df2), q = 1 - p,
  # d = S + 1/6 - (n + 1/3) * p + 0.04 * ((q/df2) - (p/df1) + (q - 0.5)/(df1 + df2)),
  # g(y) = (1 - y^2 + 2 * y * log(y)) / (1 - y^2) if y != 1, y > 0, and
  # g(0) = 1, g(1) = 0.

  # again we can simply replace x by region and E

  # some helper functions

  g <- function(x, tol = 1e-6) {
    if (x < -tol) { # x < 0
      stop("log(x) while x < 0")
    }
    if (abs(x) <= tol) { # x == 0
      return(1)
    }
    if (abs(x - 1) <= tol) { # x == 1
      return(0)
    }

    # now we know tol < x < 1-tol, or x > 1+tol
    if (x < 1-tol) {
      return((1 - x^2 + 2 * x * log(x)) / (1 - x)^2)
    }

    # x > 1 + tol
    return(-g(1/x))
  }


  # P(X >= x) = P(Z >= F2N(x)), where X ~ F(df1, df2), Z ~ N(0, 1)
  # this function can take scalor, vector or matrix
  F2N <- function(x, df1, df2, tol = 1e-6) {

    # when x is a scaler
    if (is.numeric(x) && length(x) == 1) {
      if (x <= tol) { # x <= 0
        return(-Inf)
      }
      if (x == Inf) {
        return(Inf)
      }
      # we know x > 0 and x is finite
      S <- (df2 - 1)/2
      R <- (df1 - 1)/2
      n <- (df1 + df2)/2 - 1
      p <- df2/(df1 * x + df2)
      q <- 1 - p
      d1 <- S + 1/6 - (n+1/3) * p
      d2 <- d1 + 0.02 * (q/(S+0.5) - p/(R+0.5) + (q-0.5)/(n+1))
      z <- d2 * sqrt((1 + q * g(S/(n*p)) + p * g(R/(n*q))) / ((n + 1/6) * p * q))
      return(z)
    }

    if (is.vector(x)) {
      return(vapply(X = x, FUN = F2N, FUN.VALUE = 0.1, df1 = df1, df2 = df2))
    }

    if (is.matrix(x)) {
      return(structure(vapply(X = x, FUN = F2N, FUN.VALUE = 0.1, df1 = df1, df2 = df2), dim = dim(x)))
    }

    stop("x must be a scaler, vector or matrix")
  }


  E1 <- F2N(E1, df1, df2)
  E1 <- sortE(E1) # notice that F2N can be negative
  E2 <- F2N(E2, df1, df2)
  E2 <- sortE(E2)

  # now we want P(Z in region) / P(Z in E), Z ~ N(0, 1)
  return(TNRatioApprox(E1, E2))

}





#' Probability of an F random variable in a single interval
#'
#' This function returns \eqn{P(lo \le X \le up)}, where \eqn{X} is a
#'     central F random variable with \code{df1, df2} degrees of freedom.
#'
#' @keywords internal
#'
#' @param df1,df2 degrees of freedom of the central F random variable.
#' @param lo,up quantiles.
#'
#' @return This function returns the desired probability.
TFProbEachInt <- function(df1, df2, lo, up) {

  if (up == Inf) { # P(X > lo)
    return(stats::pf(q = lo, df1 = df1, df2 = df2, lower.tail = FALSE))
  }

  # we know up < Inf, want P(lo <= X <= up).
  # we want small - small (big - big will mask small numbers),
  # since f distr is not centered at zero, we try two kinds of calculations
  try1 <- stats::pf(q = lo, df1 = df1, df2 = df2, lower.tail = FALSE) - stats::pf(q = up, df1 = df1, df2 = df2, lower.tail = FALSE)
  if (try1 != 0) return(try1)

  # we know try1 gives 0
  try2 <- stats::pf(q = up, df1 = df1, df2 = df2, lower.tail = TRUE) - stats::pf(q = lo, df1 = df1, df2 = df2, lower.tail = TRUE)
  return(try2)

}



#' Probability of an F random variable in a union of intervals
#'
#' This function returns \eqn{P(X \in E)}, where \eqn{X} is a
#'     central F random variable with \code{df1, df2} degrees of freedom.
#'
#' @keywords internal
#'
#' @param df1,df2 degrees of freedom of the central F random variable.
#' @param E an "Intervals" object or a matrix where rows represents
#'     a union of disjoint intervals.
#'
#' @return This function returns the desired probability.
TFProb <- function(df1, df2, E) {

  # sum cdf over each disjoint interval of E
  res <- sum(sapply(1:nrow(E), function(i) {
    return(TFProbEachInt(df1 = df1, df2 = df2, lo = E[i, 1], up = E[i, 2]))
  }))
  return(res)

}


#' Survival function of truncated central F distribution
#'
#' This function returns the upper tail probability of a truncated central F distribution
#'     at quantile \code{q}.
#'
#' Let \eqn{X} be a central F random variable with \code{df1, df2} degrees of freedom. Truncating
#'     \eqn{X} to the set \eqn{E} is equivalent to conditioning on \eqn{{X \in E}}. So this function
#'     returns \eqn{P(X \ge q | X \in E)}.
#'
#' @keywords internal
#'
#' @param q the quantile.
#' @param df1,df2 the degrees of freedom.
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
#' @references Peizer, David B., and John W. Pratt. "A normal approximation for binomial,
#'     F, beta, and other common, related tail probabilities, I." Journal of the American
#'     Statistical Association 63.324 (1968): 1416-1456.
TFSurv <- function(q, df1, df2, E, approx = FALSE) {
  # check if truncation is empty (i.e. truncated to the empty set)
  if (nrow(E) == 0) {
    stop("truncation set is empty")
  }

  # check if truncation is the whole support: [0, Inf]
  if (isSameIntervals(E, intervals::Intervals(c(0, Inf)))) {
    return(stats::pf(q = q, df1 = df1, df2 = df2, lower.tail = FALSE))
  }

  # E is not empty and is not the whole support
  # i.e. 0 < P(X in E) < 1
  region <- suppressWarnings(intervals::interval_intersection(E, intervals::Intervals(c(q, Inf))))

  # check if the result is 0 or 1
  if(nrow(region) == 0) return(0)
  if (isSameIntervals(E, region)) return(1)

  # now we are in a nontrivial situation

  if (approx) {
    res <- TFRatioApprox(df1 = df1, df2 = df2, E1 = region, E2 = E)
    return(max(0, min(1, res)))

  }


  # first try exact calculation
  denom <- TFProb(df1 = df1, df2 = df2, E = E)
  num <- TFProb(df1 = df1, df2 = df2, E = region)

  if (denom < 1e-100 || num < 1e-100) {

    res <- TFRatioApprox(df1 = df1, df2 = df2, E1 = region, E2 = E)
    return(max(0, min(1, res)))

  }

  # we know denom and num are reasonably > 0
  return(max(0, min(1, num/denom)))

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

# return selective ci for v^T mu (i.e. a single parameter)

#' compute selective confidence intervals
#'
#' This function computes the selective confidence intervals.
#'
#' @keywords internal
#'
#' @param v, the contrast vector.
#' @param y, the response.
#' @param sigma, the noise level \eqn{\sigma}.
#' @param truncation, the truncation set for the \eqn{Z}-statistic.
#' @param alpha, the significance level.
#'
#' @return This function returns a vector of lower and upper confidence limits.
#'
computeCI <- function(v, y, sigma, truncation, alpha) {
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
