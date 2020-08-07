TTSurv <- function(q, NC, DF, S) {
    # check if truncation is empty (i.e. truncated to the empty set)
    if (nrow(S) == 0) {
      stop("truncation set is empty")
    }
    
    # check if truncation is the whole real line
    if (isSameIntervals(S, intervals::Intervals(c(-Inf, Inf)))) {
      return(stats::pt(q, df=DF, ncp=NC, lower.tail = FALSE))
    }
    
    # E is not empty and is not the whole real line,
    # i.e. 0 < P(X in E) < 1
    
    
   
    R <- intervals::Intervals(c(q, Inf))
    region <- suppressWarnings(intervals::interval_intersection(S, R))
    if(nrow(region) == 0) return(0)
    if (isSameIntervals(S, region)) return(1)

    # transform region and E so that intervals have positive endpoints
    #region <- sortE(region)
    #S <- sortE(S)
  
    # try exact calculation
    areaNum <- sum(apply(region, 1, function(u) pt(u[2], df=DF, ncp=NC)-pt(u[1], df=DF, ncp=NC)))
    areaDenom <- sum(apply(S, 1, function(u) pt(u[2],df=DF, ncp=NC)-pt(u[1], df=DF, ncp=NC)))
    
    res <- areaNum / areaDenom
    return(max(0, min(1, res)))
}




T_conf_int <- function(q, DF, S, alpha) {
  fun <- function(x) {
    return(TTSurv(q, x, DF, S))
  }
  
  #scale <- 1
  
  # find x1 s.t. fun(x1) <= alpha/2
  # what we know:
  # fun is monotone incresing;
  # fun(x) = NaN if x too small;
  # fun(x) > alpha/2 if x too big.
  # so we can do a modified bisection search to find x1.
  #step <- 0
  #x1.up <- q * scale + scale
  #x1 <- q * scale - 10 * scale
  #f1 <- fun(x1)
  #while(step <= 20) {
  #3  if (is.na(f1)) { # x1 is too small
  #    x1 <- (x1 + x1.up) / 2
  #    f1 <- fun(x1)
  #  }
   # else if (f1 > alpha/2) { # x1 is too big
  #    x1.up <- x1
  #    x1 <- x1 - 10 * 1.4^step
  #    f1 <- fun(x1)
  #    step <- step + 1
  #  }
   # else { # fun(x1) <= alpha/2, excited!
  #    break
  #  }
 # }
  
  # find x2 s.t. fun(x2) <= 1 - alpha/2
  # what we know:
  # fun is monotone incresing;
  # fun(x) = NaN if x too big;
  # fun(x) < 1 - alpha/2 if x too small.
  # again can do a modified bisection search to find x2.
 # step <- 0
  #x2 = q * scale + 10 * scale
  #x2.lo = q * scale - scale
  #f2 = fun(x2)
  #while(step <= 20) {
  #  if (is.na(f2)) { # x2 is too big
    #  x2 <- (x2 + x2.lo) / 2
    #  f2 <- fun(x2)
    #}
    #else if (f2 < 1 - alpha/2) { # x2 is too small
   #   x2.lo <- x2
   #   x2 <- x2 + 10 * 1.4^step
    #  f2 <- fun(x2)
    #  step <- step + 1
  #  }
  #  else { # fun(x2) >= 1 - alpha/2, excited!
  #    break
  #  }
  #}

  fLOW <- function(x) {
    tryCatch({
      TTSurv(q,NC=x, DF=DF, S=S) - alpha/2
    }, error = function(e) {
      5
    })
  }
  
  fHIGH <- function(x) {
    TTSurv(q,NC=x, DF=DF, S=S) - (1-alpha/2)
  }
  
 
  
  L <- Inf
  U <- Inf
  tries <- 1
  maxWidth <- 100
  
  x11 <- x12 <- q-1
  x21 <- x22 <- q+1
  
  while (tries < maxWidth) {
    if (L==Inf) {
      L <- tryCatch({
        stats::uniroot(fLOW, c(x11, x21), extendInt = "upX", tol = 1e-3, maxiter=10000)$root
      }, error = function(e) {
        -Inf
      })
      x11 <- x11 - 1
      x21 <- x21 + 1
    }
    if (U==Inf) {
      U <- tryCatch({
        stats::uniroot(fHIGH, c(x12,x22), extendInt = "upX", tol = 1e-3, maxiter=10000)$root
      }, error = function(e) {
        Inf
      })
      x12 <- x12 - 1
      x22 <- x22 + 1
    }
    tries <- tries + 1
  }
  
  return(c(L, U))

}

