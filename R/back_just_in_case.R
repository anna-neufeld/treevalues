mysortE <- function(E) {
  E.mysortEd <- lapply(1:nrow(E), function(i){
        temp <- as.numeric(E[i, ])
    if (temp[1] <= 0 & temp[2] <= 0) {
      return(sort(-temp))
          }
    if (temp[1] >= 0 & temp[2] >= 0) {
      return(sort(temp))
    }
 #we know temp[1] < 0, temp[2] > 0 OR temp[1] > 0, temp[2] < 0
    temp <- abs(temp)
    return(rbind(c(0, temp[1]), c(0, temp[2])))
  })
  E.mysortEd <- do.call(rbind, E.mysortEd)
 #in order to use the approximation, we translate Inf to a large number
 return((E.mysortEd))
}

myFiniteE <- function(E) {
  ind.inf <- which(E == Inf)
  if (length(ind.inf) == 0) return(E)
  # we know there are some infinite entries
  E.max <- max(E[-ind.inf])
  E[which(E == Inf)] <- max(10000, E.max * 2)
  return(E)
}