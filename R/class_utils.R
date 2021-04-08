## Class utils 

#' Summary of inference for a branch in a tree. 
#' 
#' @param  x branch_inference object
#'
#' @export
print.branch_inference <- function(x, ...) {
  cat("\n ", x$alpha*100, "%confidence interval:", x$confint,  "\n")
  cat("Type: ", x$type, "\n")
  cat("\n ", "p-value for test that param=", x$c, ": ", x$pval)
}


## Class utils 

#' Summary of inference for a branch in a tree. 
#' 
#' @param  x branch_inference object
#'
#' @export
summary.branch_inference <- function(x, ...) {
  cat("\n ", x$alpha*100, "%confidence interval:", x$confint,  "\n")
  cat("\n ", "p-value for test that param=", x$c, ": ", x$pval)
}


