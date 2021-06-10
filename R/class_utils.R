## Class utils

#' Summary of inference for a branch in a tree.
#'
#' @param  x branch_inference object
#'
#' @export
#' @noRd
print.branch_inference <- function(x,...) {
  CIlev <- 1-x$alpha
  CI <- paste(round(x$confint,4), collapse=", ")
  cat("\nSample statistic: ", x$samplemean)
  cat("\n", CIlev*100, "% confidence interval: ", CI,sep="")
  cat("\nType: ", x$type,sep="")
  cat("\n", "p-value for test that param=", x$c, ": ", x$pval,sep="")
  cat("\n", "Conditioning Set: \n",sep="")
  print(x$condset)
}


## Class utils

#' Summary of inference for a branch in a tree.
#'
#' @param  object branch_inference object
#'
#' @export
#' @noRd
summary.branch_inference <- function(object,...) {
  CIlev <- 1-object$alpha
  CI <- paste(round(object$confint,4), collapse=", ")
  cat("\nSample statistic: ", object$samplemean)
  cat("\n", CIlev*100, "% confidence interval: ", CI,sep="")
  cat("\nType: ", object$type,sep="")
  cat("\n", "p-value for test that param=", object$c, ": ", object$pval,sep="")
  cat("\n", "Conditioning Set: \n",sep="")
  print(object$condset)
}


