## Class utils

#' Summary of inference for a branch in a tree.
#'
#' @param  x branch_inference object
#'
#' @export
print.branch_inference <- function(x, ...) {
  CIlev <- 1-x$alpha
  CI <- paste(round(x$confint,4), collapse=", ")
  cat("\n", CIlev*100, "% confidence interval: ", CI,sep="")
  cat("\nType: ", x$type,sep="")
  cat("\n", "p-value for test that param=", x$c, ": ", x$pval,sep="")
  cat("\n", "Conditioning Set: \n",sep="")
  print(x$condset)
}


## Class utils

#' Summary of inference for a branch in a tree.
#'
#' @param  x branch_inference object
#'
#' @export
summary.branch_inference <- function(x, ...) {
  CIlev <- 1-x$alpha
  CI <- paste(round(x$confint,4), collapse=", ")
  cat("\n", CIlev*100, "% confidence interval: ", CI,sep="")
  cat("\nType: ", x$type,sep="")
  cat("\n", "p-value for test that param=", x$c, ": ", x$pval,sep="")
  cat("\n", "Conditioning Set: \n",sep="")
  print(x$condset)
}


