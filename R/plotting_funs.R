#' This function allows users to plot their rpart object with all splits labeled with p-values and all nodes
#' labeled with confidence intervals.
#'
#' @export
#'
#' @param tree An rpart tree. Note that the tree must have been built with arguements maxcompete=0 and maxsurrogates=0.
#' @param inferenceMatrix Optionally provide an inference matrix, built with ``fullTreeInference()``. If not provided, ``fullTreeInference()``
#' will be called internally (and is relatively slow).
#' @param sigma_y Provide the standard deviation of y, if known. If not provided, the sample standard deviation of y will be used
#' as a conservative estimate.
#' @param inferenceType An integer specifying which pieces of inference information should be added to the plot. The options
#' currently available are:
#' @param nn boolean- would you like node numbers to be printed?
#' @param extra Passed on to rpart.plot. Most relevant are the following options.
#' \describe{
#' \item{0}{Do not print sample size in each node. }
#' \item{1}{Print sample size in each node. }
#' \item{2}{Label each internal node with a confidence interval and label each split with a p-value.}
#' \item{3}{Slightly less verbose version of "2". Omit the phrase "95% CI" from node labels.}
#' }
#' \describe{
#' \item{0}{No confidence intervals, p-values, or "fitted mean" label. Just calls rpart.plot().}
#' \item{1}{No confidence intervals. Each split labeled with a p-value.}
#' \item{2}{Label each internal node with a confidence interval and label each split with a p-value.}
#' \item{3}{Slightly less verbose version of "2". Omit the phrase "95% CI" from node labels.}
#' }
#' @param digits Integer- how many digits would you like the text in the plot rounded to.
#' @param ... Additional arguements are passed on to rpart.plot(). Examples include "cex".
#' @importFrom intervals interval_union
#' @importFrom intervals interval_complement
#' @importFrom rpart rpart
#' @importFrom stats update
treeval.plot <- function(tree, inferenceMatrix = NULL, sigma_y=NULL,nn=TRUE, printn=1,
                         inferenceType=2, digits=3, ...) {



  if (NROW(tree$frame) != (2*NROW(tree$splits)+1)) {
    tree$call$maxcompete <- 0
    tree$call$maxsurrogate <- 0
    tree <- update(tree)
  }
  if (inferenceType != 0) {
  if (is.null(sigma_y)) {
    y <- tree$model[,1]
    sigma_y <- sd(y)
  }
  if (is.null(inferenceMatrix)) {
    inferenceMatrix <- fullTreeInference(tree, sigma_y)
  }

  inferenceMatrix$pval <- as.numeric(inferenceMatrix$pval)

  tree$splits <- cbind(tree$splits, inferenceMatrix$pval[-1])
  colnames(tree$splits)[6] <- "pvalue"

  tree$frame <- cbind(tree$frame, NA,NA)
  names(tree$frame)[9:10] <- c("CI","pval")


  rootLower <- inferenceMatrix[1,4]
  rootUpper <- inferenceMatrix[1,5]

  tree$frame[1,]$CI <-paste("(",  round(  rootLower, digits), ", ",
                            round(  rootUpper,digits), ")",
                            sep="")

  indices <- which(inferenceMatrix$pval < 1e-3)
  if (length(indices)>0) {
  inferenceMatrix[-indices,]$pval <- paste("=", round(inferenceMatrix[-indices,]$pval, digits))
  inferenceMatrix[indices,]$pval <- "<1e-3"
  } else {
    inferenceMatrix$pval <- paste("=", round(inferenceMatrix$pval, digits))
  }

  ### Anna note to self. Now that you have figured out node-numbering-- there really must
  ### be a better way to do this!!1 Work in progress, but we don't need to guess and check anymore.
  for (row in 2:NROW(tree$frame)) {
    y <- tree$frame[row,]$yval
    child1 <- which(abs(inferenceMatrix$child1mean-y) < 1e-6)
    child2 <- which(abs(inferenceMatrix$child2mean-y) < 1e-6)

    if (length(child1) > 0) {
      tree$frame[row,]$CI <- paste("(", round(inferenceMatrix[child1,]$child1lower,3), ", ", round(inferenceMatrix[child1,]$child1upper,3), ")",
                                   sep="")
      tree$frame[row,]$pval <- inferenceMatrix[child1,]$pval
    }
    if (length(child2) > 0) {
      tree$frame[row,]$CI <- paste("(", round(inferenceMatrix[child2,]$child2lower,4), ", ", round(inferenceMatrix[child2,]$child2upper,4), ")",
                                   sep="")
      tree$frame[row,]$pval <- inferenceMatrix[child2,]$pval
    }
  }
  inner.plot(tree,inferenceType, roundint=FALSE, nn=nn, extra=(printn==TRUE),digits=digits,...)
  } else{
    rpart.plot::rpart.plot(tree, nn=nn, extra=(printn==TRUE), roundint=FALSE, digits=digits, ...)
  }
}
