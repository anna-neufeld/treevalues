#' This function allows users to plot their rpart object with all splits labeled with p-values and all nodes
#' labeled with confidence intervals.
#'
#' @export
#'
#' @param tree An rpart tree. Note that the tree must have been built with arguements maxcompete=0 and maxsurrogates=0.
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
treeval.plot <- function(tree, sigma_y=NULL,nn=TRUE, printn=1,
                         inferenceType=2, digits=3,alpha=0.05, permute=FALSE, ...) {



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
    if (names(tree$frame)[length(names(tree$frame))] != "pval") {
      if (inferenceType > 1) {
        tree <- inferenceFrame(tree, sigma_y=sigma_y, CI=TRUE, alpha=alpha, digits=digits, permute=permute)
      } else {
        tree <- inferenceFrame(tree, sigma_y=sigma_y, CI=FALSE, digits=digits)
      }

    }
    inner.plot(tree,inferenceType, roundint=FALSE, nn=nn, extra=(printn==TRUE),digits=digits,...)
  } else{
    rpart.plot::rpart.plot(tree, nn=nn, extra=(printn==TRUE), roundint=FALSE, digits=digits, ...)
  }
}




inferenceFrame <- function(tree, sigma_y = sd(tree$model[,1]), CI=TRUE, alpha=0.05,digits=3,
                           permute=FALSE) {

  if (is.null(tree$model)) {
    stop('Must build rpart object with parameter model=TRUE')
  }
  if (NROW(tree$frame) != (2*NROW(tree$splits)+1)) {
    tree$call$maxcompete <- 0
    tree$call$maxsurrogate <- 0
    tree <- update(tree)
  }

  dat <- tree$model
  X <- dat[,-1]
  y <- dat[,1]
  p <- NCOL(X)
  n <- NROW(X)

  allBranches <- getAllBranches(tree)

  tree$frame <- cbind(tree$frame, NA,NA)
  names(tree$frame)[9:10] <- c("CI","pval")

  ### FILL IN ROOT NODE
  fullMean <- mean(y)
  zstar <- qnorm(1-alpha/2)
  tree$frame[1,]$CI <-paste("(",  round(  fullMean - zstar*sigma_y/sqrt(n), digits), ", ",
                            round(  fullMean + zstar*sigma_y/sqrt(n),digits), ")",
                            sep="")
  if (NROW(tree$frame) > 1) {
    rootnodepval <- branchInference(tree, allBranches[["2"]], type="sib", computeCI=FALSE)$pval
    child1CI <- branchInference(tree, allBranches[["2"]], type="reg", computeCI=TRUE, alpha=alpha, permute=permute)$confint
    child2CI <- branchInference(tree, allBranches[["3"]], type="reg", computeCI=TRUE,alpha=alpha, permute=permute)$confint
    tree$frame["2",9] <- paste("(", round(child1CI[1],digits), ", ", round(child2CI[2],digits), ")", sep="")
    tree$frame["3",9] <- paste("(", round(child2CI[1],digits), ", ", round(child2CI[2],digits), ")", sep="")
    tree$frame["2",10] <- ifelse(rootnodepval < 1e-3, "<1e-3", paste(" = ", round(rootnodepval, digits)))
    tree$frame["3",10] <- ifelse(rootnodepval < 1e-3, "<1e-3", paste(" = ", round(rootnodepval, digits)))
  }

  for (i in 2:NROW(tree$frame)) {
    region <- row.names(tree$frame)[i]
    if (tree$frame[i,1] != "<leaf>") {
      child1 <- as.character(as.numeric(region)*2+1)
      child2 <- as.character(as.numeric(region)*2)
      splitpval <- branchInference(tree, allBranches[[child1]], type="sib", computeCI=FALSE)$pval
      child1CI <- branchInference(tree, allBranches[[child1]], type="reg", computeCI=TRUE, alpha=alpha, permute=permute)$confint
      child2CI <- branchInference(tree, allBranches[[child2]], type="reg", computeCI=TRUE,alpha=alpha, permute=permute)$confint

      if (splitpval < 1e-3) {
        tree$frame[child1,10] <-  "<1e-3"
        tree$frame[child2,10] <-  "<1e-3"
      } else {
        tree$frame[child1,10] <-  paste("=", round(splitpval, digits), sep=" ")
        tree$frame[child2,10] <-  paste("=", round(splitpval, digits), sep=" ")

      }
      tree$frame[child1,9] <- paste("(", round(child1CI[1],digits), ", ", round(child1CI[2],digits), ")", sep="")
      tree$frame[child2,9] <- paste("(", round(child2CI[1],digits), ", ", round(child2CI[2],digits), ")", sep="")
    }
 }
  return(tree)
}
