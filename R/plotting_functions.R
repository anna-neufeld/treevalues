#' Plot an rpart object with all splits labeled with p-values and all nodes
#' labeled with confidence intervals.
#'
#' Essentially a wrapper function for \code{rpart.plot()} from the \code{rpart.plot} package, but
#' additional arguments allow user to add p-values and confidence intervals to plots.
#'
#' @export
#'
#' @param tree An rpart tree. The tree must have been build with parameter \code{model=TRUE}.
#' @param sigma_y Provide the standard deviation of y, if known. If not provided, the sample standard deviation of y will be used
#' as a conservative estimate.
#' @param nn boolean- would you like node numbers to be printed? Nodes are numbered using the same methodology as the \code{rpart} package. If node \code{n} has children,
#' its children are numered \code{2n} and \code{2n+1}.
#' @param printn boolean - would you like the number of observations to be printed in each node?
#' @param inferenceType An integer specifying which pieces of inference information should be added to the plot. The options
#' currently available are
#' (0) No confidence intervals, p-values, or "fitted mean" label. Just calls \code{rpart.plot()}.
#' (1) No confidence intervals. Each split labeled with a p-value.
#' (2)  Label each internal node with a confidence interval and label each split with a p-value. This is the default, but
#' can also be a little messy/hard to read. Options 3 and 4 print the same information but with small
#' formatting tweaks.
#' @param digits Integer- how many digits would you like the text in the plot rounded to.
#' @param alpha If inferenceType is such that confidence intervals will be printed, \code{(1-alpha)} confidence intervals will be printed.
#' @param permute If inferenceType is such that confidence intervals will be printed, should the conditioning set for the confidence intervals include
#' all permutations of the relevant branch? Setting this to TRUE will lead to slightly narrower confidence intervals, but will make computations more expensive.
#' See paper for more details.
#' @param ... Additional arguments are passed on to rpart.plot(). Examples include \code{"cex"}.
#' @importFrom intervals interval_union
#' @importFrom intervals interval_complement
#' @importFrom rpart rpart
#' @importFrom stats update
#' @importFrom stats qnorm
#' @importFrom stats sd
#' @examples
#' bls.tree <-rpart::rpart(kcal24h0~hunger+disinhibition+resteating,
#'     model = TRUE, data = blsdata, maxdepth=1)
#' treeval.plot(bls.tree, inferenceType=0)
#' treeval.plot(bls.tree, inferenceType=1)
#' treeval.plot(bls.tree, inferenceType=2)
treeval.plot <- function(tree, sigma_y=NULL,nn=TRUE, printn=TRUE,
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
    inner.plot(tree,inferenceType=0, roundint=FALSE, nn=nn, extra=(printn==TRUE),digits=digits,...)
  }
}


#' This function can optionally be used prior to running \code{treeval.plot} to make \code{treeval.plot} run more efficiently.
#'
#' This function is computationally expensive, especially if \code{CI=TRUE} and/or \code{permute=TRUE}. This function is called internally by \code{treeval.plot()},
#' as it updates \code{tree$frame} to store information (pvalues and confidence intervals) that will be printed by \code{treeval.plot()}. If you will be
#' making several plots while playing around with font size and formatting, it is a good idea to call this function first so that it need not be called
#' repeatedly by different calls of treeval.plot
#'
#' @param tree The tree that you will be plotting.
#' @param sigma_y The standard deviation of the response. If known, should be provided. Otherwise, a conservative estimate (the sample
#' standard deviation of the response) is used.
#' @param CI Boolean. Should confidence intervals be computed? As confidence intervals are inefficient to compute, this should be set to
#' \code{FALSE} if you intend to make a plot that does not display confidence intervals.
#' @param alpha If \code{CI=TRUE}, the confidence intervals that are computed will be \code{(1-alpha)} confidence intervals.
#' @param digits Integer. The number of digits that the p-values and confidence intervals will be rounded to in the later plot.
#' @param permute If \code{CI=TRUE}, this boolean says whether or not the
#' @return An rpart object. Identical to \code{tree} expect that now \code{tree$frame} has two extra columns; one storing p-values for splits and the other
#' storing confidence intervals for regions. If this object is passed in to \code{treeval.plot}, the plots will be made more efficiently.
#' @importFrom stats qnorm
#' @importFrom stats sd
#' @export
#' @examples
#' \dontrun{
#' library(rpart)
#' bls.tree <-rpart(
#'   kcal24h0~hunger+disinhibition+resteating+rrvfood+liking+wanting,
#'   model = TRUE, data = blsdata, cp=0.02
#' )
#' bls.tree2 <- inferenceFrame(bls.tree)
#' treeval.plot(bls.tree2, inferenceType=1)
#' treeval.plot(bls.tree2, inferenceType=2)
#' treeval.plot(bls.tree2, inferenceType=2, nn=FALSE)
#' }
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
    tree$frame["2",9] <- paste("(", round(child1CI[1],digits), ", ", round(child1CI[2],digits), ")", sep="")
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
