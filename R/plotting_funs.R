#' This function allows users to plot their rpart object with all splits labeled with p-values and all nodes
#' labeled with confidence intervals.
#'
#' @export
#'
#' @param tree An rpart tree.
#' @param inferenceMatrix The corrsponding inference matrix. This matrix should have been build with the "getFullTreeInference" function.
treeval.plot <- function(tree, inferenceMatrix = NULL, sigma_y=NULL) {
  if (is.null(inferenceMatrix) & is.null(sigma_y)) {
    stop('Must provide either inference matrix or noise variance')
  }
  if (is.null(inferenceMatrix)) {
    inferenceMatrix = fullTreeInference(tree, sigma_y)
  }

  inferenceMatrix$pval <- as.numeric(inferenceMatrix$pval)

  tree$splits <- cbind(tree$splits, inferenceMatrix$pval[-1])
  colnames(tree$splits)[6] <- "pvalue"

  tree$frame <- cbind(tree$frame, NA,NA)
  names(tree$frame)[9:10] <- c("CI","pval")

  #fullMean <- tree$frame[1,]$yval
  #fulln <- tree$frame[1,]$n

  rootLower <- inferenceMatrix[1,4]
  rootUpper <- inferenceMatrix[1,5]

  tree$frame[1,]$CI <-paste("(",  round(  rootLower, 4), ", ",
                            round(  rootUpper,4), ")",
                            sep="")

  indices <- which(inferenceMatrix$pval < 1e-6)
  if (length(indices)>0) {
  inferenceMatrix[-indices,]$pval <- paste("=", round(inferenceMatrix[-indices,]$pval, 6))
  inferenceMatrix[indices,]$pval <- "<1e-6"
  } else {
    inferenceMatrix$pval <- paste("=", round(inferenceMatrix$pval, 6))
  }

  for (row in 2:NROW(tree$frame)) {
    y <- tree$frame[row,]$yval
    child1 <- which(abs(inferenceMatrix$branch1mean-y) < 1e-6)
    child2 <- which(abs(inferenceMatrix$branch2mean-y) < 1e-6)

    if (length(child1) > 0) {
      tree$frame[row,]$CI <- paste("(", round(inferenceMatrix[child1,]$branch1lower,4), ", ", round(inferenceMatrix[child1,]$branch1upper,4), ")",
                                   sep="")
      tree$frame[row,]$pval <- inferenceMatrix[child1,]$pval
    }
    if (length(child2) > 0) {
      tree$frame[row,]$CI <- paste("(", round(inferenceMatrix[child2,]$branch2lower,4), ", ", round(inferenceMatrix[child2,]$branch2upper,4), ")",
                                   sep="")
      tree$frame[row,]$pval <- inferenceMatrix[child2,]$pval
    }
  }
  inner.plot(tree,roundint=FALSE)
}
