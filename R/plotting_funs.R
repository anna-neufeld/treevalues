treeval.plot <- function(tree, inferenceMatrix) {
  tree$frame$yval <- round(tree$frame$yval, 6)
  inferenceMatrix$pval <- as.numeric(inferenceMatrix$pval)

  tree$splits <- cbind(tree$splits, inferenceMatrix$pval)
  colnames(tree$splits)[6] <- "pvalue"

  tree$frame <- cbind(tree$frame, NA,NA)
  names(tree$frame)[9:10] <- c("CI","pval")

  fullMean <- tree$frame[1,]$yval
  fulln <- tree$frame[1,]$n
  tree$frame[1,]$CI <-paste("(",  round(fullMean - 1.96*sigma_y/sqrt(fulln), 4), ", ",
                            round(fullMean + 1.96*sigma_y/sqrt(fulln),4), ")",
                            sep="")


  for (row in 2:NROW(tree$frame)) {
    y <- tree$frame[row,]$yval
    child1 <- which(inferenceMatrix$branch1mean==y)
    child2 <- which(inferenceMatrix$branch2mean==y)

    if (length(child1) > 0) {
      tree$frame[row,]$CI <- paste("(", inferenceMatrix[child1,]$branch1lower, ", ", inferenceMatrix[child1,]$branch1upper, ")",
                                   sep="")
      tree$frame[row,]$pval <- inferenceMatrix[child1,]$pval
    }
    if (length(child2) > 0) {
      tree$frame[row,]$CI <- paste("(", inferenceMatrix[child2,]$branch2lower, ", ", inferenceMatrix[child2,]$branch2upper, ")",
                                   sep="")
      tree$frame[row,]$pval <- inferenceMatrix[child2,]$pval
    }
  }
  inner.plot(tree,roundint=FALSE)
}
