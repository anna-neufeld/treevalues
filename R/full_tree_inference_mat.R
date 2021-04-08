#' Constructs a large matrix storing a pvalue for every split and a CI for every node.
#' Used as a precursor to plotting.
#'
#' @param tree An rpart object. Must have been built with rpart arguement "model=TRUE".
#' @param sigma_y The known error variance. If not provided, is estimated with a conservative guess.
#'
#' @return A large matrix storing lots of pvalues and confidence intervals.
#' @importFrom stats sd
#' @importFrom stats var
fullTreeInference <- function(tree, sigma_y =
                                sd(tree$model[,1])) {

  if (is.null(tree$model)) {
    stop('Must build rpart object with parameter model=TRUE')
  }

  ## This is unfortunately necessary to make sure I can read the splits correctly
  ## out of the R output
  dat <- tree$model
  X <- dat[,-1]
  y <- dat[,1]
  p <- NCOL(X)
  n <- NROW(X)

  terminalNodes <- sort(unique(tree$where))
  splitResults <- data.frame(split = NA, pval = NA,
                             effectSize = NA, lower = NA, upper=NA,
                             child1mean=NA,  child1lower=NA, child1upper=NA,
                             child2mean=NA,child2lower=NA, child2upper=NA)

  fullMean <- mean(y)
  splitResults[1,] <- c(" ", NA, NA, fullMean - 1.96*sigma_y/sqrt(n),
                        fullMean + 1.96*sigma_y/sqrt(n), NA, NA,NA,NA, NA, NA)

  j=2
  if (length(terminalNodes) > 1) {
    splitList <- getAllBranches(tree)
    i=1
    while (i < length(splitList)) {
      #print(i)
      splits <- splitList[i][[1]]
      splitText1 <- paste(splitList[i][[1]], collapse=" & ")
      splitText2 <- paste(splitList[i+1][[1]], collapse=" & ")
      i <- i+2
      node1 <- eval(parse(text =splitText1))
      node2 <- eval(parse(text =splitText2))*2
      where <- node1+node2
      y1 <- y[where==1]
      y2 <- y[where==2]
      nu <- (where==1)/sum(where==1) - (where==2)/sum(where==2)

      sample_signal <- t(nu)%*%y
      phi_bounds_split <- getInterval_full(tree, nu,splits)
      p_split <- correctPVal(phi_bounds_split, nu, y, sigma_y)
      CI_split <- computeCI(nu,y,sigma_y, phi_bounds_split, 0.05)

      nu1 <- (where==1)/sum(where==1)
      nu2 <- (where==2)/sum(where==2)

      phiBounds1 <- getInterval_full(tree, nu1, splits)
      phiBounds2 <- getInterval_full(tree, nu2, splits)

      CI1 <- computeCI(nu1,y,sigma_y, phiBounds1, 0.05)
      CI2 <- computeCI(nu2,y,sigma_y, phiBounds2, 0.05)


      splitResults[j,] <- c(splitText1,p_split, mean(y1)-mean(y2),
                            CI_split, mean(y1), CI1 , mean(y2),
                            CI2)
      j <- j+1
    }
    splitResults[,2:11] <- apply(splitResults[,2:11],2,as.numeric)
  }
    return(splitResults)
  }
