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
                             branch1mean=NA,  branch1lower=NA, branch1upper=NA,
                             branch2mean=NA,branch2lower=NA, branch2upper=NA)

  fullMean <- mean(y)
  splitResults[1,] <- c(" ", NA, NA, fullMean - 1.96*sigma_y/sqrt(n),
                        fullMean + 1.96*sigma_y/sqrt(n), NA, NA,NA,NA, NA, NA)

  j=2
  if (length(terminalNodes) > 1) {
    splitList <- getAllTests_rpart(tree)
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

getAllTests_rpart <- function(base_tree) {
  allNodes <- sort(as.numeric(as.character(unique(rpart.utils::rpart.rules.table(base_tree)$Rule)[-1])))
  ruleTab <- rpart.utils::rpart.rules.table(base_tree)[-1,]
  ruleTab$Subrule <- factor(ruleTab$Subrule) #levels=levels(ruleTab$Subrule)[-1])
  #ruleTab$Subrule <- factor(ruleTab$Subrule, levels=levels(ruleTab$Subrule)[-1])
  subruleTab <-  rpart.utils::rpart.subrules.table(base_tree)

  allSplits <- list()
  index <- 1
  for (node in allNodes) {
    rules <- ruleTab[ruleTab$Rule==node,]
    subRules <- c()
    for (i in 1:NROW(rules)) {
      rule=rules[i,]
      subRules <- rbind(subRules,
                        subruleTab[subruleTab$Subrule == rule$Subrule,])
    }



    subRules$Variable <- as.factor(subRules$Variable)
    splits <- c()
    if (NROW(subRules)==1) {
      split <- subRules
      if (is.na(split[4])) {
        splits <- c(splits, paste("dat$", split$Variable, " >=", split$Greater))
      } else {
        splits <- c(splits, paste("dat$", split$Variable, " <", split$Less))
      }
    } else{
      for (j in 1:NROW(subRules)) {
        split <- subRules[j,]
        if (is.na(split[4])) {
          splits <- c(splits, paste("dat$", split$Variable, " >=", split$Greater))
        } else {
          splits <- c(splits, paste("dat$", split$Variable, " <", split$Less))
        }
      }
      #splits <- paste(splits, collapse=" & ")
    }

    allSplits[[index]] <- splits
    index <- index+1
  }
  return(allSplits)
}


















