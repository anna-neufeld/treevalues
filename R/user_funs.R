#' Pass in a branch!!!
#'
#' @param tree An rpart object. Must have been built with model=TRUE
#' @param branch A vector of splits.
#' @param type Pass in either "reg" or "sib"
#' @param alpha The significance level forthe confidence interval
#' @param sigma_y the error variance.Assumed known but if not we'll use the sample variance.
#'
#' @return A fitted object!! With summary info
#' @export
#' @importFrom intervals Intervals
branchInference <- function(tree, branch, type="reg", alpha=0.05,sigma_y=NULL,c=0) {
  dat <- tree$model
  y <- dat[,1]
  if (is.null(sigma_y)) {
    sigma_y <- sd(dat[,1])
  }

  branch <- paste("dat$", branch)

  if (type=="reg") {
    splitText <- paste(branch, collapse=" & ")
    node1 <- eval(parse(text =splitText))
    nu <- (node1==1)/sum(node1==1)

    sample_signal <- t(nu)%*%y
    phiBounds <- getInterval_full(tree, nu, branch)
    CI <- computeCI(nu,y,sigma_y, phiBounds, alpha)
    if (c!=0) {
      print("error: not yet implemented")
    }
    pval <- correctPVal(phiBounds,nu,y,sigma_y)
  }

  if (type=="sib") {
    branchsib <- branch
    branchsib[length(branchsib)] <- paste("!", branchsib[length(branchsib)])
    splitText1 <- paste(branch, collapse=" & ")
    splitText2 <- paste(branchsib, collapse=" & ")
    node1 <- eval(parse(text =splitText1))
    node2 <- eval(parse(text =splitText2))*2
    where <- node1+node2
    nu <- (where==1)/sum(where==1) - (where==2)/sum(where==2)
    sample_signal <- t(nu)%*%y
    phiBounds <- getInterval_full(tree, nu,branch)
    pval <- correctPVal(phiBounds, nu, y, sigma_y)
    CI <- computeCI(nu,y,sigma_y, phiBounds, alpha)
  }

  out <- list(
    confint = CI, pval = pval, samplemean = t(nu)%*%y, condset = phiBounds, type=type,
    branch=branch,c=c,alpha=alpha
  )
  class(out) <- "branch_inference"
  return(out)
}


getAllBranches <- function(tree) {
  allNodes <- sort(as.numeric(as.character(unique(rpart.utils::rpart.rules.table(tree)$Rule)[-1])))
  ruleTab <- rpart.utils::rpart.rules.table(tree)[-1,]
  ruleTab$Subrule <- factor(ruleTab$Subrule)

  subruleTab <-  rpart.utils::rpart.subrules.table(tree)

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
        splits <- c(splits, paste(split$Variable, " >=", split$Greater))
      } else {
        splits <- c(splits, paste(split$Variable, " <", split$Less))
      }
    } else{
      for (j in 1:NROW(subRules)) {
        split <- subRules[j,]
        if (is.na(split[4])) {
          splits <- c(splits, paste(split$Variable, " >=", split$Greater))
        } else {
          splits <- c(splits, paste(split$Variable, " <", split$Less))
        }
      }
    }

    allSplits[[index]] <- splits
    index <- index+1
  }

  possibleNames <- rownames(tree$frame)[-1]

  names(allSplits) <- as.character(sort(as.numeric(possibleNames)))

  return(allSplits)
}


getBranch <- function(tree, nn) {
  branches <- getAllBranches(tree)
  if (as.numeric(nn)) {nn <- as.character(nn)}
  return(branches$nn)
}





