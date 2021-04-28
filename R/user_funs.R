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
branchInference <- function(tree, branch, type="reg", alpha=0.05,sigma_y=NULL,c=0, computeCI=TRUE) {
  dat <- tree$model
  y <- dat[,1]
  if (is.null(sigma_y)) {
    sigma_y <- sd(dat[,1])
  }

  branch2 <- paste("dat$", branch)

  if (type=="reg") {
    splitText <- paste(branch2, collapse=" & ")
    node1 <- eval(parse(text =splitText))
    nu <- (node1==1)/sum(node1==1)

    sample_signal <- t(nu)%*%y
    phiBounds <- getInterval_full(tree, nu, branch)
    if (computeCI) {
    CI <- computeCI(nu,y,sigma_y, phiBounds, alpha)
    } else {CI=NA}
    if (c!=0) {
      print("error: not yet implemented")
    }
    pval <- correctPVal(phiBounds,nu,y,sigma_y)
  }

  if (type=="sib") {
    branchsib <- branch2
    branchsib[length(branchsib)] <- paste("!", branchsib[length(branchsib)])
    splitText1 <- paste(branch2, collapse=" & ")
    splitText2 <- paste(branchsib, collapse=" & ")
    node1 <- eval(parse(text =splitText1))
    node2 <- eval(parse(text =splitText2))*2
    where <- node1+node2
    nu <- (where==1)/sum(where==1) - (where==2)/sum(where==2)
    sample_signal <- t(nu)%*%y
    ### SWITCH TO SIB=TRUE for computation later!!!
    phiBounds <- getInterval_full(tree, nu,branch,sib=FALSE)
    pval <- correctPVal(phiBounds, nu, y, sigma_y)
    if (computeCI) {
      CI <- computeCI(nu,y,sigma_y, phiBounds, alpha)
    } else {CI=NA}
  }

  out <- list(
    confint = CI, pval = pval, samplemean = t(nu)%*%y, condset = phiBounds, type=type,
    branch=branch,c=c,alpha=alpha
  )
  class(out) <- "branch_inference"
  return(out)
}

#' Get a list of branches in tree to easily find out what you want to do inferenceon
#'
#' @param tree An rpart object. Must have been built with model=TRUE
#' @export
getAllBranches <- function(tree) {
  if (length(unique(tree$where))==1) {return(NA)}
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

#' Pass in a tree and a node number (node number like the type that comes up in rpart.plot!!!)
#'
#' @param tree An rpart object. Must have been built with model=TRUE
#' @param nn A node number!! As a string or as an integer is fine
#' @export
getBranch <- function(tree, nn) {
  branches <- getAllBranches(tree)
  if (as.numeric(nn)) {nn <- as.character(nn)}
  return(branches[[nn]])
}


#' Pass in a tree and a node number (node number like the type that comes up in rpart.plot!!!)
#'
#' @param tree An rpart object. Must have been built with model=TRUE
#' @param nn A node number!! As a string or as an integer is fine
#' @export
getRegion <- function(tree, nn){
  rule <- path.rpart(tree, nn,print.it=FALSE)
  dat <- tree$model
  rule_2 <- sapply(rule[[1]][-1], function(x) strsplit(x, '(?<=[><=])(?=[^><=])|(?<=[^><=])(?=[><=])', perl = TRUE))
  ind <- apply(do.call(cbind, lapply(rule_2, function(x) eval(call(x[2], dat[,x[1]], as.numeric(x[3]))))), 1, all)
  return(ind)
}





