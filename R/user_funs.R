#' This is the main function for carrying out inference on an ``rpart`` regression tree. This function can be used to
#' carry out the "inference on a pair of sibling regions" or the "efficient alternative inference on a single region" described in our paper.
#'
#' @param tree An ``rpart`` object corresponding to the tree that you wish to do inference on. This tree must have been built
#' with ``rpart`` parameters ``model=TRUE``, ``maxcompete=0``,``maxsurrogate=``.
#' @param branch A vector of splits describing the branch that you wish to do inference on. You should obtain this using the function
#' ``getBranch()``. Must actually correspond to a branch in ``tree``; otherwise, errors will occur.
#' @param type A string that should be set to either "reg" (default) or "sib". This specifies whether you are doing
#' inference on the mean of single region defined by the end of this branch ("reg"), or doing inference on the difference
#' between this region and its subling.
#' @param alpha The significance level for the confidence interval
#' @param sigma_y The true error variance. If known, this should be passed in. Otherwise, the sample variance will be computed
#' as a conservative estimate.
#' @param c The p-value c
#' @param computeCI Would you like a confidence interval to be computed? Confidence intervals are much slower to compute than
#' p-values, and so if you are performing simulations it may be wise to set this to false.
#' @param permute If you are doing type="reg", do you want to consider all possible permutations of the branch? This leads to the highest power
#' for region inference, but is computationally inefficient. If you are doing type="sib", there is no distinction.
#' @return An object of class ``branch_inference`` that contains a confidence interval, a p-value,
#' the sample statistic, the conditioning set, and a flag reminding the user if type="reg" or type="sib".
#' @export
#' @importFrom intervals Intervals
#' @importFrom intervals size
branchInference <- function(tree, branch, type="reg", alpha=0.05,sigma_y=NULL,c=0, computeCI=TRUE,
                            permute=FALSE) {

  if (c!=0) {stop("At this time, this package only supports testing hypotheses of the form param = 0")}
  dat <- tree$model
  n <- NROW(dat)
  X <- dat[,-1]
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


    #####
    if (!permute) {
      phiBounds <- getInterval(tree, nu, branch)
      if (sum(size(phiBounds))==0) {
        pval <- 1
        CI <- Intervals(c(-Inf,Inf))
      } else {
        if (computeCI) {
          CI <- computeCI(nu,y,sigma_y, phiBounds, alpha)
        } else {CI=NA}
        if (c!=0) {
          print("error: not yet implemented")
        }
        pval <- correctPVal(phiBounds,nu,y,sigma_y)
      }
    } else {
      list1 <- combinat::permn(branch)
      branch2 <- paste("dat$", branch)
      splitText <- paste(branch2, collapse=" & ")
      node1 <- eval(parse(text =splitText))
      nu <- (node1==1)/sum(node1==1)

      ## START US OFF
      phiBounds <- getInterval(tree, nu, branch,sib=FALSE)
      if (length(branch) > 1) {
        list1 <- list1[-1]

        while (length(list1) > 0) {
          perm = list1[[1]]
          sgrow <- getInterval(tree,nu,perm,sib=FALSE, grow=TRUE)
          #print(sgrow)
          if (sum(size(sgrow) > 0)) {
            if (sgrow[1,1] == -Inf) {phiSample <- sgrow[1,2]-5
            } else {
              if (sgrow[1,2] == Inf) {phiSample <- sgrow[1,1]+5
              } else {phiSample <- stats::runif(1,sgrow[1,1],sgrow[1,2])}
            }
            Pi_perp <- diag(rep(1,n)) - nu%*%t(nu)/sum(nu^2)

            yphi <- Pi_perp%*%y + nu/sum(nu^2)*phiSample
            dat2 <- data.frame(cbind(yphi,X))
            names(dat2) = names(dat)
            lambda <- tree$control$cp * (sum((y-mean(y))^2))
            newcontrol <- tree$control
            newcontrol$cp <- 0

            tree0 <-  rpart::rpart(tree$call$formula, data=dat2, control=newcontrol, model=TRUE)
            if (lambda > 0) {
              nodenew <- row.names(tree0$frame)[tree0$where[nu != 0][1]]
              possibleNodes <- getAnc(tree0, nodenew)
              nn <- possibleNodes[length(possibleNodes)-length(branch)]
              treeKL <- tree_KL(tree0,lambda,nn,y=yphi)
            } else {treeKL <- tree0}

            sprune <- getInterval(treeKL,nu,perm,sib=FALSE,prune=TRUE)
            theseBounds <- suppressWarnings(interval_intersection(sgrow, sprune))
            phiBounds <- interval_union(phiBounds, theseBounds)
            list1 <- list1[-1]
          } else {list1 <- list1[-1]}
        }

        #phiBounds <-reduce(phiBounds)
      }

      ### FINALLY WE ARE ALL DONE GETTING THE GIANT CONDITIONING SET
      ### and we can simply do inference!! Yay.
      sample_signal <- t(nu)%*%y
      if (computeCI) {
        CI <- computeCI(nu,y,sigma_y, phiBounds, alpha)
      } else {CI=NA}
      if (c!=0) {
        print("error: not yet implemented")
      }
      pval <- correctPVal(phiBounds,nu,y,sigma_y)
    }
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
    phiBounds <- getInterval(tree, nu,branch,sib=FALSE)
    if (sum(size(phiBounds))==0) {
      pval <- 1
      CI <- Intervals(c(-Inf,Inf))
    } else {
    pval <- correctPVal(phiBounds, nu, y, sigma_y)
    if (computeCI) {
      CI <- computeCI(nu,y,sigma_y, phiBounds, alpha)
    } else {CI=NA}
    }
  }

  out <- list(
    confint = CI, pval = pval, samplemean = t(nu)%*%y, condset = phiBounds, type=type,
    branch=branch,c=c,alpha=alpha
  )
  class(out) <- "branch_inference"
  return(out)
}

#' This function takes an rpart tree and
#'
#' @param tree An rpart object. Must have been built with model=TRUE
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

#' Pass in an rpart tree and also a node number. The nodes are numbered strangely in rpart.
#' If no node numbers are passed in, a list of all branches in the tree will be returned.
#' Pass in a tree and a node number (node number like the type that comes up in rpart.plot!!!)
#'
#' @param tree An rpart object.
#' @param nn A node number that corresponds to a valid node in ``tree``. The list of valid node numbers can be obtained with
#' ``row.names(tree$frame)`` or by plotting ``tree`` with ``treeval.plot`` and ``nn=TRUE``. The node number can be passed in as
#' either a character string or an integer.
#' @return Either a single branch (which is a vector of splits) or (if nn=NULL), a list of all branches in the tree.
#' @export
getBranch <- function(tree, nn=NULL) {
  branches <- getAllBranches(tree)
  if (is.null(nn)) {return(branches)}
  if (as.numeric(nn)) {nn <- as.character(nn)}
  return(branches[[nn]])
}


#' Pass in a tree and a node number (node number like the type that comes up in rpart.plot!!!)
#'
#' @param tree An rpart object. Must have been built with model=TRUE
#' @param nn A node number!! As a string or as an integer is fine
#' @return The indices of data that belong to this region in the training set. The training set is stored in tree$model.
#' @export
#' @importFrom rpart path.rpart
getRegion <- function(tree, nn){
  rule <- path.rpart(tree, nn,print.it=FALSE)
  dat <- tree$model
  rule_2 <- sapply(rule[[1]][-1], function(x) strsplit(x, '(?<=[><=])(?=[^><=])|(?<=[^><=])(?=[><=])', perl = TRUE))
  ind <- apply(do.call(cbind, lapply(rule_2, function(x) eval(call(x[2], dat[,x[1]], as.numeric(x[3]))))), 1, all)
  return(ind)
}





