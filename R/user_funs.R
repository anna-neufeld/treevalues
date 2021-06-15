#' Main function for carrying out inference on an \code{rpart} regression tree.
#'
#' Can be used to carry out inference on a pair of sibling regions or on a single region.
#'
#' @param tree An \code{rpart} object corresponding to the tree that you wish to do inference on. This tree must have been built
#' with \code{rpart} parameter \code{model=TRUE}.
#' @param branch A vector of splits describing the branch that you wish to do inference on. You should obtain this using the function
#' `getBranch()`. Must actually correspond to a branch in \code{tree}.
#' @param type A string that should be set to either \code{"reg"} (default) or \code{"sib"}. This specifies whether you are doing
#' inference on the mean of single region at the end of this branch ("reg"), or doing inference on the difference
#' between this region and its sibling.
#' @param alpha Function will compute an equi-tailed \code{(1-alpha)} confidence interval. Default is \code{0.05}.
#' @param sigma_y The true standard deviation of the response. If known, this should be passed in. Otherwise, the sample standard deviation
#' will be used as a conservative estimate.
#' @param c The p-value returned will test the null hypothesis that the parameter of interest is equal to c. Currently, only c=0 is a valid
#' input.
#' @param computeCI Boolean that specifies if you would like a confidence interval to be computed. Confidence intervals are much slower to compute than
#' p-values, and so if you only wish to see a p-value you may want to set this to be false.
#' @param permute Boolean. Only relevant if type="reg". If \code{FALSE} (default), inference is conducted conditional on the event that this exact branch
#' appeared in the tree. If \code{TRUE}, inference is conducted conditional on the event that some permutation of this branch appeared in the tree. While the
#' latter achieves higher power, it can be computationally prohibitive in large trees.
#' @return An object of class \code{branch_inference} that contains a confidence interval, a p-value,
#' the sample statistic, the conditioning set, and a flag reminding the user if \code{type="reg"} or \code{type="sib"}.
#' @export
#' @importFrom intervals Intervals
#' @importFrom intervals size
#' @examples
#' bls.tree <- rpart::rpart(kcal24h0~hunger+disinhibition+resteating+rrvfood+liking+wanting,
#'     model = TRUE, data = blsdata, cp=0.02)
#' branch <- getBranch(bls.tree, 8)
#' branchInference (bls.tree, branch, type="sib")
#' branchInference (bls.tree, branch, type="reg", permute=TRUE)
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
    confint = CI, pval = pval, samplemean = as.numeric(t(nu)%*%y), condset = phiBounds, type=type,
    branch=branch,c=c,alpha=alpha
  )
  class(out) <- "branch_inference"
  return(out)
}

#' Given an rpart tree, returns a list describing every branch in the tree. Useful for extracting individual branches, which are necessary inputs to \code{branchInference()}.
#'
#' @param tree An rpart object. Must have been built with model=TRUE
#' @return A list. Each entry corresponds to a branch in the tree. The entries are named with the node-number of the node that falls at the end of the branch.
#' Each entry is a vector of strings.
#' @keywords internal
#' @noRd
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

#' Get branch or list of branches from tree.
#'
#' Given an rpart tree and a node number, returns a vector of strings that describes the branch which defines the node.
#' If no node number is provided, returns a list describing every branch in the tree.
#' Useful for extracting individual branches, which are necessary inputs to \code{branchInference()}.
#'
#' @param tree An rpart object.
#' @param nn A node number that corresponds to a valid node in \code{tree}. The list of valid node numbers can be obtained with
#' \code{row.names(tree$frame)} or by plotting \code{tree} with \code{treeval.plot()}. The node number can be passed in as
#' either a character string or an integer. If no node number is provided, a list of all branches in the tree will be returned.
#' @return Either a single branch (which is a vector of splits) or (if nn=NULL), a list of all branches in the tree.
#' @export
#' @examples
#' bls.tree <- rpart::rpart(kcal24h0~hunger+disinhibition+resteating+rrvfood+liking+wanting,
#'     model = TRUE, data = blsdata, cp=0.02)
#' branch <- getBranch(bls.tree, 8)
#' branchInference (bls.tree, branch, type="sib")
getBranch <- function(tree, nn=NULL) {
  branches <- getAllBranches(tree)
  if (is.null(nn)) {return(branches)}
  if (as.numeric(nn)) {nn <- as.character(nn)}
  return(branches[[nn]])
}


#' Get observations belonging to a region in \code{tree}.
#'
#' Pass in a tree and a node number. This returns a vector of booleans identifying which members of the training set
#' belong to the given region.
#' Mainly used to form vectors nu to define parameters from a given tree. Since this is called internally by \code{branchInference()},
#' will rarely be needed directly by users.
#'
#' @param tree An rpart object. Must have been built with model=TRUE
#' @param nn A node number. Can be a string or an integer.
#' @return The indices of data that belong to this region in the training set. The training set is stored in tree$model.
#' @export
#' @importFrom rpart path.rpart
#' @examples
#'data(blsdata, package="treevalues")
#' bls.tree <-rpart::rpart(kcal24h0~hunger+disinhibition+resteating+rrvfood+liking+wanting,
#'     model = TRUE,  data = blsdata, cp=0.02)
#' branch <- getBranch(bls.tree, 2)
#' left_child <- getRegion(bls.tree,2)
#' right_child <- getRegion(bls.tree,3)
#' nu_sib <- left_child/sum(left_child) -  right_child/sum(right_child)
#' S_sib <- getInterval(bls.tree, nu_sib,branch)
#' correctPVal(S_sib, nu_sib, blsdata$kcal24h0, sd(blsdata$kcal24h0))
getRegion <- function(tree, nn){
  rule <- path.rpart(tree, nn,print.it=FALSE)
  dat <- tree$model
  rule_2 <- sapply(rule[[1]][-1], function(x) strsplit(x, '(?<=[><=])(?=[^><=])|(?<=[^><=])(?=[><=])', perl = TRUE))
  ind <- apply(do.call(cbind, lapply(rule_2, function(x) eval(call(x[2], dat[,x[1]], as.numeric(x[3]))))), 1, all)
  return(ind)
}





