#' #' Perform inference on a single region in a tree.
#' #' @param tree An rpart object containing the region that you want to do inference on.
#' #' @param nn The node number for the region that you want to do inference on.
#' #' @param permutations A boolean. If TRUE, inference will be performed conditional on the event that the
#' #' region was in the tree. This leads to higher power, but is computationally expensive. If FALSE, inference will
#' #' condition on both the region and its ancestors being in the tree. This is slightly less powerful but much more efficient.
#' #' @param alpha significance level for confidence intervals.
#' #' @param sigma_y The standard deviation of the response. If not provided, we will use the sample sd, which is conservative.
#' #' @param c Defaults to 0. The p-value provided tests the null that the mean in the region is equal to c.
#' #' @param computeCI Confidence intervals are more computaitonally difficult than p-values. Can set to false if you only want a p-value and you want the code
#' #' to run faster.
#' regInference <- function(tree, nn, permutations=TRUE, alpha=0.05,sigma_y=NULL,c=0, computeCI=TRUE) {
#'   dat <- tree$model
#'   y <- dat[,1]
#'   if (is.null(sigma_y)) {
#'     sigma_y <- sd(dat[,1])
#'   }
#'
#'   branch <- getBranch(tree, nn)
#'   list1 <- combinat::permn(branch)
#'   branch2 <- paste("dat$", branch)
#'   splitText <- paste(branch2, collapse=" & ")
#'   node1 <- eval(parse(text =splitText))
#'   nu <- (node1==1)/sum(node1==1)
#'
#'   ## START US OFF
#'   phiBounds <- getInterval(tree, nu, branch,sib=FALSE)
#'   if (permutations==TRUE & length(branch) > 1) {
#'   list1 <- list1[-1]
#'
#'   while (length(list1) > 0) {
#'     perm = list1[[1]]
#'     sgrow <- getInterval(tree,nu,perm,sib=FALSE, grow=TRUE)
#'     #print(sgrow)
#'     if (sum(size(sgrow) > 0)) {
#'       if (sgrow[1,1] == -Inf) {phiSample <- sgrow[1,2]-5
#'       } else {
#'         if (sgrow[1,2] == Inf) {phiSample <- sgrow[1,1]+5
#'         } else {phiSample <- runif(1,sgrow[1,1],sgrow[1,2])}
#'       }
#'       Pi_perp <- diag(rep(1,n)) - nu%*%t(nu)/sum(nu^2)
#'
#'       yphi <- Pi_perp%*%y + nu/sum(nu^2)*phiSample
#'       dat2 <- data.frame(cbind(yphi,X))
#'       names(dat2) = names(dat)
#'       lambda <- tree$control$cp * (sum((y-mean(y))^2))
#'       newcontrol <- tree$control
#'       newcontrol$cp <- 0
#'
#'       tree0 <-  rpart::rpart(tree$call$formula, data=dat2, control=newcontrol, model=TRUE)
#'       if (lambda > 0) {
#'       nodenew <- row.names(tree0$frame)[tree0$where[nu != 0][1]]
#'       possibleNodes <- getAnc(tree0, nodenew)
#'       nn <- possibleNodes[length(possibleNodes)-length(branch)]
#'       treeKL <- tree_KL(tree0,lambda,nn,y=yphi)
#'       } else {treeKL <- tree0}
#'
#'       sprune <- getInterval(treeKL,nu,perm,sib=FALSE,prune=TRUE)
#'       theseBounds <- suppressWarnings(interval_intersection(sgrow, sprune))
#'       phiBounds <- interval_union(phiBounds, theseBounds)
#'       list1 <- list1[-1]
#'       } else {list1 <- list1[-1]}
#'   }
#'
#'   #phiBounds <-reduce(phiBounds)
#'   }
#'
#'   ### FINALLY WE ARE ALL DONE GETTING THE GIANT CONDITIONING SET
#'   ### and we can simply do inference!! Yay.
#'   sample_signal <- t(nu)%*%y
#'   if (computeCI) {
#'     CI <- computeCI(nu,y,sigma_y, phiBounds, alpha)
#'   } else {CI=NA}
#'   if (c!=0) {
#'     print("error: not yet implemented")
#'   }
#'   pval <- correctPVal(phiBounds,nu,y,sigma_y)
#'
#'   out <- list(
#'     confint = CI, pval = pval, samplemean = t(nu)%*%y, condset = phiBounds, type="reg",
#'     branch=branch,c=c,alpha=alpha
#'   )
#'   class(out) <- "branch_inference"
#'   return(out)
#' }
