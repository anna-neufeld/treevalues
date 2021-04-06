### ANNA TO DO:: MAKE THESE WORK FOR ANY ROWS IN FRAME?? NOT JUST TERMINAL??


#' Get a confidence interval for the mean response in a terminal node.
#'
#' @export
#'
#' @param tree An rpart object. The tree must have been built with arguement model=TRUE,
#' @param  node identify the node!!. Need to explain this. Must be terminal node.
#' @param sigma_y enter the known noise SD (for now)
#' @param alpha The significance level. We build a 1-alpha CI.
#' @return a list with 3 things.
#' @importFrom intervals Intervals
#' @importFrom intervals interval_complement
#' @importFrom intervals interval_union
#' @importFrom intervals interval_intersection
#' @importFrom intervals size
nodeInference <- function(tree, node, sigma_y=NULL, alpha=0.05) {
  splits <- getBranch(tree, node)
  nu <- (tree$where==node)/sum((tree$where==node))
  phi_bounds <- getInterval_full(tree,nu, splits)

  y <- tree$model[,1]
  if (is.null(sigma_y)) {
    sigma_y <- sqrt(sum((y-mean(y))^2)/(length(y)-1))
  }
  CI <- computeCI(nu,y,sigma_y, phi_bounds, alpha)
  pval <- correctPVal(phi_bounds, nu, y, sigma_y)
  results <- list(confint = CI, pval = pval, samplemean = t(nu)%*%y, condset = phi_bounds)
  return(results)
}



#' Get inference for a difference in means between two terminal nodes
#'
#' @export
#'
#' @param tree an rpart object. Must have been built with model=TRUE.
#' @param locTest identify the pair of nodes that you are testing for. this should be a vector of length 2. the numbers should correspond to entries in tree$where
#' @param sigma_y enter the known SD of y. If none is entered, the default is to use the total SST of y,
#' which is conservative.
#' @importFrom intervals Intervals
#' @importFrom intervals interval_complement
#' @importFrom intervals interval_union
#' @importFrom intervals interval_intersection
splitInference <- function(tree, locTest, sigma_y = NULL, alpha=0.05) {
  splits <- getBranch(tree, locTest[1])
  nu <- (tree$where==locTest[1])/sum((tree$where==locTest[1])) - (tree$where==locTest[2])/sum(tree$where==locTest[2])
  phi_bounds <- getInterval_full(tree,nu, splits)
  y <- tree$model[,1]
  if (is.null(sigma_y)) {
    sigma_y <- sqrt(sum((y-mean(y))^2)/(length(y)-1))
  }
  pval <- correctPVal(phi_bounds, nu, y, sigma_y)
  CI <- computeCI(nu,y,sigma_y, phi_bounds, alpha)
  results <- list(confint = CI, pval = pval, samplemean = t(nu)%*%y, condset = phi_bounds)
  return(results)
}


#' Obtain the list of splits that refine a certain region.
#' Necessary precursor to calling "getInterval"
#'
#' @param tree An rpart object. Must have been built with model=TRUE
#' @param node A number identifying the region that you want the ancestors of. Numbers correspond to elements in tree$where.
#'
#' @return A vector of strings describing the splits that define the node in the tree.
#' @export
getBranch <- function(tree, node)
{
  try1 <- try({
    object <- partykit::as.party(tree)
    rules <- partykit:::.list.rules.party(object)
    relevantRules <- rules[[as.character(node)]]
    relevantRules <- strsplit(relevantRules, '&')[[1]]
    relevantRules <- sapply(relevantRules, trimws)
    relevantRules <- sapply(relevantRules, function(u) paste0("dat$", u))
  })
  if (class(try1)=="try-error") {
    return(errorCondition("Could not obtain splits from tree. Make sure that the rpart tree
                          was built with model=TRUE parameter."))
  } else {
    return(relevantRules)
  }
}


branchInference <- function(tree, branch, type="reg", alpha=0.05,sigma_y=NULL) {

  dat <- tree$model
  y <- dat[,1]
  if (is.null(sigma_y)) {
    sigma_y <- sd(dat[,1])
  }

  if (type=="reg") {
    splitText <- paste(branch, collapse=" & ")
    node1 <- eval(parse(text =splitText))
    nu <- (node1==1)/sum(node1==1)

    sample_signal <- t(nu)%*%y
    phiBounds <- getInterval_full(tree, nu, branch)
    CI <- computeCI(nu,y,sigma_y, phiBounds, alpha)
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
    phi_bounds <- getInterval_full(tree, nu,branch)
    pval <- correctPVal(phi_bounds, nu, y, sigma_y)
    CI <- computeCI(nu,y,sigma_y, phi_bounds, alpha)
  }


  results <- list(confint = CI, pval = pval, samplemean = t(nu)%*%y, condset = phi_bounds)
  return(results)
}






