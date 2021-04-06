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





