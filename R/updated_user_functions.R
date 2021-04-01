getSreg <- function(tree, node=NULL, nu=NULL, splits=NULL, permutations=FALSE) {
  if (is.null(nu) | is.null(splits)) {
    if (is.null(node)) {stop("Must provide valid node OR provide vector nu and split vector")}
    check <- which(tree$where==node)
    if (length(check)==0) {stop("Must provide valid node OR provide vector nu and split vector")}
    nu <- (tree$where==node)/sum((tree$where==node))
    splits <- getAncestors(tree, node)
  }
  
  if (!permutations) {
    return(getInterval_full(tree,nu,splits))
  } else{
    return(getInterval_permutations(tree, nu,splits))
  }
}

getSbranch <- function(tree, nu, splits) {
    return(getInterval_full(tree,nu,splits))
}

#' 
getSsib <- function(tree, locTest=NULL, nu=NULL, splits=NULL) {
  return(getInterval_full_sibs(tree,nu,splits))
}

getNodeInterval <- function(tree, node, sigma_y=NULL, alpha=0.05, permutations=FALSE) {
  phi_bounds <- getSreg(tree,node,permutations)
  nu <- (tree$where==node)/sum((tree$where==node))
  y <- tree$model[,1]
  if (is.null(sigma_y)) {
    sigma_y <- sum((y-mean(y))^2)
  }
  CI <- computeCI(nu,y,sigma_y, phi_bounds, alpha)
  return(CI)
}

