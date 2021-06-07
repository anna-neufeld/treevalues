#' Returns vector of descendents in ``tree`` for node ``nn``.
#' @param tree an rpart object.
#' @param nn the number of a node in tree. Note that node numbers correspond to the row names of tree$frame.
#' @return a vector of descendents.
#' @export
getDesc <- function(tree, nn) {
  allNodes <- as.numeric(row.names(tree$frame))
  if (nn %in% allNodes) return(c(nn, getDesc(tree,2*nn), getDesc(tree,2*nn+1)))
  return(c())
}


#' Returns vector of ancestors in ``tree`` for node ``nn``.
#' @param tree an rpart object.
#' @param nn the number of a node in tree. Note that node numbers correspond to the row names of tree$frame.
#' @return a vector of ancestors.
#' @export
getAnc <- function(tree, nn) {
  anc <- c()
  allNodes <- as.numeric(row.names(tree$frame))
  for (node in allNodes) {
    if (nn %in% getDesc(tree,node)) {
      anc <- c(anc, node)
    }
  }
  return(sort(anc, decreasing=TRUE))
}

#' Returns the tree that remains after all nodes in "tree" BESIDES the ancestors of node 'nn' have been visited
#' for bottom-up pruning with parameter lambda.
#' @param tree An rpart object. If y is not provided, must have been built with model=TRUE
#' @param lambda The cost complexity parameter for pruning.
#' @param nn the node number of the region whose ancestors you want to leave in the tree.
#' @param y the response vector to evaluate the pruning on.
#' @return An rpart object; intermediately pruned.
#' @noRd
#' @internal
tree_KL <- function(tree,lambda,nn, y=NULL) {
  if (is.null(y)) {
    y <- tree$model[,1]
  }
  regions <- sort(as.numeric(row.names(tree$frame)), decreasing = TRUE)
  K <- length(regions)
  BRANCH <- getAnc(tree,nn)
  ### AKA JUST DONT VISIT THE ONES IN BRANCH??????
  regions <- c(regions[(regions %in% BRANCH)==FALSE], BRANCH[1])

  terminal <- as.numeric(row.names(tree$frame))[sort(unique(tree$where))]

  #max_non_terminal <- (max(regions)-1)/2
  to_visit <- regions[(regions %in% terminal)==FALSE]


  prev_tree <- tree
  prev_SSE <- sum((y-predict(prev_tree))^2)
  for (region in to_visit) {
    temp_tree <-rpart::snip.rpart(prev_tree, region)
    temp_SSE <- sum((y-predict(temp_tree))^2)
    if ((length(unique(prev_tree$where))-length(unique(temp_tree$where))) > 0) {
      improvement <- (temp_SSE - prev_SSE)/(length(unique(prev_tree$where))-length(unique(temp_tree$where)))
      if (improvement < lambda) {
        prev_tree <- temp_tree
        prev_SSE <- sum((y-predict(prev_tree))^2)
      }
    }
  }
  return(prev_tree)
}
