getDesc <- function(tree, nn) {
  allNodes <- as.numeric(row.names(tree$frame))
  if (nn %in% allNodes) return(c(nn, getDesc(tree,2*nn), getDesc(tree,2*nn+1)))
  return(c())
}


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


bup_branch <- function(tree,lambda,nn, y=NULL) {
  if (is.null(y)) {
    y <- tree$model[,1]
  }
  regions <- sort(as.numeric(row.names(tree$frame)), decreasing = TRUE)
  BRANCH <- getAnc(tree,nn)
  regions <- c(regions[(regions %in% BRANCH)==FALSE], BRANCH)

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
