#' Get a confidence interval for the mean response in a terminal node.
#'
#' @export
#'
#' @param tree An rpart object. The tree must have been built with arguement model=TRUE,
#' @param X the X data; as a matrix. names must match those
#' @param y the y data
#' @param  node identify the node!!
#' @param sigma_y enter the known noise SD (for now)
#' @param alpha The significance level. We build a 1-alpha CI.
#' @return a vector of two numbers; the lower and the upper bound of the CI
getNodeInterval <- function(tree, node, sigma_y=NULL, alpha=0.05) {
  splits <- getAncestors(tree, node)
  nu <- (tree$where==node)/sum((tree$where==node))
  phi_bounds <- getInterval_permutation(tree,nu, splits)
  y <- tree$model[,1]
  CI <- computeCI(nu,y,sigma_y, phi_bounds, alpha)
  return(CI)
}

#' Compute the conditioning set for a region
#'
#' @param tree the tree
#' @param nu the contrast vector
#' @param splits the original set of splits describing the region.
#'
#' @return An object of class Intervals describing the conditioning set
#' @export
getInterval_permutation <- function(tree, nu, splits) {
  list1 <- combinat::permn(splits)
  allbounds1 <- sapply(list1, function(u) getInterval(tree, nu, u), simplify=FALSE)
  phi_bounds <- allbounds1[[1]]
  if (length(allbounds1) > 1) {
    for (l in 2:length(allbounds1)) {
      phi_bounds <- interval_union(phi_bounds, allbounds1[[l]])
    }
  }
  return(phi_bounds)
}


#' Get a pvalue for a difference in means between two terminal nodes
#'
#' @export
#'
#' @param tree an rpart object. Must have been built with model=TRUE.
#' @param locTest identify the pair of nodes that you are testing for
#' @param sigma_y enter the known noise SD (for now)
getSplitPval <- function(tree, locTest, sigma_y) {
  splits <- getAncestors(tree, locTest[1])
  nu <- (tree$where==locTest[1])/sum((tree$where==locTest[1])) - (tree$where==locTest[2])/sum(tree$where==locTest[2])
  phi_bounds <- getInterval(tree,nu, splits)
  y <- tree$model[,1]
  return(correctPVal(phi_bounds, nu, y, sigma_y))
}

#' Get a pvalue for a difference in means between two terminal nodes that are not siblings
#'
#' @export
#'
#' @param tree an rpart object
#' @param X the X data; as a matrix. names must match those
#' @param y the y data
#' @param locTest identify the pair of nodes that you are testing for
#' @param sigma_y enter the known noise SD (for now)
getSplitPval_nonSibs <- function(tree, locTest, sigma_y) {
  nu <- (tree$where==locTest[1])/sum((tree$where==locTest[1])) - (tree$where==locTest[2])/sum(tree$where==locTest[2])

  splits1 <- getAncestors(tree, locTest[1])
  splits2 <- getAncestors(tree, locTest[2])

  list1 <- combinat::permn(splits1)
  list2 <-  combinat::permn(splits2)

  allbounds1 <- sapply(list1, function(u) getInterval(tree, nu, u), simplify=FALSE)
  allbounds2 <- sapply(list2, function(u) getInterval(tree, nu, u), simplify=FALSE)
  bounds1 <- allbounds1[[1]]
  bounds2 <- allbounds2[[1]]
  if (length(allbounds1) > 1) {
  for (l in 2:length(allbounds1)) {
    bounds1 <- interval_union(bounds1, allbounds1[[l]])
  }
  }
  if (length(allbounds2) > 1) {
  for (l in 2:length(allbounds2)) {
    bounds2 <- interval_union(bounds2, allbounds2[[l]])
  }
  }
  phi_bounds <- interval_intersection(bounds1, bounds2)
  y <- tree$model[,1]
  return(correctPVal(phi_bounds, nu, y, sigma_y))
}



#' Constructs a large matrix storing a pvalue for every split and a CI for every node.
#' Used as a precursor to plotting.
#'
#' @param tree An rpart object. Must have been built with rpart arguement "model=TRUE".
#' @param sigma_y The known error variance. If not provided, is estimated with a conservative guess.
#'
#' @return A large matrix storing lots of pvalues and confidence intervals.
fullTreeInference <- function(tree, sigma_y =
                               sd(var(tree$model[,1]))) {

  if (is.null(tree$model)) {
    stop('Must build rpart object with parameter model=TRUE')
  }

  ## This is unfortunately necessary to make sure I can read the splits correctly
  ## out of the R output
  dat <- tree$model
  X <- dat[,-1]
  y <- dat[,1]
  p <- NCOL(X)
  n <- NROW(X)

  terminalNodes <- sort(unique(tree$where))
  splitResults <- data.frame(split = NA, pval = NA,
                             effectSize = NA, lower = NA, upper=NA,
                             branch1mean=NA,  branch1lower=NA, branch1upper=NA,
                             branch2mean=NA,branch2lower=NA, branch2upper=NA)

  fullMean <- mean(y)
  splitResults[1,] <- c(" ", NA, NA, fullMean - 1.96*sigma_y/sqrt(n),
                        fullMean + 1.96*sigma_y/sqrt(n), NA, NA,NA,NA, NA, NA)

  j <- 1
  k <- 2

  #### DO THE TOP LAYER (only once).
  splits <- getAncestors(tree, terminalNodes[1])
  splitText <- splits[1]

  this_split <- eval(parse(text =splitText))
  temp_tree <- tree
  temp_tree$where <- as.numeric(this_split)
  locTest <- c(0,1)
  y1 <- y[temp_tree$where==locTest[1]]
  y2 <- y[temp_tree$where==locTest[2]]
  nu <- (temp_tree$where==locTest[1])/sum((temp_tree$where==locTest[1])) - (temp_tree$where==locTest[2])/sum(temp_tree$where==locTest[2])
  sample_signal <- t(nu)%*%y
  phi_bounds <- getInterval(temp_tree, nu, splits[1])



  p_split <- correctPVal(phi_bounds, nu, y, sigma_y)
  CI_split <- computeCI(nu,y,sigma_y, phi_bounds, 0.05)


  nu1 <- (temp_tree$where==locTest[1])/sum((temp_tree$where==locTest[1]))



  phiBounds1 <- getInterval(temp_tree, nu1, splits[1])
  CI1 <- computeCI(nu1,y,sigma_y, phiBounds1, 0.05)

  nu2 <- (temp_tree$where==locTest[2])/sum((temp_tree$where==locTest[2]))
  phiBounds2 <- getInterval(temp_tree, nu2, splits[1])
  CI2 <- computeCI(nu2,y,sigma_y, phiBounds2, 0.05)

  ### NOTE TO SELF. MAKE SURE I UNDERSTAND THAT THE "TRUE" side of the branch
  ## is actually labeld as 2 here...
  ## MIGHT BE BACKWARD.
  splitResults[k,] <- c(splitText,p_split, mean(y1)-mean(y2), CI_split, mean(y1), CI1 , mean(y2), CI2)
  k <- k+1


  while(j < length(terminalNodes)) {
    #print(paste("Working on Branch",j))
    node <- terminalNodes[j]
    if (node+1==terminalNodes[j+1]) {
      j <- j+2
    } else {j <- j+1}

    splits <- getAncestors(tree, node)

    #### LOOP THROUGH THE LAYERS BELOW
    for (i in 1:(length(splits)-1)) {
      splitText1 <- paste(splits[1:(i+1)], collapse=" & ")
      splitText2 <- paste(paste(splits[1:i], collapse=" & "), paste("& !",splits[i+1]))
      ### CHECK THAT WE DIDN"T ALREADY DO THIS WORK.
      if (!(splitText1 %in% splitResults$split)) {
        #print("Found a new node!!")
        node1 <- eval(parse(text =splitText1))
        node2 <- eval(parse(text =splitText2))*2


        temp_tree$where <- node1+node2

        locTest <- c(1,2)
        y1 <- y[temp_tree$where==locTest[1]]
        y2 <- y[temp_tree$where==locTest[2]]

        repeated <- which(splitResults$branch2mean==mean(y1))
        if (length(repeated) > 0) {
          if (splitResults[repeated,]$branch1mean == mean(y2)) {
            next
          }
        }

        nu <- (temp_tree$where==locTest[1])/sum((temp_tree$where==locTest[1])) - (temp_tree$where==locTest[2])/sum(temp_tree$where==locTest[2])
        sample_signal <- t(nu)%*%y
        phi_bounds <- getInterval(temp_tree, nu,splits[1:(i+1)])
        p_split <- correctPVal(phi_bounds, nu, y, sigma_y)
        CI_split <- computeCI(nu,y,sigma_y, phi_bounds, 0.05)

        #### I THINK THAT ONE OF THESE SPLITS THINGS COULD BE BACKWARD
        nu1 <- (temp_tree$where==locTest[1])/sum((temp_tree$where==locTest[1]))
        phiBounds1 <- getInterval(temp_tree,nu1, splits[1:(i+1)])
        CI1 <- computeCI(nu1,y,sigma_y, phiBounds1, 0.05)

        nu2 <- (temp_tree$where==locTest[2])/sum((temp_tree$where==locTest[2]))
        phiBounds2 <- getInterval(temp_tree, nu2, splits[1:(i+1)])
        CI2 <- computeCI(nu2,y,sigma_y, phiBounds2, 0.05)

        #if (mean(y1) < CI1[1] | mean(y1) > CI1[2]) {
        #  print(splitText1)
        #}
        splitResults[k,] <- c(splitText1,p_split, mean(y1)-mean(y2), CI_split, mean(y1), CI1 , mean(y2), CI2)
        k <- k+1
      }
    }
  }
  splitResults[,2:11] <- apply(splitResults[,2:11],2,as.numeric)

  return(splitResults)
}



