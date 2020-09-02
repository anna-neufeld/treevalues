fullTreeInference_CAT <- function(tree, sigma_y) {

  ### Add a warning to make sure they build their tree right!!!

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

  j <- 1
  k <- 1

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
  phi_bounds <- getInterval_minBUCKET(temp_tree, nu, splits[1])



  p_split <- correctPVal(phi_bounds, nu, y, sigma_y)
  CI_split <- computeCI(nu,y,sigma_y, phi_bounds, 0.05)


  nu1 <- (temp_tree$where==locTest[1])/sum((temp_tree$where==locTest[1]))



  phiBounds1 <- getInterval_minBUCKET(temp_tree, nu1, splits[1])
  CI1 <- computeCI(nu1,y,sigma_y, phiBounds1, 0.05)

  nu2 <- (temp_tree$where==locTest[2])/sum((temp_tree$where==locTest[2]))
  phiBounds2 <- getInterval_minBUCKET(temp_tree, nu2, splits[1])
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
        phi_bounds <- getInterval_minBUCKET(temp_tree, nu,splits[1:(i+1)])
        p_split <- correctPVal(phi_bounds, nu, y, sigma_y)
        CI_split <- computeCI(nu,y,sigma_y, phi_bounds, 0.05)

        #### I THINK THAT ONE OF THESE SPLITS THINGS COULD BE BACKWARD
        nu1 <- (temp_tree$where==locTest[1])/sum((temp_tree$where==locTest[1]))
        phiBounds1 <- getInterval_minBUCKET(temp_tree,nu1, splits[1:(i+1)])
        CI1 <- computeCI(nu1,y,sigma_y, phiBounds1, 0.05)

        nu2 <- (temp_tree$where==locTest[2])/sum((temp_tree$where==locTest[2]))
        phiBounds2 <- getInterval_minBUCKET(temp_tree, nu2, splits[1:(i+1)])
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



