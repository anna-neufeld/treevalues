library(treevalues)
library(rpart)
devtools::load_all()
bls.tree2 <-rpart(kcal24h0~hunger+disinhibition+resteating+rrvfood+liking+wanting, model = TRUE, data = blsdata, cp=0.02)

treeval.plot(bls.tree2)
treeval.plot2(bls.tree2)



### No CI, No Pval
treeval.plot(bls.tree2, inferenceType=0,inferenceMatrix = mat)
treeval.plot(bls.tree2, inferenceType=0,inferenceMatrix = mat,printn=FALSE)


# No fitted mean, yes pval
treeval.plot(bls.tree2, inferenceType=1,inferenceMatrix = mat)
treeval.plot(bls.tree2, inferenceType=1,inferenceMatrix = mat, printn=FALSE)


treeval.plot(bls.tree2, inferenceType=2, inferenceMatrix = mat)
treeval.plot(bls.tree2, inferenceType=2,inferenceMatrix = mat, printn=FALSE)


 treeval.plot2(bls.tree2)
 treeval.plot(bls.tree2, inferenceType=3,inferenceMatrix = mat, printn=TRUE)
 treeval.plot(bls.tree2, inferenceType=3,inferenceMatrix = mat, printn=FALSE)

