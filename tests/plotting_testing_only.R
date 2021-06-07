# library(treevalues)
# library(rpart)
# devtools::load_all()
# bls.tree2 <-rpart(kcal24h0~hunger+disinhibition+resteating+rrvfood+liking+wanting, model = TRUE, data = blsdata, cp=0.02)
#
# ## Fast
# treeval.plot(bls.tree2, inferenceType=0)
#
# ## Slow
# treeval.plot(bls.tree2, inferenceType=1)
#
# ## A trick for speeding up!!!!!
# bls.tree2 <- inferenceFrame(bls.tree2)
# treeval.plot(bls.tree2, inferenceType=2)
# treeval.plot(bls.tree2, inferenceType=3)
# treeval.plot(bls.tree2, inferenceType=4)
# treeval.plot2(bls.tree2, inferenceType=4,printn=FALSE)
#
# treeval.plot(bls.tree2, inferenceType=2, space=1.6)
