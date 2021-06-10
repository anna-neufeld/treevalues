# treevalues: inference for CART regression trees. 

What is treevalues?
-----

The ``treevalues`` R package computes confidence intervals and p-values for the mean response within a region or the difference in mean response between two regions in a CART regression tree (built using the package ``rpart``).  As the regions in a regression tree are selected using the data, the same data cannot naively be used to do inference on these means. The ``treevalues`` package implements a selective inference approach to conduct inference without *double dipping* in the data. 


How can I get treevalues?
-----

Make sure that ``remotes`` is installed by running ``install.packages("remotes")``, then type

```R
remotes::install_github("anna-neufeld/treevalues")
```

Where can I learn more? 
-----

See [https://anna-neufeld.github.io/treevalues/articles/inference_tutorial.html](https://anna-neufeld.github.io/treevalues/articles/inference_tutorial.html) for instructions on how to use this package on real data. See [???](???) for the preprint.

See [https://github.com/anna-neufeld/treevalues-simulations](https://github.com/anna-neufeld/treevalues-simulations) for code to reproduce the experiements and figures in the preprint. 




