% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_functions.R
\name{treeval.plot}
\alias{treeval.plot}
\title{Plot an rpart object with all splits labeled with p-values and all nodes
labeled with confidence intervals.}
\usage{
treeval.plot(
  tree,
  sigma_y = NULL,
  nn = TRUE,
  printn = TRUE,
  inferenceType = 2,
  digits = 3,
  alpha = 0.05,
  permute = FALSE,
  ...
)
}
\arguments{
\item{tree}{An rpart tree. The tree must have been build with parameter \code{model=TRUE}.}

\item{sigma_y}{Provide the standard deviation of y, if known. If not provided, the sample standard deviation of y will be used
as a conservative estimate.}

\item{nn}{boolean- would you like node numbers to be printed? Nodes are numbered using the same methodology as the \code{rpart} package. If node \code{n} has children,
its children are numbered \code{2n} and \code{2n+1}.}

\item{printn}{boolean - would you like the number of observations to be printed in each node?}

\item{inferenceType}{An integer specifying which pieces of inference information should be added to the plot. The options
currently available are
(0) No confidence intervals, p-values, or "fitted mean" label. Just calls \code{rpart.plot()}.
(1) No confidence intervals. Each split labeled with a p-value.
(2)  Label each internal node with a confidence interval and label each split with a p-value. This is the default, but
can also be a little messy/hard to read. Options 3 and 4 print the same information but with small
formatting tweaks.}

\item{digits}{Integer- how many digits would you like the text in the plot rounded to.}

\item{alpha}{If inferenceType is such that confidence intervals will be printed, \code{(1-alpha)} confidence intervals will be printed.}

\item{permute}{If inferenceType is such that confidence intervals will be printed, should the conditioning set for the confidence intervals include
all permutations of the relevant branch? Setting this to TRUE will lead to slightly narrower confidence intervals, but will make computations more expensive.
See paper for more details.}

\item{...}{Additional arguments are passed on to rpart.plot(). Examples include \code{"cex"}.}
}
\description{
Essentially a wrapper function for \code{rpart.plot()} from the \code{rpart.plot} package, but
additional arguments allow user to add p-values and confidence intervals to plots.
}
\examples{
bls.tree <-rpart::rpart(kcal24h0~hunger+disinhibition+resteating,
    model = TRUE, data = blsdata, maxdepth=1)
treeval.plot(bls.tree, inferenceType=0)
treeval.plot(bls.tree, inferenceType=1)
treeval.plot(bls.tree, inferenceType=2)
}
