% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/user_funs.R
\name{getBranch}
\alias{getBranch}
\title{Get branch or list of branches from tree.}
\usage{
getBranch(tree, nn = NULL)
}
\arguments{
\item{tree}{An rpart object.}

\item{nn}{A node number that corresponds to a valid node in \code{tree}. The list of valid node numbers can be obtained with
\code{row.names(tree$frame)} or by plotting \code{tree} with \code{treeval.plot()}. The node number can be passed in as
either a character string or an integer. If no node number is provided, a list of all branches in the tree will be returned.}
}
\value{
Either a single branch (which is a vector of splits) or (if nn=NULL), a list of all branches in the tree.
}
\description{
Given an rpart tree and a node number, returns a vector of strings that describes the branch which defines the node.
If no node number is provided, returns a list describing every branch in the tree.
Useful for extracting individual branches, which are necessary inputs to \code{branchInference()}.
}
\examples{
bls.tree <- rpart::rpart(kcal24h0~hunger+disinhibition+resteating+rrvfood+liking+wanting,
    model = TRUE, data = blsdata, cp=0.02)
branch <- getBranch(bls.tree, 8)
branchInference (bls.tree, branch, type="sib")
}
