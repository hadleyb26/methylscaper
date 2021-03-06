% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/refinement.R
\name{forceReverse}
\alias{forceReverse}
\title{Force reversal of a subset of the ordering}
\usage{
forceReverse(
  orderObject,
  reverseStart = 1,
  reverseEnd = length(orderObject$order1)
)
}
\arguments{
\item{orderObject}{An object of class \code{orderObject}, 
generated with the \code{initialOrder} function.}

\item{reverseStart}{The first index to be included in the reversal.}

\item{reverseEnd}{The last index to be included in the reversal.}
}
\value{
The new complete ordering, with the reversal applied.
}
\description{
This reverses a subset of the ordering, as determined by the user. 
By default, the entire ordering is reversed.
}
