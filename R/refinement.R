#' Refinement step
#' 
#' Reorders a subset of the methylation data. 
#' 
#' @param orderObject An object of class \code{orderObject}, 
#'                    generated with the \code{initialOrder} function.
#' @param refineStart The index of the first sample (row) used in the 
#'                    refinement.
#' @param refineEnd The index of the last sample used in the refinement.
#' @param Method The seriation method used to perform the refinement.
#' 
#' @return The new complete ordering with the refinement applied.
#' @export

refineFunction <- function(orderObject, refineStart, refineEnd, 
                           Method="HC_average") {
  
    toClust <- orderObject$toClust
    order1 <- orderObject$order1
    toRefineOrder <- order1[refineStart:refineEnd]
  
  
    toRefineClust <- toClust[toRefineOrder, ]
  
    if (Method == "PCA") {
    colCentered <- apply(toRefineClust, 2, function(x) x - mean(x))
    try1 <- svd(colCentered, nu=1, nv=0)
    orderNew <- order(try1$u[,1])
    } else { # Methods available for refining are: ARSA, HC_complete, 
             ## HC_average, HC_ward.
    ## Shouldn't need to recalculate the distance each time!!
    if (is.null(orderObject$distMat)) distMat <- dist(toRefineClust,
                                                      method="euclidean")
    else distMat <- as.dist(as.matrix(orderObject$distMat)
                            [toRefineOrder,toRefineOrder])
    orderNew <- seriation::seriate(distMat, method=Method, verbose=FALSE)
    orderNew <- seriation::get_order(orderNew)
  }
  
  
    ## New order:
    orderNew <- order1[refineStart:refineEnd][orderNew]
    orderFinal <- order1
    orderFinal[refineStart:refineEnd] <- orderNew
  
    return(orderFinal)
}

#' Force reversal of a subset of the ordering
#' 
#' This reverses a subset of the ordering, as determined by the user. 
#' By default, the entire ordering is reversed.
#' 
#' @param orderObject An object of class \code{orderObject}, 
#'                    generated with the \code{initialOrder} function.
#' @param reverseStart The first index to be included in the reversal.
#' @param reverseEnd The last index to be included in the reversal.
#' 
#' @return The new complete ordering, with the reversal applied.
#' @export

forceReverse <- function(orderObject, reverseStart=1, 
                         reverseEnd=length(orderObject$order1))
{
    order1 <- orderObject$order1
    order1[reverseStart:reverseEnd] <- rev(order1[reverseStart:reverseEnd])
    return(order1)
}
