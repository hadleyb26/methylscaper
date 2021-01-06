#' Ordering of methylation data
#'
#' Orders methylation data using a given seriation method. This function 
#' can also perform a weighted seriation if the method is set to "PCA".
#'
#' @param inputGCH The GCH input data.
#' @param inputHCG The HCG input data.
#' @param Method Indicates the seriation method to use. 
#'               The default option is "PCA", which orders the data using 
#'               the first principal component. Any seriation method provided 
#'               in the \code{seriation} package is valid.
#' @param weightStart Index of the first column used in the weighted seriation.
#' @param weightEnd Index of the last column used in the weighted seriation.
#' @param weightFeature Indicates whether to weight the GCH or HCG data.
#' @param updateProgress A function to handle the progress bar for the 
#'                       Shiny app. Should not be used when using the function 
#'                       independently.
#' @return An object of class \code{orderObject}, which contains the generated 
#'         ordering and the cleaned data matrix.
#' @importFrom seriation seriate get_order
#' @importFrom stats as.dist dist
#' @importFrom Rfast Dist
#' @export
initialOrder <- function(inputGCH, inputHCG, Method="PCA", weightStart=NULL, 
                         weightEnd=NULL, weightFeature="red", 
                         updateProgress=NULL){

    ## File checks:
    if (nrow(inputHCG) != nrow(inputGCH)) 
       {stop("Input files have different numbers of rows.")}

    if (all(rownames(inputGCH) == inputGCH[, 1])) inputGCH <- inputGCH[, -1]
    if (all(rownames(inputHCG) == inputHCG[, 1])) inputHCG <- inputHCG[, -1]


    if (is.function(updateProgress)) 
        updateProgress(message="Recoding input data", value=0.1)
    recoded <- recode(inputGCH, inputHCG)
    inputGCH <- recoded$inputGCH
    inputHCG <- recoded$inputHCG


    ## Clustering:
    toClust <- cbind(inputGCH, inputHCG)
    weighted=FALSE

    ## (Optional) Weighting: Adds a variable indicating the number of 
    ## red or yellow patches at specific DNA location
    if (!is.null(weightStart) & !is.null(weightEnd)) {
        if (is.function(updateProgress)) 
            updateProgress(message="Weighting selected columns", value=0.2)
        weighted=TRUE
        if (weightFeature == "red") {
            FEATURE=3
            weightVector <- apply(inputHCG[, weightStart:weightEnd], 1, 
                                  function(x) 
                sum(x == FEATURE))
        }
        if (weightFeature == "yellow") {
            FEATURE=-3
            weightVector <- apply(inputGCH[, weightStart:weightEnd], 1, 
                                  function(x) 
                sum(x == FEATURE))
        }
        weightVector[weightVector == 0] <- 1 # We dont want to have 0 weights
    }


    if (is.function(updateProgress)) 
        updateProgress(message=paste("Ordering with", Method), value=0.35)
    ## PCA should be the default method:
    if (Method == "PCA") {
        if (weighted)
        {
            w <- weightVector / sum(weightVector)
            wSqrt <- sqrt(w)
            toClustWeighted <- diag(wSqrt) %*% toClust

            colCentered <- apply(toClustWeighted, 2, function(x) x - mean(x))
            try1 <- svd(colCentered, nu=1, nv=0)
            order1 <- order(try1$u[, 1])
        }
        else
        {
            colCentered <- apply(toClust, 2, function(x) x - mean(x))
            try1 <- svd(colCentered, nu=1, nv=0)
            order1 <- order(try1$u[, 1])
        }

    } 
    else {

        if (weighted)
        {
            w <- weightVector / sum(weightVector)
            wSqrt <- sqrt(w)
            toClustWeighted <- diag(wSqrt) %*% toClust

            distMat <- asDist(Rfast::Dist(toClustWeighted, 
                                          method="euclidean"))
            order1 <- seriation::seriate(distMat, method=Method)
            order1 <- seriation::get_order(order1)
        }
        else
        {
            distMat <- as.dist(Rfast::Dist(toClust, method="euclidean"))
            order1 <- seriation::seriate(distMat, method=Method)
            order1 <- seriation::get_order(order1)

        }
    }
    orderObject <- list(toClust=toClust, order1=order1)
    if (Method != "PCA") orderObject$distMat <- distMat
    if (weighted) orderObject$weights <- weightVector
    if (is.function(updateProgress)) updateProgress(message="Done", value=1)
    return(orderObject)
}


recode <- function(inputGCH, inputHCG)
{
    inputGCH[inputGCH == "."] <- 99
    inputGCH <- apply(inputGCH, 2, as.numeric)
    inputGCH[inputGCH == 2] <- -4
    inputGCH[inputGCH == 1] <- -3
    inputGCH[inputGCH == 0] <- -2.5
    inputGCH[inputGCH == -1] <- -99
    inputGCH[inputGCH == -2] <- -1
    inputGCH[inputGCH == -99] <- -2
    inputGCH[inputGCH == 99] <- 0

    inputHCG[inputHCG == "."] <- 99
    inputHCG <- apply(inputHCG, 2, as.numeric)
    inputHCG[inputHCG == 2] <- 4
    inputHCG[inputHCG == 1] <- 3
    inputHCG[inputHCG == 0] <- 2.5
    inputHCG[inputHCG == -1] <- 2
    inputHCG[inputHCG == -2] <- 1
    inputHCG[inputHCG == 99] <- 0

    return(list(inputGCH=inputGCH, inputHCG=inputHCG))

}

Delete 
