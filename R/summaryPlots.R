#' Calculates the percentage of C at each methylation site
#' 
#' @param orderObject An object of class \code{orderObject}
#' @param plotPercents Logical, indicates whether to 
#'                     generate the percentage plot
#' @param ... Additional parameters used by the \code{plot} function.
#' 
#' @importFrom graphics hist lines plot points
#' @export
percentC <- function(orderObject, plotPercents=FALSE, ...){
    dat <- orderObject$toClust
    redSites <- which(dat[1, 1:ncol(dat)] == 4 |
                      dat[1, 1:ncol(dat)] == 1)

    yellowSites <- which(dat[1,1:ncol(dat)] == -4 |
                         dat[1, 1:ncol(dat)] == -1)
    cRed <- sapply(redSites, function(i){
        sum(dat[, i] == 4) / nrow(dat)
    })
    cYellow <- sapply(yellowSites, function(i){
        sum(dat[, i] == -4) / nrow(dat)
    })
    if (plotPercents)
    {
        plot(x=redSites - ncol(dat)/2, y=cRed,
            col="brown1", pch=19, ylim=c(0, 1), xlab="Region (Base Pair)", 
            ylab="%C",
            bty='n', cex.lab=1.3, xaxt='n', yaxt='n', ...)
        axis(side=1, lwd=2, cex.axis=1.2)
        axis(side=2, lwd=2, cex.axis=1.2)
        lines(x=redSites - ncol(dat)/2, y=cRed, col="brown1")
        points(x=yellowSites, y=cYellow, col="gold2", pch=19)
        lines(x=yellowSites, y=cYellow, col="gold2")

        nSites <- length(union(redSites, yellowSites))
        labs <- union(redSites, yellowSites)[seq(1, nSites, by=nSites/12)]
    }
    final <- list(cRed, cYellow)
    names(final) <- c("red", "yellow")
    return(final)

}

#' Calculate the proportion of methylated bases for the GCH and HCG data sets.
#'
#' @param orderObject An object of class \code{orderObject}
#' @param color Indicates which data set to compute proportions for
#' @param plotHistogram Indicates whether to plot a histogram of 
#'                      the proportions across all reads.
#' @param ... Additional parameters used by the \code{hist} function.
#'
#' @importFrom graphics hist
#' @export
proportionColor <- function(orderObject, color="YELLOW", 
                            plotHistogram=FALSE, ...){
    colorIndicator <- ifelse(color == "YELLOW", -1, 1)
    proportion <- apply(orderObject$toClust, 1, function(x){
    sum(x == colorIndicator * 3 | x == colorIndicator * 4) / (length(x) / 2)
  })
    if (plotHistogram) {
    opar <- par(lwd=4)
    h <- hist(proportion, plot=F, breaks=15)
    plot(h, xlim=c(0, 1), border=ifelse(color == "YELLOW", "gold2", "brown1"),
         col="gray75",
         lwd=2, ...)
    par(opar)
  }
    return(proportion)
}

#' Calculate the average methylation/accessibility status across all reads.
#'
#' @param orderObject An object of class \code{orderObject}
#' @param window_length Length of the window to be used to compute 
#'                      a moving average.
#' @param plotAverages Logical, indicates whether to generate a line plot 
#'                     of average status.
#' @param ... Addition parameters used by the \code{plot} function.
#'
#' @importFrom stats filter
#' @importFrom graphics legend
#' @export
averageStatus <- function(orderObject, windowLength=1, plotAverages=FALSE, ...)
{
    gchNum <- orderObject$toClust[,1:(ncol(orderObject$toClust) / 2)]
    hcgNum <- orderObject$toClust[,(ncol(orderObject$toClust) / 2 + 1)
                                  :ncol(orderObject$toClust)]

    accSum <- colSums(gchNum == -3)
    methSum <- colSums(hcgNum == 3)
    accDenom <- colSums(gchNum != 0)
    methDenom <- colSums(hcgNum != 0)
    accAvg <- accSum / accDenom
    methAvg <- methSum / methDenom

    width <- windowLength
    movingAccAvg <- filter(x=accAvg, filter=rep(1, width)) / width
    movingMethAvg <- filter(x=methAvg, filter=rep(1, width)) / width

    if (plotAverages)
    {
        plot(movingAccAvg, type="l", col="gold2",
             xlab="Position along read", ylab="Population-averaged status", 
             ylim=c(0,1))
        lines(movingMethAvg, col="brown1")
        legend("topright", legend=c("Methylation", "Accessibility"), 
               fill=c("brown1", "gold2"))
    }
    return(list(methAvg=movingMethAvg, accAvg=movingAccAvg))
}
