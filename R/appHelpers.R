buildOrderObjectShiny <- function(gch, hcg, method,
    coordinatesObject, updateProgress)
{
    if (coordinatesObject$weightStart == 0 | 
        coordinatesObject$weightStop == 0)
    {
        orderObject <- initialOrder(gch, hcg, Method=method)
    }
    else orderObject <- initialOrder(gch, hcg, Method=method,
        weightStart=coordinatesObject$weightStart,
        weightEnd=coordinatesObject$weightStop,
        weightFeature=coordinatesObject$weightColor,
        updateProgress=updateProgress)
    return(orderObject)
}

refineOrderShiny <- function(orderObject, refineMethod, coordinatesObject)
{
    refineFunction(orderObject, coordinatesObject$refineStart, 
        coordinatesObject$refineStop, Method=refineMethod)
}

makePlot <- function(orderObject, coordinatesObject, drawLines=TRUE, ...)
{

    plotSequence(orderObject, ...)
    if (coordinatesObject$refineStart != 0 & coordinatesObject$refineStop != 0)
    ## Draw the horizontal lines
    {
        n <- nrow(orderObject$toClust)
        ymin <- (((n:1)[coordinatesObject$refineStart] / 
                    n * (n - 10)) + 10) / n 
        ## Convert back to raw coordinates
        ymax <- (((n:1)[coordinatesObject$refineStop] / 
                    n * (n - 10)) + 10) / n
        if (drawLines)
        {
            abline(b=0, a=ymax, col="blue", lwd=2.5)
            abline(b=0, a=ymin, col="blue", lwd=2.5)
        }

    }
    if (coordinatesObject$weightStart != 0 & coordinatesObject$weightStop != 0) 
    ## Draw the vertical lines
    {
        m <- ncol(orderObject$toClust) / 2 # Convert back to raw coordinates
        xmin <- (coordinatesObject$weightStart / m) * 0.45
        xmax <- (coordinatesObject$weightStop / m) * 0.45
        if (coordinatesObject$weightColor == "yellow")
        {
            xmin <- xmin + 0.55
            xmax <- xmax + 0.55
        }
        if (drawLines)
        {
            abline(v=xmin, col="green", lwd=2.5)
            abline(v=xmax, col="green", lwd=2.5)
        }
    }
}

handleBrushCoordinates <- function(plotBrush, n, m){
    weightColor <- "red"

    firstRowRaw <- round(plotBrush$ymin * n) - 10
    lastRowRaw <- round(plotBrush$ymax * n) - 10
    firstRow <- round((firstRowRaw / (n - 10)) * n)
    lastRow <- round((lastRowRaw / (n - 10)) * n)

    if (firstRow <= 2) firstRow <- 1
    if (lastRow >= n - 1) lastRow <- n

    if (firstRow >= n - 1 | lastRow <= 2)
    {
        firstRow <- 0
        lastRow <- 0
    }


    firstCol <- round(plotBrush$xmin, 2)
    lastCol <- round(plotBrush$xmax, 2)

    if (firstCol <= 0.45) # Red weighting
    {
        if (lastCol >= 0.45) lastCol <- 0.45 
        ## Force the last column to be in the red
        firstCol <- firstCol / 0.45
        firstCol <- round(firstCol * m)
        lastCol <- lastCol / 0.45
        lastCol <- round(lastCol * m)

    }
    else if (firstCol >= 0.55) # Yellow weighting
    {
        weightColor <- "yellow"
        firstCol <- firstCol - 0.55
        lastCol <- lastCol - 0.55

        firstCol <- firstCol / 0.45
        firstCol <- round(firstCol * m)
        lastCol <- lastCol / 0.45
        lastCol <- round(lastCol * m)
    }
    else # In the middle, just set them to 0
    {
        firstCol <- 0
        lastCol <- 0
    }

    if (firstCol <= 2) firstCol <- 1
    if (lastCol >= (m - 2)) lastCol <- m

    return(list(firstRow=ifelse(firstRow == 0, 0, (n:1)[firstRow]),
                lastRow=ifelse(lastRow == 0, 0, (n:1)[lastRow]),
                firstCol=firstCol,
                lastCol=lastCol,
                weightColor=weightColor))

}

