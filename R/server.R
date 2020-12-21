server <- function(input, output) {

    actionsLog <- reactiveValues(log=c("")) 
    ## logs the actions taken wrt the plot
  
##Single-cell data 

    scSeqData <- reactiveValues(gch=NULL, hcg=NULL) # for raw data
    scInputData <- reactiveValues(gch=NULL, hcg=NULL) # for state matrices

    observe({
    if (!is.null(input$gchSeqFile) & !is.null(input$hcgSeqFile))
    {
        isolate({
        progress <- Progress$new()
        progress$set(message="Loading GC data", value=0)
        on.exit(progress$close())

        scSeqData$gch <- readRDS(input$gchSeqFile$datapath)
        progress$set(message="Loading CG data", value=0.5)
        scSeqData$hcg <- readRDS(input$hcgSeqFile$datapath)
        actionsLog$log <- c(actionsLog$log, paste("Loading GCH RDS file:"
                                                  , input$gchSeqFile$name))
        actionsLog$log <- c(actionsLog$log, paste("Loading HCG RDS file:"
                                                  , input$hcgSeqFile$name))
        })
        }
    })


    output$positionSlider <- renderUI({
        if (!is.null(scSeqData$gch) & !is.null(scSeqData$hcg) &
            input$startPos!=0 & input$endPos!=0)
        {
            cgMaxPos <- max(sapply(scSeqData$hcg, FUN=function(x) max(x$pos)))
            cgMinPos <- min(sapply(scSeqData$hcg, FUN=function(x) min(x$pos)))
            gcMaxPos <- max(sapply(scSeqData$gch, FUN=function(x) max(x$pos)))
            gcMinPos <- min(sapply(scSeqData$gch, FUN=function(x) min(x$pos)))
            start <- input$startPos
            end <- input$endPos

            if (start < cgMinPos | start < gcMinPos | end > cgMaxPos | 
                end > gcMaxPos)
            {
                showNotification("Selected range is out of bounds. 
                                 Please choose a valid starting and 
                                 end position to generate the plot.",
                                 type="error")
                return(NULL)
            }
            if (end -  start > 10000)
            {
                showNotification("Selected range is longer than 10k bp,
                                 reducing length for stability.",
                                 type="warning")
                end <- start + 10000
            }
            len <- end - start

            sliderInput(inputId="positionSliderInput",
                        label="Position adjustment slider",
                        min=start - len, max=end + len,
                        value=c(start, end))
        }

    })

    observe({
    if (!is.null(input$positionSliderInput))
    {
        progress <- Progress$new()
        progress$set(message="Beginning single-cell processing", value=0)
        on.exit(progress$close())

        updateProgress <- function(value=NULL, message=NULL, detail=NULL) {
            progress$set(value=value, message=message, detail=detail)}

        prepOut <- prepSC(scSeqData$gch, scSeqData$hcg, 
                          input$positionSliderInput[1],
                           input$positionSliderInput[2],
                           updateProgress=updateProgress)

        tempGch <- prepOut$gch
        tempHcg <- prepOut$hcg
        if (nrow(tempGch) == nrow(tempHcg))
        {
            scCoordinatesObject$refineStart <- 0
            scCoordinatesObject$refineStop <- 0
            scCoordinatesObject$weightStart <- 0
            scCoordinatesObject$weightStop <- 0
            scInputData$gch <- tempGch
            scInputData$hcg <- tempHcg
            isolate({
                actionsLog$log <- c(actionsLog$log, paste("Beginning 
                                    single-cell data analysis"))
                actionsLog$log <- c(actionsLog$log, paste("From position",
                                    input$positionSliderInput[1],
                                    "to", input$positionSliderInput[2]))
        })
        }
        }
    })

    ## this object keeps track of the coordinates for refinement and weighting
    scCoordinatesObject <- reactiveValues(refineStart=0, refineStop=0,
                                         weightStart=0, weightStop=0, 
                                         weight.color = "red")
    ## now construct the scOrderObject
    scOrderObject <- reactiveValues(toClust=0, order1=0)
    observe({ if (!is.null(scInputData$gch) & !is.null(scInputData$hcg))
    {
    progress <- Progress$new()
    progress$set(message="Beginning seriation", value=0)
    on.exit(progress$close())

    updateProgress <- function(value=NULL, message=NULL, detail=NULL) {
        progress$set(value=value, message=message, detail=detail)}

    tempObj <- buildOrderObjectShiny(scInputData$gch, scInputData$hcg,
                                     input$scSerMethod,
                                     scCoordinatesObject, updateProgress)
    scOrderObject$order1 <- tempObj$order1
    scOrderObject$toClust <- tempObj$toClust
    isolate({
        actionsLog$log <- c(actionsLog$log,
                          paste("Ordering with", input$scSerMethod))
        })
    }

    })

    ## this handles updates to scCoordinatesObject
    observeEvent(input$scPlotBrush, {
    n <- nrow(scInputData$gch)
    m <- ncol(scInputData$hcg)
    processed.brush <- handleBrushCoordinates(input$scPlotBrush, n, m)

    if (isolate(input$scBrushChoice) == "Weighting")
    {
        scCoordinatesObject$refineStart <- 0
        scCoordinatesObject$refineStop <- 0
        scCoordinatesObject$weightStart <- processed.brush$first.col
        scCoordinatesObject$weightStop <- processed.brush$last.col
        scCoordinatesObject$weight.color <- processed.brush$weight.color
        isolate({
        actionsLog$log <- c(actionsLog$log,
                            paste("Weighting", processed.brush$weight.color,
                                  "columns",
                                  processed.brush$first.col, "to",
                                  processed.brush$last.col))
        })
    }
    if (isolate(input$scBrushChoice) == "Refinement")
    {
        s <- processed.brush$first.row
        f <- processed.brush$last.row
        if (s!=f)
        {
        scCoordinatesObject$refineStart <- s
        scCoordinatesObject$refineStop <- f
        scOrderObject$order1 <- refineOrderShiny(isolate(scOrderObject),
                                  refine.method=isolate(input$scRefineMethod),
                                  scCoordinatesObject)
        isolate({
            actionsLog$log <- c(actionsLog$log,
                              paste("Refining rows",
                                    processed.brush$first.row, "to",
                                    processed.brush$last.row))
            actionsLog$log <- c(actionsLog$log,
                              paste("Applying refinement with",
                                    input$scRefineMethod))
            })
        }
    }    
    })

    observeEvent( input$scForceReverse, {
    isolate({
        if (scCoordinatesObject$refineStart == scCoordinatesObject$refineStop)
        {
        scOrderObject$order1 <- rev(scOrderObject$order1)
        actionsLog$log <- c(actionsLog$log, paste("Reversing rows 1 to",
                                                  nrow(scInputData$gch)))
        }
        else
        {
        scOrderObject$order1[scCoordinatesObject$refineStart : 
                               scCoordinatesObject$refineStop] <-
            scOrderObject$order1[scCoordinatesObject$refineStop : 
                                   scCoordinatesObject$refineStart]
        actionsLog$log <- c(actionsLog$log,
                            paste("Reversing rows", 
                                  scCoordinatesObject$refineStart,
                                  "to", scCoordinatesObject$refineStop))
            }
        })
    })

    output$scSeqPlot <- renderPlot({
        obj <- scOrderObject
        if (sum(obj$toClust) == 0) {showNotification("Select methylation
                                    data files to generate the plot.",
                                                     type="message");NULL}
       else makePlot(obj, isolate(scCoordinatesObject))
    }, height=600, width=600)

    output$scPlotDown <- downloadHandler(
    filename = function(){
        if (input$scFileType == "PNG") return("plot.png")
        if (input$scFileType == "SVG") return("plot.svg")
        if (input$scFileType == "PDF") return("plot.pdf")
    },
    content=function(file){
        if (input$scFileType == "PNG") png(file)
        if (input$scFileType == "SVG") svglite::svglite(file)
        if (input$scFileType == "PDF") pdf(file)

        makePlot(scOrderObject, scCoordinatesObject, drawLines=FALSE,
                 plotFAST=FALSE)
        dev.off()
        }
    )

    output$scLogDown <- downloadHandler(
        filename=function(){
        "changes.txt"
    },
    content=function(file){
        fileConn <- file(file)
        writeLines(actionsLog$log, fileConn)
        close(fileConn)
        }
    )

    output$scInfo <- renderText({
      paste0("Refinement selection: ", scCoordinatesObject$refineStart,
             " ", scCoordinatesObject$refineStop, "\n",
          "Weighting selection: ", scCoordinatesObject$weightStart, " ",
          scCoordinatesObject$weightStop)
    })

    output$scProportionColorHistogram <- renderPlot({
    obj <- scOrderObject
    if (sum(obj$toClust) == 0)
    {showNotification("Select methylation data files to generate the plot.",
                      type="message");NULL}
    else proportionColor(obj, plotHistogram=TRUE,
                         color=toupper(input$scProportionChoice))
    })

    output$scProportionHistDownload <- downloadHandler(
    filename=function(){
        if (input$filetype == "PNG") return("hist.png")
        if (input$filetype == "SVG") return("hist.svg")
        if (input$filetype == "PDF") return("hist.pdf")
    },
    content=function(file){
        if (input$filetype == "PNG") png(file)
        if (input$filetype == "SVG") svglite::svglite(file)
        if (input$filetype == "PDF") pdf(file)

         proportionColor(scOrderObject, plotHistogram=TRUE,
                       color=toupper(input$scProportionChoice))
        dev.off()
        }
        )
    output$scProportionDataDownload <- downloadHandler(
    filename=function(){
        return("proportion_data.csv")
    },
    content=function(file){
        dat <-  proportionColor(scOrderObject, plotHistogram=FALSE,
                               color=toupper(input$scProportionChoice))
        write.csv(dat, file=file)
        }
    )

    output$scPercentC <- renderPlot({
        obj <- scOrderObject
        if (sum(obj$toClust) == 0)
        {showNotification("Select methylation data files to generate the plot."
                          , type="message");NULL}
        else percent_C(obj, plotPercents=TRUE)
    })

    output$scPercentCPlotDownload <- downloadHandler(
        filename=function(){
        if (input$filetype == "PNG") return("percentC.png")
        if (input$filetype == "SVG") return("percentC.svg")
        if (input$filetype == "PDF") return("percentC.pdf")
    },
    content=function(file){
        if (input$filetype == "PNG") png(file)
        if (input$filetype == "SVG") svglite::svglite(file)
        if (input$filetype == "PDF") pdf(file)
  
      percent_C(scOrderObject, plotPercents=TRUE)
      dev.off()
      }
    )

    output$scPercentCDataDownload <- downloadHandler(
      filename=function(){
      return("proportion_data.RData")
    },
    content=function(file){
        dat <-  percent_C(scOrderObject, plotPercents=FALSE)
        save(dat, file=file)
        }
    )

##Single-molecule data 
    
    smInputData <- reactiveValues(gch=NULL, hcg=NULL)

    ##alignment handling
    observeEvent(input$run.align, {
        ref <- read.fasta(input$ref.file$datapath)
        fasta <- read.fasta(input$fasta.file$datapath)

    progress <- Progress$new()
    progress$set(message="Beginning alignment", value=0)
    on.exit(progress$close())

    updateProgress <- function(value=NULL, message=NULL, detail=NULL) {
        progress$set(value=value, message=message, detail=detail)}

    alignOut <- runAlign(ref, fasta, updateProgress=updateProgress,
                          log.file=input$processing.log.name)

    hcgFileName <- input$hcgFileName
    gchFileName <- input$gchFileName

    writeMethylationData(dat=alignOut$hcg, filepath=hcgFileName)
    writeMethylationData(dat=alignOut$gch, filepath=gchFileName)

    })

    observe({if (!is.null(input$smGchFile) & !is.null(input$smHcgFile))
    {
    tempGch <- readMethylationData(filepath=input$smGchFile$datapath)
    tempHcg <- readMethylationData(filepath=input$smHcgFile$datapath)
    if (nrow(tempGch) == nrow(tempHcg))
    {
        smCoordinatesObject$refineStart <- 0
        smCoordinatesObject$refineStop <- 0
        smCoordinatesObject$weightStart <- 0
        smCoordinatesObject$weightStop <- 0
        smInputData$gch <- tempGch
        smInputData$hcg <- tempHcg
        smInputData$datatype <- "sm"
        isolate({
            actionsLog$log <- c(actionsLog$log, paste("Beginning 
                                          single-molecule data analysis"))
            actionsLog$log <- c(actionsLog$log, paste("Loading GCH file:",
                                                      input$gch.file$name))
            actionsLog$log <- c(actionsLog$log, paste("Loading HCG file:",
                                                      input$hcg.file$name))
            })
        }

    }})

    ##this object keeps track of the coordinates for refinement and weighting
    smCoordinatesObject <- reactiveValues(refineStart=0, refineStop=0,
                                      weightStart=0, 
                                      weightStop=0, weight.color="red")
    # now construct the smOrderObject
    smOrderObject <- reactiveValues(toClust=0, order1=0)
    observe({ if (!is.null(smInputData$gch) & !is.null(smInputData$hcg))
    {
        progress <- Progress$new()
         progress$set(message="Beginning seriation", value=0)
    on.exit(progress$close())

    updateProgress <- function(value=NULL, message=NULL, detail=NULL) {
        progress$set(value=value, message=message, detail=detail)}

    tempObj <- buildOrderObjectShiny(smInputData$gch, smInputData$hcg,
                                     input$smSerMethod, smCoordinatesObject,
                                     updateProgress)
    smOrderObject$order1 <- tempObj$order1
    smOrderObject$toClust <- tempObj$toClust
    isolate({
        actionsLog$log <- c(actionsLog$log,
                          paste("Ordering with", input$smSerMethod))
        })
    }

  })

    ##this handles updates to smCoordinatesObject
    observeEvent(input$smPlotBrush, {
        n <- nrow(smInputData$gch)
        m <- ncol(smInputData$hcg)
        processed.brush <- handleBrushCoordinates(input$smPlotBrush, n, m)

    if (isolate(input$smBrushChoice) == "Weighting")
    {
        smCoordinatesObject$refineStart <- 0
        smCoordinatesObject$refineStop <- 0
        smCoordinatesObject$weightStart <- processed.brush$first.col
        smCoordinatesObject$weightStop <- processed.brush$last.col
        smCoordinatesObject$weight.color <- processed.brush$weight.color
        isolate({
        actionsLog$log <- c(actionsLog$log,
                            paste("Weighting", processed.brush$weight.color,
                                  "columns",
                                  processed.brush$first.col, "to",
                                  processed.brush$last.col))
         })
    }
    if (isolate(input$smBrushChoice) == "Refinement")
    {
        s <- processed.brush$first.row
        f <- processed.brush$last.row
        if (s!=f)
        {
        smCoordinatesObject$refineStart <- s
        smCoordinatesObject$refineStop <- f
        smOrderObject$order1 <- refineOrderShiny(isolate(smOrderObject),
                                  refine.method=isolate(input$smRefineMethod),
                                  smCoordinatesObject)
        isolate({
            actionsLog$log <- c(actionsLog$log,
                              paste("Refining rows",
                                    processed.brush$first.row, "to",
                                    processed.brush$last.row))
            actionsLog$log <- c(actionsLog$log,
                              paste("Applying refinement with", 
                                    input$smRefineMethod))
            })
        }
    }


    })

    observeEvent( input$smForceReverse, {
        isolate({
        if (smCoordinatesObject$refineStart == smCoordinatesObject$refineStop)
        {
        smOrderObject$order1 <- rev(smOrderObject$order1)
        actionsLog$log <- c(actionsLog$log, paste("Reversing rows 1 to",
                                                  nrow(smInputData$gch)))
        }
        else
        {
        smOrderObject$order1[smCoordinatesObject$refineStart : 
                               smCoordinatesObject$refineStop] <-
            smOrderObject$order1[smCoordinatesObject$refineStop : 
                                   smCoordinatesObject$refineStart]
        actionsLog$log <- c(actionsLog$log,
                            paste("Reversing rows", 
                                  smCoordinatesObject$refineStart,
                                  "to", smCoordinatesObject$refineStop))

            }
        })
    })



    output$smSeqPlot <- renderPlot({
        obj <- smOrderObject
        if (sum(obj$toClust) == 0) {showNotification("Select methylation data
                                                     files to generate the 
                                                     plot.", type="message")
          ;NULL}
        else makePlot(obj, isolate(smCoordinatesObject))
    }, height=600, width=600)

    output$smPlotDown <- downloadHandler(
        filename=function(){
            if (input$smFiletype == "PNG") return("plot.png")
            if (input$smFiletype == "SVG") return("plot.svg")
            if (input$smFiletype == "PDF") return("plot.pdf")
        },
    content=function(file){
        if (input$smFiletype == "PNG") png(file)
        if (input$smFiletype == "SVG") svglite::svglite(file)
        if (input$smFiletype == "PDF") pdf(file)

        makePlot(smOrderObject, smCoordinatesObject, drawLines=FALSE,
                 plotFAST=FALSE)
      dev.off()
        }
    )

    output$smLogDown <- downloadHandler(
        filename=function(){
            "changes.txt"
        },
    content=function(file){
        fileConn <- file(file)
        writeLines(actionsLog$log, fileConn)
        close(fileConn)
        }
    )

    output$smInfo <- renderText({
        paste0("Refinement selection: ", smCoordinatesObject$refineStart,
               " ", smCoordinatesObject$refineStop, "\n",
           "Weighting selection: ", smCoordinatesObject$weightStart, " ",
           smCoordinatesObject$weightStop)
    })

    output$smProportionColorHistogram <- renderPlot({
        obj <- smOrderObject
        if (sum(obj$toClust) == 0)
        {showNotification("Select methylation data files to generate the plot."
                          , type="message");NULL}
        else proportionColor(obj, plotHistogram=TRUE, 
                             color=toupper(input$smProportionChoice))
     })

    output$smProportionHistDownload <- downloadHandler(
        filename=function(){
            if (input$filetype == "PNG") return("hist.png")
            if (input$filetype == "SVG") return("hist.svg")
            if (input$filetype == "PDF") return("hist.pdf")
        },
    content=function(file){
        if (input$filetype == "PNG") png(file)
        if (input$filetype == "SVG") svglite::svglite(file)
        if (input$filetype == "PDF") pdf(file)

        proportionColor(smOrderObject, plotHistogram=TRUE,
                       color=toupper(input$smProportionChoice))
        dev.off()
        }
    )
    output$smProportionDataDownload <- downloadHandler(
        filename=function(){
        return("proportion_data.csv")
        },
    content=function(file){
        dat <-  proportionColor(smOrderObject, plotHistogram=FALSE,
                               color=toupper(input$smProportionChoice))
        write.csv(dat, file=file)
        }
    )

    output$smPercentC <- renderPlot({
      obj <- smOrderObject
      if (sum(obj$toClust) == 0)
      {showNotification("Select methylation data files to generate the plot.",
                        type="message");NULL}
      else percent_C(obj, plotPercents=TRUE)
    })

    output$smPercentCPlotDownload <- downloadHandler(
        filename=function(){
            if (input$filetype == "PNG") return("percentC.png")
            if (input$filetype == "SVG") return("percentC.svg")
            if (input$filetype == "PDF") return("percentC.pdf")
        },
    content=function(file){
        if (input$filetype == "PNG") png(file)
        if (input$filetype == "SVG") svglite::svglite(file)
        if (input$filetype == "PDF") pdf(file)

        percent_C(smOrderObject, plotPercents=TRUE)
        dev.off()
        }
    )

    output$smPercentCDataDownload <- downloadHandler(
        filename=function(){
        return("proportion_data.RData")
        },
    content=function(file){
        dat <-  percent_C(smOrderObject, plotPercents=FALSE)
        save(dat, file=file)
        }
     )
    }
