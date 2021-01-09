#' runAlign
#'
#' Runs the preprocessing methods for single-molecule data on sequences.
#'
#' @param ref A reference sequence
#' @param fasta A list of sequences
#' @param fastaSubset A vector of indices indicating which 
#'                      sequences to process.
#' @param multicoreParam A MulticoreParam object, 
#'                          used to align sequences in parallel.
#' @param updateProgress Used to add a progress bar to the Shiny app. 
#'                          Should not be used otherwise.
#' @param logFile String indicating where to save a log of the alignment 
#'                  process. If left NULL, no log is saved.
#'
#' @importFrom Biostrings DNAString DNA_ALPHABET reverseComplement
#'              pairwiseAlignment score alignedPattern alignedSubject
#' @importFrom seqinr c2s s2c read.fasta
#' @importFrom BiocParallel bplapply
#' @export
runAlign <- function(ref, fasta, fastaSubset=(seq_len(fasta)),
    multicoreParam=NULL, updateProgress=NULL, logFile=NULL)
{
    fasta <- fasta[fastaSubset]
    refString <- DNAString(toupper(c2s(ref[[1]])))

    logVector <- c("Beginning preprocessing")

    if (is.function(updateProgress)) updateProgress(message="Aligning 
        sequences", value=0.1)
    alignmentOut <- alignSequences(fasta, refString, logVector,
        multicoreParam, updateProgress)

    alignedSeq <- alignmentOut$alignedSeq
    logVector <- alignmentOut$logVector


    if (is.function(updateProgress)) 
        updateProgress(message="Identifying sites", value=0.75)
    ## We want to avoid GCG sites:
    gCSites <- gregexpr("GC",c2s(refString),fixed=TRUE)[[1]] + 1
    cGSites <- gregexpr("CG",c2s(refString),fixed=TRUE)[[1]]

    logVector <- c(logVector, 
        paste("Throwing out", 
            length(which(s2c(paste(refString))[gCSites+1] == "G")) +
            length(which(s2c(paste(refString))[cGSites-1] == "G"))
            , "GCG sites"))

    cGSites <- cGSites[which(s2c(paste(refString))[cGSites-1] != "G")]
    gCSites <- gCSites[which(s2c(paste(refString))[gCSites+1] != "G")]

    if (is.function(updateProgress)) 
        updateProgress(message="Mapping sites", value=0.8)


    if (is.null(multicoreParam))
    {
        gCMap <- lapply(alignedSeq, mapSeq, sites=gCSites)
        cGMap <- lapply(alignedSeq, mapSeq, sites=cGSites)
    } else {
    gCMap <- bplapply(alignedSeq, mapSeq, sites=gCSites, 
        BPPARAM=multicoreParam)
    cGMap <- bplapply(alignedSeq, mapSeq, sites=cGSites, 
        BPPARAM=multicoreParam)
    }

    if (is.function(updateProgress)) 
        updateProgress(message="Preparing matrices", value=0.95)

    saveCG <- data.matrix(do.call(rbind, lapply(cGMap, function(x) (x))))
    saveCG <- cbind(rownames(saveCG), saveCG)

    saveGC <- data.matrix(do.call(rbind, lapply(gCMap, function(x) (x))))
    saveGC <- cbind(rownames(saveGC), saveGC)
    if (is.function(updateProgress)) updateProgress(message="Done", value=1)

    if (!is.null(logFile))
    {
        writeLines(logVector, con=logFile)
    }

    return(list(hcg=saveCG, gch=saveGC))
}


## this handles the alignment of ALL the sequences,
## and returns the alignedSeq object used in the runAlign function

## this needs the logVector, multicoreParam, and updateProgress 
## so that we can continue keeping track of these things
alignSequences <- function(fasta, refString, logVector, multicoreParam=NULL,
    updateProgress=NULL)
{
    ## this creates the substitution matrix for use in alignment
    penaltyMat <- matrix(0,length(DNA_ALPHABET[seq_len(4)]), 
        length(DNA_ALPHABET[seq_len(4)]))
    penaltyMat[seq_len(4), seq_len(4)] <- c(1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1,
        0, 0, 1, 0, 1)
    penaltyMat[penaltyMat == 0] <- -5
    penaltyMat <- cbind(penaltyMat, c(0, 0, 0, 0))
    penaltyMat <- rbind(penaltyMat, c(0, 0, 0, 0, 1))
    rownames(penaltyMat) <- colnames(penaltyMat) <- 
        c(DNA_ALPHABET[seq_len(4)], "N")

    if (is.null(multicoreParam)) seqAlignOut <- lapply(seq_len(fasta), 
        function(i) {

        if (is.function(updateProgress)) 
            updateProgress(message="Aligning seqences",
                detail=paste(i, "/", length(fasta)),
                value=(0.1+ 0.65/length(fasta) * i))
        seqAlign(fasta[[i]], refString, substitutionMatrix=penaltyMat)
        }) else seqAlignOut <- bplapply(seq_len(length(fasta)),
            function(i) 
                seqAlign(fasta[[i]], refString, 
                    substitutionMatrix=penaltyMat)
                , BPPARAM=multicoreParam)
    useSeqs <- sapply(seqAlignOut, function(i) i$u)
    scores <- sapply(seqAlignOut, function (i) i$score)
    maxAligns <- sapply(seqAlignOut, function (i) i$maxAlign)

    scoreCutoffIdx <- which.max(diff(sort(scores))) + 1
    scoreCutoff <- sort(scores)[scoreCutoffIdx]

    goodAlignmentIdxs <- which(scores > scoreCutoff)

    alignedSeq <- lapply(goodAlignmentIdxs, function(i){
        SEQ1=s2c(paste(alignedPattern(useSeqs[[i]])))
        SEQ2=s2c(paste(alignedSubject(useSeqs[[i]])))

        toReplace <- SEQ1[which(SEQ2 == "-")]
        toReplace[toReplace != "C"] <- "."
        toReplace[toReplace == "C"] <- "." #or T?. Leave as "." for now.

        SEQ2[which(SEQ2 == "-")] <- toReplace
        SEQ2 <- SEQ2[which(SEQ1 != "-")]
        if (maxAligns[i] == 1) 
            SEQ2 <- s2c(paste(reverseComplement(DNAString(c2s(SEQ2)))))

        SEQ2
    })
    logVector <- c(logVector, paste("Throwing out", 
        length(useSeqs) - length(goodAlignmentIdxs), "alignments"))

    names(alignedSeq) <- names(fasta)[goodAlignmentIdxs]
    return(list(alignedSeq=alignedSeq, logVector=logVector))


}


## aligns a single read to the reference, returns the useSeq string. 
## Alignment is finished in the alignSequences fn
seqAlign <- function(read, refString, substitutionMatrix) {

    fastaString <- DNAString(toupper(c2s(read)))

    alignBB <- pairwiseAlignment(reverseComplement(refString),
        reverseComplement(fastaString),
        type="global-local", gapOpening=8, 
        substitutionMatrix=substitutionMatrix)
    alignAB <- pairwiseAlignment(refString, reverseComplement(fastaString),
        type="global-local", gapOpening=8, 
        substitutionMatrix=substitutionMatrix)
    alignAA <- pairwiseAlignment(refString, fastaString,type="global-local", 
        gapOpening=8, 
        substitutionMatrix=substitutionMatrix)

    maxAlign <- which.max(c(score(alignBB), score(alignAB), score(alignAA)))
    allSeq <- list(alignBB, alignAB, alignAA)
    useSeq <- allSeq[[maxAlign]]

    return(list(u=useSeq, score=score(useSeq), maxAlign=maxAlign))
}

mapSeq <- function(i, sites) {
    editSeq <- i
    editSeq[sites][editSeq[sites] == "."] <- "T"
    editSeq[sites][editSeq[sites] == "T"] <- "-2"
    editSeq[sites][editSeq[sites] == "C"] <- "2"
    editSeq[sites][editSeq[sites] == "G"] <- "."
    editSeq[sites][editSeq[sites] == "A"] <- "."
    editSeq[sites][editSeq[sites] == "N"] <- "." 
    ## we need to make sure that the N sites stay marked with a "."
    missingBP <- which(editSeq == ".")

    sitesTemp <- c(0, sites, length(editSeq)+1)
    for (j in seq_len(length(sitesTemp) - 1)) {
        toFill <- seq(sitesTemp[j] + 1, (sitesTemp[j + 1] - 1))
        s1 <- editSeq[pmax(1, sitesTemp[j])]
        s2 <- editSeq[pmin(length(i), sitesTemp[j + 1])]
        isFirst <- (pmax(1, sitesTemp[j]) == 1)

    if (s1 == "2" & s2 == "2") { fillVec <- 1 }
    else if (s1 == "2" & s2 == "-2") {fillVec <- 0}
    else if (s1 == "-2" & s2 == "2") {fillVec <- 0}
    else if (s1 == "-2" & s2 == "-2") {fillVec <- -1}
    else if (s1 == "." & s2 == "."){fillVec <- "."}
    else if (isFirst & s2 == ".") {fillVec <- "."}
    else fillVec <- 0
    fillVec <- rep(fillVec, length(toFill))
    editSeq[toFill] <- fillVec
    editSeq[intersect(toFill, missingBP)] <- "."
    }

    substringTable <- getContigSubstrings(editSeq)
    longMissing <- which(substringTable$char == "." & substringTable$count > 3)
    shortNonMissing <- which(substringTable$char == "-" 
        & substringTable$count <= 20)
    shortNonMissingSurrounded <- shortNonMissing[(shortNonMissing + 1)
                                    %in% longMissing | 
                                    (shortNonMissing - 1) 
                                    %in% longMissing]
    counts <- as.numeric(substringTable$count)
    for (idx in shortNonMissingSurrounded) 
    ## this part is tricky... we want to change the short non-missing sections
    ## to missing
    {
        if (idx == 1) editSeq[seq_len(counts[idx])] <- "."
        else
        {
            first <- sum(counts[seq_len(idx-1)]) + 1
            last <- first + counts[idx] - 1
            editSeq[first:last] <- "."
        }
    }

    return(editSeq)
}


## we want to be able to get all contiguous substrings of a certain string... 
## in particular one of the editseq strings used above

## i want to return a table
getContigSubstrings <- function(s)
{
    ## we want some sort of table to keep track of the substrings
    substringTable <- data.frame(char="a", count=0)

    s <- ifelse(s == ".", ".", "-") 
    ## use "-" to denote the non-missing ones... 
    ## so we only keep track of missing and non missing
    i <- 1
    while (i <= length(s))
    {
        j <- i
        while(s[j] == s[i] & j <= length(s)) j <- j + 1
        substringTable <- rbind(substringTable, c(s[i], j - i))
        i <- j
    }
    substringTable <- data.frame(char=substringTable$char[-1], 
                                count=as.numeric(substringTable$count[-1]))
    return(substringTable)
}
