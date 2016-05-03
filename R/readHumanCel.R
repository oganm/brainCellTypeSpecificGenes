library(ogbox)
library(oligo)
readHumanCel = function(GSMs, fileOut=NULL, humanDir){
    cels = oligoClasses::list.celfiles(humanDir,listGzipped=T)
    whichCels = sapply(GSMs, function(x){which(grepl(x,cels))})
    sampleCels = cels[whichCels]
    # browser()
    affyRaw = oligo::read.celfiles(paste0(humanDir,'/',sampleCels))
    exonTS <- oligo::rma(affyRaw, target = "core")
    rm(affyRaw)
    featureData(exonTS) <- getNetAffx(exonTS, "transcript")
    #View the features of the obtained data
    # exonTS
    #Extract the expression data
    
    aned = gemmaAnnotOligo(normalized=exonTS,chipFile='data/GemmaAnnots/GPL5175')
    
    #Create the expression file
    if (!is.null(fileOut)){
        write.csv(aned, fileOut, row.names=FALSE)
    }
    invisible(aned)
}


readOligoCel = function(GSMs, platform,fileOut=NULL, celdir){
    cels = oligoClasses::list.celfiles(celdir,listGzipped=T)
    whichCels = sapply(GSMs, function(x){which(grepl(x,cels))})
    sampleCels = cels[whichCels]
    affyRaw = oligo::read.celfiles(paste0(celdir,'/',sampleCels))
    exonTS <- oligo::rma(affyRaw, target = "core")
    rm(affyRaw)
    featureData(exonTS) <- getNetAffx(exonTS, "transcript")
    #View the features of the obtained data
    # exonTS
    #Extract the expression data
    
    aned = gemmaAnnotOligo(normalized=exonTS,chipFile=paste0('data/GemmaAnnots/',platform))
    
    #Create the expression file
    if (!is.null(fileOut)){
        write.csv(aned, fileOut, row.names=FALSE)
    }
    invisible(aned)