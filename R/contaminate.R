contamination = function(desFile, exprLoc, defMarkers, outDes){
    design = read.table(desFile,header=T,sep='\t')
    allDataPre = read.csv(exprLoc, header = T)
    geneData = allDataPre[,1:3]
    exprData = allDataPre[,4:ncol(allDataPre)]
    exprData = exprData
    design = design[match(colnames(exprData),make.names(design$sampleName),),]
    defMarkers = read.table(defMarkers,header=T,sep='\t')
    defMarkers = lapply(as.list(defMarkers),trimElement,'')

    cindexes = vector(mode = 'list', length= len(defMarkers))
    names(cindexes) = names(defMarkers)

    newExprData = t(scale(t(exprData)))

    for (i in 1:len(defMarkers)){
        mi = tryCatch({apply(newExprData[which(geneData$Gene.Symbol %in% defMarkers[[i]]), !design[,names(cindexes)[i]]],1,min)},
                      error = function(e){min(newExprData[which(geneData$Gene.Symbol %in% defMarkers[[i]]), !design[,names(cindexes)[i]]])})
        ma = tryCatch({apply(newExprData[which(geneData$Gene.Symbol %in% defMarkers[[i]]),!!design[,names(cindexes)[i]]],1,max)},
                      error = function(e){max(newExprData[which(geneData$Gene.Symbol %in% defMarkers[[i]]),!!design[,names(cindexes)[i]]])})

        contaminations = tryCatch({apply((newExprData[which(geneData$Gene.Symbol %in% defMarkers[[i]]),]-mi)/(ma-mi),2,mean)},
                                  error = function(e){(newExprData[which(geneData$Gene.Symbol %in% defMarkers[[i]]),]-mi)/(ma-mi)})
        cindexes[[i]] = contaminations

    }
    cindexes = as.data.frame(cindexes)

    cindexes = cbind(cindexes,apply(cindexes,1,mean))
    names(cindexes)[len(cindexes)] = 'MeanCont'

    newDesign = cbind(design,cindexes)
    newDesign = newDesign[order(as.numeric(rownames(newDesign))),]

    write.table(newDesign, row.names=FALSE,sep = '\t', quote=F ,file = outDes)

}

