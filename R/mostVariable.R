# this version merges cell types into a single sample then looks for a variable.
# specific to our cell type data. not ideal but oh well...
mostVariableCT = function(whichFile,outFile=NULL,selectionNaming, design='data/meltedDesign.tsv'){

    allDataPre = read.csv(whichFile, header = T)
    design = read.design(design)
    
    list[,exprData]= sepExpr(allDataPre)
    
    cellTypes = trimNAs(unique(design[,selectionNaming]))
    
    cellTypeExpr = lapply(cellTypes,function(x){
        apply(exprData[,design[,selectionNaming] %in% x,drop=F],1,mean)
    })
    exprData = as.data.frame(cellTypeExpr)
    # remove the ones with highest expression below 6
    rowmax = apply(exprData, 1, max)
    discludeGenes = (rowmax<6)
    allDataPre = allDataPre[!discludeGenes,]
    exprData = exprData[!discludeGenes,]
    
    # ignore multiple matching probesets while mostVariable selection
    allDataMulti = allDataPre[grepl('[|]',allDataPre$Gene.Symbol),]
    exprData = exprData[!grepl('[|]',allDataPre$Gene.Symbol),]
    allDataPre = allDataPre[!grepl('[|]',allDataPre$Gene.Symbol),]
    
    # you bloody idiot... taken from lila
    decreasingVar = order(apply(exprData,1,var), decreasing = T)
    allDataPre = allDataPre[decreasingVar,]
    allDataPre = allDataPre[!duplicated(allDataPre$Gene.Symbol),]
    
    # add the multiple matching probesets back
    allDataPre = rbind(allDataPre,allDataMulti)
    allDataPre = allDataPre[!allDataPre$Gene.Symbol=='',]
    
    write.csv(allDataPre, file = outFile, row.names=FALSE)
}

# this function is a generic function that looks for the most variable probeset
# of a gene. unlike the previous one, it takes in objects and outputs objects 
mostVariable = function(allDataPre,genes = 'Gene.Symbol'){
    list[,exprData]= sepExpr(allDataPre)
    rowmax = apply(exprData, 1, max)
    discludeGenes = (rowmax<6)
    allDataPre = allDataPre[!discludeGenes,]
    exprData = exprData[!discludeGenes,]
    
    decreasingVar = order(apply(exprData,1,var), decreasing = T)
    allDataPre = allDataPre[decreasingVar,]
    allDataPre = allDataPre[!duplicated(allDataPre[,genes]),]
    allDataPre = allDataPre[!allDataPre[,genes]=='',]
    return(allDataPre)
}
