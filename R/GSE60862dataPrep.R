GSE60862dataPrep = function(){
    
    library(ogbox)
    library(data.table)
    library(stringr)
    library(sva)
    source('R/mostVariable.R')
    source('R/puristOut.R')
    
    
    humanExp = fread('data/GSE60862_expression')
    humanExp = humanExp[!is.na(humanExp$Gene.Symbol),]
    humanExp = humanExp[humanExp$Gene.Symbol!='',]
    
    # use median expression as elimination treshold
    list[geneData,exprData] = sepExpr(humanExp)
    medExp = median(unlist(exprData))
    humanExp = mostVariable(humanExp,medExp) 
    list[geneData,exprData] = sepExpr(humanExp)
    
    # get rid of file extensions
    colnames(exprData) = str_extract(string=colnames(exprData), pattern='GSM.*?(?=\\.)')
    
    
    softFile = read.design('data/GSE60862_meta')
    softFile = softFile[match(colnames(exprData) ,softFile$GSM),]
    
    # remove outliers detected in the previous step
    outliers = unlist(read.table('analysis//04.Human Coexpression/outliers'))
    softFile = softFile[!softFile$GSM %in% outliers,]
    exprData = exprData[,!colnames(exprData) %in% outliers,with=F]
    
    # remove white matter from the mix
    # exprData = exprData[,!softFile$brainRegion %in% 'white matter',with=F]
    # softFile = softFile[!softFile$brainRegion %in% 'white matter',]
    
    exprData = as.data.frame(exprData)
    rownames(exprData) = geneData$Gene.Symbol
    
    #exprData = exprData[,str_extract(softFile$scanDate,'....-..-..') %in% names(table(str_extract(softFile$scanDate,'....-..-..')))[table(str_extract(softFile$scanDate,'....-..-..'))>1]]
    #softFile = softFile[str_extract(softFile$scanDate,'....-..-..') %in% names(table(str_extract(softFile$scanDate,'....-..-..')))[table(str_extract(softFile$scanDate,'....-..-..'))>1],]
    set.seed(1)
    exprData = ComBat(exprData,batch = kmeans(as.Date(str_extract(softFile$scanDate,'....-..-..')),centers=4)$cluster, mod = model.matrix(~brainRegion,softFile))
    return(list(exprData,softFile))
}
