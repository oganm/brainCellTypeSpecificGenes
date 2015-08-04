require(RCurl)
eval( expr = parse( text = getURL(
    "https://raw.githubusercontent.com/oganm/toSource/master/ogbox.R",
    ssl.verifypeer=FALSE) ))

sexFind = function(input, output, expr){

    allDataPre = read.csv(expr, header = T)
    design = read.table(input,header=T,sep='\t',quote='"')

    list[geneData, exprData] = sepExpr(allDataPre)

    design = design[match(colnames(exprData),make.names(design$sampleName)),]


    Xist = which(geneData$Gene.Symbol %in% 'Xist')

    sex = rep('varys', ncol(exprData))

    sex[exprData[Xist,]>=7] = 'female'
    sex[exprData[Xist,]<7] = 'male'


    newDesign = cbind(design,sex)
    newDesign = newDesign[order(as.numeric(rownames(newDesign))),]

    write.table(newDesign,file=output, quote=F, sep= '\t',row.names=F)

}