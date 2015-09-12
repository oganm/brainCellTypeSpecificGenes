library(corpcor)
library(igraph)
library(RBGL)
library(ogbox)
library(data.table)
library(stringr)
source('R/mostVariable.R')
source('R/puristOut.R')
sourceGithub(oganm, toSource, homologene)
homoloGeneTarget = 'data/homologene.tsv'

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
exprData = exprData[,!softFile$brainRegion %in% 'white matter',with=F]
softFile = softFile[!softFile$brainRegion %in% 'white matter',]

rownames(exprData) = geneData$Gene.Symbol


# handcrafting sets for coexpression
sets = list(all = rep(T,ncol(exprData)),
            cortex = softFile$brainRegion %in% c('frontal cortex', 'occipital cortex', 'temporal cortex'),
            substantiaNigra = softFile$brainRegion %in% 'substantia nigra',
            cerebellum = softFile$brainRegion %in% 'cerebellar cortex')


# get all genes from all regions
allGenes = allPuristOut('analysis/01.Gene Selection/FinalGenes/PyramidalDeep')
allGenes = lapply(allGenes,lapply,function(x){
    mouse2human(x)$humanGene
})
# allGenes = unlist(allGenes,recursive=F)

# handcrafting gene sets
geneSets = list(all = allGenes[[1]][c('Astrocyte', 'Microglia','Oligo')],
                cortex = allGenes[['Cortex']],
                substantiaNigra = allGenes$Midbrain,
                cerebellum = allGenes$Cerebellum)




for (i in 1:len(sets)){
    setExpr = exprData[,sets[[i]],with=F]
    rownames(setExpr) = rownames(exprData)
    coexpNetwork(exprData=setExpr,genes=geneSets[[i]],
                 plotOut=paste0('analysis//04.Human Coexpression/coexpPlots/',names(sets)[i]),
                 dataOut=paste0('analysis//04.Human Coexpression/coexpData/',names(sets)[i]))
}

purge()

files = list.files('analysis/04.Human Coexpression/coexpData',full.names=T,recursive=T)
files = files[grepl(pattern='RData',files)]

coexpData = lapply(files,function(x){
    load(x)
    return(list(coexpressions=lapply(coexpressions,sort),
                daCorvec=sort(daCorvec),
                subFolderName=subFolderName,
                plotOut=plotOut,
                regionsIn=regionsIn))
})

coexpBoxViol(coexpData)

