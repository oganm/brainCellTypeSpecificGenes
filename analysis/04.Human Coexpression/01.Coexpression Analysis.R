
library(ogbox)
library(data.table)
library(stringr)
source('R/mostVariable.R')
source('R/puristOut.R')
source('R/coexpBoxViol.R')
source('R/coexpNetwork.R')


library(homologene)

source('R/GSE60862dataPrep.R')
list[exprData,softFile] = GSE60862dataPrep()


# remove white matter from the mix
exprData = exprData[,!softFile$brainRegion %in% 'white matter']
softFile = softFile[!softFile$brainRegion %in% 'white matter',]


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

# remove genes found by Lilah
allGenes$Cortex$GabaPV = allGenes$Cortex$GabaPV [!allGenes$Cortex$GabaPV  %in% 'LPL']
allGenes$Cortex$GabaVIPReln = allGenes$Cortex$GabaVIPReln[!allGenes$Cortex$GabaVIPReln %in% 'HTR3A']



# allGenes = unlist(allGenes,recursive=F)

# handcrafting gene sets
geneSets = list(all = allGenes[[1]][c('Astrocyte', 'Microglia','Oligo')],
                #all = allGenes[[1]],
                cortex = allGenes[['Cortex']],
                substantiaNigra = allGenes$Midbrain,
                cerebellum = allGenes$Cerebellum)


for (i in 1:len(sets)){
    setExpr = exprData[,sets[[i]]]
    rownames(setExpr) = rownames(exprData)
    coexpNetwork(exprData=setExpr,genes=geneSets[[i]],
                 plotOut=paste0('analysis//04.Human Coexpression/coexpPlots/',names(sets)[i]),
                 dataOut=paste0('analysis//04.Human Coexpression/coexpData/',names(sets)[i]))
}


files = list.files('analysis/04.Human Coexpression/coexpData',full.names=T,recursive=T)
files = files[grepl(pattern='RData',files)]
files = files[grepl(pattern='spear',files)]

coexpData = lapply(files,function(x){
    load(x)
    return(list(coexpressions=lapply(coexpressions,sort),
                daCorvec=sort(daCorvec),
                subFolderName=subFolderName,
                plotOut=plotOut,
                regionsIn=regionsIn))
})

coexpBoxViol(coexpData)
