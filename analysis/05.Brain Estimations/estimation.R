library(ogbox)
source('R/estimate.R')
source('R/puristOut.R')
source('runVars.R')
require('reshape2')
require(ggplot2)
source('R/mostVariable.R')
# for human regions ------
puristList = puristOut('analysis/01.Gene Selection/FinalGenes/PyramidalDeep/Cortex')
# filtering shit for human data
softFile = read.design('data/GSE46706_meta.tsv')
softExpr = read.exp('data/HumanRegionExpr/Mixed/cortex-white')

names(softExpr)[4:len(softExpr)] = sub('_.*','',names(softExpr)[4:len(softExpr)] )
# determine groups based on sample names from soft file
list[softGenes,softExp] = sepExpr(softExpr)
groups = softFile$Region[match(names(softExp), softFile$GSM)]


medExpr = median(unlist(softExp))
keep = apply(softExp,1,function(row){
    return(mean(row[groups %in% unique(groups)[1]])>medExpr & mean(row[groups %in% unique(groups)[2]])>medExpr)
})
softExpr=softExpr[keep,]
mostVarSoft = mostVariable(softExpr,'Gene.Symbol')




fullEstimate(mostVarSoft,
             genes=puristList,
             geneColName="Gene.Symbol",
             groups=groups,
             outDir='Data/Estimates/Cortex-White/',
             seekConsensus=T,
             groupRotations=T,
             outlierSampleRemove=T,
             controlBased=NA)


