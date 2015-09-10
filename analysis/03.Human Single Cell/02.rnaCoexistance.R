library(gplots)
library(ogbox)
source('R/puristOut.R')
source('R/rnaCoexist.R')
sourceGithub(OganM,toSource,homologene)
homoloGeneTarget = 'data/homologene.tsv'

rnaExp = read.table('data/humanRNASeq.csv', sep= ',', comment.char= "",stringsAsFactors=F, row.names=1, header=T)
rnaExp = rnaExp[1:(nrow(rnaExp)-3),]
maximExp = apply(rnaExp,1,max)
rnaExp = rnaExp[maximExp >= 1,]
genes = rn(rnaExp)
rnaExp = apply(rnaExp,2,as.numeric)
rownames(rnaExp) = genes
markerGenes = puristOut('analysis//01.Gene Selection//FinalGenes/PyramidalDeep/Cortex/')

markerGenes = sapply(markerGenes,function(x){ as.character(mouse2human(x)$humanGene)})


# no treshold run ----------------
print('beginning thresholding')
tresholds = read.table('analysis/03.Human Single Cell/noTresh')
tresholds = tresholds[maximExp >= 1,]
for (i in c(0.3333, 0.40, 0.5, 0.6, 0.7, 0.8, 0.9 ,1)){
    print(i)
    rnaCoexist(rnaExp,
               tresholds,
               markerGenes,
               i,
               T,
               paste0('analysis//03.Human Single Cell//Output/NoTreshold',i),
               paste0('analysis//03.Human Single Cell//Plots/NoTreshold',i))
}


# treshold run  ---------
tresholds = read.table('analysis/03.Human Single Cell//tresholds')
tresholds = tresholds[maximExp >= 1,]
for (i in c(0.3333, 0.40, 0.5, 0.6, 0.7, 0.8, 0.9 ,1)){
rnaCoexist(rnaExp,
           tresholds,
           markerGenes,
           i,
           T,
           paste0('analysis//03.Human Single Cell/Output/treshold',i),
           paste0('analysis//03.Human Single Cell/Plots/treshold',i))
}