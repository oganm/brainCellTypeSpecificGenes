library(gplots)
library(ogbox)
library(homologene)
source('R/puristOut.R')
source('R/rnaCoexist.R')

rnaExp = read.table('data/humanRNASeq.csv', sep= ',', comment.char= "",stringsAsFactors=F, row.names=1, header=T)
rnaExp = rnaExp[1:(nrow(rnaExp)-3),]
maximExp = apply(rnaExp,1,max)
rnaExp = rnaExp[maximExp >= 1,]
genes = rn(rnaExp)
rnaExp = apply(rnaExp,2,as.numeric)
rownames(rnaExp) = genes
markerGenes = puristOut('analysis//01.Gene Selection//FinalGenes/PyramidalDeep/Cortex/')

markerGenes = sapply(markerGenes,function(x){ as.character(mouse2human(x)$humanGene)})

# removal of genes from lilah's findings
markerGenes$GabaPV = markerGenes$GabaPV[!markerGenes$GabaPV %in% 'LPL']
markerGenes$GabaVIPReln = markerGenes$GabaVIPReln[!markerGenes$GabaVIPReln %in% 'HTR3A']

# no treshold run ----------------
print('beginning thresholding')
tresholds = read.table('analysis/03.Human Single Cell/noTresh')
tresholds = tresholds[maximExp >= 1,]


rnaCoexist(rnaExp,
           tresholds,
           markerGenes,
           T,
           paste0('analysis//03.Human Single Cell//Output/NoTreshold'),
           paste0('analysis//03.Human Single Cell//Plots/NoTreshold'),
           cores  = 16)


# treshold run  ---------
tresholds = read.table('analysis/03.Human Single Cell//tresholds')
tresholds = tresholds[maximExp >= 1,]

rnaCoexist(rnaExp,
           tresholds,
           markerGenes,
           T,
           paste0('analysis//03.Human Single Cell//Output/Threshold'),
           paste0('analysis//03.Human Single Cell//Plots/Threshold'),
           cores  = 16)