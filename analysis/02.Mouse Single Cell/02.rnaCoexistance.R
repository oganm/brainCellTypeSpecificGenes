library(gplots)
library(ogbox)
source('R/puristOut.R')
source('R/rnaCoexist.R')

# get the file from Linarrson lab's website if this line doesn't work in the future, contact us or the Linnarson lab
rnaSeq = read.table('data/mouseRNASeq_Zeisel 2015.txt', sep= '\t', comment.char= "",stringsAsFactors=F)
rnaMeta = rnaSeq[1:10,3:ncol(rnaSeq)]
rnaMeta = as.data.frame(t(rnaMeta))
colnames(rnaMeta) = rnaSeq[1:10,2]

rnaExp = rnaSeq[12:nrow(rnaSeq),3:ncol(rnaSeq)]
rnaExp = apply(rnaExp,2,as.numeric)
rnaExp = matrix(unlist(rnaExp), nrow = nrow(rnaExp))
rownames(rnaExp) = rnaSeq[12:nrow(rnaSeq),1]
# rnaCelIDs = as.numeric(as.character(rnaSeq[12:nrow(rnaSeq), 2]))
rnaExp = rnaExp[,rnaMeta$tissue %in% 'sscortex']

# remove low expressed ones -----
maximExp = apply(rnaExp,1,max)
rnaExp = rnaExp[maximExp >= 1,]

markerGenes = puristOut('analysis//01.Gene Selection//FinalGenes/PyramidalDeep/Cortex/')

# no treshold run
tresholds = read.table('analysis/02.Mouse Single Cell/noTresh')
tresholds = tresholds[maximExp >= 1,]
for (i in c(0.3333, 0.40, 0.5, 0.6, 0.7, 0.8, 0.9 ,1)){
    rnaCoexist(rnaExp,
               tresholds,
               markerGenes,
               i,
               F,
               paste0('analysis//02.Mouse Single Cell/Output/NoTreshold',i),
               paste0('analysis//02.Mouse Single Cell/Plots/NoTreshold',i))
}
# treshold run 
tresholds = read.table('analysis/02.Mouse Single Cell/tresholds')
tresholds = tresholds[maximExp >= 1,]
rnaCoexist(rnaExp,
           tresholds,
           markerGenes,
           0.33,
           F,
           'analysis//02.Mouse Single Cell/Output/treshold',
           'analysis//02.Mouse Single Cell/Plots/treshold')


# with duplicated genes removed
tresholds = read.table('analysis/02.Mouse Single Cell/noTresh')
tresholds = tresholds[maximExp >= 1,]
for (i in c(0.3333, 0.40, 0.5, 0.6, 0.7, 0.8, 0.9 ,1)){
    rnaCoexist(rnaExp,
               tresholds,
               markerGenes,
               i,
               T,
               paste0('analysis//02.Mouse Single Cell/Output/NoTresholdNoDub',i),
               paste0('analysis//02.Mouse Single Cell/Plots/NoTresholdNoDup',i))
}

tresholds = read.table('analysis/02.Mouse Single Cell/tresholds')
tresholds = tresholds[maximExp >= 1,]
rnaCoexist(rnaExp,
           tresholds,
           markerGenes,
           1/3,
           T,
           'analysis//02.Mouse Single Cell/Output/tresholdNoDup',
           'analysis//02.Mouse Single Cell/Plots/tresholdNoDup')

