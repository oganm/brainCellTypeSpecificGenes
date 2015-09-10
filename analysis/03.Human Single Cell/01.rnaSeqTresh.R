library(ogbox)
source('R/rnaSeqTresh.R')

rnaExp = read.csv('data/humanRNASeq.csv', row.names=1)
rnaExp = as.matrix(rnaExp)
rnaExp = rnaExp[1:(nrow(rnaExp)-3),] # remove info lines at the bottom

rnaSeqTresh(rnaExp,'analysis/03.Human Single Cell/tresholds',cores=16) -> tresh
tresholds = matrix(rep(1,len(rn(rnaExp))))
rownames(tresholds) = rn(rnaExp)
write.table(tresholds,file='analysis/03.Human Single Cell/noTresh',col.names=F,row.names=T, quote=F)
