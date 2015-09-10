
library(ogbox)
source('R/rnaSeqTresh.R')

# load rna seq data--------------

rnaSeq = read.table('data/mouseRNASeq_Zeisel 2015.txt', sep= '\t', comment.char= "",stringsAsFactors=F)
rnaMeta = rnaSeq[1:10,3:ncol(rnaSeq)]
rnaMeta = as.data.frame(t(rnaMeta))
colnames(rnaMeta) = rnaSeq[1:10,2]

rnaExp = rnaSeq[12:nrow(rnaSeq),3:ncol(rnaSeq)]
rnaExp = apply(rnaExp,2,as.numeric)
rnaExp = matrix(unlist(rnaExp), nrow = nrow(rnaExp))
rownames(rnaExp) = rnaSeq[12:nrow(rnaSeq),1]
rnaCelIDs = as.numeric(as.character(rnaSeq[12:nrow(rnaSeq), 2]))
rnaExp = rnaExp[,rnaMeta$tissue %in% 'sscortex']

# calculate regular tresholds
rnaSeqTresh(rnaExp,'analysis/02.Mouse Single Cell/tresholds',cores=16) -> tresh

# trial treshold of 1 ------------
rn(rnaExp)

tresholds = matrix(rep(1,len(rn(rnaExp))))
rownames(tresholds) = rn(rnaExp)
write.table(tresholds,file='analysis/02.Mouse Single Cell/noTresh',col.names=F,row.names=T, quote=F)

