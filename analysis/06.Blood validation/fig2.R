library(ogbox)
source('R/puristOut.R')
source('R/estimate.R')
source('R/superImpose.R')
library(gtable)
library(grid)
library(parallel)

theirPred = read.table('data/bloodCellType/lymphomaCounts.tsv',sep='\t',header=T)
gsms = sapply(theirPred$Sample.ID, function(x){gsmFind(GSE='GSE65135', regex=x)})
theirPred$Sample.ID = gsms
names(theirPred) = c('Sample.ID', 'Cell.type', 'RealCounts', 'CibersortPred')

theirPred = theirPred %>% arrange(Cell.type, Sample.ID )


# l11Genes = puristOut('analysis//06.Blood validation//humanGenes//Fold/Relax/Leuk11/')
l11Genes = puristOut('analysis//06.Blood validation//humanGenes//rotSel/Relax/Leuk11//')

l11Genes = l11Genes[c("B Cells",
                      'CD4',
                      'CD8 T cells')]
names(l11Genes) = unique(theirPred$Cell.type)
l11Genes -> l11GenesHuman

expression = read.exp('data//bloodCellType//lymphoma.csv')
names(expression)[2] = 'GeneSym'
list[gene,exp] = sepExpr(expression)
rownames(exp) = gene$GeneSym


humanGenePredictions = superImpose(l11GenesHuman, paste0('analysis//06.Blood validation/lymphomaSuperImpose','/'))


l11Genes = puristOut('analysis//06.Blood validation//mouseGenes//RotSel//Relax/lm11/') %>% lapply(function(x){mouse2human(x)$humanGene})
l11Genes = l11Genes[c("B Cell",
                      'CD4',
                      'CD8')]
names(l11Genes) = c("B Cells",
               'CD4',
               'CD8 T cells')
names(l11Genes) = unique(theirPred$Cell.type)

mouseGenePredictions = superImpose(l11Genes, paste0('analysis//06.Blood validation/lymphomaSuperImposeMouseGenes','/'))

c('/home/omancarci/brainCellTypeSpecificGenes/analysis/06.Blood validation/mouseGenes/Fold/Relax/lm11/B Cell',
  '/home/omancarci/brainCellTypeSpecificGenes/analysis/06.Blood validation/mouseGenes/Fold/Relax/lm11/CD4',
  '/home/omancarci/brainCellTypeSpecificGenes/analysis/06.Blood validation/mouseGenes/Fold/Relax/lm11/CD8') %>% lapply(function(x){
     read.table(x, header=F) %>% arrange(desc(V3))
  }) ->geneQuality


accuracyDistribution = mclapply(1:1000, function(j){
    print(j)
    genes = lapply(names(l11Genes), function(i){
        l11Genes[[i]] %>% sample(l11GenesHuman[[i]] %>% len)
    })
    names(genes) = names(l11Genes)
    superImpose(genes, NULL)
}, mc.cores = 10
)

accuracyDistribution %<>% as.data.frame 
sapply(1:3, function(i){
    plot(density(accuracyDistribution[i,] %>% unlist), main = names(l11Genes)[i])
    abline(v=humanGenePredictions[i])
    abline(v = mouseGenePredictions[i],col='red')
    abline(v= accuracyDistribution[i,] %>% unlist %>% median , col = 'green',lty=2 )
})

apply(1,function(x){plot(density(x))})

accuracyDistribution %>% as.data.frame