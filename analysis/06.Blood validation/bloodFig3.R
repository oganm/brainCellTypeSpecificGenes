library(ogbox)
source('R/puristOut.R')
source('R/estimate.R')
source('R/superImpose.R')
library(ggplot2)
library(gtable)
library(grid)
library(dplyr)
library(grDevices)
library(magrittr)

realCounts = read.table('data//bloodCellType//PBMCcounts.tsv',sep='\t',header=T,row.names=1)
theirPred = read.table('data/bloodCellType/PBMCcibersort.tsv', sep='\t', header=T)
names(theirPred) = c('Sample.ID', 'Cell.type', 'RealCounts', 'CibersortPred')

# l22Genes = puristOut('analysis//06.Blood validation//Fold/Relax/Abreviated.name/')
# l22Genes2 = puristOut('analysis//06.Blood validation//Fold/Relax//')

l22Genes = puristOut('analysis//06.Blood validation//humanGenes/rotSel/Relax/Abreviated.name/')
l22GenesjustPBMC = puristOut('analysis//06.Blood validation//rotSel/Relax/justPBMC/')
l22GenesRelax = puristOut('analysis//06.Blood validation//Fold/Relax/Abreviated.name/')
l22GenesjustPBMCRelax = puristOut('analysis//06.Blood validation/Fold/Relax/justPBMC/')
# deal with mouse gene --------------
l22GenesMouse = puristOut('analysis//06.Blood validation//mouseGenes/RotSel//Relax/lm22/') %>% lapply(function(x){mouse2human(x)$humanGene}) 
l22GenesMouse = l22GenesMouse[c('B cell naive',
                                'B cell memory',
                                'T_CD8',
                                'T_CD4_naive',
                                'T_CD4_memory',
                                'T_CD4_activated',
                                'NK resting',
                                'Monoctye')]
names(l22GenesMouse) = c('Naïve B cells',
                         'Memory B cells',
                         'CD8 T cells',
                         'CD4 naïve T cells',
                         'CD4 memory T cells-', 
                         'CD4 memory T cells+',
                         'NK cells-', # subject to change be careful
                         'Monos')
# back to reality ------------
geneLists = list(l22Genes, l22GenesjustPBMC, l22GenesRelax, l22GenesjustPBMCRelax,l22GenesMouse)
names(geneLists) =  c('l22Genes', 'l22GenesjustPBMC', 'l22GenesRelax', 'l22GenesjustPBMCRelax','l22GenesMouse')

# renaming for matching purposes
theirPred$Cell.type = make.names(theirPred$Cell.type)
theirPred$Cell.type[theirPred$Cell.type == 'Naïve.CD4.T.cells'] = 'CD4.naïve.T.cells'
theirPred$Cell.type[theirPred$Cell.type == 'Resting.memory.CD4.T.cells'] = 'CD4.memory.T.cells.resting'

geneLists = lapply(geneLists,function(x){
    x=x[c('Naïve B cells',
               'Memory B cells',
               'CD8 T cells',
               'CD4 naïve T cells',
               'CD4 memory T cells-', 
               'CD4 memory T cells+',
               'NK cells-', # subject to change be careful
               'Monos')]
    names(x) =  unique(theirPred$Cell.type)
    return(x)
})


expression = read.exp('data//bloodCellType//PBMCs.csv')
list[gene,exp] = sepExpr(expression)
rownames(exp) = gene$GeneSym

lapply(1:len(geneLists), function(i){
    superImpose(geneLists[[i]], paste0('analysis//06.Blood validation/superImpose',names(geneLists)[i],'/'))
})



# try to improve correlation by adding correlated genes. see what happens -------------

estimates = cellTypeEstimate(exprData = expression, 
                             genes= geneLists$l22GenesRelax,
                             geneColName = 'GeneSym',
                             outlierSampleRemove = F,
                             synonymTaxID = NULL,
                             geneTransform = NULL, 
                             groups = NA,
                             controlBased = NULL, 
                             tableOut = NULL,
                             indivGenePlot = NULL,
                             seekConsensus = F,
                             plotType = NULL,
                             PC = 1)
list[genes,exp]=sepExpr(expression)
estimCors = cor(estimates$estimates %>% as.data.frame, exp %>% t,method='spearman') %>% 
    t %>%
    as.data.frame %>% 
    mutate(genes = genes$GeneSym)

newGeneList = lapply(1:len(geneLists$l22GenesRelax),function(i){
    estimCors %>% 
        arrange_(names(estimCors[,i,drop=F])) %>% 
        filter(!genes %in% geneLists$l22GenesRelax[[i]]) %>%
        tail(n=len(geneLists$l22GenesRelax[[i]])/2 %>% select_('genes',names(estimCors[,i,drop=F])) %T>% print %>% 
        select(genes) %>%
        unlist %>% c(geneLists$l22GenesRelax[[i]])
})


names(newGeneList) = names(geneLists$l22GenesRelax)

newEstimates = cellTypeEstimate(exprData = expression, 
                             genes= newGeneList,
                             geneColName = 'GeneSym',
                             outlierSampleRemove = F,
                             synonymTaxID = NULL,
                             geneTransform = NULL, 
                             groups = NA,
                             controlBased = NULL, 
                             tableOut = NULL,
                             indivGenePlot = NULL,
                             seekConsensus = F,
                             plotType = NULL,
                             PC = 1)

countCors = data.frame(old=c(0),new=c(0))
for (i in 1:len(newGeneList)){
    print(names(newGeneList)[i])
hede = theirPred %>%
    filter(Cell.type == names(estimates$estimates[i])) %>%
    select(RealCounts) %>% cor ( cbind(estimates$estimates[[i]], newEstimates$estimates[[i]]),method='spearman') %>% print
hede = hede %>% as.data.frame
names(hede) = c('old','new')
countCors = rbind(countCors,hede %>% as.data.frame)
}
countCors = countCors[2:nrow(countCors),]
rownames(countCors) = names(newGeneList)
knitr::kable(countCors[2:nrow(countCors),])

# merge cell types and retry
realCounts = read.table('data//bloodCellType//PBMCcounts.tsv',sep='\t',header=T,row.names=1)

simpleCounts = data.frame('B Cells' = realCounts$Naïve.B.cells + realCounts$Memory.B.cells,
                          'CD4' = realCounts$CD4.naïve.T.cells + realCounts$CD4.memory.T.cells.resting + realCounts$CD4.memory.T.cells.activated)

bloodDes = read.design('analysis/06.Blood validation/BloodCells.tsv')

translation = bloodDes %>% select(Abreviated.name, Leuk11) %>% unique

theirPred %>% mutate(Cell.Type.)

