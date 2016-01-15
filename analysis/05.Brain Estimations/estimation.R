library(ogbox)
source('R/estimate.R')
source('R/puristOut.R')
require('reshape2')
library(ggplot2)
source('R/mostVariable.R')
library(data.table)
library(sva)
library(dplyr)
# for human regions ------

source('R/GSE60862dataPrep.R')
library(homologene)

list[exprData,softFile] = GSE60862dataPrep()

exprData  = data.frame(Gene.Symbol = rownames(exprData), exprData)

groups = softFile$brainRegion

# overall --------
genes = puristOut('analysis/01.Gene Selection/FinalGenes/PyramidalDeep/All/')
genes = genes[sapply(genes,len)>2]


fullEstimate(exprData,
             genes=genes,
             geneColName="Gene.Symbol",
             groups=groups,
             outDir='analysis//05.Brain Estimations/estimations/overall/',
             seekConsensus=T,
             groupRotations=T,
             outlierSampleRemove=F,
             controlBased=NULL,
             comparisons = NULL,
             estimateFile = 'analysis//05.Brain Estimations/estimations/overall/estimations')

# cortex white -------
cortex_white = exprData[, c(T, groups %in% c('frontal cortex', 'white matter'))]
groupsCortex_white = groups[groups %in% c('frontal cortex', 'white matter')]

genes = puristOut('analysis/01.Gene Selection/FinalGenes/PyramidalDeep/Cortex/')
genes = genes[sapply(genes,len)>2]

fullEstimate(cortex_white,
             genes=genes,
             geneColName="Gene.Symbol",
             groups=groupsCortex_white,
             outDir='analysis//05.Brain Estimations/estimations/cortex_white/',
             seekConsensus=T,
             groupRotations=T,
             outlierSampleRemove=F,
             controlBased=NULL,
             comparisons = 'all',
             estimateFile = 'analysis//05.Brain Estimations/estimations/cortex_white/estimations')

# cortex general -----
cortexGeneral = exprData[, c(T, grepl(pattern='(?<!cerebellar )cortex',groups,perl=T))]
groupsCortexGeneral = groups[grepl(pattern='(?<!cerebellar )cortex',groups,perl=T)]

genes = puristOut('analysis/01.Gene Selection/FinalGenes/PyramidalDeep/Cortex/')
genes = genes[sapply(genes,len)>2]

fullEstimate(cortexGeneral,
             genes=genes,
             geneColName="Gene.Symbol",
             groups=groupsCortexGeneral,
             outDir='analysis//05.Brain Estimations/estimations/cortexGeneral/',
             seekConsensus=T,
             groupRotations=T,
             outlierSampleRemove=F,
             controlBased=NULL,
             comparisons = 'all',
             estimateFile = 'analysis//05.Brain Estimations/estimations/cortexGeneral/estimations')



list[gene,exp] = 'data/finalExp.csv' %>% read.exp %>% sepExpr

exp[gene$Gene.Symbol %in% 'Tac2',c('GSM1372754','GSM1372755','GSM1372756')]


# huntington disease cortex
exp = read.exp('data/huntington/Frontal Cortex_exp.csv')
exp = mostVariable(exp)
des = read.design('data/huntington/Frontal Cortex_des.tsv')

list[geneDat,expDat]=  exp %>% sepExpr 

expDat = ComBat(expDat, des$scanDate,model.matrix(data = des, ~huntington))

exp = cbind(geneDat,expDat)

genes = puristOut('analysis/01.Gene Selection/FinalGenes/PyramidalDeep/Cortex/')
genes = genes[sapply(genes,len)>2]

fullEstimate(exp,
             genes=genes,
             geneColName="Gene.Symbol",
             groups=des$huntington,
             outDir='analysis//05.Brain Estimations/estimations/huntingonCortex/',
             seekConsensus=T,
             groupRotations=T,
             outlierSampleRemove=F,
             controlBased=NULL,
             comparisons = 'all',
             estimateFile = 'analysis//05.Brain Estimations/estimations/huntingonCortex/estimations')

# huntington cerebellum -----------------------
exp = read.exp('data/huntington/Cerebellum_exp.csv')
exp = mostVariable(exp)
des = read.design('data/huntington/Cerebellum_des.tsv')

list[geneDat,expDat]=  exp %>% sepExpr 

expDat = ComBat(expDat, des$scanDate,model.matrix(data = des, ~huntington))

exp = cbind(geneDat,expDat)

genes = puristOut('analysis/01.Gene Selection/FinalGenes/PyramidalDeep/Cerebellum//')
genes = genes[sapply(genes,len)>2]

fullEstimate(exp,
             genes=genes,
             geneColName="Gene.Symbol",
             groups=des$huntington,
             outDir='analysis//05.Brain Estimations/estimations/huntingonCereb/',
             seekConsensus=T,
             groupRotations=T,
             outlierSampleRemove=F,
             controlBased=NULL,
             comparisons = 'all',
             estimateFile = 'analysis//05.Brain Estimations/estimations/huntingonCereb/estimations')


# parkinson substantia nigra -----
exp = read.exp('data/GSE7621_parkinsonsExp.csv')
list[geneDat, expDat] = sepExpr(exp)
des = read.design('data/GSE7621_parkinsonsMeta.tsv')
des %>% mutate(parkinson = grepl('PD',title)) %>% 
    mutate(female = grepl('female',Characteristic)) -> des

expDatFemale = expDat[des$female==T]
groupsFemale = des$parkinson[des$female==T]

expDatMale = expDat[des$female==F]
groupsMale = des$parkinson[des$female==F]

expDatHealty = expDat[des$parkinson==F]
groupsHealty = des$female[des$parkinson==F]

rownames(expDatHealty) = geneDat$Gene.Symbol
expDatHealty %>% as.matrix -> expDatHealtyHeat
colnames(expDatHealtyHeat) = groupsHealty
heatmap.2(expDatHealtyHeat[geneDat$Gene.Symbol %in% (genes$Dopaminergic %>% mouse2human)$humanGene,] %>% as.matrix,scale='row',trace='none')

#exp %>% filter(Gene.Symbol %in% c(mouse2human(genes$Dopaminergic)$humanGene,'XIST')) %>% 
#    select_( .dots  = c(des %>% filter(parkinson) %>% select(GSM) %>% unlist %>% paste0('.cel')),'Gene.Symbol') ->bok




genes = puristOut('analysis/01.Gene Selection/FinalGenes/PyramidalDeep/Midbrain/')
genes = genes[sapply(genes,len)>2]


fullEstimate(exp,
             genes=genes,
             geneColName="Gene.Symbol",
             groups=des$parkinson,
             outDir='analysis//05.Brain Estimations/estimations/parkinson/',
             seekConsensus=T,
             groupRotations=T,
             outlierSampleRemove=F,
             controlBased=FALSE,
             comparisons = 'all',
             estimateFile = 'analysis//05.Brain Estimations/estimations/parkinson/estimations')

fullEstimate(cbind(geneDat,expDatFemale),
             genes=genes,
             geneColName="Gene.Symbol",
             groups=groupsFemale,
             outDir='analysis//05.Brain Estimations/estimations/parkinsonFemale/',
             seekConsensus=T,
             groupRotations=T,
             outlierSampleRemove=F,
             controlBased=FALSE,
             comparisons = 'all',
             estimateFile = 'analysis//05.Brain Estimations/estimations/parkinsonFemale/estimations')


fullEstimate(cbind(geneDat,expDatHealty),
             genes=genes,
             geneColName="Gene.Symbol",
             groups=groupsHealty,
             outDir='analysis//05.Brain Estimations/estimations/parkinsonMalevsF/',
             seekConsensus=T,
             groupRotations=T,
             outlierSampleRemove=F,
             controlBased=FALSE,
             comparisons = 'all',
             estimateFile = 'analysis//05.Brain Estimations/estimations/parkinsonMalevsF/estimations')


fullEstimate(cbind(geneDat,expDatMale),
             genes=genes,
             geneColName="Gene.Symbol",
             groups=groupsFemale,
             outDir='analysis//05.Brain Estimations/estimations/parkinsonFemale/',
             seekConsensus=T,
             groupRotations=T,
             outlierSampleRemove=F,
             controlBased=NULL,
             comparisons = 'all',
             estimateFile = 'analysis//05.Brain Estimations/estimations/parkinsonFemale/estimations')



hede = cellTypeEstimate(exprData=exp,
                 genes=genes,
                 geneColName='Gene.Symbol',
                 outlierSampleRemove=F,
                 groups=des$parkinson,
                 controlBased= NULL,
                 seekConsensus = T,
                 PC = 1)

sapply(hede$estimates, function(x){
    wilcox.test(x[!hede$groups$Dopaminergic] %>% unlist,x[hede$groups$Dopaminergic] %>% unlist)$p.value
}) %>% p.adjust(method='fdr')


# female male substantia nigra in big dataset -------------
source('R/GSE60862dataPrep.R')
sourceGithub(oganm, toSource, homologene)
homoloGeneTarget = 'data/homologene.tsv'

list[exprData,softFile] = GSE60862dataPrep()

exprData  = data.frame(Gene.Symbol = rownames(exprData), exprData)

list[gene,exp] = exprData %>% sepExpr
exp = exp[,softFile$brainRegion %in% 'substantia nigra']
softFile = softFile[softFile$brainRegion %in% 'substantia nigra',]

groups = softFile$sex

genes = puristOut('analysis/01.Gene Selection/FinalGenes/PyramidalDeep/Midbrain/')
genes = genes[sapply(genes,len)>2]


fullEstimate(cbind(gene,exp),
             genes=genes,
             geneColName="Gene.Symbol",
             groups=groups,
             outDir='analysis//05.Brain Estimations/estimations/substantiaNigraMvsF/',
             seekConsensus=T,
             groupRotations=T,
             outlierSampleRemove=F,
             controlBased=NULL,
             comparisons = 'all',
             estimateFile = 'analysis//05.Brain Estimations/estimations/substantiaNigraMvsF/estimations')

