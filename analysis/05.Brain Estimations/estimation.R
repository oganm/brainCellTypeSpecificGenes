library(ogbox)
source('R/estimate.R')
source('R/puristOut.R')
require('reshape2')
library(ggplot2)
source('R/mostVariable.R')
library(data.table)
library(sva)
library(dplyr)
library(lme4)
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
exp = read.exp('data/parkinsons/GSE7621_parkinsonsExp.csv')
list[geneDat, expDat] = sepExpr(exp)
des = read.design('data/parkinsons/GSE7621_parkinsonsMeta.tsv')
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
             groups=setNames(c('parkinson\'s','control'), c(T,F))[des$parkinson %>% as.character],
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
             groups=setNames(c('parkinson\'s','control'), c(T,F))[groupsFemale %>% as.character],
             outDir='analysis//05.Brain Estimations/estimations/parkinsonFemale/',
             seekConsensus=F,
             groupRotations=F,
             outlierSampleRemove=F,
             controlBased=FALSE,
             comparisons = 'all',
             estimateFile = 'analysis//05.Brain Estimations/estimations/parkinsonFemale/estimations')


fullEstimate(cbind(geneDat,expDatHealty),
             genes=genes,
             geneColName="Gene.Symbol",
             groups=setNames(c('parkinson\'s','control'), c(T,F))[groupsHealty %>% as.char], 
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
             groups=setNames(c('parkinson\'s','control'), c(T,F))[groupsMale %>% as.char],
             outDir='analysis//05.Brain Estimations/estimations/parkinsonMale/',
             seekConsensus=T,
             groupRotations=T,
             outlierSampleRemove=F,
             controlBased=NULL,
             comparisons = 'all',
             estimateFile = 'analysis//05.Brain Estimations/estimations/parkinsonMale/estimations')



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

# alzheimer's estimations --------------
genes = puristOut('analysis/01.Gene Selection/FinalGenes/PyramidalDeep/Cortex/')
genes = genes[sapply(genes,len)>2]


list[gene,exp] = 'data/alzheimers/GSE36980_alzheimersExp.csv' %>% read.exp %>% sepExpr

colnames(exp) = gsub(pattern= '[.]cel',replacement='', x=colnames(exp))

design = read.design('data/alzheimers/GSE36980_alzheimersMeta.tsv')

expCortex = exp[,design$Tissue %in% 'Frontal cortex']
designCortex = design %>%  filter(Tissue %in% 'Frontal cortex')


# remove outliers
expCortex= expCortex[,-c(24,30)]
designCortex = designCortex[-c(24,30),]

fullEstimate(cbind(gene,expCortex),
             genes=genes,
             geneColName="Gene.Symbol",
             groups=designCortex$alzheimers,
             outDir='analysis//05.Brain Estimations/estimations/alzheimers/frontalCortex',
             seekConsensus=T,
             groupRotations=T,
             outlierSampleRemove=F,
             controlBased=NULL,
             comparisons = 'all',
             estimateFile = 'analysis//05.Brain Estimations/estimations/alzheimers/frontalCortex/estimations')

expTempCortex = exp[,design$Tissue %in% 'Temporal cortex']
designTempCortex = design %>%  filter(Tissue %in% 'Temporal cortex')

fullEstimate(cbind(gene,expTempCortex),
             genes=genes,
             geneColName="Gene.Symbol",
             groups=designTempCortex$alzheimers,
             outDir='analysis//05.Brain Estimations/estimations/alzheimers/temporalCortex',
             seekConsensus=T,
             groupRotations=T,
             outlierSampleRemove=T,
             controlBased=NULL,
             comparisons = 'all',
             estimateFile = 'analysis//05.Brain Estimations/estimations/alzheimers/temporalCortex/estimations')

genesHippo = puristOut('analysis/01.Gene Selection/FinalGenes/PyramidalDeep/Hippocampus//')
genesHippo = genesHippo[sapply(genesHippo,len)>2]
expHippo = exp[,design$Tissue %in% 'Hippocampus']
designHippo = design %>%  filter(Tissue %in% 'Hippocampus')

fullEstimate(cbind(gene,expHippo),
             genes=genesHippo,
             geneColName="Gene.Symbol",
             groups=designHippo$alzheimers,
             outDir='analysis//05.Brain Estimations/estimations/alzheimers/hippocampus',
             seekConsensus=T,
             groupRotations=T,
             outlierSampleRemove=F,
             controlBased=NULL,
             comparisons = 'all',
             estimateFile = 'analysis//05.Brain Estimations/estimations/alzheimers/hippocampus/estimations')

# advanced analysis with linear models. factoring out sex and shit

fcEstimation = read.table('analysis/05.Brain Estimations/estimations/alzheimers/frontalCortex/estimations')

fcFit = lapply(1:(ncol(fcEstimation)-1), function(i){
    data = data.frame(estimation = fcEstimation[,i], 
                      disease = designCortex$alzheimers, 
                      age = designCortex$Age,
                      sex = designCortex$Sex)
    lm(estimation~disease+age+sex,data=data) %>% summary %$% coefficients
    
    anova(lm(estimation~disease,data=data),
          lm(estimation~disease+age,data=data)#,
         # lm(estimation~disease+age+sex,data=data)
         )
    
    
})

names(fcFit) = colnames(fcEstimation)[1:(ncol(fcEstimation)-1)]


tcEstimation = read.table('analysis/05.Brain Estimations/estimations/alzheimers/temporalCortex/estimations')

tcFit = lapply(1:(ncol(tcEstimation)-1), function(i){
    data = data.frame(estimation = tcEstimation[,i], 
                      disease = designTempCortex$alzheimers, 
                      age = designTempCortex$Age,
                      sex = designTempCortex$Sex)
    lm(estimation~disease+age+sex,data=data) %>% summary %$% coefficients
    
})

names(tcFit) = colnames(tcEstimation)[1:(ncol(tcEstimation)-1)]

hipEstimation = read.table('analysis/05.Brain Estimations/estimations/alzheimers/hippocampus//estimations')

hipFit = lapply(1:(ncol(hipEstimation)-1), function(i){
    data = data.frame(estimation = hipEstimation[,i], 
                      disease = designHippo$alzheimers, 
                      age = designHippo$Age,
                      sex = designHippo$Sex)
    lm(estimation~disease+age+sex,data=data) %>% summary %$% coefficients
    
})

names(hipFit) = colnames(hipEstimation)[1:(ncol(hipEstimation)-1)]


# SIV monkey head --------------
SIVexp = read.exp('data/macacca_HIVdementia/GSE2377_exp.csv')
SIVdes = read.design('data/macacca_HIVdementia/GSE2377_meta.tsv')
list[gene,exp] = SIVexp %>% sepExpr

SIVfrontal = SIVdes %>% filter(Tissue == 'frontal lobe')
SIVexpFrontal = exp[,SIVfrontal$GSM]

genes = puristOut('analysis/01.Gene Selection/FinalGenes/PyramidalDeep/Cortex/')

fullEstimate(cbind(gene,SIVexpFrontal),
             genes=genes,
             geneColName="Gene.Symbol",
             groups=SIVfrontal$SIV,
             outDir='analysis//05.Brain Estimations/estimations/SIVDementedMonkey/frontalLobe',
             seekConsensus=T,
             groupRotations=T,
             outlierSampleRemove=F,
             controlBased=NULL,
             comparisons = 'all',
             estimateFile = 'analysis//05.Brain Estimations/estimations/SIVDementedMonkey/frontalLobe/estimations')

SIVfrontal = SIVdes %>% filter(Tissue == 'frontal lobe')
SIVexpFrontal = exp[,SIVfrontal$GSM]

genes = puristOut('analysis/01.Gene Selection/FinalGenes/PyramidalDeep/Cortex/')

fullEstimate(cbind(gene,SIVexpFrontal),
             genes=genes,
             geneColName="Gene.Symbol",
             groups=SIVfrontal$SIV,
             outDir='analysis//05.Brain Estimations/estimations/SIVDementedMonkey/frontalLobe',
             seekConsensus=T,
             groupRotations=T,
             outlierSampleRemove=F,
             controlBased=NULL,
             comparisons = 'all',
             estimateFile = 'analysis//05.Brain Estimations/estimations/SIVDementedMonkey/frontalLobe/estimations')
