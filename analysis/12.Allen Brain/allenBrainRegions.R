library(data.table)
library(dplyr)
library(ogbox)
library(homologene)
source('R/mostVariable.R')
source('R/puristOut.R')
source('R/cellColors.R')
source('R/heatmap.3.R')
source('analysis//12.Allen Brain//regionAnalysis.R')
library(tidyr)
library(cowplot)
library(memoise)

memoMedian = memoise(median)
memoMatrix = memoise(as.matrix)
memoMostVariable = memoise(mostVariable)
# check for variation of estimation within and between human brains.

# load ontology from any one of the samples get cortex samples
ontology = fread('data/allenBrain/H0351.1009/Ontology.csv', data.table=F)
relOntology = ontology %>% filter(grepl('4006',structure_id_path))

# load samples from all and get the common samples
commonSamples = intersectMult( list = 
                   lapply(dirs, function(x){
    samples = fread(paste0(x,'/SampleAnnot.csv'), 
                    data.table=F)
    samples %<>% filter(structure_id %in% relOntology$id)
    return(samples$structure_id)
}))

# load data
dirs = list.dirs('data/allenBrain/')[-1]

humanBrains = lapply(dirs, function(x){
    expression = fread(paste0(x,'/MicroarrayExpression.csv'), 
                       data.table=F)
    list[probe, expression] = sepExpr(expression)

    
    
    # read sample names and probe names
    samples = fread(paste0(x,'/SampleAnnot.csv'), 
                    data.table=F)
    probes = fread(paste0(x,'/Probes.csv'),
                   data.table=F)
    
    
    # only take the cortical samples
    expression =  expression[,samples$structure_id %in% commonSamples]
    samples %<>% filter(structure_id %in% commonSamples)
    
   # expression = lapply(commonSamples,function(z){
#        apply(expression[,samples$structure_id %in% z, drop = F],1,mean)
 #   }) %>% as.data.frame
    
    # colnames(expression) = commonSamples
    
    return(list(expression = expression,samples = samples))
})

probes = fread('data/allenBrain/H0351.1009/Probes.csv',
               data.table=F)

flatBrains = cbind(humanBrains[[1]]$expression,
                   humanBrains[[2]]$expression,
                   humanBrains[[3]]$expression,
                   humanBrains[[4]]$expression,
                   humanBrains[[5]]$expression,
                   humanBrains[[6]]$expression)

design = rbind(humanBrains[[1]]$samples,
               humanBrains[[2]]$samples,
               humanBrains[[3]]$samples,
               humanBrains[[4]]$samples,
               humanBrains[[5]]$samples,
               humanBrains[[6]]$samples)

design$donor = repIndiv(1:6, sapply(humanBrains,function(x){nrow(x$samples)}))

genes = puristOut('analysis/01.Gene Selection/FinalGenes/PyramidalDeep/All//') %>%
    lapply(function(x){mouse2human(x)$humanGene})
genes = genes[sapply(genes, len)>2]


regionAnalysis(flatBrains, design, genes,'analysis/12.Allen Brain/all')


# just cortex ======
corticalRegions = ontology %>% filter(grepl('4009',structure_id_path))

flatBrainsSub = flatBrains[,design$structure_id %in% corticalRegions$id]
designSub  = design[design$structure_id %in% corticalRegions$id,]

genes = puristOut('analysis/01.Gene Selection/FinalGenes/PyramidalDeep/Cortex/') %>% 
    lapply(function(x){mouse2human(x)$humanGene})
genes = genes[sapply(genes, len)>2]

regionAnalysis(flatBrainsSub, designSub, genes,'analysis/12.Allen Brain/cortex')



