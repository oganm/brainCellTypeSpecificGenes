
# tempGenes = c('Ank1',
#               'Cox6a2',
#               'Gabrd',
#               'Kcnh2',
#               'Kcnmb2',
#               'Lpl',
#               'Mme',
#               'Pvalb',
#               'Tac1') %>% mouse2human %$% humanGene
# tempGenes = tempGenes[tempGenes %in% rownames(setExpr)]
# trial part
#######################################

library(ogbox)
library(data.table)
library(stringr)
library(corpcor)
library(magrittr)
source('R/mostVariable.R')
source('R/puristOut.R')
source('R/selectRandom.R')

library(homologene)

source('R/GSE60862dataPrep.R')
list[exprData,softFile] = GSE60862dataPrep()


etienneCortex = read.exp('data/etienneCortex/GSE71620_exp.csv')
etienneDesign = read.design('data/etienneCortex/GSE71620_design', comment.char='')

list[etienneGenes,etienneCortex] = etienneCortex %>% sepExpr
rownames(etienneCortex) = etienneGenes %>% unlist

load("data/stanleyData/StanleyData.RData")

stanleyOut = function(stanleyStud){
    stud = stanleyStud$aned_good
    meta = stanleyStud$Metadata
    list[stanleyGenes,stanleyStud] = stud %>%
        mostVariable(.,genes='GeneSymbol',treshold=.[,-1] %>% unlist %>% median) %>% 
        sepExpr
    
    rownames(stanleyStud) = stanleyGenes$GeneSymbol
    
    return(list(stanleyStud,meta))
    
}
list[stanleyStud1,stanleyMeta1] = stanleyOut(studyFinal$study1)
list[stanleyStud3,stanleyMeta3] = stanleyOut(studyFinal$study3)
list[stanleyStud5,stanleyMeta5] = stanleyOut(studyFinal$study5)
list[stanleyStud7,stanleyMeta7] = stanleyOut(studyFinal$study7)

# handcrafting sets for coexpression
expObject = list(all = 'exprData',
                 cortex = 'exprData',
                 substantiaNigra= 'exprData',
                 cerebellum = 'exprData',
                 hippocampus = 'exprData',
                 thalamus = 'exprData',
                 etienneCortex = 'etienneCortex',
                 stanley1Cortex = 'stanleyStud1',
                 stanley3Cortex = 'stanleyStud3',
                 stanleyStud5 = 'stanleyStud5',
                 stanleyStud7 = 'stanleyStud7')


sets = list(all = rep(T,ncol(exprData)),
            cortex = softFile$brainRegion %in% c('frontal cortex', 
                                                 'occipital cortex',
                                                 'temporal cortex'
            ),
            substantiaNigra = softFile$brainRegion %in% 'substantia nigra',
            cerebellum = softFile$brainRegion %in% 'cerebellar cortex',
            hippocampus = softFile$brainRegion %in% 'hippocampus',
            thalamus = softFile$brainRegion %in% 'thalamus',
            etienneCortex = rep(T,ncol(etienneCortex)),
            stanley1Cortex = stanleyMeta1$Profile %in% 'Cont',
            stanley3Cortex = stanleyMeta3$Profile %in% 'Cont',
            stanley5Cortex = stanleyMeta5$Profile %in% 'Cont',
            stanley7Cortex = stanleyMeta7$Profile %in% 'Cont')


# get all genes from all regions
allGenes = allPuristOut('analysis/01.Gene Selection/FinalGenes/PyramidalDeep')
allGenes = lapply(allGenes,lapply,function(x){
    mouse2human(x)$humanGene
})

# remove genes found by Lilah
allGenes$Cortex$GabaPV = allGenes$Cortex$GabaPV [!allGenes$Cortex$GabaPV  %in% 'LPL']
allGenes$Cortex$GabaVIPReln = allGenes$Cortex$GabaVIPReln[!allGenes$Cortex$GabaVIPReln %in% 'HTR3A']



# allGenes = unlist(allGenes,recursive=F)

# handcrafting gene sets
geneSets = list(all = allGenes[[1]][c('Astrocyte', 'Microglia','Oligo')],
                #all = allGenes[[1]],
                cortex = allGenes[['Cortex']][(allGenes[['Cortex']] %>% sapply(len))>3] %>% {.[!names(.) %in% c('Astrocyte', 'Microglia','Oligo','Microglia_activation','Microglia_deactivation','Microglia_regionCorrected')]},
                substantiaNigra = allGenes$Midbrain %>% {.[!names(.) %in% c('Astrocyte', 'Microglia','Oligo','Microglia_activation','Microglia_deactivation','Microglia_regionCorrected')]},
                cerebellum = allGenes$Cerebellum %>% {.[!names(.) %in% c('Astrocyte', 'Microglia','Oligo','Microglia_activation','Microglia_deactivation','Microglia_regionCorrected')]},
                hippocampus = allGenes$Hippocampus %>% {.[!names(.) %in% c('Astrocyte', 'Microglia','Oligo','Microglia_activation','Microglia_deactivation','Microglia_regionCorrected')]},
                thalamus = allGenes$Thalamus %>% {.[!names(.) %in% c('Astrocyte', 'Microglia','Oligo','Microglia_activation','Microglia_deactivation','Microglia_regionCorrected')]},
                etienneCortex = allGenes[['Cortex']][(allGenes[['Cortex']] %>% sapply(len))>3] %>% {.[!names(.) %in% c('Astrocyte', 'Microglia','Oligo','Microglia_activation','Microglia_deactivation','Microglia_regionCorrected')]},
                stanleyCortex1 = allGenes[['Cortex']][(allGenes[['Cortex']] %>% sapply(len))>3] %>% {.[!names(.) %in% c('Astrocyte', 'Microglia','Oligo','Microglia_activation','Microglia_deactivation','Microglia_regionCorrected')]},
                stanleyCortex3 = allGenes[['Cortex']][(allGenes[['Cortex']] %>% sapply(len))>3] %>% {.[!names(.) %in% c('Astrocyte', 'Microglia','Oligo','Microglia_activation','Microglia_deactivation','Microglia_regionCorrected')]},
                stanleyCortex5 = allGenes[['Cortex']][(allGenes[['Cortex']] %>% sapply(len))>3] %>% {.[!names(.) %in% c('Astrocyte', 'Microglia','Oligo','Microglia_activation','Microglia_deactivation','Microglia_regionCorrected')]},
                stanleyCortex7 = allGenes[['Cortex']][(allGenes[['Cortex']] %>% sapply(len))>3] %>% {.[!names(.) %in% c('Astrocyte', 'Microglia','Oligo','Microglia_activation','Microglia_deactivation','Microglia_regionCorrected')]}
                )

# actual analysis

dupResolve = T
ps = lapply(1:len(sets), function(i){
    print(i)
    # take the relevant expression profile data
    expression = teval(expObject[[i]])
    # subset the human expression data to take in only samples from specific regions  
    setExpr =  expression[,sets[[i]]]
    rownames(setExpr) = rownames(expression)
    genes = geneSets[[i]]
    genes %<>% lapply(function(geneSub){geneSub = geneSub[geneSub %in% rn(setExpr)]})
    #expression[genes$GabaPV,] %>%t %>% cor
    #pearCor = setExpr[geneSets[[i]]$GabaPV,] %>% t %>% cor(method='spearman')
    # pearCor[geneSets[[i]]$GabaPV,geneSets[[i]]$GabaPV] %>% sm2vec %>% density %>% plot   
    #pearCor[tempGenes,tempGenes] %>% sm2vec %>% density %>% plot   
    medianExps = setExpr %>% apply(1,median)
    simuGenes = genes %>% lapply(function(geneSub){sapply(geneSub,selectRandom,500, medianExps)})
    
    # resolve duplicates
    if (dupResolve==T){
        simuGenes = lapply(1:len(simuGenes),function(j){
            simuGenes = simuGenes[[j]]
            while (any(apply(apply(simuGenes,1,duplicated),2,any))){
                print(paste("had to resolve equality",names(genes)[j],'in',
                            sum(apply(apply(simuGenes,1,duplicated),2,any))))
                simuGenes[apply(apply(simuGenes,1,duplicated),2,any),] = 
                    t( apply(simuGenes[apply(apply(simuGenes,1,duplicated),2,any),,drop=F],1, function(x){
                        x = sample(x,len(x), replace = F)
                        x[duplicated(x)] = sapply(1:sum(duplicated(x)), function(y){
                            selectRandom(x[duplicated(x)][y],n=1, criteriaValue = medianExps,
                                         invalids=x[!x %in% x[duplicated(x)][y]])
                        })
                        return(x)
                    }))
            }
            return(simuGenes)
        })
        names(simuGenes) = names(genes)
    }
    
    # create the corr matrix from only the necesarry genes. hopefully will be slightly faster (it indeed is faster)
    realCors = genes %>% lapply(function(gene){
        tempCor = setExpr[gene,] %>% t %>% cor %>% sm2vec #%>% median
        #corMat[gene,gene] %>% sm2vec
    })
    
    simuCoexp = simuGenes %>% lapply(function(simuGene){
        simuGene %>% apply(1, function(x){
            tempCor = setExpr[x,] %>% t %>% cor %>% sm2vec # %>% median
        })
    })
    
    ps = sapply(1:len(simuCoexp), function(j){
        print('dunnit')
       wilcox.test(simuCoexp[[j]] %>% as.vector, realCors[[j]], alternative = 'less')$p.value
    })
    
    ps[!sapply(genes,len)>2] = NA
    
#     
#     ps = sapply(1:len(realCors), function(j){
#         1-ecdf(simuCoexp[[j]])(realCors[j])
#     })
    names(ps) = names(realCors)
    return(data.frame(ps,geneCounts = genes %>% sapply(len)))
})

names(ps) = names(geneSets)

pAdjusted  = relist(flesh=p.adjust(ps %>% lapply(function(x){x$ps} %>% unlist),
                                   method='fdr'),
                    skeleton=ps %>% lapply(function(x){x$ps}))

ps = mapply(function(ps,pAdjusted){
    ps$ps = pAdjusted
    return(ps)},
    ps,pAdjusted,SIMPLIFY =FALSE)


lapply(1:len(ps),function(i){    
    write.table(ps[[i]] %>% as.data.frame,col.names=F, sep= '\t',quote=F,file=paste0('analysis//04.Human Coexpression/bootstrapPs/',names(ps)[i]))
})


