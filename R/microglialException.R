library(ogbox)
library(stringr)
library(GOADquery)

microglialException = function(restDir=NULL, genelist = NULL ,updateList = F, cores = 1){
    #applies to all files inside a directory recursively
    require(foreach)
    require(doMC)
    require(parallel)# so that I wont fry my laptop
    if (detectCores()<cores){ 
        cores = detectCores()
        print('max cores exceeded')
        print(paste('set core no to',cores))
    }
    registerDoMC(cores)
    # genes effected by old age and LPS stimulation
    
    effectedGenes = read.table('data/GOAD/LPS-STIMULATED MICROGLIA (4HR) vs. CONTROL MICROGLIA.tsv',
                                  header=T,sep='\t',stringsAsFactors = FALSE, quote = '') %>%
        filter(Adjusted.p.value<=0.05) %$% Gene.symbol
    
    activationGenes = read.table('data/GOAD/LPS-STIMULATED MICROGLIA (4HR) vs. CONTROL MICROGLIA.tsv',
                               header=T,sep='\t',stringsAsFactors = FALSE, quote = '') %>%
        filter(Adjusted.p.value<=0.05, Fold.Change >1 ) %$% Gene.symbol
    
    deactivationGenes = read.table('data/GOAD/LPS-STIMULATED MICROGLIA (4HR) vs. CONTROL MICROGLIA.tsv',
                                 header=T,sep='\t',stringsAsFactors = FALSE, quote = '') %>%
        filter(Adjusted.p.value<=0.05, Fold.Change <1 ) %$% Gene.symbol
    
    oldAgeGenes = read.table('data/GOAD/OLD (24M) MICROGLIA vs. YOUNG (5M) MICROGLIA.tsv',
                             header=T,sep='\t',stringsAsFactors = FALSE, quote = '') %>% 
        filter(Adjusted.p.value<=0.05, Fold.Change >1 ) %$% Gene.symbol
  
    cerebellumLow = read.table('data/GOAD/MICROGLIA BRAIN REGIONS: CEREBELLUM vs. OTHER REGIONS.tsv',
                                header=T,sep='\t',stringsAsFactors = FALSE, quote = '') %>% 
        filter(Adjusted.p.value<=0.05, Fold.Change <1.5 ) %$% Gene.symbol
    
    
    cortexLow = read.table('data/GOAD/MICROGLIA BRAIN REGIONS: GREY MATTER (CORTEX) vs. OTHER REGIONS.tsv',
                            header=T,sep='\t',stringsAsFactors = FALSE, quote = '') %>% 
        filter(Adjusted.p.value<=0.05, Fold.Change <1.5 ) %$% Gene.symbol
    
    
    hippocampusLow = read.table('data/GOAD/MICROGLIA BRAIN REGIONS: HIPPOCAMPUS vs. OTHER REGIONS.tsv',
                                 header=T,sep='\t',stringsAsFactors = FALSE, quote = '') %>% 
        filter(Adjusted.p.value<=0.05, Fold.Change <1.5 ) %$% Gene.symbol
    
    
    if (!is.null(restDir)){
        fileNames = list.files(restDir, recursive =T )
        fileNames = fileNames[grepl('Microglia$',fileNames)]
        #for(i in fileNames){
        foreach (i = fileNames) %dopar% {
            micro = read.table(paste0(restDir,'/',i))
            microAll = micro[!toupper(micro$V1) %in% effectedGenes,]
            write.table(microAll, quote = F, row.names = F, col.names = F, paste0(restDir,'/',i))
            
            actiMicro = micro[toupper(micro$V1) %in% activationGenes,]
            write.table(actiMicro, quote = F, row.names = F, col.names = F, paste0(restDir,'/',i,'_activation'))
            
            deactiMicro = micro[toupper(micro$V1) %in% deactivationGenes,]
            write.table(deactiMicro, quote = F, row.names = F, col.names = F, paste0(restDir,'/',i,'_deactivation'))
            
            # create corrected version based on region specific microarray expression data
            if(grepl(pattern='Cortex',x=i)){
                cortexMicro = microAll[!toupper(microAll) %in% cortexLow]
                write.table(cortexMicro, quote = F, row.names = F, col.names = F, paste0(restDir,'/',i,'_regionCorrected'))
            } else if (grepl(pattern='Hippocampus',x=i)){
                hippocampusMicro =  microAll[!toupper(microAll) %in% hippocampusLow]
                write.table(hippocampusMicro, quote = F, row.names = F, col.names = F, paste0(restDir,'/',i,'_regionCorrected'))
            } else if (grepl(pattern ='Cerebellum',x=i)){
                cerebellumMicro =  microAll[!toupper(microAll) %in% cerebellumLow]
                write.table(cerebellumMicro, quote = F, row.names = F, col.names = F, paste0(restDir,'/',i,'_regionCorrected'))
            }
            
        }
    }
    
    # just apply to a single microglia list
    if (!is.null(genelist)){
        return(geneList[!geneList %in% effectedGenes])
    }
}
