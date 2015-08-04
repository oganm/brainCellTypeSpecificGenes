library(ogbox)
sourceGithub(OganM,toSource,'GOAD/goadDifGenes')

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
    effectedGenes = unlist(goadDifGenes(ids=c(22,23),p=0.05,foldChange='2'))
    effectedGenes = effectedGenes[!duplicated(effectedGenes)]
        
    if (!is.null(restDir)){
        fileNames = list.files(restDir, recursive =T )
        fileNames = fileNames[grepl('Microglia',fileNames)]
        foreach (i = fileNames) %dopar% {
            micro = read.table(paste0(restDir,'/',i))
            micro = micro[!toupper(micro$V1) %in% effectedGenes,]
            write.table(micro, quote = F, row.names = F, col.names = F, paste0(restDir,'/',i))
        }
    }
    
    if (!is.null(genelist)){
        return(geneList[!geneList %in% effectedGenes])
    }
}
