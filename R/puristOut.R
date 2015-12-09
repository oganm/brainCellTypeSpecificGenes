
library(ogbox)

allPuristOut = function(genesLoc,lilah=F,regex='*'){
    allGenLocs = list.dirs(genesLoc)
    allGenLocs = allGenLocs[-1]
    allGenLocs = grep(regex,allGenLocs,value=T)
    geneLists = lapply(allGenLocs, puristOut)
    names(geneLists) = basename(allGenLocs)
    return(geneLists)
}


puristOut = function(geneLoc, lilah = F){
    filenames = list.files(geneLoc,include.dirs = FALSE)
    fileContents = lapply(paste0(geneLoc,'/', filenames), function(x){
        tryCatch(
            read.table(x),
            error = function(e){
                NULL
            })
    })
    names(fileContents) = filenames
    
    # empty first file protection
    lengths = sapply(fileContents,len)
    destinedToBeFirst = which.max(lengths>0)
    
    theFirst = fileContents[1]
    fileContents[1] = fileContents[destinedToBeFirst]
    fileContents[destinedToBeFirst] = theFirst
    
    geneList = vector(mode = 'list', length = length(fileContents))
    names(geneList) = names(fileContents)
    
    if (ncol(fileContents[[1]])==3 & lilah == F){
        # this if for a combination of fold change and silhouette coefficient
        for (i in 1:length(fileContents)){
            tempContent = fileContents[[i]][!fileContents[[i]][,1] %in% 
                                                      unlist(sapply((1:len(fileContents))[-i], function(x){
                                                          fileContents[[x]][,1]
                                                          })),]
                
            geneList[[i]] = as.character(tempContent$V1[as.numeric(as.character(tempContent$V3))>0.5
                                                              & as.numeric(as.character(tempContent$V2))>log(10,base=2)])
        }
    }else if (ncol(fileContents[[1]])==3 & lilah == T){
        # this if for lilah's selection method
        for (i in 1:length(fileContents)){
            tempContent = fileContents[[i]][!fileContents[[i]][,1] %in% 
                                                unlist(sapply((1:len(fileContents))[-i], function(x){
                                                    fileContents[[x]][,1]
                                                })),]
            geneList[[i]] = as.character(tempContent$V1[as.numeric(as.character(tempContent$V3))*
                                                              as.numeric(as.character(tempContent$V2))>2])
        }
    } else if (ncol(fileContents[[1]])==1){
        # this is for a mere gene list
        for (i in 1:length(fileContents)){
            geneList[[i]] = as.character(fileContents[[i]]$V1)
        }
    } else if(ncol(fileContents[[1]])==2){
        # this is for selection of percentages from confidence output
        for (i in 1:length(fileContents)){
        geneList[[i]] = as.character(fileContents[[i]]$V1[as.numeric(as.character(fileContents[[i]]$V2))>0.95])
        }
    } else if(ncol(fileContents[[1]])==4){
        # this is for 10 fold changed markers. none of the other collumns matter. just get the genes dammit
        for (i in 1:length(fileContents)){
            geneList[[i]] = as.character(fileContents[[i]]$V1)
        }
    } else {
        stop('What kind of gibberish is this')
    }

    puristList = vector(mode = 'list', length = length(geneList))
    for (i in 1:length(geneList)){
        puristList[[i]] = trimElement(geneList[[i]], unlist(geneList[-i]))
    }
    names(puristList) = names(geneList)
    puristList = lapply(puristList, as.character)
    
    theFirst = fileContents[1]
    fileContents[1] = fileContents[destinedToBeFirst]
    fileContents[destinedToBeFirst] = theFirst
    
    
    return(puristList)
}

# lapply(c('Eno2','Mog'), findGene, list)
findGene = function(gene,list){
    out = lapply(list, function(x){
        findInList(gene,x)
    })
    matches = out[lapply(out,len)>0]
    matches = sapply(1:len(matches), function(i){
        paste0(names(matches[i]),'_', names(list[[names(matches[i])]][matches[