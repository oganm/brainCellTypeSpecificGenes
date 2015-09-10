# requires internet connection to get n if it is not provided
# a general puprose soft parser that will break if people placed columns in weird places
library(RCurl)
library(ogbox)
softParser = function(softFile, # file to read
                       mergeFrame = c('intersect', 'union'), # union not implemented
                       n=NULL # number of samples in the file
                       ){
    con  = file(softFile, open = "r")
    oneLine = readLines(con, n = 1, warn = FALSE)
    
    # if n isn't provided grab it online
    if (is.null(n)){
        while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0){
             if (grepl('\\^SERIES', oneLine)){
                 GSE = strsplit(oneLine,' = ')[[1]][2]
                 n = len(gsmFind(GSE))
                 break
             }
        }
    }
    
    i=0
    sampleData = vector(mode ='list', length = n)
    # get relevant information
    while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
        if (grepl('\\^SAMPLE', oneLine)){
            sampLines = vector(mode = 'character', length=0)
            while (oneLine != '!sample_table_begin'){
                sampLines = c(sampLines, oneLine)
                oneLine = readLines(con, n = 1, warn = FALSE)
            }
            i = i+1
            sampleData[[i]] = sampLines
            print(i)
        }
    }
    
    close(con)
    
    names(sampleData) = sapply(sampleData,function(x){
        strsplit(x[1],' = ')[[1]][2]
    })
    
    sampleData = sapply(sampleData,function(x){
        x[grepl('^\\!',x)]
    })
    
    samples = lapply(sampleData,function(x){
        singleSample = sapply(x, function(y){
            out = strsplit(y, '(\ =\ (?!.*?:\ ))|(:\ )', perl=T)[[1]]
            if (len(out)==1){
                out[2] = "NULL"
            }
            return(out[2])
        })
        
        names(singleSample) = sapply(x,function(y){
            out = strsplit(y, '(\ =\ (?!.*?:\ ))|(:\ )', perl=T)[[1]][1]
            return(out)
        })
        
        # some fields occur more than once. merge'em
        dups = unique(names(singleSample)[duplicated(names(singleSample))])       
        
        for (i in 1:len(dups)){
            temp  = paste0(singleSample[names(singleSample) %in% dups[i]], collapse = ' ')
            singleSample = singleSample[!names(singleSample) %in% dups[i]]
        }
        
        
        # print(names(singleSample))
        return(singleSample)
    })
    
    
    
    fields = table(unlist(lapply(samples,names)))
    
    if (mergeFrame[1] =='intersect'){
        fields = names(fields[fields==max(fields)])
        samples = lapply(samples, function(x){
            x[fields]
        })
    }
    samples = as.data.frame(t(as.data.frame(samples)))
    return(samples)
}