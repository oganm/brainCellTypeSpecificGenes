#treshold detection for absance presence call in RNA seq.

library(ogbox)
library(foreach)
library(doMC)
library(parallel)
library(mixtools)

rnaSeqTresh = function(rnaExp, # matrix of expressio values with row names = gene names
                       outfile=NULL,
                       cores=16){    

    if (detectCores()<cores){ 
        cores = detectCores()
        print('max cores exceeded')
        print(paste('set core no to',cores))
    }
    registerDoMC(cores)
    # to every gene, fit two gaussians, lower one is counted as non expressed,
    # higher is expressed
    gaus = foreach(x = 1:nrow(rnaExp)) %dopar% {
        print(x)
        tryCatch({
            normalmixEM(rnaExp[x,],maxrestarts=5, verb = F)
        }, error = function(e){
            return(NA)
        })
    }
    
    # find the probability of each count. set the treshold to the point where the 
    # probability of a being in the expressed group is higher than probability of 
    # being in the non expressed group.
    tresholds = sapply(gaus,function(x){
        if (is.na(x)){
            return(1)
        }
        prob1=pnorm(0:max(x$x),mean = x$mu[1], x$sigma[1],lower.tail=F)
        prob2=pnorm(0:max(x$x),mean = x$mu[2], x$sigma[2],lower.tail=T)
        return(min(which(prob2>prob1)))
    })
    
    tresholds = as.matrix(tresholds)
    rownames(tresholds) = rn(rnaExp)
    
    if(!is.null(outfile)){
        write.table(tresholds,file=outfile,col.names=F,quote=F)
    }
    invisible(tresholds)
}
