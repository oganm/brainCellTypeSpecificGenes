library(ogbox)
library(dplyr)
library(limma)
library(cluster)
source('R/regionize.R')
source('R/puristOut.R')

exp = read.exp('data/finalExp.csv')
des = read.design('data/meltedDesign.tsv')

list[genes,exp]=sepExpr(exp) 

markerGenes = allPuristOut('analysis//01.Gene Selection//FinalGenes/forContanim/')

exp  = exp %>% select(which(!is.na(des$forContanim)))

des = des %>% filter(!is.na(forContanim))


regionGroups = regionize(des,'Region','forContanim')
regionGroupsUnique = lapply(regionGroups,function(x){
    x %>% unique %>% trimNAs
})



contaminations = lapply(1:len(regionGroups), function(i){
    markers = markerGenes[[names(markerGenes) %in% 
                              (regionGroups[i] %>% names %>% strsplit('_') %>% .[[1]] %>% .[1]) %>% which
                           ]]
    cindexes = vector(mode = 'list', length= len(markers))
    for (j in 1:len(markers)){
        mi = apply(exp[which(genes$Gene.Symbol %in% markers[[j]]), (!des$forContanim %in% names(markers)[j]) & (!is.na(regionGroups[[i]])) , drop=F],
                   1,min)
        ma = apply(exp[which(genes$Gene.Symbol %in% markers[[j]]),(des$forContanim %in% names(markers)[j]) & (!is.na(regionGroups[[i]])), drop=F],
                   1,max)
        
        contaminations = apply((exp[which(genes$Gene.Symbol %in% markers[[j]]),,drop=F]-mi)/(ma-mi),2,mean)
        contaminations[des$forContanim %in% names(markers)[j] | is.na(regionGroups[[i]])] = NA
    
        cindexes[[j]] = contaminations
        
    } 
    names(cindexes) = names(markers)
    return(cindexes)
})

names(contaminations) = (regionGroups %>% names %>% strsplit('_') %>% sapply(function(x){x[1]}))

contaminations %>% unlist(recursive=F) %>% as.data.frame %>% cbind(des,.) %>% write.design('analysis//11.Contamination/contaminDes.tsv')
    