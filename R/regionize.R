regionize = function(design,regionNames,groupNames){
    regions =
        trimNAs(
            trimElement(
                unique(
                    unlist(
                        strsplit(as.character(design[,regionNames]),',')))
                ,c('ALL','All','all','Cerebrum'))) #S pecial names
    regionBased = expand.grid(groupNames, regions)
    regionGroups = vector(mode = 'list', length = nrow(regionBased))
    names(regionGroups) = paste0(regionBased$Var2,'_',regionBased$Var1)
    for (i in 1:nrow(regionBased)){
        regionGroups[[i]] = design[,as.character(regionBased$Var1[i])]
        
        # remove everything except the region and ALL labeled ones. for anything but cerebellum, add Cerebrum labelled ones as well
        if (regionBased$Var2[i] == 'Cerebellum'){
            regionGroups[[i]][!grepl(paste0('(^|,)((',regionBased$Var2[i],')|((A|a)(L|l)(l|l)))($|,)'),design[,regionNames])] = NA
        } else {
            # look for cerebrums
            cerebrums = unique(regionGroups[[i]][grepl('(Cerebrum)',design[,regionNames])])
            
            # find which cerebrums are not represented in the region
            cerebString = paste(cerebrums[!cerebrums %in% regionGroups[[i]][grepl(paste0('(^|,)((',regionBased$Var2[i],')|((A|a)(L|l)(l|l)))($|,)'),design[,regionNames])]],
                                collapse = ')|(')
            
            # add them as well (or not remove them as well) with all the rest of the region samples
            regionGroups[[i]][(!grepl(paste0('(^|,)((',regionBased$Var2[i],')|((A|a)(L|l)(l|l)))($|,)'),design[,regionNames])
                               & !(grepl(paste0('(',cerebString,')'),design[,as.character(regionBased$Var1[i])]) & grepl('Cerebrum',design[,regionNames])))] =  NA
            
        } 
    }
    return(regionGroups)
}