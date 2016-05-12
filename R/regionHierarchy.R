regionHierarchy = list(All = list(Cerebrum = 
                                      list(Cortex = '',
                                           BasalForebrain ='',
                                           Striatum = '',
                                           Amygdala ='',
                                           Hippocampus = ''),
                                  Subependymal = '',
                                  Thalamus = '',
                                  Brainstem = 
                                      list(Midbrain = '',
                                           LocusCereuleus=''),
                                  Cerebellum = '',
                                  SpinalCord ='')
)


library(data.tree)
# to translate this into a data tree
listToDataTree = function(list){
    members = list %>% unlist %>% names
    # they have to have a commont root
    rootName = members[[1]] %>% strsplit(split='[.]') %>% {.[[1]][1]}
    teval(paste0(rootName, 
                 " = Node$new(members[[1]] %>% strsplit(split='[.]') %>% {.[[1]][1]})"),
          envir = parent.frame())
    for (x in members){
        newNodes = x %>% strsplit(split='[.]') %>% {.[[1]]}
        for (i in 2:len(newNodes)){
            if(!is.null(teval(rootName)$FindNode(newNodes[i]))){
                next
            }
            teval(paste0(newNodes[i], ' = ',newNodes[i-1],'$AddChild("',newNodes[i],'")'),
                  envir = parent.frame())
            
        }
    }
    
    return(teval(rootName))
}