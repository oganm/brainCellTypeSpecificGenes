# 02.figure mouse single cell ----------
library(dplyr)
library(viridis)
library(ggplot2)
library(reshape2)
library(cowplot)
library(ogbox)

simuProbs = read.table('analysis/02.Mouse Single Cell/Output/NoTresholdNoDub0.3333/simuProbs', header= TRUE)
pVals = read.table("analysis/02.Mouse Single Cell/Output/NoTresholdNoDub0.3333/realProbs", header=TRUE)

sigMarks = signifMarker(pVals$ps)

allFiles = list.files('analysis/02.Mouse Single Cell/Output/NoTresholdNoDub0.3333/',full.names=TRUE)


allFiles = allFiles[!grepl('(simuProbs)|(realProbs)', allFiles)]

existanceMatri = lapply(allFiles, function(x){
    read.table(x,skip=1)
})
names(existanceMatri) = basename(allFiles)


#largeExistance = do.call(rbind, existanceMatri)
# cellOrder = largeExistance[2:ncol(largeExistance)] %>% t %>% dist %>% hclust %>% as.dendrogram %>% labels

existanceMatri = existanceMatri[c('Astrocyte', 'Oligo', 'Microglia',
                                  names(existanceMatri)[!names(existanceMatri) %in% c('Astrocyte', 'Oligo', 'Microglia')])]

pVals = pVals[c('Astrocyte', 'Oligo', 'Microglia',
                names(existanceMatri)[!names(existanceMatri) %in% c('Astrocyte', 'Oligo', 'Microglia')]),]

sigMarks = signifMarker(pVals$ps)


plots = lapply(1:len(existanceMatri), function(i){
    
    GeneOrder = existanceMatri[[i]][existanceMatri[[i]][2:ncol(existanceMatri[[i]])] %>% dist %>% hclust %>% as.dendrogram %>% labels,1]
    cellOrder = existanceMatri[[i]][2:ncol(existanceMatri[[i]])] %>% t %>% dist %>% hclust %>% as.dendrogram %>% labels
    
    frame = existanceMatri[[i]] %>% (reshape2::melt)
    frame$V1 = factor(frame$V1, levels = GeneOrder)
    frame$variable = factor(frame$variable, levels = cellOrder)
    
    
    p=frame %>% 
        ggplot(aes(y = V1, x = variable)) + 
        geom_tile(aes(fill = value)) + 
        scale_fill_gradient(low = 'white', high = muted('blue')) + theme_bw() + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              legend.position= 'none',
              plot.title = element_text(size = 22)) + 
        ylab('Marker genes') + 
        xlab('Cells') +
        ggtitle(bquote(~ .(names(existanceMatri)[[i]])^.(sigMarks[i]) ))
    return(p)
})


plot = plot_grid(plotlist = plots, ncol=5)
save_plot('analysis//09.Plots/singleMouseCellPlot.png',plot, base_height=10,base_aspect_ratio=2)


# 03. figure human single cell -------------
simuProbs = read.table('analysis/03.Human Single Cell/Output/NoTreshold0.3333/simuProbs', header= TRUE)
pVals = read.table("analysis/03.Human Single Cell/Output/NoTreshold0.3333/realProbs", header=TRUE)

allFiles = list.files('analysis/03.Human Single Cell/Output/NoTreshold0.3333/',full.names=TRUE)


allFiles = allFiles[!grepl('(simuProbs)|(realProbs)', allFiles)]

existanceMatri = lapply(allFiles, function(x){
    read.table(x,skip=1)
})
names(existanceMatri) = basename(allFiles)
    

existanceMatri = existanceMatri[c('Astrocyte', 'Oligo', 'Microglia',
                                  names(existanceMatri)[!names(existanceMatri) %in% c('Astrocyte', 'Oligo', 'Microglia')])]

pVals = pVals[c('Astrocyte', 'Oligo', 'Microglia',
                names(existanceMatri)[!names(existanceMatri) %in% c('Astrocyte', 'Oligo', 'Microglia')]),]

sigMarks = signifMarker(pVals$ps)


plots = lapply(1:len(existanceMatri), function(i){
    
    GeneOrder = existanceMatri[[i]][existanceMatri[[i]][2:ncol(existanceMatri[[i]])] %>% dist %>% hclust %>% as.dendrogram %>% labels,1]
    cellOrder = existanceMatri[[i]][2:ncol(existanceMatri[[i]])] %>% t %>% dist %>% hclust %>% as.dendrogram %>% labels
    if(names(existanceMatri)[i] %in% c('Astrocyte', 'Oligo', 'Microglia')){
       dendro =  existanceMatri[[i]][2:ncol(existanceMatri[[i]])] %>% dist %>% hclust %>% as.dendrogram 
       apply(existanceMatri[[i]][labels(dendro[[1]]),2:ncol(existanceMatri[[i]])],1,sum)
       apply(existanceMatri[[i]][labels(dendro[[2]]),2:ncol(existanceMatri[[i]])],1,sum)
       
    }
    
    frame = existanceMatri[[i]] %>% (reshape2::melt)
    frame$V1 = factor(frame$V1, levels = GeneOrder)
    frame$variable = factor(frame$variable, levels = cellOrder)
        
    p=frame %>% 
        ggplot(aes(y = V1, x = variable)) + 
        geom_tile(aes(fill = value)) + 
        scale_fill_gradient(low = 'white', high = muted('blue')) + theme_bw() + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              legend.position= 'none',
              plot.title = element_text(size = 22)) + 
        ylab('Marker genes') + 
        xlab('Cells') +
        ggtitle(bquote(~ .(names(existanceMatri)[[i]])^.(sigMarks[i]) ))
    return(p)
})


plot = plot_grid(plotlist = plots, ncol=5)
save_plot('analysis//09.Plots/singleHumanCellPlot.png',plot, base_height=10,base_aspect_ratio=2)

# 04.human coexpression ------

