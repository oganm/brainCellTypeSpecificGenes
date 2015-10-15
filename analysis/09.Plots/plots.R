library(dplyr)
library(viridis)
library(ggplot2)
library(reshape2)
library(cowplot)
library(ogbox)
library(scales)
library(gplots)
library(stringr)
source('R/cellColors.R')
source('R/puristOut.R')
source('R/heatmap.3.R')

# 01. expression plot of the selected genes -------------
list[geneDat, exp] = read.exp('data/finalExp.csv')  %>% sepExpr
design = read.design('data/meltedDesign.tsv')

dpylrFriendly = cbind(design, t(exp))
cellTypes = sort(trimNAs(unique(dpylrFriendly$PyramidalDeep)))
dpylrFriendly$JustPyra = factor(dpylrFriendly$PyramidalDeep, 
                                levels = c('Astrocyte', 'Microglia', 'Oligo', 
                                           cellTypes[!cellTypes %in% c('Astrocyte', 'Microglia', 'Oligo')]))


dpylrFriendly = dpylrFriendly %>% filter(!is.na(JustPyra)) %>% arrange(JustPyra)

list[design,exp] = dpylrFriendly %>% sepExpr
cellTypes = unique(design$PyramidalDeep)

# get all genes for a given cell type
genes  = allPuristOut('analysis/01.Gene Selection/FinalGenes/PyramidalDeep/')

colors = cellColors()
              
colors = colors[cellTypes]

geneCellTypes = sapply(names(unlist(genes)), function(x){
    str_split(x,'\\.')[[1]][2]
})

geneCellTypes = str_extract(geneCellTypes, regexMerge(cellTypes))

cellTypeGenes = lapply(cellTypes,function(x){
    a = unique(unlist(genes)[ geneCellTypes %in% x])
})
names(cellTypeGenes) = cellTypes
geneCellTypes = str_extract(names(unlist(cellTypeGenes)) , regexMerge(cellTypes))

relExp = exp[,match(unlist(cellTypeGenes), geneDat$Gene.Symbol)]
relGene = geneDat[match(unlist(cellTypeGenes), geneDat$Gene.Symbol),]


relExp = apply(relExp,2,scale01)


heatCols = toColor(design$PyramidalDeep, colors)
geneCols = toColor(geneCellTypes, colors)

png('analysis/09.Plots/markerGenes.png',height=1000,width=1000)
heatmap.2(t(relExp),Rowv=F,Colv=F,trace='none',col=viridis(20),ColSideColors=heatCols$cols,RowSideColors=geneCols$cols,labRow='',labCol='')
legend(title= 'Cell Types','bottomleft' ,legend= names(heatCols$palette), fill = heatCols$palette )
dev.off()
w# 02.figure mouse single cell ----------

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

files = list.files('analysis/04.Human Coexpression/coexpData',full.names=T,recursive=T)
files = files[grepl(pattern='RData',files)]
files = files[grepl(pattern='spear', files)]
coexpData = lapply(files,function(x){
    load(x)
    return(list(coexpressions=lapply(coexpressions,sort),
                daCorvec=sort(daCorvec),
                subFolderName=subFolderName,
                plotOut=plotOut,
                regionsIn=regionsIn))
})

files = list.files('analysis/04.Human Coexpression/coexpData',full.names=T,recursive=T)
files = files[grepl(pattern='pVals',files)]
files = files[grepl(pattern='spear', files)]

pVals = lapply(files,function(x){
    tb = read.table(x)
    vec = tb$V2
    names(vec) = tb$V1
    return(vec)
})

# remove glial ones from others to prevent clutter
for (i in 2:len(pVals)){
    pVals[[i]] = pVals[[i]][!names(pVals[[i]]) %in% c('Astrocyte', 'Microglia','Oligo','GabaRelnCalb', 'Pyramidal_Glt_25d2', 'Pyramidal_Thy1')] 
    coexpData[[i]]$coexpressions = coexpData[[i]]$coexpressions[!names(coexpData[[i]]$coexpressions) %in%  c('Astrocyte', 'Microglia','Oligo','GabaRelnCalb', 'Pyramidal_Glt_25d2','Pyramidal_Thy1' ) ]
}

pAdjusted  = relist(flesh=p.adjust(unlist(pVals),method='fdr'),
                    skeleton=pVals)

plotVioBox = function(relevants, ps, all, name){
    relevantsFull = rbind(all,relevants)
    # values are presorted so just take the last
    yMax = sapply(unique(relevantsFull$L1), function(x){
        tail(relevantsFull$value[relevantsFull$L1 %in% x],n=1)
    })
    # yMax = max(relevantsFull$value)
    cellTypes = unique(relevants$L1)
    
    plot = ggplot(relevantsFull, aes(x=L1, y = value, group = L1)) + 
        list(geom_violin(color = "#C4C4C4", fill = "#C4C4C4"),
             geom_boxplot(width = 0.2, fill = "lightblue")) + theme_bw() + 
        scale_y_continuous(name = 'Gene to gene correlations') +
        scale_x_discrete(name = '') +
        theme(axis.text.x  = element_text(size=25, angle = 45,hjust=1),
              axis.title.y = element_text(vjust=0.5, size=25),
              title = element_text(vjust=0.5, size=25),
              axis.text.y = element_text(size = 13))+
        geom_signif(c(1,ps),maxY=yMax) + ggtitle(name)
    return(plot)
}

set.seed(1)
names(pAdjusted) = sapply(coexpData,function(x){basename(x$plotOut)})
names(pAdjusted)[1] = 'All Regions'

plots = mclapply(1:len(coexpData), function(i){
    relevants = coexpData[[i]]$coexpressions
    ps = pAdjusted[[i]]
    relevants = melt(relevants)
    
    allTemp = sample(coexpData[[i]]$daCorvec,size = 2000)
    all = data.frame(value = allTemp, L1 = 'all')
    # all = data.frame(value = coexpData[[i]]$daCorvec, L1 = 'all')
    plot = plotVioBox(relevants, ps, all, names(pAdjusted)[i])    
},mc.cores = len(coexpData))

plot = plot_grid(plotlist = plots, ncol=2)
save_plot('analysis//09.Plots/humanCoexpressionPlot.png',plot, base_height=20,base_aspect_ratio=1.3)



# 05. Estimation plots --------
estims = read.table('analysis/05.Brain Estimations/estimations/cortex_white/estimations', header=T,sep='\t')
estims[1:(ncol(estims)-1)] = apply(estims[1:(ncol(estims)-1)],2,scale01)
frame =melt(estims)
names(frame) = c('brainRegions','cellType','estimation')
frame$cellType = factor(frame$cellType, levels = c('Astrocyte', 'Oligo', 'Microglia',
                         as.char(unique(frame$cellType)[!unique(frame$cellType) %in% c('Astrocyte', 'Oligo', 'Microglia')])))
frame = frame %>% group_by(cellType)
maxY = frame %>% summary(max())
ps= by(frame, frame$cellType, function(x){
    a1 = x %>% filter(brainRegions == 'frontal cortex') %>% select(estimation) %>% unlist
    a2 = x %>% filter(brainRegions == 'white matter') %>% select(estimation) %>% unlist
    wilcox.test(a1,a2)$p.value
})

ps = p.adjust(ps)

markers = signifMarker(ps)
signifFrame = data.frame(markers, x= 1.5,y= 1,cellType = names(ps))

p  = frame %>%
    ggplot(aes(x=brainRegions, y = estimation)) + facet_wrap(~cellType) + 
    geom_violin( color="#C4C4C4", fill="#C4C4C4") +
    geom_boxplot(width=0.1,fill = 'lightblue') + 
    theme(axis.text.x = element_text(angle=45, hjust = 1),
          strip.text.x = element_text(size = 16)) +
    geom_text(data=signifFrame , aes(x = x, y=y, label = markers))
ggsave(filename='analysis//09.Plots/cortex_WhiteMatterEstimations.png',p,width=11,height=10)
