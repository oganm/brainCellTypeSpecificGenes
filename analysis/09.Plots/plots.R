library(dplyr)
library(viridis)
library(ggplot2)
library(reshape2)
library(cowplot)
library(ogbox)
library(homologene)
library(scales)
library(gplots)
library(stringr)
library(magrittr)
source('R/cellColors.R')
source('R/puristOut.R')
source('R/heatmap.3.R')
source('R/regionize.R')
order = read.design('data/meltedDesign.tsv') %>%
    arrange(MajorType,Neurotransmitter1,PyramidalDeep) %>% 
    filter(!is.na(PyramidalDeep)) %>% .$PyramidalDeep %>% unique

# expression plot of selected genes by region -----------------
list[geneDat, exp] = read.exp('data/finalExp.csv')  %>% sepExpr
design = read.design('data/meltedDesign.tsv')

dpylrFriendly = cbind(design, t(exp))

dpylrFriendly %<>% arrange(MajorType,Neurotransmitter1,PyramidalDeep) %>% filter(!is.na(PyramidalDeep))
regionSamples = regionize(design = dpylrFriendly,regionNames = 'Region',groupNames = 'PyramidalDeep')

regionSamples = c(regionSamples,list(All = dpylrFriendly$PyramidalDeep))
names(regionSamples) = gsub('_.*','',names(regionSamples))
list[design,exp] = dpylrFriendly %>% sepExpr

colors = cellColors()

dir.create('analysis/09.Plots/GenePlots')

for (i in 1:len(regionSamples)){
    
    genes = puristOut(paste0('analysis/01.Gene Selection/FinalGenes/PyramidalDeep/',names(regionSamples)[i]))
    genes = genes[regionSamples[[i]] %>% unique %>% trimNAs]
    genes %<>% lapply(function(x){
        x[!grepl('\\|',x)]
    })
    
    relExp = exp[!is.na(regionSamples[[i]]),match(unlist(genes), geneDat$Gene.Symbol)]
    relGene = geneDat[match(unlist(genes), geneDat$Gene.Symbol),]
    relExp = apply(relExp,2,scale01)
    
    
    geneCellTypes = str_extract(names(unlist(genes)) , regexMerge(cellTypes))
    
    heatCols = toColor(design$PyramidalDeep[!is.na(regionSamples[[i]])], colors)
    geneCols = toColor(geneCellTypes, colors)
    
    png(paste0('analysis/09.Plots/GenePlots/',names(regionSamples)[i]),height=1000,width=1000)
    heatmap.2(t(relExp),Rowv=F,Colv=F,trace='none',col=viridis(20),ColSideColors=heatCols$cols,RowSideColors=geneCols$cols,labRow='',labCol='', main = names(regionSamples[[i]]))
    legend(title= 'Cell Types','bottomleft' ,legend= names(heatCols$palette), fill = heatCols$palette )
    dev.off()
}



# 01. expression plot of the selected genes -------------
list[geneDat, exp] = read.exp('data/finalExp.csv')  %>% sepExpr
design = read.design('data/meltedDesign.tsv')

dpylrFriendly = cbind(design, t(exp))
cellTypes = sort(trimNAs(unique(dpylrFriendly$PyramidalDeep)))
dpylrFriendly$JustPyra = factor(dpylrFriendly$PyramidalDeep, 
                                levels = c('Astrocyte', 'Microglia', 'Oligo', 
                                           cellTypes[!cellTypes %in% c('Astrocyte', 'Microglia', 'Oligo')]))


dpylrFriendly = dpylrFriendly %>% filter(!is.na(JustPyra)) %>% arrange(JustPyra)
.
list[design,exp] = dpylrFriendly %>% sepExpr
cellTypes = unique(design$PyramidalDeep)

# get all genes for a given cell type
genes  = allPuristOut('analysis/01.Gene Selection/FinalGenes/PyramidalDeep/')

files = list.files('analysis/01.Gene Selection//RotSel/Marker/',recursive=T,full.names=T)
files = files[grepl('PyramidalDeep',files)]
genes = lapply(files,function(x){
    try({read.table(x)$V1})
})

names(genes) = basename(files)

lapply(unique(names(genes)), function(x){
    
})


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





# 02.figure mouse single cell ----------

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
    
    set.seed(1)
    GeneOrder = existanceMatri[[i]][existanceMatri[[i]][2:ncol(existanceMatri[[i]])] %>% 
                                        dist %>%
                                        hclust %>% 
                                        as.dendrogram %>% 
                                        reorder(rowMeans(existanceMatri[[i]][2:ncol(existanceMatri[[i]])])) %>%
                                        labels,1]
    set.seed(1)
    cellOrder = existanceMatri[[i]][2:ncol(existanceMatri[[i]])] %>%
        t %>% 
        dist %>% 
        hclust %>%
        as.dendrogram %>%
        reorder(colMeans(existanceMatri[[i]][2:ncol(existanceMatri[[i]])])) %>%
        labels
    
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


plot = plot_grid(plotlist = plots, ncol=5, labels = letters[1:10])
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
    set.seed(1)
    GeneOrder = existanceMatri[[i]][existanceMatri[[i]][2:ncol(existanceMatri[[i]])] %>% 
                                        dist %>%
                                        hclust %>% 
                                        as.dendrogram %>% 
                                        reorder(rowMeans(existanceMatri[[i]][2:ncol(existanceMatri[[i]])])) %>%
                                        labels,1]
    set.seed(1)
    cellOrder = existanceMatri[[i]][2:ncol(existanceMatri[[i]])] %>%
        t %>% 
        dist %>% 
        hclust %>%
        as.dendrogram %>%
        reorder(colMeans(existanceMatri[[i]][2:ncol(existanceMatri[[i]])])) %>%
        labels
    
#      if(names(existanceMatri)[i] %in% c('Astrocyte', 'Oligo', 'Microglia')){
#         dendro =  existanceMatri[[i]][2:ncol(existanceMatri[[i]])] %>% dist %>% hclust %>% as.dendrogram 
#         apply(existanceMatri[[i]][labels(dendro[[1]]),2:ncol(existanceMatri[[i]])],1,sum)
#         apply(existanceMatri[[i]][labels(dendro[[2]]),2:ncol(existanceMatri[[i]])],1,sum)
#      }
    
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


plot = plot_grid(plotlist = plots, ncol=5, labels = letters[11:20])
save_plot('analysis//09.Plots/singleHumanCellPlot.png',plot, base_height=10,base_aspect_ratio=2)

# get crappy genes for glial cells
set.seed(1)
crappyGenesAstro = existanceMatri$Astrocyte[,2:ncol(existanceMatri$Astrocyte)] %>% 
    as.matrix %>%
    heatmap.2(trace='none') %$% 
    rowDendrogram[[2]][[2]][[1]] %>% labels %>% existanceMatri$Astrocyte[.,1]
set.seed(1)
crappyGenesOligo = existanceMatri$Oligo[,2:ncol(existanceMatri$Oligo)] %>% 
    as.matrix %>%
    heatmap.2(trace='none') %$% 
    rowDendrogram[[2]][[2]][[1]] %>% labels %>% existanceMatri$Oligo[.,1]
set.seed(1)
crappyGenesMicroglia = existanceMatri$Microglia[,2:ncol(existanceMatri$Microglia)] %>% 
    as.matrix %>%
    heatmap.2(trace='none') %$% 
    rowDendrogram[[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]][[2]] %>% labels %>% existanceMatri$Microglia[.,1]

crappyGenes = list(Astrocyte = crappyGenesAstro,
                   Oligo = crappyGenesOligo,
                   Microglia = crappyGenesMicroglia)

# human coexpression crappy genes -------
dir.create('analysis//09.Plots/crappyGenesHumanCoexp')
for (i in 1:len(crappyGenes)){
    coexpression = read.csv(paste0('analysis/04.Human Coexpression/intraGroup/singleData/cortex/', names(crappyGenes)[i]))
    coexpression = coexpression[,2:ncol(coexpression)]
    isCrappy = colnames(coexpression) %in% crappyGenes[[i]]
    colors = toColor(isCrappy,c('FALSE'='white','TRUE' = muted('red')))$cols
    png(paste0('analysis//09.Plots/crappyGenesHumanCoexp/', names(crappyGenes)[i],'.png'),width=800,height=800)
    heatmap.2(coexpression %>% as.matrix,trace = 'none',col=viridis(20), dendrogram='column', ColSideColors=colors, main = names(crappyGenes)[i],labRow='')
    dev.off()
}

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
# brain region
estims = read.table('analysis/05.Brain Estimations/estimations/cortex_white/estimations', header=T,sep='\t')
estims[1:(ncol(estims)-1)] = apply(estims[1:(ncol(estims)-1)],2,scale01)
frame =melt(estims)
names(frame) = c('brainRegions','cellType','estimation')
frame$cellType = factor(frame$cellType, levels =order) %>% droplevels
frame = frame %>% group_by(cellType)
maxY = frame %>% summarise(max(estimation))
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
    # geom_point()+
    theme(axis.text.x = element_text(angle=45, hjust = 1),
          strip.text.x = element_text(size = 16)) +
    geom_text(data=signifFrame , aes(x = x, y=y, label = markers),size=10) + 
    coord_cartesian(ylim = c(-0.10, 1.10)) 

ggsave(filename='analysis//09.Plots/cortex_WhiteMatterEstimations.png',p,width=11,height=10)

# parkinsons
estims = read.table('analysis/05.Brain Estimations/estimations/parkinsonMale/estimations', header=T,sep='\t')
estims[1:(ncol(estims)-1)] = apply(estims[1:(ncol(estims)-1)],2,scale01)
frame =melt(estims)
names(frame) = c('parkinsons','cellType','estimation')
frame$parkinsons[frame$parkinsons] = "Parkinson's"
frame$parkinsons[frame$parkinsons=='FALSE'] ='Control'
frame$cellType = factor(frame$cellType, levels = order) %>% droplevels
maxY = frame %>% summarise(max(estimation))
ps= by(frame, frame$cellType, function(x){
    a1 = x %>% filter(parkinsons == 'Parkinson\'s') %>% select(estimation) %>% unlist
    a2 = x %>% filter(parkinsons == 'Control') %>% select(estimation) %>% unlist
    wilcox.test(a1,a2)$p.value
})
ps = p.adjust(ps)
markers = signifMarker(ps)
signifFrame = data.frame(markers, x= 1.5, y= 1,cellType = names(ps))

p  = frame %>%
    ggplot(aes(x=parkinsons, y = estimation)) + facet_wrap(~cellType) + 
    geom_violin( color="#C4C4C4", fill="#C4C4C4") +
    geom_boxplot(width=0.1,fill = 'lightblue') + 
    geom_point()+
    theme(axis.text.x = element_text(angle=45, hjust = 1),
          strip.text.x = element_text(size = 16)) +
    geom_text(data=signifFrame , aes(x = x, y=y, label = markers),size=10) + 
    coord_cartesian(ylim = c(-0.10, 1.10)) 

ggsave(filename='analysis//09.Plots/parkinsonEstimations.png',p,width=11,height=10)
    
frame %<>% filter(cellType == 'Dopaminergic')
signifFrame %<>% filter(cellType == 'Dopaminergic')
p  = frame %>%
    ggplot(aes(x=parkinsons, y = estimation)) + facet_wrap(~cellType) + 
    geom_violin( color="#C4C4C4", fill="#C4C4C4") +
    geom_boxplot(width=0.1,fill = 'lightblue') + 
    geom_point()+
    theme(axis.text.x = element_text(angle=45, hjust = 1),
          strip.text.x = element_text(size = 16)) +
    geom_text(data=signifFrame , aes(x = x, y=y, label = markers),size=10) + 
    coord_cartesian(ylim = c(-0.10, 1.10)) 
ggsave(filename='analysis//09.Plots/parkinsonDopaminergicAlone.png',p,width=4,height=5)



# 06. blood run plots-------------


# blood marker gene plots lm22 ------
markerPlot = function(expression,
                      design, 
                      geneList,
                      cellTypeNaming,
                      sampleNaming = 'sampleName',
                      order = cellTypeNaming, 
                      colors, 
                      geneSymbol = 'Gene.Symbol',
                      fileName, 
                      main=NULL){
    list[geneDat, exp] = expression  %>% sepExpr
    exp = exp[match(design[,sampleNaming], colnames(exp))]
    dpylrFriendly = cbind(design, t(exp))
    dpylrFriendly %<>% arrange_(.dots = order) %>% filter_(!is.na(cellTypeNaming))
    list[design,exp] = dpylrFriendly %>% sepExpr
    geneList %<>% lapply(function(x){
        x[!grepl('\\|',x)]
    })
    
    genes = geneList
    
    relExp = exp[,match(unlist(genes), geneDat[,geneSymbol]) %>% trimNAs]
    relGene = geneDat[match(unlist(genes), geneDat[,geneSymbol]) %>% trimNAs,]
    relExp = apply(relExp,2,scale01)
    
    cellTypes = names(colors)
    
    genes = unlist(genes)[unlist(genes) %in% relGene[,geneSymbol]]
    
    geneCellTypes = str_extract(names(genes) , regexMerge(cellTypes,exact=TRUE))
    
    heatCols = toColor(design[,cellTypeNaming], colors)
    geneCols = hede(geneCellTypes, colors)
    
    png(fileName,height=1000,width=1000)
    heatmap.2(t(relExp),Rowv=F,Colv=F,trace='none',col=viridis(20),ColSideColors=heatCols$cols,RowSideColors=geneCols$cols,labRow='',labCol='', main = main)
    legend(title= 'Cell Types','bottomleft' ,legend= names(heatCols$palette), fill = heatCols$palette )
    dev.off()
    
}

expr= read.exp('data/bloodCellType/finalBlood.csv') 
design = read.design('analysis/06.Blood validation/BloodCells.tsv')
genes = puristOut('analysis/06.Blood validation/humanGenes/rotSel/Relax/Abreviated.name/')
colors = c('CD4 memory T cells-' = 'darkgreen',
           'CD4 memory T cells+' = 'palegreen',
           'CD4 naïve T cells' = 'chartreuse4',
           'CD8 T cells' = 'blue4',
           'DCs-' = 'coral2',
           'DCs+' = 'coral4',
           'Eos' = 'darkmagenta',
           'M0-MΦs' = 'gold4',
           'M1-MΦs' = 'gold',
           'M2-MΦs' = 'goldenrod',
           'MCs-' = 'mediumorchid1',
           'MCs+' = 'mediumorchid4',
           'Memory B cells' = 'lightsteelblue1',
           'Monos' = 'powderblue',
           'Naïve B cells' = 'lightsteelblue4',
           'NK cells-' = 'red1',
           'NK cells+' = 'red4',
           'PCs' = 'violetred1',
           'PMNs' = 'yellow',
           'Tfh cells' = 'springgreen4',
           'Tregs' = 'springgreen',
           'γδ T cells' = 'tan1')
markerPlot(expression = expr,
           design = design,
           geneList = genes,
           cellTypeNaming = 'Abreviated.name',
           colors = colors,
           fileName = 'analysis/09.Plots/blood22plot',
           main = 'lm22 Genes')


# lm11------
genes = puristOut('analysis/06.Blood validation/humanGenes/rotSel/Relax/Leuk11//')
colors = c('B Cell' = 'lightsteelblue1',
           'CD4' = 'darkgreen',
           'CD8 T cells' = 'blue4',
           'Dendritic' = 'coral4',
           'Eos' = 'darkmagenta',
           'GammaDeltaT' = 'tan1',
           'Mast' = 'mediumorchid4',
           'MonoMacro' = 'gold',
           'Neutrophils' = 'yellow',
           'NK' = 'red1',
           'PCs' = 'violetred1')
markerPlot(expression = expr,
           design = design,
           geneList = genes,
           cellTypeNaming = 'Leuk11',
           order = c('Leuk11', 'Abreviated.name'),
           colors = colors,
           fileName = 'analysis/09.Plots/blood11plot',
           main = 'lm11 Genes')

genes = puristOut('analysis/06.Blood validation/mouseGenes/Fold/Relax/lm11/') %>% lapply(function(x){mouse2human(x)$humanGene})
markerPlot(expression = expr,
           design = design,
           geneList = genes,
           cellTypeNaming = 'Leuk11',
           order = c('Leuk11', 'Abreviated.name'),
           colors = colors,
           fileName = 'analysis/09.Plots/blood11plot',
           main = 'lm11 Genes')


# misc plots ------------------------
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


astroDendro = existanceMatri[[1]][2:ncol(existanceMatri[[1]])] %>% dist %>% hclust %>% as.dendrogram

existanceMatri[[1]][astroDendro[[2]] %>% labels,1]

existanceMatri[[1]][astroDendro[[2]][[1]][[1]] %>% labels, 2:ncol(existanceMatri[[1]])] %>% as.matrix %>% heatmap.2(trace='none')

geneNames = existanceMatri[[1]][astroDendro[[2]][[1]][[1]] %>% labels,1]

astroGroupRots = 'analysis/05.Brain Estimations/estimations/cortex_white/Astrocyte groupRots' %>% 
    read.table(header=T, sep='\t')
astroGroupRots$Gene = rownames(astroGroupRots)

astroGroupRots = astroGroupRots %>% arrange(frontal.cortex) %>% mutate(isBad =Gene %in% geneNames ) 


astroGroupRots %>% 
    ggplot(aes(y=frontal.cortex, x= 1, color = isBad)) + geom_point(size = 15,shape=95) +  scale_colour_manual(values = c(muted('blue'),'red')) + ylab('Gene rotation for PC1') +
    theme(legend.position="none",
          axis.title.y = element_text(size = 20))



pal(toColor(astroGroupRots$Gene %in% geneNames,c(muted('blue'),muted('red')))$col)
