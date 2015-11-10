library(ogbox)
library(magrittr)
library(dplyr)

geneAnalysisPipeline = function(gene){
    regionNames = read.table('data/development/Region_codes', stringsAsFactors=F)
    files = list.files('data/development//DevelopData',full.names=T)
    files = files[!grepl('full', files)]
    require(data.table)
    
    dir.create(paste0('analysis//99.GeneAnalysisPipeline/',gene,'/developPlots'),recursive=T, showWarnings = F)
    
    for (i in files){
        orbitalDat = fread(i)
        orbitalDat = orbitalDat[!is.na(orbitalDat$Gene.Symbol),]
        orbitalDat = mostVariable(orbitalDat)
        #orbitalDat = read.exp('Data/DevelopData/OFC')
        
        list[orbGenes, orbExp] = sepExpr(orbitalDat)
        meta = read.table('data//development//Metadata_GSE25219', sep='\t', stringsAsFactors=F)
        relMeta = meta[match(cn(orbExp),rn(meta)),]
        relMeta = relMeta[!duplicated(relMeta$brain.code),]
        orbExp = orbExp[,cn(orbExp) %in% rn(relMeta), with=F]
        orbExp = orbExp[,relMeta$Age_years<=20,with=F]
        relMeta = relMeta[relMeta$Age_years<=20,]
        png(paste0('analysis//99.GeneAnalysisPipeline/',gene,'/developPlots/',regionNames$V2[regionNames$V1 %in% basename(i)],'.png'))
        plot(unlist(relMeta$Age_years),unlist(orbExp[orbGenes$Gene.Symbol %in% gene,]),
             main=regionNames$V2[regionNames$V1 %in% basename(i)],xlab='Age', 
             ylab=paste(gene,'log2 Expression'))
        dev.off()
    }
    
    dir.create(paste0('analysis//99.GeneAnalysisPipeline/',gene,'/developPlotsEarly'),recursive=T, showWarnings = F)
    
    for (i in files){
        orbitalDat = fread(i)
        orbitalDat = orbitalDat[!is.na(orbitalDat$Gene.Symbol),]
        orbitalDat = mostVariable(orbitalDat)
        #orbitalDat = read.exp('Data/DevelopData/OFC')
        
        list[orbGenes, orbExp] = sepExpr(orbitalDat)
        meta = read.table('data//development//Metadata_GSE25219', sep='\t', stringsAsFactors=F)
        relMeta = meta[match(cn(orbExp),rn(meta)),]
        relMeta = relMeta[!duplicated(relMeta$brain.code),]
        orbExp = orbExp[,cn(orbExp) %in% rn(relMeta), with=F]
        orbExp = orbExp[,relMeta$Age_years<=1,with=F]
        relMeta = relMeta[relMeta$Age_years<=1,]
        png(paste0('analysis//99.GeneAnalysisPipeline/',gene,'/developPlotsEarly/',regionNames$V2[regionNames$V1 %in% basename(i)],'.png'))
        plot(unlist(relMeta$Age_years),unlist(orbExp[orbGenes$Gene.Symbol %in% gene,]),
             main=regionNames$V2[regionNames$V1 %in% basename(i)],xlab='Age', 
             ylab=paste(gene,'log2 Expression'))
        dev.off()
    }
    
    # region from brainspan data. all normalized -------------
    allExpress =sapply(files, function(i){
        orbitalDat = fread(i)
        orbitalDat = orbitalDat[!is.na(orbitalDat$Gene.Symbol),]
        list[orbGenes, orbExp] = sepExpr(orbitalDat)
        #return(mean(apply(orbExp,2,mean)))
        unlist(orbExp)
    })
    allExpress = melt(allExpress)
    names(allExpress) = c('log2Exp', 'BrainReg')
    
    allExpress[,2] = regionNames[,2][match(basename(allExpress[,2]),regionNames[,1])]
    
    allExpress[allExpress[,2] %in% 'Mediodorsal nucleus of the thalamus',2] = 'Thalamus'
    allExpress[allExpress[,2] %in%'Primary visual cortex',2] = 'Occipital Cortex'
    allExpress[allExpress[,2] %in% 'Superior temporal cortex',2] = 'S. Temporal Cortex'
    allExpress[allExpress[,2] %in% 'Posterior inferior parietal cortex',2] = 'Parietal Cortex'
    allExpress[allExpress[,2] %in%'Primary auditory cortex',2] = 'Auditory Cortex'
    allExpress[allExpress[,2] %in% 'Primary somatosensory cortex',2] = 'Somatosensory Cortex'
    allExpress[allExpress[,2] %in%'Inferior temporal cortex',2] = 'I. Temporal Cortex'
    allExpress[allExpress[,2] %in%'Primary motor cortex',2] = 'Motor Cortex'
    
    oldData = sapply(files,
                     function(i){
                         orbitalDat = fread(i)
                         orbitalDat = orbitalDat[!is.na(orbitalDat$Gene.Symbol),]
                         orbitalDat %<>% mostVariable
                         #orbitalDat = read.exp('Data/DevelopData/OFC')
                         
                         list[orbGenes, orbExp] = sepExpr(orbitalDat)
                         meta = read.table('data//development//Metadata_GSE25219', sep='\t', stringsAsFactors=F)
                         relMeta = meta[match(cn(orbExp),rn(meta)),]
                         relMeta = relMeta[!duplicated(relMeta$brain.code),]
                         orbExp = orbExp[,cn(orbExp) %in% rn(relMeta), with=F]
                         orbExp = orbExp[,relMeta$Age_years>20, with=F]
                         relMeta = relMeta[relMeta$Age_years>20,]
                         unlist(orbExp[orbGenes$Gene.Symbol %in% gene,])
                     })
    
    
    require(reshape2)
    require(ggplot2)
    oldData = melt(oldData)

    oldData$L1 = basename(oldData$L1)
    oldData = oldData[!grepl('full',oldData$L1),]
    oldData$L1 = regionNames$V2[match( oldData$L1,regionNames$V1)]
    
    oldData -> x
    x[x[,2] %in% 'Mediodorsal nucleus of the thalamus',2] = 'Thalamus'
    x[x[,2] %in% 'Primary visual cortex',2] = 'Occipital Cortex'
    x[x[,2] %in% 'Superior temporal cortex',2] = 'S. Temporal Cortex'
    x[x[,2] %in% 'Posterior inferior parietal cortex',2] = 'Parietal Cortex'
    x[x[,2] %in% 'Primary auditory cortex',2] = 'Auditory Cortex'
    x[x[,2] %in% 'Primary somatosensory cortex',2] = 'Somatosensory Cortex'
    x[x[,2] %in% 'Inferior temporal cortex',2] = 'I. Temporal Cortex'
    x[x[,2] %in% 'Primary motor cortex',2] = 'Motor Cortex'
    
    x -> oldData
    cols <-teval( paste0('c("Expression of All Genes"="lightgray",',gene,'="white")'))
    
    names(oldData) = names(allExpress)
    png(paste0('analysis//99.GeneAnalysisPipeline/',gene,'/Across regions.png'),height = 600,width=800)
    p = ggplot(oldData, aes(y = log2Exp, x = BrainReg)) +  
        geom_boxplot(data=allExpress, color ='gray',aes(fill ='Expression of All Genes'),outlier.size=0, width= 0.3)+ 
        geom_boxplot(aes(fill = gene), outlier.size=0) + 
        geom_point(position='jitter') + 
        theme_bw() +
        theme(axis.text.x = element_text(angle=90, size = 14)) +
        ylab(paste(gene,'expression')) + xlab('') + 
        scale_fill_manual(name="",values=cols)
    
    ylim1 = boxplot.stats(allExpress$log2Exp)$stats[c(1, 5)]
    ylim1[2] = ylim1[2]*1.05
    ylim1[1] = ylim1[1]*0.95
    
    p= p + coord_cartesian(ylim = ylim1)
    (p)
    dev.off()
    
    
    # just old ones ---------------
    files = list.files('data//development/DevelopDataOldens/',full.names=T)
    
    
    allExpress =sapply(files, function(i){
        orbitalDat = fread(i)
        orbitalDat = orbitalDat[!is.na(orbitalDat$Gene.Symbol),]
        list[orbGenes, orbExp] = sepExpr(orbitalDat)
        #return(mean(apply(orbExp,2,mean)))
        unlist(orbExp)
    })
    allExpress = melt(allExpress)
    names(allExpress) = c('log2Exp', 'BrainReg')
    
    allExpress[,2] = regionNames[,2][match(basename(allExpress[,2]),regionNames[,1])]
    
    allExpress[allExpress[,2] %in% 'Mediodorsal nucleus of the thalamus',2] = 'Thalamus'
    allExpress[allExpress[,2] %in% 'Primary visual cortex',2] = 'Occipital Cortex'
    allExpress[allExpress[,2] %in% 'Superior temporal cortex',2] = 'Temporal Cortex'
    
    allExpress = allExpress[allExpress[,2] %in% c('Cerebellar cortex',
                                                  'Hippocampus',
                                                  'Thalamus',
                                                  'Occipital Cortex',
                                                  'Temporal Cortex',
                                                  'Striatum'),]
    
    allExpress[,2] = factor(allExpress[,2],levels=c('Cerebellar cortex',
                                                    'Hippocampus',
                                                    'Thalamus',
                                                    'Occipital Cortex',
                                                    'Temporal Cortex',
                                                    'Striatum'))
    
    
        
        oldData = sapply(files,
                         function(i){
                             orbitalDat = fread(i)
                             orbitalDat = orbitalDat[!is.na(orbitalDat$Gene.Symbol),]
                             orbitalDat %<>% mostVariable
                             #orbitalDat = read.exp('Data/DevelopData/OFC')
                             
                             list[orbGenes, orbExp] = sepExpr(orbitalDat)
                             meta = read.table('data//development//Metadata_GSE25219', sep='\t', stringsAsFactors=F)
                             relMeta = meta[match(cn(orbExp),rn(meta)),]
                             relMeta = relMeta[!duplicated(relMeta$brain.code),]
                             orbExp = orbExp[,cn(orbExp) %in% rn(relMeta), with=F]
                             unlist(orbExp[grep(gene, orbGenes$Gene.Symbol),])
                         })  
    
    oldData %>% melt -> oldData
    
    x = oldData
    x[,2] = regionNames[,2][match(basename(x[,2]),regionNames[,1])]
    
    x[x[,2] %in% 'Mediodorsal nucleus of the thalamus',2] = 'Thalamus'
    x[x[,2] %in% 'Primary visual cortex',2] = 'Occipital Cortex'
    x[x[,2] %in% 'Superior temporal cortex',2] = 'Temporal Cortex'
    
    x = x[x[,2] %in% c('Cerebellar cortex',
                       'Hippocampus',
                       'Thalamus',
                       'Occipital Cortex',
                       'Temporal Cortex',
                       'Striatum'),]
    x[,2] = factor(x[,2],levels=c('Cerebellar cortex',
                                  'Hippocampus',
                                  'Thalamus',
                                  'Occipital Cortex',
                                  'Temporal Cortex',
                                  'Striatum'))
    
    colnames(x) = c('log2Exp', 'BrainReg')
    cols <-teval( paste0('c("Expression of All Genes"="lightgray",',gene,'="white")'))
    
    p = ggplot(x , aes(x= BrainReg, y= log2Exp))  + 
        geom_boxplot(data=allExpress, color ='gray',aes(fill ='Expression of All Genes'),outlier.size=0, width= 0.3) +
        geom_boxplot(aes(fill = gene), outlier.size=0) + 
        geom_point(position='jitter') + scale_y_continuous(paste(gene,'expession')) + 
        theme_bw() + theme(axis.text.x =  element_text(angle=90, vjust=0.5, size=16),
                           axis.title.x = element_text(size=0),
                           axis.title.y = element_text(size=17)) + 
        ggtitle(paste0(gene, " GSE25219 (BrainSpan)")) + 
        scale_fill_manual(name="",values=cols)
    
    ylim1 = boxplot.stats(allExpress$log2Exp)$stats[c(1, 5)]
    ylim1[2] = ylim1[2]*1.05
    ylim1[1] = ylim1[1]*0.95
    
    p= p + coord_cartesian(ylim = ylim1)
    ggsave(paste0('analysis/99.GeneAnalysisPipeline/',gene,'/plain accross reg','.png'),plot=p,height=8,width=7)
    
    # large region dataset------
    
    
    files = list.files('analysis//99.GeneAnalysisPipeline/HumanRegionExpr/',full.names=T)
    #files = files[grepl(regexMerge(c('frontal','hippocampus')), files)]
    outliers = read.table('analysis/04.Human Coexpression/outliers')[,1]
    
    allExpress =sapply(files, function(i){
        orbitalDat = fread(i)
        orbitalDat = orbitalDat[!is.na(orbitalDat$Gene.Symbol),]
        orbitalDat %<>% mostVariable
        list[orbGenes, orbExp] = sepExpr(orbitalDat)
        orbExp = orbExp[,!grepl(regexMerge(outliers),colnames(orbExp)),with=F]
        unlist(orbExp)
    })
    allExpress = melt(allExpress)
    names(allExpress) = c('log2Exp', 'BrainReg')
    
    allExpress[,2] = basename(allExpress[,2])
    allExpress[allExpress[,2] %in% 'thalamus',2] = 'Thalamus'
    allExpress[allExpress[,2] %in% 'occipital cortex',2] = 'Occipital Cortex'
    allExpress[allExpress[,2] %in% 'temporal cortex',2] = 'Temporal Cortex'
    allExpress[allExpress[,2] %in% 'hippocampus',2] = 'Hippocampus'
    allExpress[allExpress[,2] %in% 'cerebellar cortex',2] = 'Cerebellar cortex'
    
    allExpress = allExpress[allExpress[,2] %in% c('Cerebellar cortex',
                                                  'Hippocampus',
                                                  'Thalamus',
                                                  'Occipital Cortex',
                                                  'Temporal Cortex',
                                                  'Striatum'),]
    
    
    
    allExpress[,2] = factor(allExpress[,2],levels=c('Cerebellar cortex',
                                                    'Hippocampus',
                                                    'Thalamus',
                                                    'Occipital Cortex',
                                                    'Temporal Cortex'))
    
    
    oldData = sapply(files,
                     function(i){
                         orbitalDat = fread(i)
                         orbitalDat = orbitalDat[!is.na(orbitalDat$Gene.Symbol),]
                         #orbitalDat = read.exp('Data/DevelopData/OFC')
                         
                         list[orbGenes, orbExp] = sepExpr(orbitalDat)
                         orbExp = orbExp[,!grepl(regexMerge(outliers),colnames(orbExp)),with=F]
                         unlist(orbExp[grep(gene, orbGenes$Gene.Symbol),])
                     })  
    names(oldData) = basename(files)
    oldData =  melt(oldData)
    x = oldData
    colnames(x) = c('log2Exp', 'BrainReg')
    x[,2] = basename(x[,2])
    
    
    x[x[,2] %in% 'thalamus',2] = 'Thalamus'
    x[x[,2] %in% 'occipital cortex',2] = 'Occipital Cortex'
    x[x[,2] %in% 'temporal cortex',2] = 'Temporal Cortex'
    x[x[,2] %in% 'hippocampus',2] = 'Hippocampus'
    x[x[,2] %in% 'cerebellar cortex',2] = 'Cerebellar cortex'
    
    x = x[x[,2] %in% c('Cerebellar cortex',
                       'Hippocampus',
                       'Thalamus',
                       'Occipital Cortex',
                       'Temporal Cortex',
                       'Striatum'),]
    
    x[,2] = factor(x[,2],levels=c('Cerebellar cortex',
                                  'Hippocampus',
                                  'Thalamus',
                                  'Occipital Cortex',
                                  'Temporal Cortex'))
    colnames(x) = c('log2Exp', 'BrainReg')
    oldData = x
    
    cols <-teval( paste0('c("Expression of All Genes"="lightgray",',gene,'="white")'))
    
    p = ggplot(x , aes(x= BrainReg, y= log2Exp))  + 
        geom_boxplot(data=allExpress, color ='gray',aes(fill ='Expression of All Genes'),outlier.size=0, width= 0.3) +
        geom_boxplot(aes(fill = gene), outlier.size=0) + 
        geom_point(position='jitter') + scale_y_continuous(paste(gene,'expession')) + 
        theme_bw() + theme(axis.text.x =  element_text(angle=90, vjust=0.5, size=16),
                           axis.title.x = element_text(size=0),
                           axis.title.y = element_text(size=17)) + 
        ggtitle(paste0(gene, " GSE60862")) + 
        scale_fill_manual(name="",values=cols)
    
    ylim1 = boxplot.stats(allExpress$log2Exp)$stats[c(1, 5)]
    ylim1[2] = ylim1[2]*1.05
    ylim1[1] = ylim1[1]*0.95
    
    p= p + coord_cartesian(ylim = ylim1)
    ggsave(paste0('analysis/99.GeneAnalysisPipeline/',gene,'/plain accross reg Large Dataset','.png'),plot=p,height=8,width=7)
    
    
    
}