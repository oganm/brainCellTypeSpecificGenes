
regionAnalysis = function(brains, design, genes, outDir){
    dir.create(outDir, showWarnings=F)
    med = brains %>% memoMatrix %>% memoMedian
    print('median done')
    expression = cbind(probes$gene_symbol, brains)
    expression  = memoMostVariable(expression,'probes$gene_symbol',
                                   treshold=med)
    colnames(expression)[1] = 'Gene.Symbol'
    
    estimates = cellTypeEstimate(exprData = expression, 
                                 genes= genes,
                                 geneColName = 'Gene.Symbol',
                                 outlierSampleRemove = F,
                                 synonymTaxID = NULL,
                                 geneTransform = NULL, 
                                 groups = NA,
                                 controlBased = NULL, 
                                 tableOut = NULL,
                                 indivGenePlot = NULL,
                                 seekConsensus = F,
                                 plotType = NULL)

    estimData  = cbind(design[c('structure_id', 'structure_acronym', 'structure_name', 'donor')],
                       (estimates$estimates %>% as.data.frame %>% apply(2,scale01)))
    
    # heatmap
    estimData %<>% arrange(donor,structure_id)
    list[design, estimates] = sepExpr(estimData)
    estimates %<>% apply(2,scale01) %>% round(digits=14)
    colors = as.matrix(cbind(toColor(design$structure_id)$cols, toColor(design$donor)$cols))
    png(paste0(outDir,'/donorOrder.png'),1000,1000)
    heatmap.3(as.matrix(estimates), trace='none',Rowv=F,Colv=F, col= viridis(10),RowSideColors=t(colors), 
              dendrogram='none')
    dev.off()
    
    estimData %<>% arrange(structure_id,donor)
    list[design, estimates] = sepExpr(estimData)
    estimates %<>% apply(2,scale01) %>% round(digits=14)
    colors = as.matrix(cbind(toColor(design$structure_id)$cols, toColor(design$donor)$cols))
    png(paste0(outDir,'/regionOrder.png'),1000,1000)
    heatmap.3(as.matrix(estimates), trace='none',Rowv=F,Colv=F, col= viridis(10),RowSideColors=t(colors), 
              dendrogram='none')
    dev.off()
    
    # scatterplot
    
    frame = estimData %>% melt(id.vars = c('structure_id', 'structure_acronym', 'structure_name', 'donor'),
                               value.name = 'estimate', variable.name = 'cellType')
    
    frame$donor = as.factor(frame$donor)
    frame$structure_id = as.factor(frame$structure_id)
    
    p = ggplot(frame, aes(x = structure_acronym, y = estimate, color = donor)) + facet_wrap(~cellType) + geom_point() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
    ggsave(filename=paste0(outDir, '/scatter.png'),
           plot=p,
           height = 6,
           width = 12)
    
    # what is the order of regions? is it consistent
    regionSimplify = function(x){
        list[des,estim] = sepExpr(x)
        estim %<>% apply(2,mean) %>% as.df %>% t
        cbind(des[1,,drop=F],estim)
    }
    
    
     rankFrame = frame %>%
        group_by(structure_id, donor, cellType) %>%
        do(regionSimplify(.)) %>%
        group_by(donor, cellType) %>%
        mutate(rank = rank(estimate)) %>% 
        group_by(donor) 
    

    p = rankFrame %>% ggplot(aes(y = structure_acronym, x = donor)) + 
        geom_tile(aes(fill = rank)) +
        scale_fill_viridis() + 
        facet_wrap(~cellType)
    
    ggsave(filename=paste0(outDir, '/rankHeatmap.png'),
           plot=p,
           height = 12,
           width = 12)
    
    
}
