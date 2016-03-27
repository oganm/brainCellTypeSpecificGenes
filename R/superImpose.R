superImpose = function(genes,outDir,PC=1){
    table = {if(!is.null(outDir)){
        paste0(outDir,'/rotTable_', names(genes))
    }else {
        NULL}
    }
    estimates = cellTypeEstimate(exprData = expression, 
                                 genes= genes,
                                 geneColName = 'GeneSym',
                                 outlierSampleRemove = F,
                                 synonymTaxID = NULL,
                                 geneTransform = NULL, 
                                 groups = NA,
                                 controlBased = NULL, 
                                 tableOut = table,
                                 indivGenePlot = NULL,
                                 seekConsensus = F,
                                 plotType = NULL,
                                 PC = PC)
    
    dir.create(outDir,showWarnings=FALSE)
    for (i in 1:len(estimates$estimates)){
        if(is.na(estimates$estimates[[i]][1])){
            next
        }
        counts = theirPred %>%
            filter(Cell.type == names(estimates$estimates[i])) %>%
            select(RealCounts)
        cybersort = theirPred %>%
            filter(Cell.type == names(estimates$estimates[i])) %>%
            select(CibersortPred)
        
        estimates$estimates[[i]] = scale01(estimates$estimates[[i]])
        
        scaleCyber = scaleToInt(cybersort,max(estimates$estimates[[i]]), min(estimates$estimates[[i]]))
        
        frame1 = data.frame(estimates = estimates$estimates[[i]], counts, scaleCyber)
        frame2 = data.frame(cybersort, counts)
        
        
        frame = data.frame(estimates = estimates$estimates[[i]],cybersort ,counts)
        
        p1 = ggplot(frame1, aes(y= estimates,x = RealCounts)) + geom_point(size = 4, shape = 15) +  theme_bw() +
            geom_segment(aes(x=RealCounts,
                             y = estimates,
                             xend = RealCounts,
                             yend = CibersortPred)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            annotate('text', y= max(frame1$estimates),x= min(frame1$RealCounts),vjust = 1, hjust = 0, size = 8,
                     label = paste0('spearman correlation:\n',
                                    'Marker genes: ', format(cor(frame1$estimates,frame1$RealCounts,method='spearman'),digits=3), '\n',
                                    'Cibersort: ',format(cor(frame1$CibersortPred,frame1$RealCounts,method='spearman'),digits=3))
            ) + geom_line(size=3,alpha=1/10) + ggtitle(names(estimates$estimates[i])) + theme(axis.title= element_text(size=17))
        
        p2 = ggplot(frame2, aes(y= CibersortPred ,x = RealCounts)) + geom_point(size = 4, color='red',shape=16) + geom_line(size=3,color='red',alpha=1/10) + theme_bw() %+replace% 
            theme(panel.background = element_rect(fill = NA)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title= element_text(size=17)) +  ggtitle(names(estimates$estimates[i]))
        
        
        g1 <- ggplot_gtable(ggplot_build(p1))
        g2 <- ggplot_gtable(ggplot_build(p2))
        
        pp <- c(subset(g1$layout, name == "panel", se = t:r))
        g <- gtable_add_grob(g1, g2$grobs[[which(g2$layout$name == "panel")]], pp$t, 
                             pp$l, pp$b, pp$l)
        
        ia <- which(g2$layout$name == "axis-l")
        ga <- g2$grobs[[ia]]
        ax <- ga$children[[2]]
        ax$widths <- rev(ax$widths)
        ax$grobs <- rev(ax$grobs)
        ax$grobs[[1]]$x <- ax$grobs[[1]]$x - unit(1, "npc") + unit(0.15, "cm")
        g <- gtable_add_cols(g, g2$widths[g2$layout[ia, ]$l], length(g$widths) - 1)
        g <- gtable_add_grob(g, ax, pp$t, length(g$widths) - 1, pp$b)
        png(paste0(outDir,names(estimates$estimates[i]),'.png'),width=600,height=600)
        grid.draw(g)
        dev.off()
    }
}