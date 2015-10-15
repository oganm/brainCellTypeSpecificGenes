library(parallel)
library(ogbox)
library(reshape2)
library(ggplot2)
coexpBoxViol = function(coexpData,cores = 8){
    
    pList = mclapply(coexpData, function(i){
        print('hey!')
        pVal = sapply(i$coexpressions, function(coexpressions){
            print('ho!')
            y = coexpressions
            x = i$daCorvec[!i$daCorvec %in% y]
            print('yo!')
            test = wilcox.test(x=x, y=y ,alternative='less')
            print('madafaka!')
            return(test$p.value)  
        })
        return(pVal)
    },mc.cores=cores)
    
    pAddress = sapply(coexpData,function(x){paste0(x$plotOut,'/',x$subFolderName)})
    pAddress = gsub('Plots','Data', pAddress)
    sapply(1:len(pList), function(i){
        write.table(pList[[i]],file=paste0(pAddress[i],'/','pVals'),col.names=F,quote=F)
    })
    
    
    
    pAdjusted  = relist(flesh=p.adjust(unlist(pList),method='fdr'),
                        skeleton=pList)
    
    
    
    
    # the boxviolplot and significance test
    plotVioBox = function(relevants, ps, all){
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
            geom_signif(c(1,ps),maxY=yMax)
        return(plot)
    }
    
    set.seed(1)
    mclapply(1:len(coexpData), function(i){
        relevants = coexpData[[i]]$coexpressions
        ps = pAdjusted[[i]]
        relevants = melt(relevants)
        tempSamp = sample(coexpData[[i]]$daCorvec,2000)
        all = data.frame(value = tempSamp, L1 = 'all')
        plot = plotVioBox(relevants, ps, all)
        ggsave(paste0(coexpData[[i]]$plotOut,'/',coexpData[[i]]$subFolderName,'/','boxViol.png'),plot,height=8,width=2*(1+len(coexpData[[i]]$coexpressions)))
        
    })
    
}