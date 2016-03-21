# Some human samples are mislabeled so this analysis needs an outlier removal step
library(limma)
library(gplots)

humanExp = fread('data/UCLregions/GSE60862_expression')
humanExp = humanExp[!is.na(humanExp$Gene.Symbol),]
humanExp = humanExp[humanExp$Gene.Symbol!='',]

list[genes,expr] = sepExpr(humanExp)
colnames(expr) = gsub('\\.cel\\.gz','',colnames(expr))
softFile = read.design('data/UCLregions/GSE60862_meta')

groups = unique(softFile$brainRegion)
pairwise = combn(groups,2)

difs = vector(mode = 'list', length = ncol(pairwise))

# get differentially expressed genes between any two regions
dir.create('analysis/04.Human Coexpression/regionalHeatmaps')
for (i in 1:ncol(pairwise)){
    subsetExpr = expr[,names(expr) %in% softFile$GSM[softFile$brainRegion %in% pairwise[,i]],with=F]
    subsetExpr = as.matrix(subsetExpr)
    rownames(subsetExpr) = genes$Gene.Symbol
    
    subsetGroups = softFile$brainRegion[match(colnames(subsetExpr), softFile$GSM)]  
    
    mm = model.matrix(~ subsetGroups,data.frame(subsetGroups))
    fit <- lmFit(subsetExpr, mm)
    fit <- eBayes(fit)
    dif = topTable(fit, coef=colnames(fit$design)[2],
                   lfc = log(5,base=2),
                   number = Inf, 
                   p.value = 1e-5)
    if (nrow(dif)==0){
        print(paste('no genes found for',paste(pairwise[,i],collapse='/') ,'pair'))
        next
    }
    difs[[i]] = dif$ID
    toPlot = subsetExpr[rownames(subsetExpr) %in% c(dif$ID),]
    colnames(toPlot) = paste(#softFile$IndID[softFile$GSM %in% colnames(subsetExpr)],
        substr(softFile$bank[softFile$GSM %in% colnames(subsetExpr)],1,3),
        substr(softFile$Region[match(colnames(subsetExpr), softFile$GSM)],1,3),
        substr(softFile$GSM[match(colnames(subsetExpr), softFile$GSM)],7,10))
    
    png(filename= paste0('analysis/04.Human Coexpression/regionalHeatmaps/',paste(pairwise[,i],collapse='-'),'.png'),width = 2200, height = 1500)
    par(cex.main=3)
    heatmap.2(t(scale(t(toPlot))),trace='none', distfun = function(x){dist(x)^2},cexCol=1, 
                       hclustfun=function(x){hclust(x,method='ward.D2')}, margins= c(9,5),
                       ColSideColors = toColor(subsetGroups)$cols,
                       main = paste(pairwise[,i],collapse='-'),symbreaks=T)
    legend('bottomleft', legend= names(toColor(subsetGroups)$palette),
           fill = toColor(subsetGroups)$palette,cex = 2)
    
    dev.off()
    print(i)
}

# detect outliers based on differentially expressed genes inside a single group
outliers = lapply(1:len(groups), function(i){
    relevant = unlist(difs[apply(pairwise,2,function(x){
        groups[[i]] %in%  x  
    })])
    subsetExpr = expr[,names(expr) %in% softFile$GSM[softFile$brainRegion %in% groups[i]],with=F]
    subsetExpr = subsetExpr[genes$Gene.Symbol %in% relevant,]
    pca = prcomp(t(subsetExpr))
    names(boxplot(pca$x[,1])$out)
})

write.table(unlist(outliers),'analysis//04.Human Coexpression/outliers', row.names=F,col.names=F,quote=F)
