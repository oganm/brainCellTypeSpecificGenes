library(knitr)

design = read.design('data//meltedDesign.tsv')
design = design[!is.na(design$PyramidalDeep),]
toPrint = table(design[c('PyramidalDeep','Method')])>1

toPrint = data.frame(toPrint, Studies =apply(table(design[c('PyramidalDeep','Reference')])>1,1,sum))
toPrint[,1:6] = apply(toPrint[1:6], 2, function(x){
    sapply(x,function(x){
        if (x){
            "✔"
        } else {
            ""
        }
    })
})
kable(toPrint, align= 'c')

