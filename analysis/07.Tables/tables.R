# table 1
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

fileCon = file('analysis/07.Tables//table_1.md')
writeLines(kable(toPrint, align= 'c'), fileCon)
close(fileCon)

# before our additions
before = design[1:136,]
toPrint = table(before[c('PyramidalDeep','Method')])>1

toPrint = data.frame(toPrint, Studies =apply(table(before[c('PyramidalDeep','Reference')])>1,1,sum))
toPrint[,1:6] = apply(toPrint[1:6], 2, function(x){
    sapply(x,function(x){
        if (x){
            "✔"
        } else {
            ""
        }
    })
})

fileCon = file('analysis/07.Tables//table_1.5.md')
writeLines(kable(toPrint, align= 'c'), fileCon)
close(fileCon)
