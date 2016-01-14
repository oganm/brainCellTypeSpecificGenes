
# preprocessing of the human data -----------------
library(ogbox) 

source('R/readHumanCel.R')
softFile = read.design('data/GSE60862_meta')

softFile = softFile[!is.na(softFile$pH) 
                    & !softFile$deathCause=='Cancer',]
# regions = unique(softFile$brainRegion)
# dir.create('data/GSE60862_Expression Matrices',showWarnings=F)
# 
# for (i in regions){
#     regionSet = softFile[softFile$brainRegion == i,]
#     readHumanCel(regionSet$GSM,paste0('data/GSE60862_Expression Matrices/',i),humanDir='data/cel//GPL5175')
# }


# shitloads of time and memory! 
readHumanCel(softFile$GSM,paste0('data/GSE60862_Expression Matrices/','all'),humanDir='data/cel//GPL5175')
