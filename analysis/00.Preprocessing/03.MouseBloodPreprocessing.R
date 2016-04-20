library(ogbox)
library(praise)
library(magrittr)
source('R/readDesignMergeCel.R')
source('R/quantileNormalize.R')
source('R/mostVariable.R')



outFolder  = 'data/bloodCellType/mouseBlood/'
celDir = 'data/cel/'



readDesignMergeCel('data/bloodCellType/mouseBlood/bloodDes.tsv', 'Normalize', '(GSM.*?(?=,|$))', celDir, 
                   paste0(outFolder, 'mouseBlood.csv'),
                   paste0(outFolder,'mouseBloodMelted.tsv'))

quantileNorm(paste0(outFolder,'/mouseBlood.csv'),
             paste0(outFolder,'/','qnormExp.csv'))



mostVariableCT(paste0(outFolder,'/',"qnormExp.csv"),
               paste0(outFolder,'/','finalExp.csv'),
               selectionNaming = 'lm11',
               design=paste0(outFolder,'/mouseBloodMelted.tsv'))

