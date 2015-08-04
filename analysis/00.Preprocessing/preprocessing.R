# Preprocessing of the data. The output files of this directory are also present in the data folder.
library(ogbox)
source('R/readDesignMergeCel.R')
source('R/quantileNormalize.R')
source('R/mostVariable.R')
source('R/sexFind.R')


skipNorm = F
desFile = 'data/Design.tsv'
namingCol = 'Normalize'
namingCol2 = 'Normalize2.0'
celRegex='(GSM.*?(?=,|$))|(PC\\d....)|(Y[+].*?((?=(,))|\\d+))|((?<=:)|(?<=[,]))A((9)|(10))_[0-9]{1,}_Chee_S1_M430A|(v2_(?![G,H,r]).*?((?=(,))|($)))|(SSC.*?((?=(,))|($)))|(MCx.*?((?=(,))|($)))|(Cbx.*?((?=(,))|($)))'
celDir = 'data/cel'
tinyChip = 'mouse430a2.db'
# outFolder = 'analysis/00.Preprocessing/output/'
outFolder = 'data/'


dir.create(outFolder, showWarnings=F, recursive=T)
# normalization ----
if (skipNorm == F){
    readDesignMergeCel(desFile, namingCol, celRegex, celDir,tinyChip, 
                       paste0(outFolder, 'rmaExp.csv'),
                       paste0(outFolder,'meltedDesign.tsv'))
    
    readDesignMergeCel(desFile, namingCol2, celRegex, celDir,tinyChip,
                       paste0(outFolder, 'rmaExp2.csv'),
                       paste0(outFolder,'meltedDesign2.tsv'))
    
    
    quantileNorm(paste0(outFolder,'/rmaExp.csv'),
                 paste0(outFolder,'/','qnormExp.csv'))
    quantileNorm(paste0(outFolder,'/rmaExp2.csv'),
                 paste0(outFolder,'/',"qnormExp2.csv"))
    # system(paste0('gunzip ',outFolder,'/',qnormExp))
    
    mostVariableCT(paste0(outFolder,'/','qnormExp.csv'),
                   paste0(outFolder,'/','finalExp.csv'),
                   selectionNaming = 'GabaDeep',
                   design = paste0(outFolder,'/meltedDesign.tsv'))
    
    mostVariableCT(paste0(outFolder,'/',"qnormExp2.csv"),
                   paste0(outFolder,'/','finalExp2.csv'),
                   selectionNaming = 'GabaDeep',
                   design=paste0(outFolder,'/meltedDesign2.tsv'))
    
    system(paste0('gzip -f ',outFolder,'/','qnormExp.csv'))
    system(paste0('gzip -f ',outFolder,'/rmaExp.csv'))
    system(paste0('gzip -f ',outFolder,'/','qnormExp2.csv'))
    system(paste0('gzip -f ',outFolder,'/','rmaExp2.csv'))
    
    #     system(paste0('gunzip -f ',outFolder,'/',qnormExp,'.gz'))
    #     system(paste0('gunzip -f ',outFolder,'/rmaExp.csv.gz'))
    #     system(paste0('gunzip -f ',outFolder,'/','qnormExp2.csv.gz'))
    #     system(paste0('gunzip -f ',outFolder,'/','rmaExp2.csv.gz'))
}

if (skipNorm == T){
    source('readDesignMergeCel.R')
    meltDesign(desFile, namingCol, celRegex, paste0(outFolder,'/',finalExp), paste0(outFolder,'/meltedDesign.tsv'))
    meltDesign(desFile, namingCol2, celRegex, paste0(outFolder,'/','finalExp2.csv'), paste0(outFolder,'/meltedDesign2.tsv'))
    
}


sexFind(paste0(outFolder,'/meltedDesign.tsv'),
        paste0(outFolder,'/meltedDesign.tsv'),
        paste0(outFolder,'/',finalExp))

sexFind(paste0(outFolder,'/meltedDesign2.tsv'),
        paste0(outFolder,'/meltedDesign2.tsv'),
        paste0(outFolder,'/','finalExp2.csv'))

