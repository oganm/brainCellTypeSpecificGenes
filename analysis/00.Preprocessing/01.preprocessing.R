# Preprocessing of the data. The output files of this directory are also present in the data folder.
library(ogbox)
source('R/readDesignMergeCel.R')
source('R/quantileNormalize.R')
source('R/mostVariable.R')
source('R/sexFind.R')


# this variable is there in case you just make some minor changes in the normalization process and want to create the
# meltedDesign file again. If you flip this to false, normalization won't occur and only the design file will be 
# generated
skipNorm = F
# this variable dictates if gemma annotations and homology information will be downloaded. only set this true when you 
# are running the pipeline for the first time
skipDown = F

# this is the name of your main design file
desFile = 'data/Design.tsv'
# this controls which collumn will be used to normalize. Just have a collumn that has T or F  in it depending on what
# you want to do with that particular sample. Only useful if you want to keep things in your list but you don't actually 
# want to use it
normalize = 'Normalize'
normalize2 = 'Normalize2.0'

# if all your names start with GSM all you need is (GSM.*?(?=,|$)). Chooses file identifiers from the first collumn.
# must match your filenames
celRegex='(GSM.*?(?=,|$))|(PC\\d....)|(Y[+].*?((?=(,))|\\d+))|((?<=:)|(?<=[,]))A((9)|(10))_[0-9]{1,}_Chee_S1_M430A|(v2_(?![G,H,r]).*?((?=(,))|($)))|(SSC.*?((?=(,))|($)))|(MCx.*?((?=(,))|($)))|(Cbx.*?((?=(,))|($)))'
# which directory to look for cel files
celDir = 'data/cel'
# outFolder = 'analysis/00.Preprocessing/output/'
outFolder = 'data/'



dir.create(outFolder, showWarnings=F, recursive=T)
# normalization ----
if (skipNorm == F){
    readDesignMergeCel(desFile, normalize, celRegex, celDir, 
                       paste0(outFolder, 'rmaExp.csv'),
                       paste0(outFolder,'meltedDesign.tsv'))
    
    readDesignMergeCel(desFile, normalize2, celRegex, celDir,
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
    meltDesign(desFile, normalize, celRegex, paste0(outFolder,'/','finalExp.csv'), paste0(outFolder,'/meltedDesign.tsv'))
    meltDesign(desFile, normalize2, celRegex, paste0(outFolder,'/','finalExp2.csv'), paste0(outFolder,'/meltedDesign2.tsv'))
    
}


sexFind(paste0(outFolder,'/meltedDesign.tsv'),
        paste0(outFolder,'/meltedDesign.tsv'),
        paste0(outFolder,'/',finalExp))

sexFind(paste0(outFolder,'/meltedDesign2.tsv'),
        paste0(outFolder,'/meltedDesign2.tsv'),
        paste0(outFolder,'/','finalExp2.csv'))

