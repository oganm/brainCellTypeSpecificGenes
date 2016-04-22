# no rotation gene selection
args <- commandArgs(trailingOnly = T)
start = as.numeric(args[1])
end = as.numeric(args[2])

library(ogbox)
source('R/geneSelect.R')
source('R/microglialException.R')
source('R/rotateSelect.R')

dataDir = 'data/'
geneOut = 'analysis//01.Gene Selection/Fold'
rotationOut = 'analysis//01.Gene Selection/Rotation'
rotSelOut = 'analysis/01.Gene Selection/RotSel'

#groupNames = c('GabaDeep','PyramidalDeep','MajorType')
# groupNames = c('PyramidalDeep','MajorType')
groupNames = 'DopaSelect'


regionNames = 'Region'
# many steps requires parallel processes. set to cores cores by default
cores = 15
dir.create(geneOut, showWarnings=F, recursive=T)

if (start == 1){
    geneSelect(paste0(dataDir,'/meltedDesign.tsv'),
               paste0(dataDir,'/','finalExp.csv'),
               geneOut,
               groupNames,
               regionNames,
               cores=cores#,
               #debug='Forebrain_GabaDeep'
    )
    
    geneSelect(paste0(dataDir,'/meltedDesign2.tsv'),
               paste0(dataDir,'/','finalExp2.csv'),
               paste0(geneOut,2),
               groupNames,
               regionNames,
               cores=cores#,
               #debug = 'Cerebellum_PyramidalDeep'
    )
}

#rotation gene selection

for (i in start:end){
    print(i)
    geneSelect(paste0(dataDir,'/meltedDesign.tsv'),
               paste0(dataDir,'/','finalExp.csv'),
               paste0(rotationOut,'/',i),
               groupNames,
               regionNames,
               rotate=0.33,
               cores=cores#,
               #debug="Cortex_GabaDeep"
    )
    #microglialException(paste0(rotationOut,'/',i),cores=cores)
}

for (i in start:end){
    print(i)
    geneSelect(paste0(dataDir,'/meltedDesign2.tsv'),
               paste0(dataDir,'/','finalExp2.csv'),
               paste0(rotationOut,'2','/',i),
               groupNames,
               regionNames,
               rotate=0.33,
               cores=cores)
    #microglialException(paste0(rotationOut,'2','/',i),cores=cores)
    
}


