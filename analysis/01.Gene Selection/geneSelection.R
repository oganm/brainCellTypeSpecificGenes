# no rotation gene selection
library(ogbox)
source('R/geneSelect.R')
source('R/microglialException.R')
source('R/rotateSelect.R')

dataDir = 'data/'
geneOut = 'analysis//01.Gene Selection/Fold'
rotationOut = 'analysis//01.Gene Selection/Rotation'
rotSelOut = 'analysis/01.Gene Selection/RotSel'

groupNames = c('GabaDeep','PyramidalDeep')
regionNames = 'Region'
dir.create(geneOut, showWarnings=F, recursive=T)


geneSelect(paste0(dataDir,'/meltedDesign.tsv'),
           paste0(dataDir,'/','finalExp.csv'),
           geneOut,
           groupNames,
           regionNames,
           cores=16#,
           #debug='Forebrain_GabaDeep'
           )

geneSelect(paste0(dataDir,'/meltedDesign2.tsv'),
           paste0(dataDir,'/','finalExp2.csv'),
           paste0(geneOut,2),
           groupNames,
           regionNames,
           cores=16#,
           #debug = 'Cerebellum_PyramidalDeep'
           )


#rotation gene selection

for (i in 1:500){
    print(i)
    geneSelect(paste0(dataDir,'/meltedDesign.tsv'),
               paste0(dataDir,'/','finalExp.csv'),
               paste0(rotationOut,'/',i),
               groupNames,
               regionNames,
               rotate=0.33,
               cores=15#,
               #debug="Cortex_GabaDeep"
    )
    #microglialException(paste0(rotationOut,'/',i),cores=15)
}

for (i in 1:500){
    print(i)
    geneSelect(paste0(dataDir,'/meltedDesign2.tsv'),
               paste0(dataDir,'/','finalExp2.csv'),
               paste0(rotationOut,'2','/',i),
               groupNames,
               regionNames,
               rotate=0.33,
               cores=15)
    #microglialException(paste0(rotationOut,'2','/',i),cores=16)
    
}


rotateSelect(rotationOut,rotSelOut, cores = 15)
rotateSelect(paste0(rotationOut,2),paste0(rotSelOut,2), cores = 15)

allGenes = list(genes1 = allPuristOut(paste0(rotSelOut,'/Relax')) , genes2 = allPuristOut(paste0(rotSelOut,'2','/Relax')))
for (n in 1:len(allGenes)){
    genes = allGenes[[n]]
    for (i in 1:len(genes)){
        pieces = strsplit(names(genes)[i],'_')[[1]]
        if (is.na(pieces[2])){
            pieces[2] = pieces[1]
            pieces[1] ='All'
        }
        dir.create(paste0('analysis//01.Gene Selection/FinalGenes',n,'/',
                          pieces[2] , '/' , pieces[1]), 
                   showWarnings=F, recursive=T)
        
        
        for (j in 1:len(genes[[i]])){
            write.table(genes[[i]][[j]],
                        paste0('analysis/01.Gene Selection/FinalGenes',n,'/',
                               pieces[2],'/',pieces[1],'/', 
                               names(genes[[i]])[j]),
                        row.names=F,
                        quote=F,
                        col.names=F
                        
            )
        }
    }
}

microglialException('analysis/01.Gene Selection/FinalGenes1/',cores=16)
microglialException('analysis/01.Gene Selection/FinalGenes2/',cores=16)
file.rename('analysis/01.Gene Selection/FinalGenes1/','analysis/01.Gene Selection/FinalGenes/')
