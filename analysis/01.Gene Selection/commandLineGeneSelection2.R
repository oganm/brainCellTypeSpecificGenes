



library(ogbox)
source('R/geneSelect.R')
source('R/microglialException.R')
source('R/rotateSelect.R')

dataDir = 'data/'
geneOut = 'analysis//01.Gene Selection/Fold'
rotationOut = 'analysis//01.Gene Selection/Rotation'
rotSelOut = 'analysis/01.Gene Selection/RotSel'
cores = 15

rotateSelect(rotationOut,rotSelOut, cores = cores)
rotateSelect(paste0(rotationOut,2),paste0(rotSelOut,2), cores = cores)


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

microglialException('analysis/01.Gene Selection/FinalGenes1/',cores=cores)
microglialException('analysis/01.Gene Selection/FinalGenes2/',cores=cores)
system('rm -rf "analysis/01.Gene Selection/FinalGenes/"')
file.rename('analysis/01.Gene Selection/FinalGenes1/','analysis/01.Gene Selection/FinalGenes/')
