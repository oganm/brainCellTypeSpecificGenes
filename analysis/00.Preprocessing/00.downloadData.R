# overall this will take a huge amount of time to run.

library(ogbox)
source('R/softParser.R')
library(oligoClasses)
# download gemma annotations ------
dir.create('data/GemmaAnnots/', showWarnings=F)
getGemmaAnnotGoogle('GPL1261','data/GemmaAnnots/GPL1261')
getGemmaAnnotGoogle('GPL339','data/GemmaAnnots/GPL339')
getGemmaAnnotGoogle('GPL570','data/GemmaAnnots/GPL570')
getGemmaAnnotGoogle('GPL96','data/GemmaAnnots/GPL96')
getGemmaAnnotGoogle('GPL5175','data/GemmaAnnots/GPL5175')


# download homology information ----------
# used internally in homogene functions. I know its bad design but...
sourceGithub(OganM,toSource,homologene)
homoloGeneTarget = 'data/homogene.tsv'
setUpHomologene()


# mouse rna seq data download 02 ------
download.file(url= 'http://linnarssonlab.org/blobs/cortex/expression_mRNA_17-Aug-2014.txt', destfile='data/mouseRNASeq_Zeisel 2015.txt')

# human rna seq data download 03-----
gsms = gsmFind('GSE67835')
dir.create('data/humanRNASeq')

sapply(gsms, function(gsm){
    print(gsm)
    page = getURL(paste0('http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=',gsm))
    fileURL = URLdecode(str_extract(page,'ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM.*?csv%2Egz'))
    if (len(fileURL) == 0){
        if (warnings){
            warning(paste(gsm,"doesn't have a file attached"))
        }
        return(invisible(F))
    }
    download.file(fileURL,paste0('data/humanRNASeq/',gsm,'.csv.gz'))
    system(paste0('gunzip -f "',paste0('data/humanRNASeq/',gsm,'.csv.gz'),'"'))  
})

files = list.files('data/humanRNASeq/', full.names=T)

allExpr = sapply(files,function(x){
    read.table(x, sep='\t')[,2]
})
print('files read')

singleGenes = read.table(files[1], sep='\t', stringsAsFactors=F)[,1]
singleGenes = sapply(singleGenes, trimWS)
rownames(allExpr) = singleGenes
colnames(allExpr) = gsub('[.]csv','',basename(colnames(allExpr)))
write.csv(allExpr,'data/humanRNASeq.csv',quote=F)



# get human microarray data -----------
download.file('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60862/soft/GSE60862_family.soft.gz',
              destfile='data/GSE60862_family.soft.gz')
system('gunzip data/GSE60862_family.soft.gz')

# deal with softfile
softData = softParser('data/GSE60862_family.soft')

softData = softData[,c('!Sample_characteristics_ch1 = age at death (in years)', 
                       '!Sample_characteristics_ch1 = ancestry',
                       '!Sample_characteristics_ch1 = brain bank',
                       '!Sample_characteristics_ch1 = brain ph',
                       '!Sample_characteristics_ch1 = brain region',
                       "!Sample_characteristics_ch1 = Cause of death",
                       "!Sample_characteristics_ch1 = disease status",
                       "!Sample_characteristics_ch1 = individual id",
                       "!Sample_characteristics_ch1 = post-mortem interval (in hours)",
                       "!Sample_characteristics_ch1 = rin",
                       "!Sample_characteristics_ch1 = Sex",
                       '!Sample_title',
                       '!Sample_platform_id')]
softData$GSM = rownames(softData)

colnames(softData) = c('deatAge',
                       'ancestry',
                       'brainBank',
                       'pH',
                       'brainRegion',
                       'deathCause',
                       'disease',
                       'invidID',
                       'postMortemInt',
                       'rnaInt',
                       'sex',
                       'title',
                       'platform',
                       'GSM')


# download cel files
write.design(softData,file='data/GSE60862_meta')
# softData = read.design('data/GSE60862_meta')
dir.create('data/cel/GPL5175')
sapply(softData$GSM,function(x){
    xNow<<-x
    gsmDown(gsm=x,outfile=paste0('data//cel/GPL5175/',x),unzip=F)
})

# attach scandates
softData$scanDate = sapply(softData$GSM, function(x){
    celfileDate(paste0('/home/omancarci/masterOfCellTypes/cel/GPL5175/',x, '.cel.gz'))
})
write.design(softData,file='data/GSE60862_meta')


