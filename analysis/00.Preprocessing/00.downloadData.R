# overall this will take a huge amount of time to run.

library(dplyr)
library(ogbox)
library(oligoClasses)
library(affy)
source('R/mostVariable.R')
library(XLConnect)
library(RCurl)
# download gemma annotations ------
dir.create('data/GemmaAnnots/', showWarnings=F)
getGemmaAnnotGoogle('GPL1261','data/GemmaAnnots/GPL1261',annotType='noParents')
getGemmaAnnotGoogle('GPL339','data/GemmaAnnots/GPL339',annotType='noParents')
getGemmaAnnotGoogle('GPL570','data/GemmaAnnots/GPL570',annotType='noParents')
getGemmaAnnotGoogle('GPL96','data/GemmaAnnots/GPL96',annotType='noParents')
getGemmaAnnotGoogle('GPL97','data/GemmaAnnots/GPL97',annotType='noParents')
getGemmaAnnotGoogle('GPL5175','data/GemmaAnnots/GPL5175',annotType='noParents')
getGemmaAnnotGoogle('GPL6884','data/GemmaAnnots/GPL6884',annotType='noParents')



# mouse rna seq data download 02 ------
download.file(url= 'http://linnarssonlab.org/blobs/cortex/expression_mRNA_17-Aug-2014.txt', 
              destfile='data/linnarsonSingleCell/mouseRNASeq_Zeisel 2015.txt')

# human rna seq data download from linnarson 03-----
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



# get human microarray data for healthy brain regions -----------
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

# human microarray data for huntington disease --------------------
download.file('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE3nnn/GSE3790/soft/GSE3790_family.soft.gz',
              destfile='data/GSE3790_family.soft.gz')
system('gunzip data/GSE3790_family.soft.gz')

# deal with softfile
softData  = softParser('data/GSE3790_family.soft')
softData = softData[,c('!Sample_characteristics_ch1',
                       '!Sample_geo_accession',
                       '!Sample_description',
                       '!Sample_platform_id',
                       '!Sample_source_name_ch1',
                       '!Sample_title')]
colnames(softData) = c('description',
                       'GSM',
                       'source', 
                       'platform', 
                       'source2', 
                       'name')


sapply(unique(softData$platform), function(x){
    dir.create(paste0('data//cel/',x),showWarnings=F)
})
apply(softData,1,function(x){
    xNow<<-x[colnames(softData)=='GSM']
    gsmDown(gsm=x['GSM'],outfile=paste0('data//cel/',
                                        x['platform'],'/',
                                        x['GSM'],
                                        '.cel'),
            unzip=F)
})


softData$scanDate = apply(softData,1, function(x){
    celfileDate(paste0('/home/omancarci/masterOfCellTypes/cel/',x['platform'],'/',x['GSM'], '.cel.gz'))
})


unique(softData$source2)
dir.create('data/huntington')
for (x in unique(softData$source2)){
    readA = softData %>% filter(source2 == x & platform=='GPL96') %>% select(GSM) %>% unlist
    readB = softData %>% filter(source2 == x & platform=='GPL97') %>% select(GSM) %>% unlist
    affyA = ReadAffy(filenames=paste0('data//cel//GPL96/',readA,'.cel.gz' ))
    affyB = ReadAffy(filenames=paste0('data//cel//GPL97/',readB,'.cel.gz' ))
    affyA = rma(affyA)
    affyB = rma(affyB)
    affyA = gemmaAnnot(affyA, 'data/GemmaAnnots/GPL96')
    affyB = gemmaAnnot(affyB, 'data/GemmaAnnots/GPL97')
    aName = softData %>% filter(source2 == x & platform=='GPL96') %>% select(name) %>% unlist %>% str_extract('.*?(?= [AB]$)')
    bName = softData %>% filter(source2 == x & platform=='GPL97') %>% select(name) %>% unlist %>% str_extract('.*?(?= [AB]$)')
    list[aGene,aExp] = affyA %>% sepExpr
    list[bGene,bExp] = affyB %>% sepExpr
    names(aExp) = aName
    names(bExp) = bName
    bExp = bExp[, match(aName,bName)]
    allGenes = rbind(aGene,bGene)
    allExp = rbind(aExp,bExp)
    write.csv(cbind(allGenes,allExp), paste0('data/huntington/',x,'_exp.csv'), row.names = F)
    softData %>%
        filter(source2 == x & platform=='GPL96') %>%
        select(-GSM) %>%
        mutate(huntington = !grepl('control',description)) %>% 
        write.design(file= paste0('data/huntington/',x,'_des.tsv'))
}





# human induced cell line data ------------------------
softDown('GSE40593','data/GSE40593_family.soft.gz')
system('gunzip data/GSE40593_family.soft.gz')
# files are not uploaded with this one
list[softData,expression] = softParser('data/GSE40593_family.soft',expression=T)
softData = softData[,c('!Sample_characteristics_ch1 = cell type',
                       '!Sample_geo_accession',
                       '!Sample_platform_id',
                       '!Sample_title')]
colnames(softData) = c('cellType',
                       'GSM',
                       'platform',
                       'title')


write.design(softData,file='data/GSE40593_meta')

expTable = sapply(expression, function(x){
    x$VALUE
})

genes = gemmaGeneMatch(expression[[1]]$ID_REF,chipFile='data/GemmaAnnots/GPL6884')

rownames(expTable) = genes
expTable = expTable[!rownames(expTable)=='',]
expTable = cbind(rownames(expTable), as.data.frame(expTable))
colnames(expTable)[1] = 'Gene.Symbol'
expTable = mostVariable(expTable)
write.table(expTable,'data/GSE40593_expression',sep=',',col.names=T,row.names=F)


# mouse isolated cell line data -------
softDown('GSE9566','data/GSE9566_family.soft.gz')
system('gunzip data/GSE9566_family.soft.gz')

list[softData,expression] = softParser('data/GSE9566_family.soft',expression=T)
softData = softData[,c('!Sample_geo_accession',
                       '!Sample_description',
                       '!Sample_platform_id',
                       '!Sample_title')]

colnames(softData) = c('GSM',
                       'description',
                       'platform',
                       'title')

expression = expression[softData$platform %in% 'GPL1261']
softData = softData %>% filter(platform == 'GPL1261')

write.design(softData,file='data/GSE9566_meta')

expTable = sapply(expression, function(x){
    x$VALUE
})

genes = gemmaGeneMatch(expression[[1]]$ID_REF,chipFile='data/GemmaAnnots/GPL1261')

expTable = data.frame(Gene.Symbol = genes,expTable)

expTable = expTable %>% mostVariable(treshold=0)
write.table(expTable,'data/GSE9566_expression',sep=',',col.names=T,row.names=F)


# blood cell type data (06)--------------
# download and pre-process the files, google the missing ones.
bloodDir = 'data/bloodCellType/'
dir.create(bloodDir, showWarnings=FALSE)
download.file('http://www.nature.com/nmeth/journal/v12/n5/extref/nmeth.3337-S2.xls',
              destfile=paste0(bloodDir,'supp2.xls'))
bloodDes = loadWorkbook(paste0(bloodDir,'supp2.xls'))
bloodDes = readWorksheet(bloodDes, sheet = 2, header = TRUE)

# special case because reasons....
bloodDes$Sample.ID[bloodDes$Sample.ID %in% 'A_MF_2hrEosinophils_U133A'] = 'A_MF_2hrEosinophils'

bloodDes$Leuk11 = c(rep('B Cells',15),rep('PCs', 7), rep('CD8 T cells', 4), rep('CD4', 14), rep('GammaDeltaT', 2),
                    rep('NK', 15), rep('MonoMacro', 30), rep('Dendritic', 12), rep('Mast', 4), rep('Eos', 2), rep('Neutrophils',8))

bloodDes$originalIndex = as.numeric(factor(bloodDes$Abreviated.name))
names(bloodDes)[names(bloodDes) %in% 'Sample.ID'] ='sampleName'
write.design(bloodDes,paste0(bloodDir,'BloodCells.tsv'))

gsms = bloodDes$Sample.ID[grepl('GSM',bloodDes$Sample.ID)]

sapply(gsms,function(gsm){
     gsmDown(gsm, paste0('data/cel/GPL96/',gsm,'.CEL'))
})
# other cel file links. bah complete this later. just google it damn it
links = c('http://linkage.garvan.unsw.edu.au/public/microarrays/Arthritis_Inflammation/human/eosinophils/pma/A_MF_2hrEosinophils.CEL',
          'http://linkage.garvan.unsw.edu.au/public/microarrays/Arthritis_Inflammation/human/neutrophils/LPS/A_TS_MSNeutroLPS_U133A.CEL'
          )


cels = paste0('data/cel/GPL96/',bloodDes$Sample.ID,'.CEL')
affy = ReadAffy(filenames = cels)
norm = rma(affy)
annotated = gemmaAnnot(norm, 'data/GemmaAnnots/GPL96')
names(annotated) = gsub('[.]CEL','',names(annotated))
write.csv(annotated, paste0(bloodDir,'bloodExp.csv'), row.names = F)
quantileNorm(paste0(bloodDir,'bloodExp.csv'),
             paste0(bloodDir,'qnormBlood.csv'))

mostVariableCT(paste0(bloodDir,'bloodExp.csv'),
               paste0(bloodDir,'finalBlood.csv'),
               selectionNaming = 'Abreviated.name',
               design=paste0(bloodDir,'BloodCells.tsv'))

# lymphoma stuff
gseDown(GSE='GSE65135',regex='lymph',outDir='data/cel/GPL570')
cels =celFiles('data/cel/GPL570',full.names=T) 
cels = cels[grep(regexMerge(gsmFind('GSE65135', 'lymph')), cels)]
affy = ReadAffy(filenames = cels)
norm = rma(affy)
annotated = gemmaAnnot(norm, 'data/GemmaAnnots/GPL570')
names(annotated) = gsub('[.]cel','',names(annotated))
annotated = mostVariable(allDataPre=annotated)
write.csv(annotated, paste0(bloodDir,'lymphoma.csv'), row.names = FALSE)

# pbmc stuff
PBMCs = read.table(textConnection(getURL('https://cibersort.stanford.edu/inc/inc.download.page.handler.php?file=PBMCs-Fig3a-HumanHT-12-V4.txt')),sep='\t',header=T)
write.table(PBMCs, paste0(bloodDir,'PBMCs.csv'), sep=',', quote=FALSE, row.names=F)


PBMCcounts = read.table(textConnection(getURL('https://cibersort.stanford.edu/inc/inc.download.page.handler.php?file=PBMCs-Fig3a-Flow-Cytometry.txt')),sep='\t',header=T)
write.table(PBMCcounts, paste0(bloodDir,'PBMCcounts.tsv'), sep='\t', quote=FALSE, row.names=FALSE)


download.file('http://www.nature.com/nmeth/journal/v12/n5/source_data/nmeth.3337-f3.xlsx',
              destfile=paste0(bloodDir,'fig3.xlsx'))

PBMCcibersort = loadWorkbook(paste0(bloodDir,'fig3.xlsx'))
PBMCcibersort = readWorksheet(PBMCcibersort, sheet = 1, header = TRUE)

write.table(PBMCcibersort, paste0(bloodDir,'PBMCcibersort.tsv'),sep = '\t', quote=FALSE, row.names=FALSE)

# parkinsons substantia nigra
gseDown('GSE7621',outDir='data/cel/GPL570/')
softDown('GSE7621','data/GSE7621_family.soft.gz')
system('gunzip data/GSE7621_family.soft.gz')


write.design(softData,'data/GSE7621_parkinsonsMeta.tsv')


affy = ReadAffy(filenames = cels)

norm = rma(affy)
annotated = gemmaAnnot(norm, 'data/GemmaAnnots/GPL570')
annotated = mostVariable(annotated)
write.csv(annotated, paste0('data/','GSE7621_parkinsonsExp.csv'), row.names = F)


# developmental data
softDown('GSE25219','data/development/GSE25219_family.soft.gz')
system('gunzip data/development/GSE25219_family.soft.gz')
softParser(softFile='data/development/GSE25219_family.soft',expression=F)


# allen brain human data
dir.create('data/allenBrain')
download.file('http://human.brain-map.org/api/v2/data/MicroarraySlide/query.csv?num_rows=all&criteria=well_known_files,structures,microarray_data_sets(specimen(donor),products[id$eq2])&tabular=microarray_slides.barcode,microarray_slides.id%20as%20array_id,microarray_slides.rna_in_array,microarray_slides.rna_in_labeling_reaction,microarray_slides.sample_id%20as%20tissue_sample_id,structures.name%20as%20structure,specimens.rna_integrity_number,donors.name%20as%20donor,well_known_files.download_link',
              'data/allenBrain/metadata.csv')
download.file('http://human.brain-map.org/api/v2/well_known_file_download/178238387',
              'data/allenBrain/H0351.2001.zip')
download.file('http://human.brain-map.org/api/v2/well_known_file_download/178238373',
              'data/allenBrain/H0351.2002.zip')
download.file('http://human.brain-map.org/api/v2/well_known_file_download/178238359',
              'data/allenBrain/H0351.1009.zip')
download.file('http://human.brain-map.org/api/v2/well_known_file_download/178238316',
              'data/allenBrain/H0351.1012.zip')
download.file('http://human.brain-map.org/api/v2/well_known_file_download/178238266',
              'data/allenBrain/H0351.1015.zip')
download.file('http://human.brain-map.org/api/v2/well_known_file_download/178236545',
              'data/allenBrain/H0351.1016.zip')

dir.create('data/allenBrain/H0351.2001')
system('unzip -n data/allenBrain/H0351.2001.zip -d data/allenBrain/H0351.2001')

dir.create('data/allenBrain/H0351.2002')
system('unzip -n data/allenBrain/H0351.2002.zip -d data/allenBrain/H0351.2002')

dir.create('data/allenBrain/H0351.1009')
system('unzip -n data/allenBrain/H0351.1009.zip -d data/allenBrain/H0351.1009')

dir.create('data/allenBrain/H0351.1012')
system('unzip -n data/allenBrain/H0351.1012.zip -d data/allenBrain/H0351.1012')

dir.create('data/allenBrain/H0351.1015')
system('unzip -n data/allenBrain/H0351.1015.zip -d data/allenBrain/H0351.1015')

dir.create('data/allenBrain/H0351.1016')
system('unzip -n data/allenBrain/H0351.1016.zip -d data/allenBrain/H0351.1016')


# allen institute single cell data 
download.file('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE71nnn/GSE71585/suppl/GSE71585_RefSeq_RPKM.csv.gz',
              'data/allenBrainSingleCells/GSE71585_RefSeq_RPKM.csv.gz')
system('gunzip data/allenBrainSingleCells/GSE71585_RefSeq_RPKM.csv.gz')

download.file('http://www.nature.com/neuro/journal/vaop/ncurrent/extref/nn.4216-S8.xlsx',
              'data/allenBrainSingleCells/markerGenes.xlsx')