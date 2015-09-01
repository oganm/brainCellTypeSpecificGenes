library(ogbox)
sourceGithub(OganM,toSource,gemmaAnnotate)
sourceGithub(OganM,toSource,homologene)
sourceGithub(OganM,toSource,GEO)
# download gemma annotations ------
dir.create('data/GemmaAnnots/', showWarnings=F)
getGemmaAnnotGoogle('GPL1261','data/GemmaAnnots/GPL1261')
getGemmaAnnotGoogle('GPL339','data/GemmaAnnots/GPL339')
getGemmaAnnotGoogle('GPL570','data/GemmaAnnots/GPL570')
getGemmaAnnotGoogle('GPL96','data/GemmaAnnots/GPL96')
getGemmaAnnotGoogle('GPL5175','data/GemmaAnnots/GPL5175')


# download homology information ----------
# used internally in homogene functions. I know its bad design but...
homoloGeneTarget = 'data/homogene.tsv'
setUpHomologene()

# get human data ------------------------
download.file('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46706/soft/GSE46706_family.soft.gz',
              destfile='data/GSE46706_family.soft.gz')
system('gunzip data/GSE46706_family.soft.gz')


softData = softParse('data/GSE46706_family.soft')

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
                       '!Sample_title')]

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
                       'title')

softData$GSM = rownames(softData)
write.design(softData,file='data/GSE46706_meta')

