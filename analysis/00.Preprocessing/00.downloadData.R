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

# get human data
download.file('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE46nnn/GSE46706/soft/GSE46706_family.soft.gz',destfile='data/GSE46706_family.soft.gz')
system('gunzip data/GSE46706_family.soft.gz')