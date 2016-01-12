# brainCellTypeSpecificGenes

This project aims to identify cell type specific genes of the brain, verify them in independent datasets, and use them to estimate cell type proportions. This repository currently hosts the [raw and pre-processed data](data) and [genes selected as cell type markers](/analysis/01.Gene Selection/FinalGenes/PyramidalDeep/) in different brain regions. The rest of the project is currently being migrated from the extremely hard to navigate [oganm/masterOfCellTypes](https://github.com/oganm/masterOfCellTypes) repository. The bulk of the data was compiled by Okaty et al. in an [independent study](http://dx.doi.org/10.1371/journal.pone.0016493)

### Data
Data folder used to hold raw and pre-processed data. Cell type expression profiles were removed since some of them are acquired through personal communications. We are working to make them or a limited dataset available again.

Outputs of other analysis steps can be found in their respective directories

### Dependencies
* ogbox (see below)
* geneSynonym (see below)
* homologene (see below)
* stringr
* affy
* oligo
* RCurl
* compare
* preprocessCore
* foreach
* doMC
* parallel
* cluster
* reshape2
* knitr
* ggplot2
* corpcor
* data.table
* igraph
* RBGL
* viridis
* sva
* gplots
* mixtools

ogbox is a general purpose package that is written for personal use. Hence it is not on CRAN. To get it, you need to install devtools if you don't have it already

geneSynonym is a package to find gene synonyms and homologene finds gene homologues across species. They also need to be installed from github.

```
install.packages('devtools')
library(devtools)
install_github('oganm/ogbox')
install_github('oganm/geneSynonym')
install_github('oganm/homologene')


```
