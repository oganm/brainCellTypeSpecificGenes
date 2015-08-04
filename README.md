# brainCellTypeSpecificGenes
code and data for the manuscript.

Data folder holds raw and pre-processed data. Outputs of other analysis steps can be found in their respective directories

Dependencies: ogbox, stringr, affy, compare, preprocessCore, foreach, doMC, parallel, cluster, reshape2

ogbox is a general purpose package not available in CRAN. To get it, you need to install devtools if you don't have it already

```
install.packages('devtools')
library(devtools)
install_github('oganm/ogbox')
```
