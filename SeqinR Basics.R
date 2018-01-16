# http://a-little-book-of-r-for-bioinformatics.readthedocs.io/en/latest/src/chapter3.html

library(seqinr)
library(devtools)
choosebank(infobank=T)

choosebank('genbank')
Python <- query("Python.genome", "SP=Python")

closebank()



getName(Dengue1)

install.packages('Bioconductor')
source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("Biostrings")
install_github("mhahsler/rBLAST")
? blast
??blast

library(rBLAST)
