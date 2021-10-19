## Ref: Taneera et al., 2015 (doi:10.1093/hmg/ddu610)
## GEO Accession No: GSE50397 
## Platform: [HuGene-1_0-st] Affymetrix Human Gene 1.0 ST Array  

## set path
setwd("~/Documents/Learning/Bioinformatics/MicroArray/Nordic")

#install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery")

## Step1:Query Raw Data

#load packages 
library(GEOquery)

#extract .cell files to local machine - GEO Series records (GSExxxxx)
gse <- getGEO("GSE50397",GSEMatrix=FALSE)
head(Meta(gse))

# names of all the GSM objects contained in the GSE
names(GSMList(gse))

#first GSM object on the list
class(GSMList(gse)[[1]])

head(Meta(GSMList(gse)[[1]]))

# names of the GPLs represented
names(GPLList(gse))
