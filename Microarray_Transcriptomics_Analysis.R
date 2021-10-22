## Ref: Taneera et al., 2015 (doi:10.1093/hmg/ddu610)
## GEO Accession No: GSE50397 
## Platform: [HuGene-1_0-st] Affymetrix Human Gene 1.0 ST Array  

## set path
setwd("~/Documents/Learning/Bioinformatics/MicroArray/Nordic")

#install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery")
BiocManager::install("oligo")
BiocManager::install("pd.hugene.1.0.st.v1")
BiocManager::install("hugene10sttranscriptcluster.db")
BiocManager::install("arrayQualityMetrics")
BiocManager::install("GOstats")
install.packages(geoq)
install.packages("ggplot2")


##------------------------ QUERY RAW DATA ------------------------------

#load packages 
library(GEOquery)

## GSE50397 - 
## GPL6244  [HuGene-1_0-st] Affymetrix Human Gene 1.0 ST Array [transcript (gene) version])   


#extract .cell files to local machine - GEO Series records (GSExxxxx)
gse <- getGEO("GSE50397",GSEMatrix=TRUE)
head(Meta(gse))

# names of all the GSM objects contained in the GSE
names(GSMList(gse))

#first GSM object on the list
class(GSMList(gse)[[1]])
head(Meta(GSMList(gse)[[1]]))

# names of the GPLs represented
names(GPLList(gse))

#access raw data (downloaded file paths)
file_paths = getGEOSuppFiles("GSE50397")
file_paths

#acess GSE Data Tables from GEO
# df1 <- getGSEDataTables("GSE50397")
# lapply(df1, head)



##------------------------ 1) DATA PREPROCESSING ------------------------------

#load packages 
library(oligo)

#list of all CEL files in the directory
cel.files <- list.celfiles()
cel.files

#re-specify sample names


