
## File name: Microarray_Transcriptomics_Analysis.R
## Author: Hasani Perera
## Contact: heperera826@gmail.com
## Date created: 25/10/2021
## Date last modified: 27/10/2021
## R Version: 4.0.3 

## set path
setwd("~/Documents/HazGit/Transcriptomics_Analysis")

#install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery")
BiocManager::install("oligo")
BiocManager::install("pd.hugene.1.0.st.v1")
BiocManager::install("hugene10sttranscriptcluster.db")
BiocManager::install("arrayQualityMetrics")
BiocManager::install("GOstats")
install.packages("ggplot2")


##------------------------ QUERY RAW DATA ------------------------------

#load packages 
library(GEOquery)

## GSE50397 - Taneera et al., 2015 (doi:10.1093/hmg/ddu610)
## GPL6244  [HuGene-1_0-st] Affymetrix Human Gene 1.0 ST Array [transcript (gene) version])   


#extract .cell files to local machine - GEO Series records (GSExxxxx)
gse <- getGEO("GSE50397",GSEMatrix=FALSE)
head(Meta(gse))
# show(gse)

# names of all the GSM objects contained in the GSE
names(GSMList(gse))

#first GSM object on the list
class(GSMList(gse)[[1]])
head(Meta(GSMList(gse)[[1]]))

# names of the GPLs represented
names(GPLList(gse))

#access raw data (downloaded file paths)
file_paths = getGEOSuppFiles("GSE50397")
head(file_paths)

# choose tar file
tarfile <- file.choose()

# extract tar archives
untar(tarfile, exdir="Raw_Data")

#list of all gz files in the directory
cel.files <- list.files("Raw_Data/", pattern = "[gz]")
length(cel.files)

# list/vector/dataframe ---> vector/matrix
#  extract gz archives - gunzip
sapply(paste("Raw_Data", cel.files, sep="/"), gunzip)

#list of all cel files in the directory
cel.files <- list.files("Raw_Data/", pattern = ".CEL")
head(cel.files)
length(cel.files)
cel.files

##------------------------ 1) DATA PREPROCESSING ------------------------------

#load packages 
library(oligo)
library(affy)

#re-specify sample names
sample.names = c(1:89)
sample.names

#read files to memory
affy.raw <- read.celfiles(cel.files,sampleNames=sample.names)
head(affy.raw)

#load packages 
library(pd.hugene.1.0.st.v1)

#probe set annotation 
??pd.hugene.1.0.st.v1

#perform RMA normalization (Robust Multi-Array Average)
eset <- rma(affy.raw)
nrow(eset)

#save the expression data (output - normalized and log2 transformed)
write.exprs(eset,file="rma_norm_expr.txt")

#load packages 
library(Biobase)
library(hugene10sttranscriptcluster.db)

#gene annotation 
??Biobase
??hugene10sttranscriptcluster.db

#get a list of retrievable data 
keytypes(hugene10sttranscriptcluster.db)

#retrieve data for selected objects (ENTREZID and SYMBOL) as a data frame
gns<- select(hugene10sttranscriptcluster.db,keys(hugene10sttranscriptcluster.db),
             c("ENTREZID", "SYMBOL"))
head(gns)






