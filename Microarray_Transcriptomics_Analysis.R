## File name: Microarray Transcriptomics Analysis in R
## Author: Hasani Perera
## Contact: heperera826@gmail.com
## Date created: 20/10/2021
## Date last modified: 25/11/2021
## R Version: 4.0.3
## Ref: Taneera et al., 2015 (doi:10.1093/hmg/ddu610)

## set path
setwd("~/Documents/HazGit/Transcriptomics_Analysis")

#install packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery")
nBiocManager::install("oligo")
BiocManager::install("pd.hugene.1.0.st.v1")
BiocManager::install("hugene10sttranscriptcluster.db")
BiocManager::install("arrayQualityMetrics")
BiocManager::install("GOstats")
install.packages("ggplot2")


##------------------------ QUERY RAW DATA ------------------------------

#load packages
library(GEOquery)

## GEO Accession No: GSE50397
## Platform: [HuGene-1_0-st] Affymetrix Human Gene 1.0 ST Array  - GPL6244

#extract .cell files to local machine - GEO Series records (GSExxxxx)
gse <- getGEO("GSE50397",GSEMatrix=FALSE)

#if file is already downloaded
# gse <- file.choose()
head(Meta(gse))

# names of all the GSM objects contained in the GSE
names(GSMList(gse))

#first GSM object on the list
class(GSMList(gse)[[1]])
head(Meta(GSMList(gse)[[1]]))

# names of the GPLs represented
names(GPLList(gse))

# access raw data (downloaded file paths)
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

##------------------------ 1) DATA PREPROCESSING ------------------------------

#load packages
library(oligo)

#list of all cel files in the directory
celpath <- "~/Documents/Learning/Bioinformatics/MicroArray Analysis/Nordic Islet Analysis/Raw_Data"
celfiles_list <- list.files(celpath,pattern = ".CEL", full.names=TRUE)
length(celfiles_list)
head(celfiles_list)

# import CEL files containing raw probe-level data into an R AffyBatch object
cell_files <- read.celfiles(celfiles_list)
head(cell_files)

#load packages
library(pd.hugene.1.0.st.v1)

#probe set annotation
??pd.hugene.1.0.st.v1

getClass("GeneFeatureSet")

# max expression
max(exprs(cell_files[1:1102489,1:89]))

# replace sampleNames
filename <- sampleNames(cell_files)
pData(cell_files)$filename <- filename
pData(cell_files)$filename
sampleNames <- sub("-islet.CEL$","",filename)
sampleNames <- sub("_HTL[[:digit:]]*","",sampleNames)
sampleNames(cell_files) <- sampleNames
sampleNames(cell_files)

# information on variable values/ meta-data
pData(cell_files)

#boxplot before RMA normalization
boxplot(cell_files, target="probeset")
mtext(text="log2 Intensity", side=2, line=3, las=0) 
mtext(text="Samples", side=1, line=3, las=1) 

#perform RMA normalization (Robust Multi-Array Average)
# converts an AffyBatch object into an ExpressionSet object
normData <- rma(cell_files)
nrow(normData)

#boxplot after RMA normalization
boxplot(normData, target="probeset")
mtext(text="log2 Intensity", side=2, line=3, las=0) 
mtext(text="Samples", side=1, line=3, las=1) 

#save the expression data (output - normalized and log2 transformed)
exprs(normData)[1:3,1:5]

write.exprs(normData,file="RMA_Normalised_Original.txt")

#load packages
library(Biobase)
library(hugene10sttranscriptcluster.db)

#gene annotation
??Biobase
??hugene10sttranscriptcluster.db

#get a list of retrievable data
keytypes(hugene10sttranscriptcluster.db)

#retrieve data for selected objects (ENTREZID and SYMBOL) as a data frame
anno<- select(hugene10sttranscriptcluster.db,keys(hugene10sttranscriptcluster.db),
             c("ENTREZID", "SYMBOL"))
head(anno)
tail(anno)

## optional - to keep one match per gene
anno <- anno[!duplicated(anno[,1]),]
tail(anno)

#set row names to ProbeID (for convenience)
anno = anno[,-1]
row.names(anno) = keys(hugene10sttranscriptcluster.db)
tail(anno)

#retrieve gene expression matrix from normData as a dataframe
expr <- data.frame(exprs(normData))
head(expr)

#merge gene expression and annotation according to row names (probe IDs)
expr_anno <- merge(x=anno,y=expr,by.y=0, by.x=2,all=TRUE)
head(expr_anno)

#save the annotated gene expression matrix to local file
write.table(expr_anno, file = "RMA_Norm_Expr.txt",sep = "\t",
            row.names = FALSE, col.names = TRUE, quote =FALSE)
write.csv(expr_anno, file = "RMA_Norm_Expr_Anno.csv")
