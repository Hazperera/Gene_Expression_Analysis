# Transcriptomics Analysis 
### How to perform a standard analysis of microarray transcriptomics data using R

#### Background:
- Use of microarray data for gene expression analysis is still a common practice even though the more advanced RNA-Seq data is gaining more popularity.
- The steps of microarray analysis include, quality control, normalization, differential expression analysis and downstream analysis.
- R-[Bioconductor](https://www.bioconductor.org/) is an open-source software with numerous tools for the analysis of high-throughput genomic data.

#### Objective:
This project aims to perform a microarray gene expression analysis using human pancreatic islet data generated from Affymetrix Human Gene 1.0 ST Array [transcript (gene) version] platform by following a stepwise protocol adopted from [Ming-an Sun et al.(2018)](https://pubmed.ncbi.nlm.nih.gov/29508287/)

#### Dataset Description: The expression arrays in 89 human pancreatic islet donors with different levels of blood glucose (HbA1c).
NCBI GEO Accession No: [GSE50397](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50397) 

#### References: *[Taneera, J et al. (2015). Identification of novel genes for glucose metabolism based upon expression pattern in human islets and effect on insulin secretion and glycemia. Human molecular genetics, 24(7), 1945–1955](https://doi.org/10.1093/hmg/ddu610)*
