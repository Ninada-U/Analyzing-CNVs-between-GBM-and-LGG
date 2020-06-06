#! /usr/bin/Rscript

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")

library(TCGAbiolinks)
library(SummarizedExperiment)

query <- GDCquery(project = "TCGA-GBM", 
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  experimental.strategy = "RNA-Seq",
                  workflow.type = "HTSeq - Counts"
                  )
GDCdownload(query)
BRCARnaseq_assay <- GDCprepare(query)
BRCAMatrix <- assay(BRCARnaseq_assay)
BRCARnaseq_CorOutliers <- TCGAanalyze_Preprocessing(BRCARnaseq_assay)

query <- GDCquery(project = "TCGA-LGG", 
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  experimental.strategy = "RNA-Seq",
                  workflow.type = "HTSeq - Counts"
                  )
GDCdownload(query)
BRCARnaseq_assay <- GDCprepare(query)
BRCAMatrix <- assay(BRCARnaseq_assay)
BRCARnaseq_CorOutliers <- TCGAanalyze_Preprocessing(BRCARnaseq_assay)


clin.gbm <- GDCquery_clinic("TCGA-GBM", "clinical")
TCGAanalyze_survival(clin.gbm, "gender", main = "TCGA Set\n GBM",height = 10, width=10)

clin.gbm <- GDCquery_clinic("TCGA-LGG", "clinical")
TCGAanalyze_survival(clin.gbm, "gender", main = "TCGA Set\n LGG",height = 10, width=10)
