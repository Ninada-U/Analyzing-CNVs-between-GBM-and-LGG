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