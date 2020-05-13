# You can define a list of samples to query and download providing relative TCGA barcodes.

listSamplesLGG <- c("TCGA-HT-8558-01A-21D-2391-01","TCGA-HT-7857-10A-01D-2392-01",
                 "TCGA-DU-6406-10A-01D-1704-01","TCGA-P5-A5F0-10A-01D-A288-01",
                 "TCGA-DU-A6S3-10A-01D-A328-01","TCGA-DU-6404-10A-01D-1704-01",
                 "TCGA-HT-8010-10A-01D-2392-01","TCGA-HT-8558-10A-01D-2392-01",
                 "TCGA-TM-A84Q-10A-01D-A366-01","TCGA-S9-A7R2-10A-01D-A34L-01")

# Query platform Illumina HiSeq with a list of barcode 
query <- GDCquery(project = "TCGA-GBM", 
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  experimental.strategy = "RNA-Seq",
                  workflow.type = "HTSeq - Counts"
                  )

# Download a list of barcodes with platform IlluminaHiSeq_RNASeqV2
GDCdownload(query)

# Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
# rsem.genes.results as values
BRCARnaseq_assay <- GDCprepare(query)

BRCAMatrix <- assay(BRCARnaseq_assay)

# For gene expression if you need to see a boxplot correlation and AAIC plot to define outliers you can run
BRCARnaseq_CorOutliers <- TCGAanalyze_Preprocessing(BRCARnaseq_assay)


