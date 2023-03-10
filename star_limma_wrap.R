library (openxlsx)
library (DESeq2)
library (ggplot2)
library (limma)


a <- annot <- read.xlsx ("./output/star_gene_raw_counts.xlsx")

annot <- annot[ ,c("Geneid", "gene_name", "gene_type", "mgi_id", "external_gene_name", "description")]

torm <- c("gene_name", "gene_type", "mgi_id", "external_gene_name", "description")
a <- a[ ,!colnames (a) %in% torm]
row.names (a) <- a[ ,1]
a <- a[ ,-1]
colnames (a) <- gsub("(_S[0-9].*)" , "", colnames (a))


# samples sheet
sampleTable <- read.xlsx ("samples_rnaseq.xlsx")
sampleTable$sample <- gsub("IIT_", "", sampleTable$sample)
sampleTable$sample <- gsub("(_S[0-9].*)" , "", sampleTable$sample)

idx <- match (sampleTable$sample, colnames (a))
sampleTable <- sampleTable[idx, ]
sampleTable$condition <- factor (sampleTable$condition)
sampleTable$replicate <- factor (sampleTable$replicate)

stopifnot (sampleTable$sample == colnames (a))



