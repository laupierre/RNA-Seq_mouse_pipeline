library(DESeq2)
library(tximport)
library(openxlsx)


tx2gene <- read.delim ("gencode.vM32.annotation.txt")
tx2gene <- tx2gene[ ,c("transcript_id", "gene_id")]
colnames (tx2gene) <- c("TXNAME", "GENEID")

dir <- paste (paste (getwd (), "projects", sep="/"), "salmon_results", sep="/")
files <- list.files (path=dir, pattern=".*quant.sf", recursive=TRUE)
files <- paste (dir, files, sep="/")
names(files) <- gsub (".*IIT_", "", gsub ("/quant.sf", "", files))


## txi$counts values are dependent on the countsFromAbundance option !
## the analysis is different when plugged directly into DESeq2, see: https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html#limma-voom
## Save a generic lengthScaledTPM counts re-usable later for voom-limma and DESeq2 (without direct plug to DESeq2)

txi <- tximport(files, type = "salmon", tx2gene = tx2gene, countsFromAbundance="lengthScaledTPM")
res <- txi$counts

anno <- read.delim ("gencode.vM32.annotation.txt")
anno <- anno[ ,grep ("transcript_id", colnames (anno), invert=TRUE)]
anno <- unique (anno)

res <- merge (res, anno, by.x="row.names", by.y="gene_id", all.x=TRUE)
colnames (res)[1] <- "Geneid"

write.xlsx (res, "./output/salmon_gene_lengthScaledTPM_counts.xlsx", rowNames=F)
#write.table (res, "./output/salmon_gene_lengthScaledTPM_counts.txt", sep="\t", quote=F, row.names=F)
