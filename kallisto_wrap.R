library(DESeq2)
library(tximport)
library(openxlsx)

system ("mkdir ./output")
system ("cp ./projects/log.out ./output")


tx2gene <- read.delim ("gencode.vM32.annotation.txt")
tx2gene <- tx2gene[ ,c("transcript_id", "gene_id")]
colnames (tx2gene) <- c("TXNAME", "GENEID")

dir <- paste (paste (getwd (), "projects", sep="/"), "kallisto_results", sep="/")
files <- list.files (path=dir, pattern=".*abundance.tsv", recursive=TRUE)
files <- paste (dir, files, sep="/")
names (files) <- gsub (".*IIT_", "", gsub ("/abundance.tsv", "", files))

txi <- tximport(files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE, countsFromAbundance="lengthScaledTPM")
res <- txi$counts

anno <- read.delim ("gencode.vM32.annotation.txt")
anno <- anno[ ,grep ("transcript_id", colnames (anno), invert=TRUE)]
anno <- unique (anno)

res <- merge (res, anno, by.x="row.names", by.y="gene_id", all.x=TRUE)
colnames (res)[1] <- "Geneid"
res <- res[grep ("miRNA|Mt_tRNA|Mt_rRNA|rRNA|snRNA|snoRNA|scRNA|sRNA|misc_RNA|scaRNA|ribozyme|IG_|TR_", res$gene_type, invert=TRUE), ]


write.xlsx (res, "./output/kallisto_gene_lengthScaledTPM_counts.xlsx", rowNames=F)
#write.table (res, "./output/kallisto_gene_lengthScaledTPM_counts.txt", sep="\t", quote=F, row.names=F)
