library(DESeq2)
library(tximport)
library(openxlsx)
library (ggplot2)


system ("mkdir ./output")
system ("cp ./projects/log.out ./output")


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
res <- res[grep ("miRNA|Mt_tRNA|Mt_rRNA|rRNA|snRNA|snoRNA|scRNA|sRNA|misc_RNA|scaRNA|ribozyme|IG_|TR_", res$gene_type, invert=TRUE), ]


write.xlsx (res, "./output/salmon_gene_lengthScaledTPM_counts.xlsx", rowNames=F)
#write.table (res, "./output/salmon_gene_lengthScaledTPM_counts.txt", sep="\t", quote=F, row.names=F)


##########################
## Differential expression
a <- annot <- read.xlsx ("./output/salmon_gene_lengthScaledTPM_counts.xlsx")

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


dds <- DESeqDataSetFromMatrix(countData = round (a), colData = sampleTable, design = ~ condition)
                                 
# keep <- rowSums(counts(dds)) >= 10
keep <- rowSums(counts(dds) >= 10) >= dim (a)[2]/2
dds <- dds[keep,]

# R will choose a reference level for factors based on alphabetical order
dds <- DESeq(dds)
res <- results(dds)

## MA plot
pdf ("./output/MA_plot.pdf")
plotMA(res, ylim=c(-5,5))
dev.off()

## PCA plot
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("condition", "replicate"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=condition, shape=replicate)) +
  		geom_point(size=3) +
  		xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  		ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
		  coord_fixed ()

ggsave ("./output/PCA_plot.pdf")

## save output
res <- merge (data.frame (res), counts (dds), by="row.names")
res <- merge (res, annot, by.x="Row.names", by.y="Geneid")
colnames (res)[1] <- "Geneid"
res <- res[order (res$padj), ]
write.xlsx (res, "./output/salmon_deseq2_differential_expression.xlsx")

   

