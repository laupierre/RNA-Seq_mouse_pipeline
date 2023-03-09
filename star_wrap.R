library (openxlsx)
library (DESeq2)
library (ggplot2)

system ("mkdir ./output")
system ("cp ./projects/log.out ./output")

anno <- read.delim ("gencode.vM32.annotation.txt")
anno <- anno[ ,grep ("transcript_id", colnames (anno), invert=TRUE)]
anno <- unique (anno)

a <- read.delim ("./projects/star_results/subread.counts.txt", skip=1)
a <- a[ ,grep ("ene|bam", colnames (a))]
a <- a[grep ("miRNA|Mt_tRNA|Mt_rRNA|rRNA|snRNA|snoRNA|scRNA|sRNA|misc_RNA|scaRNA|ribozyme|IG_|TR_", a$gene_type, invert=TRUE), ]
colnames (a) <- gsub ("_S[0-9]+.*", "", colnames (a)) 
colnames (a) <- gsub (".*IIT_", "", colnames (a))

a <- merge (a, anno, by.x="Geneid", by.y="gene_id", all.x=TRUE) 
a <- a[ ,grep ("gene_type.y|gene_name.y", colnames (a), invert=TRUE)]
colnames (a) [colnames (a) == "gene_name.x"] <- "gene_name"
colnames (a) [colnames (a) == "gene_type.x"] <- "gene_type"
  
write.xlsx (a, "./output/star_gene_raw_counts.xlsx", rowNames=F)
#write.table (a, "./output/star_gene_raw_counts.txt", sep="\t", quote=F, row.names=F)


##########################
## Differential expression
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


dds <- DESeqDataSetFromMatrix(countData = a, colData = sampleTable, design = ~ condition)
                                 
# keep <- rowSums(counts(dds)) >= 10
keep <- rowSums(counts(dds) >= 10) >= dim (a)[2]/2
dds <- dds[keep,]

# R will choose a reference level for factors based on alphabetical order
dds <- DESeq(dds)
res <- results(dds)
res <- merge (data.frame (res), counts (dds), by="row.names")
res <- merge (res, annot, by.x="Row.names", by.y="Geneid")
colnames (res)[1] <- "Geneid"
res <- res[order (res$padj), ]
write.xlsx (res, "./output/star_deseq2_differential_expression.xlsx")


## MA plot
pdf ("MA_plot.pdf")
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
ggsave ("PCA_plot.pdf")
