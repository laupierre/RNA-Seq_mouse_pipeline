library (openxlsx)
library (limma)
library (edgeR)


anno <- read.delim ("gencode.vM32.annotation.txt")
anno <- anno[ ,grep ("transcript_id", colnames (anno), invert=TRUE)]
anno <- unique (anno)


a <- annot <- read.xlsx ("./output/kallisto_gene_lengthScaledTPM_counts.xlsx")

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



## without contrast
y <- DGEList(a)

# filtering using the design information:
design <- model.matrix(~condition, data = sampleTable)
colnames(design) <- gsub("condition", "", colnames(design))
keep <- filterByExpr(y, design)
y <- y[keep, ]


y <- calcNormFactors(y)
v <- voom(y, design)

vfit <- lmFit(v, design)
efit <- eBayes(vfit)
res <- topTable(efit, n=Inf)

res <- merge (res, v$E, by="row.names")
res <- merge (res, annot, by.x="Row.names", by.y="Geneid")
colnames (res)[1] <- "Geneid"
res <- res[order(res$adj.P.Val), ]
write.xlsx (res, "./output/kallisto_limma_differential_expression.xlsx", rowNames=F)





