library (openxlsx)
library (limma)
library (edgeR)

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



y <- DGEList(a)

# filtering using the design information:
design <- model.matrix(~0+condition, data = sampleTable)
colnames(design) <- gsub("condition", "", colnames(design))
keep <- filterByExpr(y, design)
y <- y[keep, ]

contr.matrix <- makeContrasts(contrast1= treated - control, levels = colnames(design))
  
# normalize and run voom transformation
y <- calcNormFactors(y)
v <- voom(y, design)

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix) 
efit <- eBayes(vfit)
res <- topTable(efit, coef=1, n=Inf)




