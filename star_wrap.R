library (openxlsx)

system ("mkdir ./output")
system ("cp ./projects/log.out ./output")
system ("cp /projects/ncrrbt_share_la/dev_pipe/gencode.vM32.annotation.txt .")

anno <- read.delim ("gencode.vM32.annotation.txt")
anno <- anno[ ,grep ("transcript_id", colnames (anno), invert=TRUE)]
anno <- unique (anno)
head (anno)

a <- read.delim ("./projects/star_results/subread.counts.txt", skip=1)
a <- a[ ,grep ("ene|bam", colnames (a))]
a <- a[grep ("miRNA|Mt_tRNA|Mt_rRNA|rRNA|snRNA|snoRNA|scRNA|sRNA|misc_RNA|scaRNA|ribozyme|IG_|TR_", a$gene_type, invert=TRUE), ]
colnames (a) <- gsub ("_S[0-9]+.*", "", colnames (a)) 
colnames (a) <- gsub (".*IIT_", "", colnames (a))

a <- merge (a, anno, by.x="Geneid", by.y="gene_id", all.x=TRUE) 
a <- a[ ,grep ("gene_type.y|gene_name.y", colnames (a), invert=TRUE)]
colnames (a) [colnames (a) == "gene_name.x"] <- "gene_name"

head (a)
  
write.xlsx (a, "./output/star_gene_raw_counts.xlsx", rowNames=F)
#write.table (a, "./output/star_gene_raw_counts.txt", sep="\t", quote=F, row.names=F)
