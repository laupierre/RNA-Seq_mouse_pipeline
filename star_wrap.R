library (openxlsx)

system ("mkdir ./output")
system ("cp ./projects/log.out ./output")
system ("cp /projects/ncrrbt_share_la/dev_pipe/ .")

a <- read.delim ("./projects/star_results/subread.counts.txt", skip=1)
a <- a[ ,grep ("ene|bam", colnames (a))]
a <- a[grep ("miRNA|Mt_tRNA|Mt_rRNA|rRNA|snRNA|snoRNA|scRNA|sRNA|misc_RNA|scaRNA|ribozyme|IG_|TR_", a$gene_type, invert=TRUE), ]
colnames (a) <- gsub ("_S[0-9]+.*", "", colnames (a)) 
colnames (a) <- gsub (".*IIT_", "", colnames (a))

a <- merge (a, 

write.xlsx (a, "./output/star_gene_raw_counts.xlsx", rowNames=F)
#write.table (a, "./output/star_gene_raw_counts.txt", sep="\t", quote=F, row.names=F)
