library(rtracklayer)
library(biomaRt)

system ("wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/gencode.vM32.annotation.gtf.gz")

my_obj <- import("gencode.vM32.annotation.gtf.gz")
my_obj <- as.data.frame (my_obj)
my_obj <- my_obj[my_obj$type == "transcript", ]
my_obj <- my_obj[ ,c("gene_id", "transcript_id", "gene_type", "gene_name", "mgi_id")]


ensembl <- useEnsembl(biomart = 'genes', 
                       dataset = 'mmusculus_gene_ensembl',
                       version = 109)

values <- unique (my_obj$gene_id)
length (values)
# 56953

res <- getBM(attributes = c('ensembl_gene_id_version', 'external_gene_name', 'description'),
             filters = 'ensembl_gene_id_version',
             values = values, mart = ensembl)      
dim (res)
# 56953

my_obj <- merge (my_obj, res, by.x="gene_id", by.y="ensembl_gene_id_version", all.x=TRUE)

write.table (my_obj, "gencode.vM32.annotation.txt", sep="\t", quote=F, row.names=F)


