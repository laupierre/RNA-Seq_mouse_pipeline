library(rtracklayer)

system ("wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/gencode.vM32.annotation.gtf.gz")

my_obj <- import("gencode.vM32.annotation.gtf.gz")
my_obj <- as.data.frame (my_obj)
my_obj <- my_obj[my_obj$type == "transcript", ]
my_obj <- my_obj[ ,c("gene_id", "transcript_id", "gene_type", "gene_name", "mgi_id")]
write.table (my_obj, "gencode.vM32.annotation.txt", sep="\t", quote=F, row.names =F)
