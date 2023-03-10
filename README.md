### This is the working place for developing a RNA-Seq pipeline based on singularity

### 1- Mouse RNA-Seq: PE Illumina sequencing, in reverse strand orientation

This is the development version deposited for the production team: the current version is v0.0.1_dev

The singularity images are available for:  
BBMap version 39.01  
FastQC version 0.11.9  
featureCounts version 2.0.4  
kallisto version 0.40.0  
multiqc version 1.14  
rseqc version 5.0.1  
salmon version 1.10.0  
sambamba version 0.8.2  
STAR version 2.7.10b  
R version 4.2.2 Patched     


The Singularity containers and mouse indexes are located in:
/projects/ncrrbt_share_la/dev_pipe/

The containers are built following this example used for STAR:  
apptainer build star.sif singularity_star.def

See mouse_indexes for details on how the different indexes are made.

   
In a new folder, put 1- the submit.pbs file below and 2- the samples_rnaseq.xlsx containing the samples information.
This excel file has 3 mandatory columns: sample   condition   replicate   

Finally, add 3- all the PE fastq files you want to analyze that are described in the samples_rnaseq.xlsx file.

In this pipeline, there are three methods available for the RNA quantification: star, kallisto, salmon. Select the method you want to use (eg --method=salmon, if you want to use salmon, or --method=star or --method=kalllisto) in the PBS file before launching the main command: qsub submit.pbs.
The differential expression is handled by DESeq2 for two-groups comparison as described in the condition column of the samples_rnaseq.xlsx file. 
The comparison is made on the condition by alphabetical order. For example, Control and Treated conditions, will be take Control as reference.

######## submit.pbs  
#!/bin/bash  
#PBS -l walltime=20:00:00   
#PBS -l select=1:ncpus=12      
#PBS -q workq  
#PBS -N rnaseq   

cd $PBS_O_WORKDIR   

cp /projects/ncrrbt_share_la/dev_pipe/rnaseq_mouse_v0.0.1.sh .  
chmod +x rnaseq_mouse_v0.0.1.sh

./rnaseq_mouse_v0.0.1.sh --method=star   
######## 

When the analysis is done, the raw counts in case of STAR, or the length scaled TPM counts of genes in case of salmon and kallisto can be found in the output folder. 
Differential expression for two-groups was performed with DESeq2 after each quantification.
 

├── kallisto_deseq2_differential_expression.xlsx  
├── kallisto_gene_lengthScaledTPM_counts.xlsx  
├── kallisto_log.out  
├── MA_plot.pdf  
├── multiqc_report_rnaseq.html  
├── PCA_plot.pdf  
├── salmon_deseq2_differential_expression.xlsx  
├── salmon_gene_lengthScaledTPM_counts.xlsx  
├── salmon.log.out  
├── star_deseq2_differential_expression.xlsx  
├── star_gene_raw_counts.xlsx  
└── star.log.out    
