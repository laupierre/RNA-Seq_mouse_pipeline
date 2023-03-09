### This is the working place for developing the RNA-Seq pipeline based on singularity

### 1- Mouse RNA-Seq: Application for PE Illumina sequencing, in reverse strand orientation

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

   
In a new folder, put 1- the submit.pbs file below and 2- the samples_rnaseq.xlsx containing the samples. Finally, add 3- all the PE fastq files you want to analyze.

The the samples_rnaseq.xlsx is an excel file containing the 3 mandatory columns
sample	condition	replicate

The RNA-Seq pipeline is found inside the rnaseq_mouse.sh file and submitted to PBS using submit.pbs. In this pipeline, there are three methods available for the RNA quantification: star, kallisto, salmon. Select the method you want to use (eg --method=salmon, if you want to use salmon, or --method=star or --method=kalllisto) in the PBS file before launching the main command: qsub submit.pbs 

######## submit.pbs  
#!/bin/bash  
#PBS -l walltime=20:00:00   
#PBS -l select=1:ncpus=12      
#PBS -q workq  
#PBS -N rnaseq   

cd $PBS_O_WORKDIR   

cp /projects/ncrrbt_share_la/dev_pipe/rnaseq_mouse_v0.0.1.sh $PBS_O_WORKDIR
chmod +x rnaseq_mouse_v0.0.1.sh

./rnaseq_mouse_v0.0.1.sh --method=star   
######## 

When the analysis is done, a log file is appearing in the output folder, together with the raw counts in case of STAR, or with the length scaled TPM counts of genes in case of salmon and kallisto. Thes counts can be further processed for differential expression.
 
