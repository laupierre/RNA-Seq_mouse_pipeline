### This is the working place for developing the RNA-Seq pipeline based on singularity

### 1- RNA-Seq: Application for PE Illumina sequencing, in reverse strand orientation

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


The Singularity containers and mouse indexes are currently located in:
/projects/ncrrbt_share_la/dev_pipe/

The containers are built following this example used for STAR:  
apptainer build star.sif singularity_star.def

The RNA-Seq pipeline is found inside the rnaseq_mouse.sh file. In the pipeline, there are three methods available for the quantification: star, kallisto, salmon.
Select the method you want to use (eg --method=salmon , if you want to use salmon  ...etc...) in the following PBS file and submit it to the PBS queue.  
In a new folder, put the 1- rnaseq_mouse_v0.0.1.sh file and 2- the submit.pbs file. Also 3- add all the PE fastq files you want to analyze.

  
######## submit.pbs  
#!/bin/bash  
#PBS -l walltime=20:00:00   
#PBS -l select=1:ncpus=12      
#PBS -q workq  
#PBS -N rnaseq   

cd $PBS_O_WORKDIR   

chmod +x rnaseq_mouse_v0.0.1.sh  
./rnaseq_mouse_v0.0.1.sh --method=star   
######## 

Submit the job using: qsub submit.pbs  
