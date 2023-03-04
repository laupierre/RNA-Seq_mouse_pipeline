### This is the working place for developing the RNA-Seq pipeline based on singularity

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


The Singularity images and mouse indexes are currently located in:
/projects/ncrrbt_share_la/dev_pipe/

The containers are built following this example for STAR:  
apptainer build star.sif singularity_star.def

The RNA-Seq pipeline is inside the rnaseq_mouse.sh file and can be submitted to qsub using the following PBS file.

###########  
#!/bin/bash  
#PBS -l walltime=20:00:00   
#PBS -l select=1:ncpus=10  
#PBS -q workq  
#PBS -N rnaseq   

cd $PBS_O_WORKDIR   

./rnaseq_mouse.sh --method=star   
###########  

In this pipeline, there are three methods available for the quantification: star, kallisto, salmon.
Just modify the file by the method you want to use (eg --method=salmon , if you want salmon  ...etc...) and submit it to the PBS queue.

