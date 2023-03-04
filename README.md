### This is the working place for developing the RNA-Seq pipeline based on singularity

This is the development version deposited for the production team: the current version is v0.0.1_dev

The singularity images are available for:
BBMap version 39.01, 
FastQC version 0.11.9,
featureCounts version 2.0.4,
kallisto version 0.40.0,
multiqc version 1.14,
rseqc version 5.0.1,
salmon version 1.10.0,
STAR version 2.7.10b.


The Singularity images and mouse indexes are currently located in:
/projects/ncrrbt_share_la/dev_pipe/

The images are built following this example:
apptainer star.sif singularity_star.def


