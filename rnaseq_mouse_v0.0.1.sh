#!/bin/bash

#### This is the RNA-Seq for the mouse species (v.0.0.1dev)

# clean the log.out file
if [ -f "log.out" ] ; then
    rm log.out
fi


# Call getopt to validate the provided input. 
options=$(getopt -o brg --long method: -- "$@")
# Check status code
[ $? -eq 0 ] || { 
    echo "Incorrect options provided"
    exit 1
}

eval set -- "$options"

while true; do
    case "$1" in
    --method)
        shift; # The arg is next in position args
        COLOR=$1
        [[ ! $COLOR =~ ^(star|kallisto|salmon)$ ]] && {
            echo "Incorrect options. Please provide a correct method: star, kallisto or salmon"
            exit 1
        }
        ;;
    --)
        shift
        break
        ;;
    esac
    shift
done


#################
#### Variables (this folder is not bound by apptainer)
CONTAINER=/projects/ncrrbt_share_la/dev_pipe
INDEX=/projects/ncrrbt_share_la/dev_pipe
CPUS=12

#### Start message
echo "The RNA-Seq pipeline (v.0.1_dev) is for the mouse species and the $COLOR method was chosen" >> log.out
echo "Starting the RNA-Seq pipeline on `date`" >> log.out
##################



##################
#### STAR analysis
if [ "$COLOR" = "star" ]; then
echo "Starting STAR ..." >> log.out

cp -r $CONTAINER/mouse_star_index .

var=(`ls *R1*.fastq.gz`)

	for i in ${var[@]}
	do
	read2=`echo ${i} | sed 's/R1/R2/g'`
	prefix=`echo ${i%%_R1*}`
	
	apptainer exec $CONTAINER/star.sif /bin/bash -c \
	"STAR --genomeDir mouse_star_index --runThreadN $CPUS \
	--readFilesIn $i $read2 \
	--outSAMtype BAM Unsorted \
	--readFilesCommand zcat \
	--outFileNamePrefix star.${prefix}"

	## Put this inside the loop
	if [ $? -eq 0 ]
	then 
    	echo "STAR processed sample ${prefix}" >> log.out
	else
	echo "STAR failed on sample ${prefix}. Pipeline terminated"  >> log.out
	exit 1
	fi
	done
	

mkdir star_results
mv star.IIT* star_results

AFTER=`date`
echo "STAR finished on ${AFTER}" >> log.out



#############################
#### featureCounts after STAR
echo "Starting featureCounts ..." >> log.out
cd star_results
files=`ls -d *bam | xargs -n1000`

apptainer exec $CONTAINER/featurecounts.sif /bin/bash -c \
"featureCounts -B -C -s 2 -p --countReadPairs -T $CPUS -t exon -g gene_id −−extraAttributes gene_name \
-a /root/gtf/gencode.vM32.annotation.gtf \
-o subread.counts.txt $files" || { 
echo "featureCounts has an error. Pipeline terminated" >> log.out
exit 1
}
	
cd ..

AFTER=`date`
echo "featureCounts finished on ${AFTER}" >> log.out



########################################
#### RSeQC analysis and sambamba markdup
echo "Starting RSeQC and sambamba ..." >> log.out

mkdir rseqc_results
cd rseqc_results

wget https://sourceforge.net/projects/rseqc/files/BED/Mouse_Mus_musculus/GRCm39_GENCODE_VM27.bed.gz
gunzip GRCm39_GENCODE_VM27.bed.gz

var=(`ls ../star_results/*bam`)

	for i in ${var[@]}
	do
	prefix=`echo ${i%%_S*}`
	prefix2=`echo ${prefix##*/}`
	apptainer exec $CONTAINER/sambamba.sif /bin/bash -c \
	"sambamba markdup -t $CPUS $i $prefix2.markdup.bam > markdup.$prefix2.log 2>&1" || { 
   	echo "Sambamba has an error. Pipeline terminated" >> log.out
    	exit 1
	}
	done

var=(`ls *.bam`)	
	
	for i in ${var[@]}
	do
	prefix=`echo ${i%%.bam}`
	apptainer exec $CONTAINER/rseqc.sif /bin/bash -c \
	"infer_experiment.py -r GRCm39_GENCODE_VM27.bed -i $i 1> rseqc.$prefix.infer_experiment.txt" || { 
   	echo "RSeQC has an error. Pipeline terminated" >> log.out
    	exit 1
	}
	done

cd ..

AFTER=`date`
echo "RSeQC and sambamba finished on ${AFTER}" >> log.out

fi
####




######################
#### kallisto analysis
if [ "$COLOR" = "kallisto" ]; then
echo "Starting kallisto ..." >> log.out

var=(`ls *R1*.fastq.gz`)

	for i in ${var[@]}
	do
	read2=`echo ${i} | sed 's/R1/R2/g'`
	prefix=`echo ${i%%_R1*}`

	apptainer exec $CONTAINER/kallisto.sif /bin/bash -c \
	"kallisto quant -i /root/kallisto/mouse_transcripts.idx -t $CPUS --rf-stranded -o kallisto.$prefix $i $read2 2> $prefix.kallisto.log"

	## Put this inside the loop
	if [ $? -eq 0 ]
	then 
    	echo "kallisto processed sample ${prefix}" >> log.out
	else
	echo "kallisto failed on sample ${prefix}. Pipeline terminated"  >> log.out
	exit 1
	fi
	done
	
mkdir kallisto_results
mv kallisto.IIT* kallisto_results
mv *kallisto.log kallisto_results

AFTER=`date`
echo "kallisto finished on ${AFTER}" >> log.out

fi
####




######################################
#### salmon analysis in automatic mode (should fall back on ISR mode)
if [ "$COLOR" = "salmon" ]; then
echo "Starting salmon ..." >> log.out

cp -r $CONTAINER/mouse_salmon_index .

var=(`ls *R1*.fastq.gz`)

	for i in ${var[@]}
	do
	read2=`echo ${i} | sed 's/R1/R2/g'`
	prefix=`echo ${i%%_R1*}`
	
	apptainer exec $CONTAINER/salmon.sif /bin/bash -c \
   	"salmon quant -i mouse_salmon_index -p $CPUS -l A --validateMappings -o salmon.${prefix} -1 $i -2 $read2"
	
	## Put this inside the loop
	if [ $? -eq 0 ]
	then
    	echo "salmon processed sample ${prefix}" >> log.out
	else
	echo "salmon failed on sample ${prefix}. Pipeline terminated"  >> log.out
	exit 1
	fi
	done

mkdir salmon_results
mv salmon.IIT* salmon_results

AFTER=`date`
echo "salmon finished on ${AFTER}" >> log.out

fi
####




####################
#### FastQC analysis
echo "Starting FastQC ..." >> log.out

files=`ls *fastq.gz | xargs -n1000`
mkdir fastqc_results

apptainer exec $CONTAINER/fastqc.sif /bin/bash -c \
	"fastqc -t $CPUS -o fastqc_results $files"

AFTER=`date`
echo "FastQC finished on ${AFTER}" >> log.out
####




####################
#### MultiQC analysis
echo "Starting MultiQC ..." >> log.out

if [ "$COLOR" = "star" ]; then
apptainer exec $CONTAINER/multiqc.sif /bin/bash -c \
"multiqc -f -n multiqc_report_rnaseq \
-m featureCounts $PBS_O_WORKDIR/star_results/*summary \
-m star $PBS_O_WORKDIR/star_results/*Log.final.out \
-m sambamba $PBS_O_WORKDIR/rseqc_results/markdup.star.IIT*.log \
-m rseqc $PBS_O_WORKDIR/rseqc_results/*infer_experiment.txt \
-m fastqc $PBS_O_WORKDIR/fastqc_results/*zip"
fi

if [ "$COLOR" = "salmon" ]; then
apptainer exec $CONTAINER/multiqc.sif /bin/bash -c \
"multiqc -f -n multiqc_report_rnaseq \
-m salmon $PBS_O_WORKDIR/salmon_results/* \
-m fastqc $PBS_O_WORKDIR/fastqc_results/*zip"
fi

if [ "$COLOR" = "kallisto" ]; then
apptainer exec $CONTAINER/multiqc.sif /bin/bash -c \
"multiqc -f -n multiqc_report_rnaseq \
-m kallisto $PBS_O_WORKDIR/kallisto_results/*.kallisto.log \
-m fastqc $PBS_O_WORKDIR/fastqc_results/*zip"
fi

AFTER=`date`
echo "MultiQC finished on ${AFTER}" >> log.out
####



#################
#### Exit message
echo "The RNA-Seq pipeline was completed on `date`. Please check the files that were produced!" >> log.out
exit 0;
#################


