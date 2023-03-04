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
#### Environments (this folder is not bound by apptainer)
CONTAINER=/projects/ncrrbt_share_la/dev_pipe
INDEX=/projects/ncrrbt_share_la/dev_pipe

#### Start message
echo "The RNA-Seq pipeline (v.0.1_dev) is for the mouse species and the $COLOR method was chosen" >> log.out
BEFORE=`date`
echo "Starting the RNA-Seq pipeline on ${BEFORE}" >> log.out
##################

