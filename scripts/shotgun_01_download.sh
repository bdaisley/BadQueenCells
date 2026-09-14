#!/bin/bash

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Step 1: Setup directories
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#Set working directory

working_dir=~/CBQC_project/shotgun
fastq_raw=${working_dir}/01_fastq

mkdir $working_dir
mkdir $fastq_raw


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Step 2: Download FASTQ files
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
cd $fastq_raw

#Download FASTA files from BioProject: PRJNA1501857
#Save to $fastq_raw directory
