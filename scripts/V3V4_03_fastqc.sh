#!/bin/bash

#This script assumes you have FastQC installed within a conda environment named "fastqc".

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Step 1: Setup directories
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
working_dir="/mnt/md1/brendan/projects/2024_09_06___BCMystery"
fastq_raw=${working_dir}/00_fastq_raw
fastq_trim=${working_dir}/01_fastq_trim
fastq_trim_QC=${working_dir}/01_fastq_trim_QC

mkdir $fastq_trim_QC

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Step 2: Inspect quality of sequences with FastQC
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

cd $fastq_trim
conda activate fastqc
for i in *_1.fastq.gz; do
  i_base=$(basename $i "_1.fastq.gz")
  mkdir -p ${fastq_trim_QC}/${i_base}_1 ${fastq_trim_QC}/${i_base}_2
  fastqc ${i_base}_1.fastq.gz -o ${fastq_trim_QC}/${i_base}_1
  fastqc ${i_base}_2.fastq.gz -o ${fastq_trim_QC}/${i_base}_2
done