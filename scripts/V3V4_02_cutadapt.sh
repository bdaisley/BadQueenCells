#!/bin/bash

#This script assumes you have Cutadapt installed within a conda environment named "cutadapt".

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Step 1: Setup directories
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#Set working directory t

working_dir=~/CBQC_project
fastq_raw=${working_dir}/00_fastq_raw
fastq_trim=${working_dir}/01_fastq_trim

mkdir $fastq_trim


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Step 2: Set primers of interest
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#v3v4 primer set
v3v4_F="CCTACGGGNGGCWGCAG"
v3v4_R="GACTACHVGGGTATCTAATCC"

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Step 3: Trim primers from demultiplexed sequencing reads
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
conda activate cutadapt

foo () {
  local fastq=$1
  base_name=$(echo $fastq | awk -F'_1' '{print $1}' | awk -F '/' '{print $NF == "" ? $(NF - 1) : $NF}')
  mkdir ${fastq_trim}
  cutadapt -e 0.2 --discard-untrimmed \
    -g $v3v4_F -G $v3v4_R \
    -o ${fastq_trim}/${base_name}_1.fastq.gz \
    -p ${fastq_trim}/${base_name}_2.fastq.gz \
    ${fastq_raw}/${base_name}_1.fastq.gz ${fastq_raw}/${base_name}_2.fastq.gz
}

for fastq in ${fastq_raw}/*_1.fastq.gz;
  do running=($(jobs -rp))
  while [ ${#running[@]} -ge 10 ] ; do
    sleep 1
    running=($(jobs -rp))
  done
  foo "$fastq" & done

wait
