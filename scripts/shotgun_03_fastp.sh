#!/bin/bash


::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Step 1: Setup directories
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
working_dir=~/CBQC_project/shotgun
fastq_nohost_dir=${working_dir}/02_fastq_nohost
fastq_nohost_trim_dir=${working_dir}/03_fastq_nohost_trim

mkdir -p $fastq_nohost_trim_dir


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Step 2: Run fastp to filter poor quality reads
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

conda activate assembly_tools
cd $fastq_nohost_dir

for i in ${fastq_nohost_dir}/*_1.fastq; do
  fastq_base=$(basename "$i" | awk -F'\\_1.fastq' '{print $1}') ;
  mkdir -p ${fastq_trim}_QCreports
  fastp -i ${fastq_base}_1.fastq -I ${fastq_base}_2.fastq \
    -o $fastq_nohost_trim_dir/${fastq_base}_trim_1.fastq  \
    -O $fastq_nohost_trim_dir/${fastq_base}_trim_2.fastq \
    --thread 16 --html ${fastq_nohost_trim_dir}_QCreports/${fastq_base}_report.html \
    --json /dev/null \
    --detect_adapter_for_pe -l 100 --cut_right --cut_right_window_size 4 --cut_right_mean_quality 15
done
