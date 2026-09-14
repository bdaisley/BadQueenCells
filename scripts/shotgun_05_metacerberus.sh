#!/bin/bash

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Step 1: Setup directories
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
working_dir=~/CBQC_project/shotgun
fastq_nohost_trim_dir=${working_dir}/03_fastq_nohost_trim
metacerberus_dir=${working_dir}/05_metacerberus

mkdir -p $working_dir $metacerberus_dir $metacerberus_dir/fastq_symlinks

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Step 2: Rename _1/_2 to _R1/_R2 FASTQ files
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
cd $fastq_nohost_trim_dir
for i in ${fastq_nohost_trim_dir}/*_nohost_1.fastq; do
  fastq_base=$(basename "$i" | awk -F'\\_1.fastq' '{print $1}') ;
  echo $fastq_base
  ln -s ${fastq_nohost_trim_dir}/${fastq_base}_1.fastq $metacerberus_dir/fastq_symlinks/${fastq_base}_R1.fastq
  ln -s ${fastq_nohost_trim_dir}/${fastq_base}_2.fastq $metacerberus_dir/fastq_symlinks/${fastq_base}_R2.fastq
done

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Step 3: Run metacerberus
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
conda activate metacerberus

metacerberus.py --prodigal $metacerberus_dir/fastq_symlinks \
  --illumina \
  --hmm KOFam_prokaryote \
  --meta --cpus 40 --skip-decon \
  --dir_out ${metacerberus_dir}
