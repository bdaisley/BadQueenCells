#!/bin/bash

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Step 1: Setup directories
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
working_dir=~/CBQC_project/shotgun
fastq_dir=${working_dir}/01_fastq
fastq_nohost_dir=${working_dir}/02_fastq_nohost

mkdir -p $working_dir $fastq_nohost_dir


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Step 2: Build bowtie2 db
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
conda activate bowtie2

bowtie2-build /mnt/md1/brendan/dbs/Genomes_reference/Amel/GCA_003254395.2_Amel_HAv3.1_genomic.fna \
 /mnt/md1/brendan/dbs/Genomes_reference/Amel/Amel_bowtie2_index


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Step 3: Map reads against Apis mellifera geneome
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
cd ${fastq_dir}

for fastq in *_1.fastq ;
   base_name=$(echo $fastq | awk -F'_1.fastq' '{print $1}')
   bowtie2 -p 16 -x /mnt/md1/brendan/dbs/Genomes_reference/Amel/Amel_bowtie2_index  \
   -1 ${base_name}_1.fastq \
   -2 ${base_name}_2.fastq \
   --very-sensitive-local \
   --un-conc \
   ${base_name} > ${base_name}_mapped_and_unmapped.sam
   mv ${base_name}.1 ${fastq_nohost_dir}/${base_name}_nohost_1.fastq
   mv ${base_name}.2 ${fastq_nohost_dir}/${base_name}_nohost_2.fastq
   mv ${base_name}_mapped_and_unmapped.sam ${fastq_nohost_dir}/${base_name}_mapped_and_unmapped.sam
done
