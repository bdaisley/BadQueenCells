#!/bin/bash

#This script assumes you have VSEARCH and USEARCH installed and available in your path, and Cutadapt installed within a conda environment named "cutadapt".

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Step 1: Setup directories & variables
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
working_dir=~/CBQC_project
fastq_dir=${working_dir}/01_fastq_trim
fastq_dir_subset=${working_dir}/01_fastq_trim_subset
output_dir=${working_dir}/02_denoised

mkdir $fastq_dir_subset
mkdir $output_dir

R1_length="250"
R2_length="250"
min_length="395"
min_qual="2"

conda activate cutadapt

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Step 2: Copy FASTQ files and rename files
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
cd ${fastq_dir_subset}
cp -r ${fastq_dir}/*[12].fastq.gz ${fastq_dir_subset}
pigz -d ${fastq_dir_subset}/*.gz --force
find . -name "*_1.fastq" -type f | rename 's/_1.fastq/_R1.fastq/g'
find . -name "*_2.fastq" -type f | rename 's/_2.fastq/_R2.fastq/g'


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Step 3: Filter FASTQ files
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#Remove in case any old files
rm -r ${fastq_dir_subset}/*.filt

foo () {
  local fastq=$1
  base_name=$(echo $fastq | awk -F'_R1.fastq' '{print $1}' | awk -F '/' '{print $NF == "" ? $(NF - 1) : $NF}')
  cutadapt -l ${R1_length} -o ${base_name}_R1.fastq.filt ${base_name}_R1.fastq
  cutadapt -l ${R2_length} -o ${base_name}_R2.fastq.filt ${base_name}_R2.fastq
  }

for fastq in *_R1.fastq* ;
  do running=($(jobs -rp))
  while [ ${#running[@]} -ge 16 ] ; do
    sleep 1
    running=($(jobs -rp))
  done
  foo "$fastq" & done

wait

#::::::::::::::::::::::::::::::::::
# Step 4: Merge FASTQ files
#::::::::::::::::::::::::::::::::::
cd ${fastq_dir_subset}
for r1_file in *_R1*.fastq.filt; do
  # Define the corresponding R2 file
  r2_file="${r1_file/_R1/_R2}"
  base_name=$(basename "${r1_file}" _R1.fastq.filt)
  vsearch --fastq_mergepairs "$r1_file" --reverse "$r2_file" \
    --fastqout "${fastq_dir_subset}/${base_name}_merged.fastq" \
    --fastq_minmergelen "${min_length}" \
    --log "${fastq_dir_subset}/${base_name}_merged_report.txt" --relabel ${base_name}.
  cat ${fastq_dir_subset}/${base_name}_merged_report.txt >> ${output_dir}/merged_report.txt
  cat ${fastq_dir_subset}/${base_name}_merged.fastq >> ${output_dir}/merged.fastq
done

#::::::::::::::::::::::::::::::::::
# Step 5: quality filter based on expected error
#::::::::::::::::::::::::::::::::::
cd ${output_dir}
vsearch -fastq_filter ${output_dir}/merged.fastq \
 -fastaout ${output_dir}/merged_filtered.fasta \
 -fastq_maxee 1.0 \
 --fastq_qmax 42

#::::::::::::::::::::::::::::::::::
# Step 6: dereplicate
#::::::::::::::::::::::::::::::::::
vsearch --fastx_uniques ${output_dir}/merged_filtered.fasta \
  -fastaout ${output_dir}/merged_filtered_derep.fasta \
  -sizeout -relabel Uniq

#::::::::::::::::::::::::::::::::::
# Step 7: denoise with UNOISE3
#::::::::::::::::::::::::::::::::::
usearch -unoise3 ${output_dir}/merged_filtered_derep.fasta \
  -zotus ${output_dir}/merged_filtered_derep_zotus.fasta \
  -ampout ${output_dir}/amplicons.fa \
  -tabbedout ${output_dir}/unoise3.txt

#::::::::::::::::::::::::::::::::::
# Step 8: make zOTU table
#::::::::::::::::::::::::::::::::::
vsearch --fastq_filter ${output_dir}/merged.fastq -fastaout ${output_dir}/merged.fasta --fastq_qmax 42

vsearch --usearch_global ${output_dir}/merged.fasta \
  --db ${output_dir}/merged_filtered_derep_zotus.fasta \
  --notrunclabels --id 0.97 --otutabout ${output_dir}/FINAL_zotutab.txt \
  --maxaccepts 1 --maxrejects 0
  