#!/bin/bash

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Step 1: Setup directories
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#Set working directory t

working_dir=~/CBQC_project/longread
fastq_dir=${working_dir}/01_fastq
fastq_trim_dir=${working_dir}/01_fastq_trim
fastq_downsample_dir=${working_dir}/01_fastq_downsample
fastq_filtlong_dir=${working_dir}/01_fastq_trim_filtlong
flye_dir=${working_dir}/02_flye
flyemeta_dir=${working_dir}/02_flye_meta
combined_dir=${working_dir}/03_manual_combined
medaka_dir=${working_dir}/04_medaka_consensus

mkdir -p $fastq_trim_dir $fastq_downsample_dir $fastq_filtlong_dir $flye_dir $flyemeta_dir $combined_dir $medaka_dir

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Step 2: Assembly with flye/metaflye
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
conda activate assembly_tools

#Try multiple assembly methods to get best genome

for fastq in ${fastq_dir}/*.fastq ; do
  THREADS_N=10
  fastq_base=$(basename "$fastq" | awk -F'\\.fastq' '{print $1}')
  filtlong_length_cut=2000
  filtlong --min_length $filtlong_length_cut --min_mean_q 20 --keep_percent 75 ${fastq_dir}/${fastq_base}.fastq > ${fastq_filtlong_dir}/${fastq_base}.fastq
  fastplong -n 0 --cut_front --cut_tail --cut_mean_quality 20 \
    --disable_quality_filtering --disable_length_filtering \
    --low_complexity_filter --complexity_threshold 30 \
    --html ${fastq_trim_dir}/${fastq_base}.html \
    --json ${fastq_trim_dir}/${fastq_base}.json \
    --verbose -i $fastq -o ${fastq_trim_dir}/${fastq_base}.fastq
  #-----------------------
  #Downsample
  #-----------------------
  coverage_set=25
  gsize=$(lrge --threads $THREADS_N $fastq)
  echo "Estimate genome size = $gsize"
  rasusa reads -g $gsize --output-type u -c $coverage_set -o ${fastq_downsample_dir}/${fastq_base}_down${coverage_set}.fastq $fastq
  #-----------------------
  #Downsample 100
  #-----------------------
  coverage_set=100
  gsize=$(lrge --threads $THREADS_N $fastq)
  echo "Estimate genome size = $gsize"
  rasusa reads -g $gsize --output-type u -c $coverage_set -o ${fastq_downsample_dir}/${fastq_base}_down${coverage_set}.fastq $fastq
  #-----------------------
  #Filtlong
  flye -t $THREADS_N --min-overlap 1000 --genome-size 2m --asm-coverage 3000 --nano-corr ${fastq_filtlong_dir}/${fastq_base}.fastq -o ${flye_dir}/${fastq_base}_flye_filtlong${filtlong_length_cut}
  flye -t $THREADS_N --min-overlap 1000 --meta --nano-corr ${fastq_filtlong_dir}/${fastq_base}.fastq -o ${flyemeta_dir}/${fastq_base}_metaflye_filtlong${filtlong_length_cut}
  #fastplong (trim length disabled)
  flye -t $THREADS_N --min-overlap 1000 --genome-size 2m --asm-coverage 3000 --nano-corr ${fastq_trim_dir}/${fastq_base}.fastq -o ${flye_dir}/${fastq_base}_flye_trimmed_lendisabled
  flye -t $THREADS_N --min-overlap 1000 --meta --nano-corr ${fastq_trim_dir}/${fastq_base}.fastq -o ${flyemeta_dir}/${fastq_base}_metaflye_trimmed_lendisabled
  #No filtering
  flye -t $THREADS_N --min-overlap 1000 --genome-size 2m --asm-coverage 3000 --nano-corr ${fastq_dir}/${fastq_base}.fastq -o ${flye_dir}/${fastq_base}_flye
  flye -t $THREADS_N --min-overlap 1000 --meta --nano-corr ${fastq_dir}/${fastq_base}.fastq -o ${flyemeta_dir}/${fastq_base}_metaflye
  #No filtering - big overlap
  overlap=7000
  flye -t $THREADS_N --min-overlap $overlap --scaffold --genome-size 2m --asm-coverage 3000 --nano-corr ${fastq_dir}/${fastq_base}.fastq -o ${flye_dir}/${fastq_base}_flye_overlap7000
  flye -t $THREADS_N --min-overlap $overlap --scaffold --meta --nano-corr ${fastq_dir}/${fastq_base}.fastq -o ${flyemeta_dir}/${fastq_base}_metaflye_overlap7000
done


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Step 3: Inspect graphs with Bandage to resolve repeats/tangles
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Bandage load assembly_graph.gfa
#Bandage image assembly_graph.gfa assembly_graph.svg

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Step 4: Polish with medaka
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#Manually add best genomes to $combined_dir

conda activate medaka
for fastq in ${combined_dir}/*.fasta; do
  #local fastq=$1
  fastq_base=$(basename "$fastq" | awk -F'\\.fasta' '{print $1}')
  medaka_consensus -i ${fastq_dir}/${fastq_base}.fastq -d ${combined_dir}/${fastq_base}.fasta \
    -o ${medaka_dir}/${fastq_base} -t 30 --bacteria -f -x
  seqkit sort -r -l ${medaka_dir}/${fastq_base}/consensus.fasta > ${medaka_dir}/${fastq_base}/consensus_sorted.fasta
done

