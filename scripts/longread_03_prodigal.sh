#!/bin/bash

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Step 1: Setup directories
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

working_dir=~/CBQC_project/longread
medaka_dir=${working_dir}/04_medaka_consensus
prodigal_dir=${working_dir}/05_prodigal

mkdir -p $prodigal_dir 

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Step 2: Predict protein sequences from DNA sequences (fna -> faa)
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
cd $medaka_dir

conda activate assembly_tools
for fna in */*.fasta;
 do base_name=$(echo $fna | awk -F'.' '{print $1}');
 prodigal -i $fna -a ${faa_dir}/${base_name}.faa -d ${faa_dir}/${base_name}.fasta
done
