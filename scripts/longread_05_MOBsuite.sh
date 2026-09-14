#!/bin/bash

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Step 1: Setup directories
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
working_dir=~/CBQC_project/longread
medaka_dir=${working_dir}/04_medaka_consensus
output_dir=${working_dir}/07_MOBsuite

mkdir $output_dir

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Step 2: Pull out plasmids from whole genome assemblies
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

cd $output_dir
conda activate mob_suite

for i in ${medaka_dir}/*.fna; do
  fna_base=$(basename "$i" | awk -F'\\.fna' '{print $1}')
  echo "Processing $fna_base"

  # Skip if any plasmid_*.gbk already exists
  if compgen -G "$output_dir/$fna_base/plasmid_*.gbk" > /dev/null; then
    echo "⚠️  Skipping $fna_base — plasmid_*.gbk already exists."
    continue
  fi

  # Run MOB-recon
  mob_recon --infile "${medaka_dir}/${fna_base}.fna" \
            --outdir "$output_dir/$fna_base" \
            --run_overhang --force \
            -d /mnt/md1/brendan/dbs/mob_database

  # Path to corresponding GBK file
  gbk_file="${working_dir}/06_metacerberus/final/genbank/prodigal_${fna_base}_template.gbk"

  # Loop over all MOB-recon plasmid FASTA outputs
  for fasta in "$output_dir/$fna_base"/plasmid_*.fasta; do
    out_gbk="${fasta%.fasta}.gbk"
    echo -n "" > "$out_gbk"

    grep "^>" "$fasta" | sed 's/^>//' | cut -d' ' -f1 | while read -r contig_name; do
      awk -v contig="$contig_name" '
        $1 == "LOCUS" && $2 == contig {in_block=1}
        in_block {print}
        /^\/\// {in_block=0}
      ' "$gbk_file"
    done >> "$out_gbk"
    echo "✓ Extracted contigs from $fasta → $out_gbk"
  done
done

