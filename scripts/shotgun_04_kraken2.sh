#!/bin/bash

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Step 1: Setup directories
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
working_dir=~/CBQC_project/shotgun
fastq_dir=${working_dir}/03_fastq_nohost_trim
kraken_output_dir=${working_dir}/04_KRAKEN2
kraken_tools="/mnt/md2/brendan/dbs/KRAKEN2/KrakenTools-master"
kraken_db_hash="/mnt/md1/brendan/dbs/KRAKEN2_GTDB_220_NCBI_merged/KRAKEN2_GTDB_220_NCBI_merged/KRAKEN2_GTDB_220_NCBI_merged"

#Input directories
export SAMPLES_DIR=${fastq_dir}
export K2_database_DIR=${kraken_db_hash}
export K2_database_DIR2=${kraken_db_hash2}
export K2_tools_DIR="${kraken_tools}"

#Output directories
export tax_out=${kraken_output_dir}
export K2_outputs_DIR=${kraken_output_dir}/k2_outputs
export K2_reports_DIR=${kraken_output_dir}/k2_reports
export B2_abundance_DIR=${kraken_output_dir}/bracken_abundance_files

mkdir -p $tax_out $tax_out/fastq_symlinks
mkdir -p $K2_outputs_DIR
mkdir -p $K2_reports_DIR ${K2_reports_DIR}_round2
mkdir -p $B2_abundance_DIR ${B2_abundance_DIR}_round2


#::::::::::::::::::::::::::::::::::::::::::::
# Step 2: Classify reads with KRAKEN2
#::::::::::::::::::::::::::::::::::::::::::::

conda activate kraken2_latest

#Navigate to folder containing fastq files:
cd $SAMPLES_DIR

for i in *_1.fastq; do
  filename=$(basename "$i")
  fname="${filename%*_*.fastq}"
  kraken2 --db $K2_database_DIR \
  --threads 16 \
  --confidence 0.1 \
  --use-names \
  --memory-mapping \
  --unclassified-out $K2_outputs_DIR/${fname}#.fastq \
  --output /dev/null  \
  --report $K2_reports_DIR/${fname}_report.txt \
  --paired ${fname}_1.fastq ${fname}_2.fastq
done


#::::::::::::::::::::::::::::::::::::::::::::
# Step 3: Run BRAKEN2 abundance estimation
#::::::::::::::::::::::::::::::::::::::::::::

# Set up directories [ ---BRACKEN2 OUTPUT DIRECTORIES--- ]
mkdir $K2_reports_DIR/bracken
mkdir $K2_reports_DIR/bracken/species
mkdir $K2_reports_DIR/bracken/species/mpa
mkdir $K2_reports_DIR/bracken/species/mpa/combined
mkdir $tax_out/bracken_abundance_files

#Abundance estimation with bracken [ ---SPECIES LEVEL--- ]
cd $K2_reports_DIR

for i in *_report.txt
do
  filename=$(basename "$i")
  fname="${filename%_report.txt}"
  bracken -d $K2_database_DIR -i $i -r 150 -t 10 -l S -o ${fname}_report_species.txt
done
rm *_report_species.txt
mv *_bracken_species.txt bracken/species/.

#Generating combined abundance tables [ ---MPA FORMAT--- ]
for i in bracken/species/*_report_bracken_species.txt
do
  filename=$(basename "$i")
  fname="${filename%_report_bracken_species.txt}"
  python $K2_tools_DIR/kreport2mpa.py -r $i -o bracken/species/mpa/${fname}_mpa.txt --display-header
done

#Combined all samples and reorganize folders [ ---CLEAN UP--- ]

python ${K2_tools_DIR}/combine_mpa.py -i bracken/species/mpa/*_mpa.txt -o bracken/species/mpa/combined/combined_species_mpa.txt
grep -E "(s__)|(#Classification)" bracken/species/mpa/combined/combined_species_mpa.txt > bracken/species/mpa/combined/bracken_abundance_species_mpa.txt
sed -i -e 's/_report_bracken_species.txt//g' bracken/species/mpa/combined/bracken_abundance_species_mpa.txt
cp bracken/species/mpa/combined/bracken_abundance_species_mpa.txt $B2_abundance_DIR/.
sed -i -e 's/_report_bracken.txt//g' bracken/species/mpa/combined/bracken_abundance_species_mpa.txt
