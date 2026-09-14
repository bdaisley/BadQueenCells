#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Step 1: Setup directories
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
working_dir=~/CBQC_project/longread
medaka_dir=${working_dir}/04_medaka_consensus
output_dir=${working_dir}/08_VIBRANT

mkdir -p $working_dir $output_dir $fna_dir_renamed ${output_dir}/vcontact3_results

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Step 2: Run VIBRANT
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

conda activate vibrant
cat ${output_dir}/*.fna > ${output_dir}/combined.fna
VIBRANT_run.py -t 12 -i ${output_dir}/combined.fna -folder $output_dir \
  -l 1000 \
  -d /mnt/md2/brendan/dbs/VIBRANT/databases \
  -m /mnt/md2/brendan/dbs/VIBRANT/files

conda activate assembly_tools
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Step 3: Convert GBK to TSV
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
python scripts/longread_06_VIBRANT.py \
  $output_dir/VIBRANT_combined/VIBRANT_phages_combined/combined.phages_combined.gbk \
  $output_dir/VIBRANT_combined/VIBRANT_phages_combined/combined.phages_combined.tsv
#Fix quotes
sed 's/"""//g' $output_dir/VIBRANT_combined/VIBRANT_phages_combined/combined.phages_combined.tsv > \
  $output_dir/VIBRANT_combined/VIBRANT_phages_combined/combined.phages_combined_fixed.tsv
