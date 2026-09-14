
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Step 1: Setup directories
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
working_dir=~/CBQC_project/shotgun
fastq_nohost_trim_dir=${working_dir}/03_fastq_nohost_trim
output_dir=${working_dir}/09_eggnog
eggnog_db="/mnt/md0/brendan/dbs/eggnog2_DB"

mkdir -p $working_dir $output_dir $output_dir/fastq_symlinks

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Step 2: Concat all FASTQ files
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
cd $fastq_nohost_trim_dir

for i in ${fastq_nohost_trim_dir}/*_trim_nohost_1.fastq; do
  fastq_base=$(basename "$i" | awk -F'_trim_nohost_1.fastq' '{print $1}')
  echo $fastq_base
  cat ${fastq_nohost_trim_dir}/${fastq_base}_trim_nohost_1.fastq \
    ${fastq_nohost_trim_dir}/${fastq_base}_trim_nohost_2.fastq \
    ${fastq_nohost_trim_dir}/${fastq_base}_trim_unpaired_nohost_1.fastq \
    ${fastq_nohost_trim_dir}/${fastq_base}_trim_unpaired_nohost_2.fastq > $output_dir/fastq_symlinks/${fastq_base}_trim_nohost.fastq
  conda activate eggnog310
  seqkit fq2fa $output_dir/fastq_symlinks/${fastq_base}_trim_nohost.fastq -o $output_dir/fastq_symlinks/${fastq_base}_trim_nohost.fasta
  conda deactivate
done

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Step 3: Run eggnog-mapper2
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
conda activate eggnog310

for i in $output_dir/fastq_symlinks/*.fasta; do
  fastq_base=$(basename "$i" | awk -F'.fasta' '{print $1}')
  rm -r $output_dir/$fastq_base
  mkdir -p $output_dir/$fastq_base
  emapper.py --cpu 24 -i $i \
    --data_dir $eggnog_db \
    --itype metagenome --genepred search \
    --cpu 30 --override \
     --output_dir $output_dir/$fastq_base --output $fastq_base
done
