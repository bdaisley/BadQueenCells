#!/bin/bash

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Step 1: Setup directories
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#Set working directory t

working_dir=~/CBQC_project
fastq_raw=${working_dir}/00_fastq_raw

mkdir $working_dir
mkdir $fastq_raw


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Step 2: Download FASTQ files
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
cd $fastq_raw

#Use https://sra-explorer.info/ to get download links, then copy and paste the curl commands as follows:

curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40570013/SRR40570013_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40570013/SRR40570013_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40570012/SRR40570012_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40570012/SRR40570012_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40570011/SRR40570011_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40570011/SRR40570011_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40570010/SRR40570010_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40570010/SRR40570010_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40570009/SRR40570009_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40570009/SRR40570009_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40570008/SRR40570008_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40570008/SRR40570008_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40570007/SRR40570007_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40570007/SRR40570007_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40570006/SRR40570006_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40570006/SRR40570006_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40570005/SRR40570005_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40570005/SRR40570005_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40570004/SRR40570004_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40570004/SRR40570004_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40570003/SRR40570003_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40570003/SRR40570003_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40570002/SRR40570002_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40570002/SRR40570002_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40570001/SRR40570001_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40570001/SRR40570001_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40570000/SRR40570000_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40570000/SRR40570000_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40569999/SRR40569999_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40569999/SRR40569999_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40569998/SRR40569998_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40569998/SRR40569998_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40569997/SRR40569997_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40569997/SRR40569997_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40569996/SRR40569996_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40569996/SRR40569996_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40569995/SRR40569995_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40569995/SRR40569995_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40569994/SRR40569994_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40569994/SRR40569994_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40569993/SRR40569993_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40569993/SRR40569993_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40569992/SRR40569992_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40569992/SRR40569992_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40569991/SRR40569991_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40569991/SRR40569991_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40569990/SRR40569990_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40569990/SRR40569990_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40569989/SRR40569989_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40569989/SRR40569989_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40569988/SRR40569988_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40569988/SRR40569988_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40569987/SRR40569987_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40569987/SRR40569987_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40569986/SRR40569986_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40569986/SRR40569986_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40569985/SRR40569985_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR405/042/SRR40569985/SRR40569985_2.fastq.gz

