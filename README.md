# Case of the Bad Queen Cells Project

This repository contains data-processing workflows, analysis scripts, and supporting code used for in the Case of the Bad Queen Cell (CBQC) project investigating an outbreak of idiopathic disease affecting developing honey bee (<i>Apis mellifera</i>) queens. 

RMarkdown document with complete walkthrough of analysis available here: https://bdaisley.github.io/BadQueenCells

### Overview:
During spring 2024, abnormal mortality was observed among developing queens from commercial queen-rearing operations across Canada. Affected queens typically died shortly after pupation, often during the white-eyed stage, and exhibited delayed development, yellow-to-black discoloration, softened cuticle, and a deflated or fluid-filled appearance. Although affected queens were frequently positive for black queen cell diseas (BQCV), the observed pathology differed from classical descriptions of BQCV. The aim of this project was characterize the microbial factors underlying this disease phenotype and to identify candidate etiological agents.

### Publication:
<i>Daisley et al. 2026. Queen-killing variants of Melissococcus plutonius exhibit phage-mediated genome diversification and are associated with black queen cell virus co-infection. In draft.</i>

### Raw data sources:
FASTQ sequence files have been uploaded to the NCBI Sequence Read Archive (SRA) under the following BioProject accessions:
- [PRJNA1525705](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1525705) - 16S rRNA gene amplicon sequencing data
- [PRJNA1509667](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1509667) - Long-read whole genome sequencing of <i>M. plutonius</i> ATCC 35311
- [PRJNA1501857](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1501857) - Shotgun metagenomic sequencing & long-read whole genome sequencing of queen-derived <i>M. plutonius</i> isolates

The assembled complete genomes of <i>M. plutonius </i> strains from this study are available under the following RefSeq accessions: 
- GCF_059945265.1 (qk204) - https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_059945265.1
- GCF_059945345.1 (qk205) - https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_059945345.1
- GCF_059945185.1 (qk206) - https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_059945185.1
- GCF_059945305.1 (qk223) - https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_059945305.1
- GCF_059945285.1 (qk224) - https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_059945285.1
- GCF_059945245.1 (qk225) - https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_059945245.1
- GCF_060436625.1 (ATCC 35311) - https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_060436625.1

### 16S rRNA microbiota profiling workflow
Microbial community composition was characterized in asymptomatic and symptomatic queen pupae using amplicon sequencing of the V3-V4 16S rRNA gene sequencing. DNA was extracted from homogenized queen pupae, followed by amplification using the established Bakt_341F (5'-CCTACGGGNGGCWGCAG-3') and Bakt_805R (5'-GACTACHVGGGTATCTAATCC-3') primer pair. Amplicon libraries were generated and sequencedusing the Illumina NextSeq 1000 platform with 2x300 bp paired-end sequencing. The following scripts were used for processing of raw sequencing data:
```
- scripts/V3V4_01_download.sh
- scripts/V3V4_02_cutadapt.sh
- scripts/V3V4_03_fastqc.sh
- scripts/V3V4_04_denoising.sh
- scripts/V3V4_05_classification.R
```

### Shotgun metagenomic sequencing workflow
```
- scripts/shotgun_01_download.sh
- scripts/shotgun_02_bowtie2.sh
- scripts/shotgun_03_fastp.sh
- scripts/shotgun_04_kraken2.sh
- scripts/shotgun_05_metacerberus.sh
- scripts/shotgun_06_eggNOG_mapper2.sh
```
### Longread genome sequencing and assembly of <i>M. plutonius</i> strains
```
- scripts/longread_01_download.sh
- scripts/longread_02_genome_assembly_pipeline.sh
- scripts/longread_03_prodigal.sh
- scripts/longread_04_metacerberus.sh
- scripts/longread_04_metacerberus_hydra.py
- scripts/longread_05_MOBsuite.sh
- scripts/longread_06_VIBRANT.py
- scripts/longread_06_VIBRANT.sh
```