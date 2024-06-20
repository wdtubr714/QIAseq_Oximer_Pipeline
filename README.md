# QIAseq_Oximer_Pipeline

#!/usr/bin/bash

### Install Anaconda3/2024.02-1 and python3/3.9.6
module load Anaconda3/2024.02-1
module load python3/3.9.6
module load java/22.0.1

### Confirm Anaconda3/2024.02-1 is up-to-date
conda update --update-all -y

### Set up conda channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

### Create conda environment
conda env create mirna_env -y
conda activate mirna_env

### Install all supporting packages for Anaconda3 to run alignment pipeline
conda install -c bioconda seqkit fastqc cutadapt bowtie samtools seqcluster mirtop multiqc wget bbmap umi_tools bwa -y

### Install R and necessary R packages
conda install -c conda-forge r-base r-essentials -y
Rscript -e "install.packages('tidyverse', repos='http://cran.rstudio.com/')"
Rscript -e "install.packages('data.table', repos='http://cran.rstudio.com/')"
Rscript -e "install.packages('readr', repos='http://cran.rstudio.com/')"

### Ensure GNU parallel is installed
conda install -c conda-forge parallel -y 

### ORDER OF RUNNING
mirna_alignment_QIAseq.sh
unmapped_to_hg38.sh
seqbuster_mirtop_multiqc.sh

#### Must provide a samples.txt document along with the collapse_mirtop.r file, 
#### downloadable hairpin.fa, and mature.fa file from miRBase
#### samples.txt file should indicate the file path of each sample (no headers)
#### Run all ".sh" in "$BASE_DIR"

### Versions used:
#### conda 24.1.2, seqkit 2.8.2, FastQC v0.12.1, Cutadapt 4.8, bowtie version 1.3.1, samtools 1.20, seqcluster 1.2.9, mirtop 0.4.25, multiqc version 1.22.2, GNU Wget 1.21.4, umi_tools 1.1.5, mirdeep2 2.0.1.3, bwa 0.7.18, bbmap 39.06
#### Rscript (R) version 4.3.3 (2024-02-29)
#### GNU Parallel 20240522 ('Tbilisi') v1
