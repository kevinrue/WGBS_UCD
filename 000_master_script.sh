l#!/bin/bash

# Check whether the md5 signature matches the original one
./001_md5.sh rawdata > out.001

# Control quality of the raw reads
./002_fastqc.sh rawdata '*.fastq.gz' fastqc 12 > out.002

# Collate the output of FastQC in a friendly CSV file
./003_fastqc_summarisation.sh rawdata > out.003

# Merge all paired Fastq files
./004_merge_fastq.sh rawdata rawdata/Merged 12 > out.004

# Control quality of the merged raw reads
./005_fastqc.sh rawdata/Merged *fastq.gz fastqc_merged 12 > out.005

# Collate the output of FastQC in a friendly CSV file
./006_fastqc_summarisation.sh fastqc_merged > out.006

# Trim adapter and low-quality bases
./007_trimgalore.sh trimgalore 12 rawdata/Merged > out.007

# Control quality of the trimmed reads
./008_fastqc.sh trimgalore/Merged *fq.gz fastqc_trimmed 12 > out.008

# Collate the output of FastQC in a friendly CSV file
./009_fastqc_summarisation.sh fastqc_trimmed > out.009

# Collate trimming statistics
./010_trimgalore_summarisation.sh trimgalore > out.010

# Align trimmed reads of bisulfite libraries
./011_bismark.sh bostaurus trimgalore bismark 4 tmp_bismark > out.011

# TODO Control quality of the aligned reads
#012

# TODO Collate the output of FastQC in a friendly CSV file
#013

# Collate alignment statistics
./014_bismark_summarisation.sh bismark > out.014

# Deduplicate aligned reads / pairs
./015_deduplicate.sh bismark 8 > out.015

# Control quality of the deduplicated reads
./016_fastqc.sh bismark/Merged *deduplicated.bam fastqc_deduplicated 12 > out.016

# Collate the output of FastQC in a friendly CSV file
./017_fastqc_summarisation.sh fastqc_deduplicated > out.017

# Extract methylation calls from deduplicated aligned reads / pairs
./018_extractor.sh bismark extract 8 > out.018

# Based on the M-bias plots
# Extract methylation calls again ignoring 7 bp from the 5' of each mate
./018_extractor.sh bismark extract 4 --ignore 7 --ignore_r2 7  > out.018.refined

# Sort, add read group, and index BAM files.
./019_sort_index.sh bismark/Merged 12 > out.019







