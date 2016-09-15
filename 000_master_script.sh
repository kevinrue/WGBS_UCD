l#!/bin/bash

# Check whether the md5 signature matches the original one
./001_md5.sh rawdata > out.001

# Control quality of the raw reads
./002_fastqc_raw.sh -t 12 rawdata '*.fastq.gz' fastqc_raw > out.002

# Merge all paired Fastq files
./003_merge_fastq.sh rawdata rawdata/Merged 12 > out.004

# Control quality of the merged raw reads
./004_fastqc_merged.sh -t 12 rawdata/Merged *fastq.gz fastqc_merged 12 > out.004

# Trim adapter and low-quality bases
./005_trimgalore.sh trimgalore 12 rawdata/Merged > out.005

# Control quality of the trimmed reads
./006_fastqc_trimmed.sh -t 12 trimgalore/Merged *fq.gz fastqc_trimmed > out.006

# Collate trimming statistics
./007_trimgalore_summarisation.sh trimgalore > out.007

# Align trimmed reads of bisulfite libraries
./008_bismark.sh bostaurus trimgalore bismark 4 tmp_bismark > out.008

# TODO Control quality of the aligned reads
#009

# Collate alignment statistics
./010_bismark_summarisation.sh bismark > out.010

# Deduplicate aligned reads / pairs
./011_deduplicate.sh bismark 8 > out.011

# Control quality of the deduplicated reads
./012_fastqc_deduplicated.sh -t 12 bismark/Merged *deduplicated.bam fastqc_deduplicated > out.012

# Collate deduplication statistics
./013_deduplicate_summarisation.sh bismark > out.013

# Extract methylation calls from deduplicated aligned reads / pairs
./014_extractor.sh bismark extract 8 > out.014

# Based on the M-bias plots
# Extract methylation calls again ignoring 7 bp from the 5' of each mate
./014_extractor.sh bismark extract_refined 4 --ignore 7 --ignore_r2 7  > out.014.refined

# TODO Collate methylation extraction statistics
#015

# Sort, add read group, and index BAM files.
./016_sort_index.sh bismark/Merged 12 > out.016

# TODO Collate HTML report for Bismark
./017_bismark_report.sh "bismark extract_refined" reports_bismark 12 > out.017





