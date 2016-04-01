#!/bin/bash

# Check whether the md5 signature matches the original one
./001_md5.sh rawdata > out.001

# TODO Control quality of the raw reads
./002_fastqc.sh rawdata '*.fastq.gz' fastqc 12 > out.002

# TODO Collate the output of FastQC in a friendly CSV file
./003_fastqc_summarisation.sh rawdata > out.003

# Merge all paired Fastq files
./004_merge_fastq.sh rawdata rawdata/Merged 12 > out.004

# TODO Control quality of the merged raw reads
#005

# TODO Collate the output of FastQC in a friendly CSV file
#006

# TODO Trim adapter and low-quality bases
./007_trimgalore.sh

# TODO Control quality of the trimmed reads
#008

# TODO Collate the output of FastQC in a friendly CSV file
#009

# TODO Collate trimming statistics
#010

# TODO Align trimmed reads of bisulfite libraries
./011_bismark.sh
