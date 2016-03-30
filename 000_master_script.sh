#!/bin/bash

# Check whether the md5 signature matches the original one
./001_md5.sh rawdata > out.001

# Control quality of the raw reads
./002_fastqc.sh trim_galore '*.fastq.gz' fastqc 12 > out.002

# Collate the output of FastQC in a friendly CSV file
./003_fastqc_summarisation.sh rawdata > out.003

# Merge all paired Fastq files
./004_merge_fastq.sh rawdata rawdata/Merged 4 > out.004
