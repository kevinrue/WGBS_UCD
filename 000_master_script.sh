#!/bin/bash

# Check whether the md5 signature matches the original one
./001_md5.sh rawdata > out.001

# Control quality of the raw reads
./002_fastqc.sh trim_galore '*.fastq.gz' fastqc 12 > out.002

# Collate the output of FastQC in a friendly CSV file
./003_fastqc_summarisation.sh rawdata > out.003

# Trim adapter and low quality nucleotides
./004_trimgalore.sh rawdata trim_galore 12 > out.004

# Control quality of the raw reads
./005_fastqc.sh trim_galore '*.fq.gz' fastqc_trimmed 12 > out.005

# Collate the output of FastQC in a friendly CSV file
./006_fastqc_trimmed_summarisation.sh fastqc_trimmed > out.006

# Collate the output of trim_galore in a friendly CSV file
./007_trimgalore_summarisation.sh fastqc_trimmed > out.007

# Align bisulfite libraries (single-end, paired-end, trimmed unpaired)
./008_bismark_align.sh <genome> <rootdir> <outdir> <theads> <tempdir>

# Collate the output of bismark in a friendly CSV file
./009_bismark_summarisation.sh bismark > out.009

# Control quality of the aligned reads
./010_fastqc.sh bismark *.bam fastqc_aligned 6 > out.010

# TODO: ./011 Collate the output of FastQC in a friendly CSV file


# Remove duplicate reads / pairs aligned
./012_deduplicate.sh bismark 6 > out.012


