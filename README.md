# WGBS_UCD
Analysis of whole-genome bisulfite sequencing (WGBS) data.

## Overview

This pipeline is being developed for the analysis of a particular data set, with the hope that it may be useful for the analysis of other data sets.

The pipeline currently includes {scripts}:

- :white_check_mark: Checksum of downloaded files ([sum](http://man7.org/linux/man-pages/man1/sum.1.html)) {001}
- :white_check_mark: Control quality of raw FASTQ files ([FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)) {002-003}
- :white_check_mark: Merge raw (paired-end) FASTQ files from multiple sequencing runs for each biological replicate ([cat](http://manpages.ubuntu.com/manpages/saucy/man1/cat.1.html)) {004}
- :white_check_mark: Control quality of merged FASTQ files ([FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)) {005-006}
- :white_check_mark: Trim adapter and low-quality bases from merged FASTQ files ([Trim Galore!](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)) {007}
- :white_check_mark: Control quality of trimmed FASTQ files ([FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)) {008-009}
- :white_check_mark: Report of trimming statistics {010}
- :white_check_mark: Alignment of trimmed read pairs to bisulfite-converted genome ([Bismark](http://www.bioinformatics.babraham.ac.uk/projects/bismark/)) {011}
- Upcoming :red_circle: Control quality of aligned BAM files ([FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)) {012-013}
- :white_check_mark: Report of alignment statistics {014}
- :white_check_mark: Deduplication of aligned read pairs ([Bismark](http://www.bioinformatics.babraham.ac.uk/projects/bismark/)) {015}
- :white_check_mark: Control quality of deduplicated BAM files ([FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)) {016-017}
- :white_check_mark: Report of deduplication statistics {018}
- :white_check_mark: Extraction of methylation calls from deduplicated aligned read pairs ([Bismark](http://www.bioinformatics.babraham.ac.uk/projects/bismark/)) {019}
- Upcoming :red_circle: Report of alignment statistics {020}
- :white_check_mark: Addition of read group, sorting and indexing of BAM files {021}
- Upcoming :red_circle: Combined Bismark report in HTML format ([Bismark](http://www.bioinformatics.babraham.ac.uk/projects/bismark/)) {022}

## Important notes

- Bismark: using v0.15.0 cloned from [GitHub](https://github.com/FelixKrueger/Bismark) to benefit of a recent bug fix allowing the proper use of the `--cytosine_report` option of the script `bismark_methylation_extractor`. The bug fix is **not yet** available in the version available on the [Bismark website](http://www.bioinformatics.babraham.ac.uk/projects/download.html#bismark)
