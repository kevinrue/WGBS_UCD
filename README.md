# WGBS_UCD
Analysis of whole-genome bisulfite sequencing (WGBS) data.

## Overview

This pipeline is being developed for the analysis of a particular data set, with the hope that it may be useful for the analysis of other data sets.

The pipeline currently includes {scripts}:

- :white_check_mark: Checksum of downloaded files ([sum](http://man7.org/linux/man-pages/man1/sum.1.html)) {001}
- :white_check_mark: Control quality of raw FASTQ files ([FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)) {002}
- :white_check_mark: Merge raw (paired-end) FASTQ files from multiple sequencing runs for each biological replicate ([cat](http://manpages.ubuntu.com/manpages/saucy/man1/cat.1.html)) {003}
- :white_check_mark: Control quality of merged FASTQ files ([FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)) {004}
- :white_check_mark: Trim adapter and low-quality bases from merged FASTQ files ([Trim Galore!](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)) {005}
- :white_check_mark: Control quality of trimmed FASTQ files ([FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)) {006}
- :white_check_mark: Report of trimming statistics {007}
- :white_check_mark: Alignment of trimmed read pairs to bisulfite-converted genome ([Bismark](http://www.bioinformatics.babraham.ac.uk/projects/bismark/)) {008}
- :white_check_mark: Control quality of aligned BAM files ([FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)) {009}
- :white_check_mark: Report of alignment statistics {010}
- :white_check_mark: Deduplication of aligned read pairs ([Bismark](http://www.bioinformatics.babraham.ac.uk/projects/bismark/)) {011}
- :white_check_mark: Control quality of deduplicated BAM files ([FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)) {012}
- :white_check_mark: Report of deduplication statistics {013}
- :white_check_mark: Extraction of methylation calls from deduplicated aligned read pairs ([Bismark](http://www.bioinformatics.babraham.ac.uk/projects/bismark/)) {014}
- Upcoming :red_circle: Report of methylation extraction statistics {015}
- :white_check_mark: Addition of read group, sorting and indexing of BAM files ([Picard](http://broadinstitute.github.io/picard/)) {016}
- :white_check_mark: Combined Bismark report in HTML format ([Bismark](http://www.bioinformatics.babraham.ac.uk/projects/bismark/)) {017}
- :white_check_mark: Import methylation calls and format GRanges for subsequent analyses ([bsseq/Bioconductor](https://bioconductor.org/packages/release/bioc/html/bsseq.html)) {018}
- :white_check_mark: Import coordinates of unmasked CpG islands (CGI) and format as GRanges ([UCSC tracks](http://genome.ucsc.edu/cgi-bin/hgTables)) {019}
- :white_check_mark: Import coordinates of genes and format as GRanges ([biomaRt/Bioconductor](https://bioconductor.org/packages/release/bioc/html/biomaRt.html)) {020}
- :white_check_mark: Overview of CGI coverage in the imported data set {021}
- :white_check_mark: Identification of differentially methylated regions (DMRs) ([bsseq/Bioconductor](https://bioconductor.org/packages/release/bioc/html/bsseq.html)) {021}


## Important notes

- Bismark: using v0.15.0 cloned from [GitHub](https://github.com/FelixKrueger/Bismark) to benefit of a recent bug fix allowing the proper use of the `--cytosine_report` option of the script `bismark_methylation_extractor`. The bug fix is **not yet** available in the version available on the [Bismark website](http://www.bioinformatics.babraham.ac.uk/projects/download.html#bismark)
