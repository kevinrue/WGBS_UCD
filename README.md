# WGBS_UCD
Analysis of whole-genome bisulfite sequencing (WGBS) data.

## Overview

This pipeline is being developed for the analysis of a particular data set, with the hope that it may be useful for the analysis of other data sets.

The pipeline currently includes:

- Control quality of raw FASTQ files ([FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
- Merge raw (paired-end) FASTQ files from multiple sequencing runs for each biological replicate ([cat](http://manpages.ubuntu.com/manpages/saucy/man1/cat.1.html))
- Trim adapter and low-quality bases from merged FASTQ files ([Trim Galore!](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
- Alignment of trimmed read pairs to bisulfite-converted genome ([Bismark](http://www.bioinformatics.babraham.ac.uk/projects/bismark/))
- Deduplication of aligned read pairs ([Bismark](http://www.bioinformatics.babraham.ac.uk/projects/bismark/))
- Extraction of methylation calls from deduplicated aligned read pairs ([Bismark](http://www.bioinformatics.babraham.ac.uk/projects/bismark/))


## Important notes

- Bismark: using v0.15.0 cloned from [GitHub](https://github.com/FelixKrueger/Bismark) to benefit of a recent bug fix allowing the proper use of the `--cytosine_report` option of the script `bismark_methylation_extractor`. The bug fix is **not yet** available in the version available on the [Bismark website](http://www.bioinformatics.babraham.ac.uk/projects/download.html#bismark)
