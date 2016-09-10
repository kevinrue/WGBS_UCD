
library(bsseq)
library(ggplot2)

# Set parameters ----------------------------------------------------------

genomeDir <- 'bostaurus'
outdir <- 'bsseq'

# Load unstranded CG methylation calls
BS <- readRDS(file.path(outdir, "BS.unstranded.rds"))

# Load CGI genomic coordinates
CGI <- readRDS(file.path(outdir, "CpG.gr.rds"))

# Match contig names of CGI coordinates (UCSC) and methylation calls (ensembl)
seqlevels(CGI) <- gsub("^chr(Un_)?", "", seqlevels(CGI))
seqlevels(BS) <- gsub("(\\.1)?", "", seqlevels(BS))

# Subset methylation calls ----

BSinCGI <- subsetByOverlaps(query = BS, subject = CGI)

length(BSinCGI)
length(BS)

# Average coverage of CGs in CGIs
meanCov.BSinCGI <- getCoverage(
  BSseq = BSinCGI, what = "perRegionAverage")
range(meanCov.BSinCGI)
