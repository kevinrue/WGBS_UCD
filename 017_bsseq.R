
# Load libraries ----------------------------------------------------------

library(bsseq)
library(data.table)
library(BiocParallel)
library(tools)
library(broom)

# Set parameters ----------------------------------------------------------

CpG.folder <- 'extract_refined/Merged'

Cpg.file.pattern <- 'CpG_report.txt.gz'

CPUs <- 2

outdir <- 'bsseq'

M.file <- 'M.txt.gz'
Cov.file <- 'Cov.txt.gz'
gr.file <- 'gr.txt.gz'

CpGislands.file <- "bostaurus/Bt_UMD31_CpG_islands_unmasked.bed"

# Create required folders -------------------------------------------------

if (file.access(outdir) != 0){
  dir.create(outdir)
}

# Assemble phenotypic information -----------------------------------------

CpG.files <- list.files(path = CpG.folder, pattern = Cpg.file.pattern)

pdata <- data.frame(
  Sample = sapply(X = strsplit(CpG.files, '_'), FUN = "[[", 1),
  Infection = factor(
    x = c("Control", "M. bovis")[
      as.numeric(factor(
        x = sapply(
          X = sapply(X = strsplit(x = CpG.files, split = '_'), FUN = "[[", 1),
          FUN = "substr", start = 1, stop = 1),
        levels = c("C", "M")))],
    levels = c("Control", "M. bovis")),
  Animal = as.factor(sapply(
    X = sapply(X = strsplit(x = CpG.files, split = '_'), FUN = "[[", 1),
    FUN = "gsub",
    pattern = "[CM]",
    replacement = "")),
  Filename = CpG.files,
  stringsAsFactors = FALSE
)
write.csv(x = pdata, file = file.path(outdir, 'phenodata.csv'), row.names = FALSE)

read.CpGreport.gz <- function(
  index, pData, folder = ".", sample = "Sample", file = "Filename"
  ){
  fread_cmd <- paste("gunzip -c", file.path(folder, pData[index, file]))
  message(fread_cmd)
  dat <- as.data.frame(fread(fread_cmd))
  gr <- GRanges(
    seqnames = dat[,1],
    ranges = IRanges(start = dat[,2], end = dat[,2]),
    strand = dat[,3])
  # Methylation calls
  M <- as.matrix(dat[,4], drop=FALSE)
  # Total calls
  Cov <- as.matrix(dat[,4] + dat[,5], drop=FALSE)
  BSseq(
    M = M,
    Cov = Cov,
    pData = pData[index,],
    sampleNames = pData[index, sample],
    gr = gr)
}

# Import methylation calls ------------------------------------------------

# Issue with read.bismark
# sampleNames
# Sent an email to peter.hickey@gmail.com

# BS.combined1 <- read.bismark(
#   files = file.path(CpG.folder, pdata$Filename[1]),
#   sampleNames = pdata$Sample[1],
#   rmZeroCov = FALSE,
#   strandCollapse = FALSE,
#   fileType = "cytosineReport",
#   verbose = TRUE)
# 
# BS.combined2 <- read.bismark(
#   files = file.path(CpG.folder, pdata$Filename[1:2]),
#   sampleNames = c("A","B"),
#   rmZeroCov = FALSE,
#   strandCollapse = FALSE,
#   fileType = "cytosineReport",
#   verbose = TRUE)

BS.list <- bplapply(
  1:nrow(pdata), 
  read.CpGreport.gz,
  pData = pdata,
  folder = CpG.folder,
  BPPARAM = MulticoreParam(workers = CPUs))

BS.combined <- combineList(BS.list)

saveRDS(object = BS.combined, file = file.path(outdir, 'BS.combined.rds'))
# BS.combined <- readRDS(file.path(outdir, 'BS.combined.rds'))

rm(BS.list)

# TODO: clean up code and apply to other slots if file is of decent size
gz1 <- gzfile(file.path(outdir, M.file), "w")
write.table(
  x = getBSseq(BSseq = BS.combined, type = "M"),
  file = gz1,
  row.names = FALSE)
close(gz1)

gz1 <- gzfile(file.path(outdir, Cov.file), "w")
write.table(
  x = getBSseq(BSseq = BS.combined, type = "Cov"),
  file = gz1,
  row.names = FALSE)
close(gz1)

gz1 <- gzfile(file.path(outdir, gr.file), "w")
write.table(
  x = data.frame(
    chr = seqnames(getBSseq(BSseq = BS.combined, type = "gr")),
    pos = start(getBSseq(BSseq = BS.combined, type = "gr")),
    strand = strand(getBSseq(BSseq = BS.combined, type = "gr"))
  ),
  file = gz1,
  row.names = FALSE)
close(gz1)

# Raw statistics ----------------------------------------------------------

meanCov.stranded <- colMeans(getCoverage(BSseq = BS.combined))
as.data.frame(meanCov.stranded)
summary(meanCov.stranded)

# Collapse information from both both strands -----------------------------

BS.unstranded <- strandCollapse(BSseq = BS.combined, shift = TRUE)
saveRDS(object = BS.unstranded, file = file.path(outdir, 'BS.unstranded.rds'))
# BS.unstranded <- readRDS(file.path(outdir, 'BS.unstranded.rds'))

rm(BS.combined)

# Statistics of strand-collapsed calls ------------------------------------

meanCov.unstranded <- colMeans(getCoverage(BSseq = BS.unstranded))
as.data.frame(meanCov.unstranded)
summary(meanCov.unstranded)

# Discard CG with zero coverage -------------------------------------------

# sum(rowSums(getCoverage(BSseq = BS.unstranded)) == 0)
# sum(rowSums(getCoverage(BSseq = BS.unstranded)) > 0) / nrow(BS.unstranded)
# 
# BS.nonEmpty <- BS.unstranded[rowSums(getCoverage(BSseq = BS.unstranded)) > 0,]
# saveRDS(object = BS.nonEmpty, file = file.path(outdir, 'BS.nonEmpty.rds'))
# BS.nonEmpty <- readRDS(file.path(outdir, 'BS.nonEmpty.rds'))

# rm(BS.unstranded)

# dim(BS.nonEmpty)

# Statistics of non-zero CpGs ------------------------------------

# meanCov.nonEmpty <- colMeans(getCoverage(BSseq = BS.nonEmpty))
# as.data.frame(meanCov.nonEmpty)
# summary(meanCov.nonEmpty)
# 
# pois.fit <- poissonGoodnessOfFit(BSseq = BS.nonEmpty)
# plot(pois.fit)
# bino.fit <- binomialGoodnessOfFit(BSseq = BS.nonEmpty)
# plot(bino.fit)

# Download CG islands -----------------------------------------------------

# Manually download CpG islands tracks from:
# http://genome.ucsc.edu/cgi-bin/hgTables
# Downloaded both masked and unmasked for testing

CpG.unmasked <- read.table(
  file = CpGislands.file,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE)

CpG.gr <- GRanges(
  seqnames = CpG.unmasked[,1],
  ranges = IRanges(
    start = CpG.unmasked[,2],
    end = CpG.unmasked[,3],
    names = CpG.unmasked[,4]))

rm(CpG.unmasked)

# Count methylations/calls per CpG island ---------------------------------

# seqlevels(CpG.gr)
# Rename BS.nonEmpty contigs to match names in CpG.gr
# seqlevels(BS.nonEmpty) <- gsub(
#   pattern = '^GJ',
#   replacement = 'Un_GJ',
#   x = seqlevels(BS.nonEmpty))
# seqlevels(BS.nonEmpty) <- gsub(
#   pattern = '^',
#   replacement = 'chr',
#   x = seqlevels(BS.nonEmpty))

# Subset both objects to overlapping contigs (save memory)
# seqlevels.intersect <- intersect(
#   seqlevels(BS.nonEmpty),
#   seqlevels(CpG.gr)
# )
# seqlevels(CpG.gr, force = TRUE) <- seqlevels.intersect
# seqlevels(BS.nonEmpty, force = TRUE) <- seqlevels.intersect

# For each CpG island, count the number of CG with at least a call
# CpG.gr$CG.covered <- countOverlaps(query = CpG.gr, subject = BS.nonEmpty)

# CpG.gr[1,]
# BS.nonEmpty[1,]
# Find/count methylation calls overlapping CpG islands
# subsetByOverlaps(query = BS.nonEmpty[1,], subject = CpG.gr[1,])
# getCoverage(subsetByOverlaps(query = BS.nonEmpty[1,], subject = CpG.gr[1,]))
# colSums(getCoverage(subsetByOverlaps(query = BS.nonEmpty[1,], subject = CpG.gr[1,])))
# getBSseq(subsetByOverlaps(query = BS.nonEmpty[1,], subject = CpG.gr[1,]), type = "M")
# getMeth(subsetByOverlaps(query = BS.nonEmpty[1,], subject = CpG.gr[1,]), type = "raw")
# colSums(getMeth(subsetByOverlaps(query = BS.nonEmpty[1,], subject = CpG.gr[1,])))
# tmp.base <- getMeth(BSseq = BS.nonEmpty, regions = CpG.gr, type = "raw", what = "perBase")
# tmp.region <- getMeth(BSseq = BS.nonEmpty, regions = CpG.gr, type = "raw", what = "perRegion")
# 
# Function to test differential methylation level in a paired design
# Pairs of samples must be ordered: group1rep1, group1rep2, ..., group2rep1, group2rep2, ...
# region.paired.t.test <- function(x, factor = rep(1:(length(x)/2), 2)){
#   values1 <- x[1:(length(x)/2)]
#   values2 <- x[(length(x)/2+1):length(x)]
#   if (sum(!is.na(values1+values2)) < 2){
#     return(rep(NA,6))
#   }
#   tidy(t.test(
#     x = values1,
#     y = values2,
#     alternative = "two.sided",
#     paired = TRUE))
# }
# 
# Run the paired t-test on all CpG islands with at least two pairs of samples
# covered by at least one read
# region.paired.t.test(x = tmp.region)
# 
# paired.t.test <- do.call(
#   rbind,
#   apply(X = tmp.region, MARGIN = 1, FUN = region.paired.t.test))
# QQ-plot
# plot(
#   x = -log10(sort(paired.t.test$p.value)),
#   y = -log10(sort(runif(n = sum(!is.na(paired.t.test$p.value)), min = 0, max = 1))))
# rm(paired.t.test)

# Comparison CpG islands vs shores ----------------------------------------

# Consider all CG, even not covered by a read (BS.unstranded)
# Calculate the coverage of each CpG island and its shores
# (Define the shores as length/2 of the island each side)
# Plot, and define a cutoff to retain a reasonable number of CG islands
# For the remaining islands, get their mean methylation level
# get the mean methylation level of the shores (combined)
# for that, combine their methylation calls and divide by their combined coverage
# Calculate the paired difference between replicates (M. bovis - control)
# Calculate mean, mean-sd, mean+sd
# Sort by increasing mean difference
# Plot mean +/-sd represented by error bars using ggplot2 

