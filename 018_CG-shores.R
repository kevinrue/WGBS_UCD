
# Load libraries ----------------------------------------------------------

library(bsseq)
library(ggplot2)
# library(data.table)
# library(BiocParallel)
# library(tools)
# library(broom)
# library(ggbio)
# library(BSgenome.Btaurus.UCSC.bosTau6)

# Set parameters ----------------------------------------------------------

CpG.gr <- readRDS("bsseq/CpG.gr.rds")

BS.nonEmpty <- readRDS("bsseq/BS.nonEmpty.rds")

# Comparison CpG islands vs shores ----------------------------------------

seqlevels(CpG.gr)
# Rename BS.nonEmpty contigs to match names in CpG.gr
seqlevels(BS.nonEmpty) <- gsub(
  pattern = '^GJ',
  replacement = 'Un_GJ',
  x = seqlevels(BS.nonEmpty))
seqlevels(BS.nonEmpty) <- gsub(
  pattern = '^',
  replacement = 'chr',
  x = seqlevels(BS.nonEmpty))

# Make BS seqlevels comparable to CpG.gr
seqlevels(BS.nonEmpty, force = TRUE) <- seqlevels(CpG.gr)

# Consider all CG, even not covered by a read (BS.unstranded)
# Calculate the average coverage of each CpG island across the 16 samples
CpG.gr$CG.avgTotalCov <- rowMeans(getCoverage(
  BSseq = BS.nonEmpty,
  regions = CpG.gr,
  type = "Cov",
  what = "perRegionTotal"))
CpG.gr$CG.sdTotalCov <- apply(
  X = log10(getCoverage(
    BSseq = BS.nonEmpty,
    regions = CpG.gr,
    type = "Cov",
    what = "perRegionTotal")),
  MARGIN = 1,
  FUN = sd)

ggplot(data = as.data.frame(mcols(CpG.gr))) +
  geom_point(mapping = aes(x = seq_along(CG.avgTotalCov), y = log10(CG.avgTotalCov)))

quantile(log10(rowMeans(CpG.gr$CG.covSum)), na.rm = TRUE)
# About 50% of CG islands cumulate 100+ calls on average across 16 samples
sum(rowMeans(CpG.gr$CG.covSum) > 500, na.rm = TRUE)
# only 2034 CG islands cumulate 500+ calls on average across 16 samples
CpG.gr$shores.avgTotalCov <- rowMeans(
  getCoverage(
    BSseq = BS.nonEmpty,
    regions = flank(x = CpG.gr, width = width(CpG.gr)/2, start = TRUE),
    type = "Cov",
    what = "perRegionTotal") +
    getCoverage(
      BSseq = BS.nonEmpty,
      regions = flank(x = CpG.gr, width = width(CpG.gr)/2, start = FALSE),
      type = "Cov",
      what = "perRegionTotal"))
CpG.gr$shores.sdTotalCov <- apply(
  X = log10(getCoverage(
    BSseq = BS.nonEmpty,
    regions = flank(x = CpG.gr, width = width(CpG.gr)/2, start = TRUE),
    type = "Cov",
    what = "perRegionTotal") +
      getCoverage(
        BSseq = BS.nonEmpty,
        regions = flank(x = CpG.gr, width = width(CpG.gr)/2, start = FALSE),
        type = "Cov",
        what = "perRegionTotal")),
  MARGIN = 1,
  FUN = sd)

# Log10(avg-coverage) of CG islands and shores
ggplot(data = as.data.frame(mcols(CpG.gr))) +
  geom_pointrange(
    mapping = aes(
      x = seq_along(CpG.gr),
      y = log10(shores.avgTotalCov),
      ymin = log10(shores.avgTotalCov) - shores.sdTotalCov,
      ymax = log10(shores.avgTotalCov) + shores.sdTotalCov), colour = "brown") +
  geom_pointrange(
    mapping = aes(
      x = seq_along(CpG.gr),
      y = log10(CG.avgTotalCov),
      ymin = log10(CG.avgTotalCov) - CG.sdTotalCov,
      ymax = log10(CG.avgTotalCov) + CG.sdTotalCov),
    colour = "blue")

# Log-ratio of coverage CGI/shores ~2.5-fold higher for CG
# (makes sense as CGI are more CG-dense)
ggplot(data = as.data.frame(mcols(CpG.gr))) +
  geom_point(mapping = aes(
    x = seq_along(CpG.gr),
    y = log2(CG.avgTotalCov/shores.avgTotalCov)))

# Differential methylation CGI/shores -------------------------------------

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