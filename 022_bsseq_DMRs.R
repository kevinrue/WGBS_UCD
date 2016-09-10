
# Load libraries ----------------------------------------------------------

library(bsseq)

# Set parameters ----------------------------------------------------------

outdir <- 'bsseq'

min.coverage <- 2

min.samples <- 2/3 # proportion

# Import previous data ----------------------------------------------------

# BS <- readRDS(file.path(outdir, "BS.unstranded.rds"))
BS <- readRDS(file.path(outdir, "BS.rmZero.rds"))
sampleNames(BS)
colnames(BS)
colnames(getCoverage(BS))
rownames(colData(BS))
colData(BS)

# Smooth methylation calls ------------------------------------------------

BS.smoothed <- BSmooth(BS, mc.cores = 6, verbose = TRUE)
BS.smoothed
saveRDS(BS.smoothed, file.path(outdir, "BS.smoothed.rds"))
# saveRDS(BS.smoothed, file.path(outdir, "BS.unstranded.smoothed.rds"))

colnames(BS.smoothed)
sampleNames(BS.smoothed)
colnames(getCoverage(BS.smoothed))
colData(BS.smoothed)

# Filter loci with enough calls -------------------------------------------

keepLoci <- which(
  rowSums(
    getCoverage(BSseq = BS.smoothed)[
      , BS.smoothed$Infection == "Control"] >= min.coverage) >= min.samples &
    rowSums(
      getCoverage(BSseq = BS.smoothed)[
        , BS.smoothed$Infection == "M. bovis"] >= min.coverage) >= min.samples)

length(keepLoci)
length(keepLoci) / length(BS.smoothed)

BS.keepLoci <- BS.smoothed[keepLoci,]

saveRDS(BS.keepLoci, file.path(outdir, "BS.keepLoci.rds"))

# Compute t-statistics ----------------------------------------------------

BS.tstat <- BSmooth.tstat(
  BS.keepLoci,
  group1 = paste0("M", levels(BS.keepLoci$Animal)),
  group2 = paste0("C", levels(BS.keepLoci$Animal)),
  estimate.var = "paired",
  local.correct = TRUE,
  verbose = TRUE)

saveRDS(BS.tstat, file.path(outdir, "BS.tstat.rds"))

# Filter NA t-stats -------------------------------------------------------

# plot(BS.tstat) does not cope with NA values
# Remove from tstat all the NA and NaN values
BS.tstat <- BS.tstat[!is.na(getStats(BS.tstat)[,"tstat"]),]

# Plot distribution of t-stats --------------------------------------------

pdf(file = file.path(outdir, "BS.tstat_density.pdf"), width = 9, height = 7)
plot(BS.tstat)
dev.off()

range(getStats(BS.tstat)[,"tstat"])

# Identify DMRs -----------------------------------------------------------

dmrs0 <- dmrFinder(BS.tstat, cutoff = c(-4.6, 4.6))
dmrs <- subset(dmrs0, n >= 3 & abs(meanDiff) >= 0.1)
nrow(dmrs)

ucsc.coord <- function(index, DMRs = dmrs){
  cat(sprintf(
    "%s:%i-%i",
    gsub("^(chr)?", "chr", dmrs[index,"chr"]),
    dmrs[index,"start"],
    dmrs[index,"end"]))
}

bed.coord <- function(index, DMRs = dmrs){
  cat(sprintf(
    "%s\t%i\t%i\tDMR_%04i",
    gsub("^(chr)?", "chr", dmrs[index,"chr"]),
    dmrs[index,"start"],
    dmrs[index,"end"],
    index))
}

write.csv(x = dmrs, file = file.path(outdir, "dmrs.csv"))

# Plotting regions --------------------------------------------------------

BS.keepLoci$col <- c("blue", "red")[as.numeric(BS.keepLoci$Infection)]
colData(BS.keepLoci)

pdf(file = file.path(outdir, "dmr_0002.pdf"), width = 9, height = 7)
plotRegion(BS.keepLoci, dmrs[2,], extend = 5000, addRegions = dmrs)
dev.off()

pdf(file = file.path(outdir, "dmr_0002_GeneCGI.pdf"), width = 9, height = 7)
plotRegion(
  BS.keepLoci, dmrs[2,], extend = 5000,
  addRegions = rbind(
    dmrs[,c("chr","start","end")],
    data.frame(
      chr = rep("20", 2),
      start = c(40967082, 41039108),
      end = c(41041629, 41043399)
      )
    )
  )
dev.off()

pdf(
  file = file.path(outdir, "dmr_0002_GeneCGI_20kb.pdf"),
  width = 9, height = 7)
plotRegion(
  BS.keepLoci, dmrs[2,], extend = 20E3,
  addRegions = rbind(
    dmrs[,c("chr","start","end")],
    data.frame(
      chr = rep("20", 2),
      start = c(40967082, 41039108),
      end = c(41041629, 41043399)
    )
  )
)
dev.off()

# T-stats restricted to autosomes -----------------------------------------

# Subset the results only to the "main" chromosomes
BS.tstat.auto <- subset(
  BS.tstat, as.logical(seqnames(BS.tstat) %in% as.character(1:22)))

pdf(
  file = file.path(outdir, "BS.tstat.auto_density.pdf"),
  width = 9, height = 7)
plot(BS.tstat.auto)
dev.off()

rm(BS.tstat.auto)

# Compare with t-stats using randomised samples ---------------------------

sample.rand <- sample(
  x = sampleNames(BS.keepLoci), size = length(sampleNames(BS.keepLoci)))

BS.tstat.rand <- BSmooth.tstat(
  BS.keepLoci,
  group1 = sample.rand[1:(length(sample.rand)/2)],
  group2 = sample.rand[(length(sample.rand)/2+1):length(sample.rand)],
  estimate.var = "paired",
  local.correct = TRUE,
  verbose = TRUE)
BS.tstat.rand <- BS.tstat.rand[!is.na(getStats(BS.tstat.rand)[,"tstat"]),]
BS.tstat.auto <- subset(
  BS.tstat.rand, as.logical(seqnames(BS.tstat.rand) %in% as.character(1:22)))

pdf(
  file = file.path(outdir, "BS.tstat.auto.RAND_density.pdf"),
  width = 9, height = 7)
plot(BS.tstat.auto)
dev.off()

range(getStats(BS.tstat.auto)[,"tstat"])

plot(density(t.stat.auto))
lines(density(t.stat.auto.rand))
