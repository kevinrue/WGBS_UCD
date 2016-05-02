
# Load libraries ----------------------------------------------------------

library(data.table)
library(BiocParallel)
library(bsseq)

# Set parameters ----------------------------------------------------------

CpG.folder <- 'extract_refined/Merged'

Cpg.file.pattern <- 'CpG_report.txt.gz'

CPUs <- 8

outdir <- 'bsseq'

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

# Import methylation calls ------------------------------------------------

# Code only working in the current devel branch
# BS.combined2 <- read.bismark(
#   files = file.path(CpG.folder, pdata$Filename[1:2]),
#   sampleNames = c("A","B"),
#   rmZeroCov = FALSE,
#   strandCollapse = FALSE,
#   fileType = "cytosineReport",
#   verbose = TRUE)

read.CpGreport.gz <- function(
  index, pData, folder = ".", sample = "Sample", file = "Filename"
){
  fread_cmd <- paste("gunzip -c", file.path(folder, pData[index, file]))
  cat(pData[index, sample], ":", fread_cmd, sep = " ", fill = TRUE)
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

BS.list <- bplapply(
  1:nrow(pdata), 
  read.CpGreport.gz,
  pData = pdata,
  folder = CpG.folder,
  BPPARAM = MulticoreParam(workers = CPUs))

BS.combined <- combineList(BS.list)

# lapply(X = BS.list, FUN = sampleNames)
colnames(BS.combined)
sampleNames(BS.combined)
colnames(getCoverage(BS.combined))
sort(colnames(BS.combined))

saveRDS(object = BS.combined, file = file.path(outdir, 'BS.combined.rds'))
# BS.combined <- readRDS(file.path(outdir, 'BS.combined.rds'))

rm(BS.list)

# gz1 <- gzfile(file.path(outdir, M.file), "w")
# write.table(
#   x = getBSseq(BSseq = BS.combined, type = "M"),
#   file = gz1,
#   row.names = FALSE)
# close(gz1)
# 
# gz1 <- gzfile(file.path(outdir, Cov.file), "w")
# write.table(
#   x = getBSseq(BSseq = BS.combined, type = "Cov"),
#   file = gz1,
#   row.names = FALSE)
# close(gz1)
# 
# gz1 <- gzfile(file.path(outdir, gr.file), "w")
# write.table(
#   x = data.frame(
#     chr = seqnames(getBSseq(BSseq = BS.combined, type = "gr")),
#     pos = start(getBSseq(BSseq = BS.combined, type = "gr")),
#     strand = strand(getBSseq(BSseq = BS.combined, type = "gr"))
#   ),
#   file = gz1,
#   row.names = FALSE)
# close(gz1)

# Raw statistics ----------------------------------------------------------

meanCov.stranded <- getCoverage(BSseq = BS.combined, what = "perRegionAverage")
covStats <- as.data.frame(meanCov.stranded)
write.csv(x = covStats, file = file.path(outdir, "covStats.csv"))
summary(meanCov.stranded)
rm(meanCov.stranded)

# Collapse information from both both strands -----------------------------

BS.unstranded <- strandCollapse(BSseq = BS.combined, shift = TRUE)

colnames(BS.unstranded)
sampleNames(BS.unstranded)
colnames(getCoverage(BS.unstranded))
sort(colnames(BS.unstranded))

saveRDS(object = BS.unstranded, file = file.path(outdir, 'BS.unstranded.rds'))
# BS.unstranded <- readRDS(file.path(outdir, 'BS.unstranded.rds'))

rm(BS.combined)

# Statistics of strand-collapsed calls ------------------------------------

meanCov.unstranded <- getCoverage(
  BSseq = BS.unstranded, what = "perRegionAverage")
covStats <- cbind(covStats, as.data.frame(meanCov.unstranded))
write.csv(x = covStats, file = file.path(outdir, "covStats.csv"))
summary(meanCov.unstranded)
rm(meanCov.unstranded)

tapply(X = getCoverage(BSseq = BS.unstranded, what = "perRegionAverage"), INDEX = BS.unstranded$Infection, FUN = sum)

# Count of non-zero CpG per sample
colSums(getCoverage(BSseq = BS.unstranded, what = "perBase") > 0) / 1E6
summary(colSums(getCoverage(BSseq = BS.unstranded, what = "perBase") > 0) / 1E6)

# Discard CG with zero coverage -------------------------------------------

# sum(rowSums(getCoverage(BSseq = BS.unstranded)) == 0)
# sum(rowSums(getCoverage(BSseq = BS.unstranded)) > 0) / nrow(BS.unstranded)
# 
BS.rmZero <- BS.unstranded[rowSums(getCoverage(BSseq = BS.unstranded)) > 0,]

colnames(BS.rmZero)
sampleNames(BS.rmZero)
colnames(getCoverage(BS.rmZero))
sort(colnames(BS.rmZero))

saveRDS(object = BS.rmZero, file = file.path(outdir, 'BS.rmZero.rds'))
# BS.rmZero <- readRDS(file.path(outdir, 'BS.rmZero.rds'))

rm(BS.unstranded)

dim(BS.rmZero)

# Statistics of non-zero CpGs ------------------------------------

meanCov.rmZero <- getCoverage(BSseq = BS.rmZero, what = "perRegionAverage")
covStats <- cbind(covStats, meanCov.rmZero <- as.data.frame(meanCov.rmZero))
write.csv(x = covStats, file = file.path(outdir, "covStats.csv"))
summary(meanCov.rmZero)
rm(meanCov.rmZero)

# Keep only CG with coverage in all samples -------------------------------

BS.allCovered <- BS.rmZero[rowSums(getCoverage(BSseq = BS.rmZero) > 0) == ncol(BS.rmZero),]

colnames(BS.allCovered)
sampleNames(BS.allCovered)
colnames(getCoverage(BS.allCovered))
sort(colnames(BS.allCovered))

saveRDS(object = BS.allCovered, file = file.path(outdir, 'BS.allCovered.rds'))
# BS.allCovered <- readRDS(file.path(outdir, 'BS.allCovered.rds'))

rm(BS.rmZero)

dim(BS.allCovered)

# Statistics of CG coverage in all samples --------------------------------

meanCov.allCovered <- getCoverage(
  BSseq = BS.allCovered, what = "perRegionAverage")
covStats <- cbind(covStats, as.data.frame(meanCov.allCovered))
write.csv(x = covStats, file = file.path(outdir, "covStats.csv"))
summary(meanCov.allCovered)
rm(meanCov.allCovered)

# Testing some functions of the bsseq package -----------------------------

# goodness of fit statistics for BSSeq objects
# For each methylation loci, the Poisson goodness of fit statistic tests 
# whether the coverage (at that loci) is independent and identically Poisson
# distributed across the samples.
pois.fit <- poissonGoodnessOfFit(BSseq = BS.rmZero)
plot(pois.fit)
bino.fit <- binomialGoodnessOfFit(BSseq = BS.rmZero)
plot(bino.fit)
