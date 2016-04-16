
# Load libraries ----------------------------------------------------------

library(bsseq)
library(data.table)
library(BiocParallel)
library(tools)

# Set parameters ----------------------------------------------------------

CpG.folder <- 'extract_refined/Merged'

Cpg.file.pattern <- 'CpG_report.txt.gz'

CPUs <- 2

outdir <- 'bsseq'

M.file <- 'M.txt.gz'
Cov.file <- 'Cov.txt.gz'
gr.file <- 'gr.txt.gz'

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

BS.combined1 <- read.bismark(
  files = file.path(CpG.folder, pdata$Filename[1]),
  sampleNames = pdata$Sample[1],
  rmZeroCov = FALSE,
  strandCollapse = FALSE,
  fileType = "cytosineReport",
  verbose = TRUE)

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

sum(rowSums(getCoverage(BSseq = BS.unstranded)) == 0)
sum(rowSums(getCoverage(BSseq = BS.unstranded)) > 0) / nrow(BS.unstranded)

BS.nonEmpty <- BS.unstranded[rowSums(getCoverage(BSseq = BS.unstranded)) > 0,]
saveRDS(object = BS.nonEmpty, file = file.path(outdir, 'BS.nonEmpty.rds'))
# BS.nonEmpty <- readRDS(file.path(outdir, 'BS.nonEmpty.rds'))

rm(BS.unstranded)

dim(BS.nonEmpty)

# Statistics of non-zero CpGs ------------------------------------

meanCov.nonEmpty <- colMeans(getCoverage(BSseq = BS.nonEmpty))
as.data.frame(meanCov.nonEmpty)
summary(meanCov.nonEmpty)

# Download CG islands -----------------------------------------------------

