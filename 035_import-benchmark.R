
# setup ----

library(data.table)
library(BiocParallel)
library(bsseq)
library(microbenchmark)

CpG.folder <- 'extract_refined/Merged'
Cpg.file.pattern <- 'CpG_report.txt.gz'
CPUs <- 2

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
  Animal = factor(sapply(
    X = sapply(X = strsplit(x = CpG.files, split = '_'), FUN = "[[", 1),
    FUN = "gsub",
    pattern = "[CM]",
    replacement = ""),
    levels = sort(as.numeric(unique(
      sapply(
        X = sapply(X = strsplit(x = CpG.files, split = '_'), FUN = "[[", 1),
        FUN = "gsub",
        pattern = "[CM]",
        replacement = ""))))),
  Filename = CpG.files,
  stringsAsFactors = FALSE
)

# A function that I copied online ~3 years ago, and used since then
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

pdata$Filename[1:2]
# [1] "C1_ATCACG_R1_merged_val_1.fq.gz_bismark_bt2_pe.deduplicated.CpG_report.txt.gz" 
# [2] "C11_CAGATC_R1_merged_val_1.fq.gz_bismark_bt2_pe.deduplicated.CpG_report.txt.gz"

# File sizes:
# C1...gz: 231.1 MB
# C2...gz: 231.4 MB

# default block.size ----

options(DelayedArray.block.size = 4.5E6)

ptime1 <- proc.time()
BS.list <- bplapply(
  1:2, 
  read.CpGreport.gz,
  pData = pdata,
  folder = CpG.folder,
  BPPARAM = MulticoreParam(workers = CPUs))
ptime2 <- proc.time()
stime <- (ptime2 - ptime1)[3L]
cat(sprintf("done in %.1f secs\n", stime)) # done in 222.9 secs, 209.5 secs (2 attempts)
BS.combined <- combineList(BS.list)
ptime3 <- proc.time()
stime <- (ptime3 - ptime1)[3L]
cat(sprintf("done in %.1f secs\n", stime)) # done in 391.4 secs (including combineList)
stime <- (ptime3 - ptime2)[3L]
cat(sprintf("done in %.1f secs\n", stime)) # done in 181.8 secs (combineList)

# 10x increased block.size ----

options(DelayedArray.block.size = 45E6)

ptime1 <- proc.time()
BS.list <- bplapply(
  1:2, 
  read.CpGreport.gz,
  pData = pdata,
  folder = CpG.folder,
  BPPARAM = MulticoreParam(workers = CPUs))
ptime2 <- proc.time()
stime <- (ptime2 - ptime1)[3L]
cat(sprintf("done in %.1f secs\n", stime)) # done in 153.4 secs, 164.2 secs (2 attempts)
BS.combined <- combineList(BS.list)
ptime3 <- proc.time()
stime <- (ptime3 - ptime1)[3L]
cat(sprintf("done in %.1f secs\n", stime)) # done in 284.5 secs (including combineList)
stime <- (ptime3 - ptime2)[3L]
cat(sprintf("done in %.1f secs\n", stime)) # done in 120.4 secs (combineList)

# using the built-in read.bismark() function (default block size) ----

options(DelayedArray.block.size = 4.5E6)

ptime1 <- proc.time()
bs1 <- read.bismark(
  files = file.path(CpG.folder, pdata$Filename[1:2]),
  sampleNames = pdata$Sample[1:2],
  rmZeroCov = FALSE,
  strandCollapse = FALSE,
  fileType = "cytosineReport",
  mc.cores = 2
)
ptime2 <- proc.time()
stime <- (ptime2 - ptime1)[3L]
cat(sprintf("done in %.1f secs\n", stime)) # done in 409.4 secs (including merge)

# using the built-in read.bismark() function (10x block size) ----

options(DelayedArray.block.size = 45E6)

ptime1 <- proc.time()
bs1 <- read.bismark(
  files = file.path(CpG.folder, pdata$Filename[1:2]),
  sampleNames = pdata$Sample[1:2],
  rmZeroCov = FALSE,
  strandCollapse = FALSE,
  fileType = "cytosineReport",
  mc.cores = 2
)
ptime2 <- proc.time()
stime <- (ptime2 - ptime1)[3L]
cat(sprintf("done in %.1f secs\n", stime)) # done in 281.0 secs (including merge)

# 10x increased block.size on 16 samples ----

# using my local function to time import and combineList separately

options(DelayedArray.block.size = 45E6)

ptime1 <- proc.time()
BS.list <- bplapply(
  1:nrow(pdata), 
  read.CpGreport.gz,
  pData = pdata,
  folder = CpG.folder,
  BPPARAM = MulticoreParam(workers = CPUs))
ptime2 <- proc.time()
stime <- (ptime2 - ptime1)[3L]
cat(sprintf("done in %.1f secs\n", stime)) # done in 1365.3 secs (22min)
BS.combined <- combineList(BS.list)
ptime3 <- proc.time()
stime <- (ptime3 - ptime1)[3L]
cat(sprintf("done in %.1f secs\n", stime)) # done in ___ secs (including combineList)
stime <- (ptime3 - ptime2)[3L]
cat(sprintf("done in %.1f secs\n", stime)) # done in ___ secs (combineList)


# Conclusion ----

# Increasing DelayedArray.block.size 10-fold (namely, to 45E6 instead of the
# current default 4.5E6) seems to have:
# a considerable impact on the import time
# a possible, smaller, impact on the combineList time

# similar performance between the function that I copied a few years ago
# and the current version of read.bismark