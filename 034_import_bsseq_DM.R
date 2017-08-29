
library(bsseq)

# Set parameters ----------------------------------------------------------

CpG.folder <- 'extract_refined/Merged'

Cpg.file.pattern <- 'CpG_report.txt.gz'

CPUs <- 2

outdir <- 'bsseq'

# Create required folders -------------------------------------------------

if (!dir.exists(outdir)){
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

# Import methylation calls ----

bs1 <- read.bismark(
  files = file.path(CpG.folder, pdata$Filename[1]),
  sampleNames = pdata$Sample[1],
  rmZeroCov = 0,
  strandCollapse = TRUE,
  fileType = "cytosineReport",
  mc.cores = 2
)
