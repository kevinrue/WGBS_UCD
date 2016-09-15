
# Set parameters ----------------------------------------------------------

CpGislands.file <- "bostaurus/Bt_UMD31_CpG_islands_unmasked.bed"

outdir <- 'bsseq'

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

saveRDS(object = CpG.gr, file = file.path(outdir, 'CpG.gr.rds'))
