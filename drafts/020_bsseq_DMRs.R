
# Set parameters ----------------------------------------------------------

outdir <- 'bsseq'

# BS <- readRDS(file.path(outdir, "BS.unstranded.rds"))
BS <- readRDS(file.path(outdir, "BS.rmZero.rds"))
sampleNames(BS)
colnames(BS)
colnames(getCoverage(BS))

BS.fit <- BSmooth(BS, mc.cores = 4, verbose = TRUE)
BS.fit