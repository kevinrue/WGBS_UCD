
library(bsseq)
library(ggplot2)

# Set parameters ----------------------------------------------------------

outdir <- 'bsseq'

# BS <- readRDS(file.path(outdir, "BS.unstranded.rds"))
BS <- readRDS(file.path(outdir, "BS.rmZero.rds"))
CGI <- readRDS(file.path(outdir, "CpG.gr.rds"))
genes <- readRDS(file.path(outdir, "biomart_ensembl_genes.rds"))

seqlevels(CGI) <- gsub("^chr(Un_)?", "", seqlevels(CGI))
seqlevels(BS) <- gsub("(\\.1)?", "", seqlevels(BS))

# Define cutoffs to keep CGIs ---------------------------------------------

mcols(CGI)[,"sumCovTotal"] <- rowSums(getCoverage(BSseq = BS, regions = CGI, type = "Cov", what = "perRegionTotal"))
mcols(CGI)[,"sumMTotal"] <- rowSums(getCoverage(BSseq = BS, regions = CGI, type = "M", what = "perRegionTotal"))
mcols(CGI)[,"MethTotal"] <- mcols(CGI)[,"sumMTotal"] / mcols(CGI)[,"sumCovTotal"]
mcols(CGI)[,"countCG"] <- countOverlaps(query = CGI, subject = BS)

ggplot(data = as.data.frame(mcols(CGI)), mapping = aes(x = seq_along(CGI), y = log2(countCG))) +
  geom_point() +
  geom_hline(yintercept = log2(50), colour = "red")

# Define cutoffs to keep CGIs
# e.g. at least 500 calls over at least 50 CG
ggplot(data = as.data.frame(mcols(CGI)), mapping = aes(x = log2(countCG), y = log2(sumCovTotal))) +
  geom_point() +
  geom_vline(xintercept = log2(50), colour = "red") +
  geom_hline(yintercept = log2(500), colour = "red")

CGI.pass <- CGI[which(
  mcols(CGI)[,"countCG"] >= 50 &
    mcols(CGI)[,"sumCovTotal"] >= 500)]

# Plot CGI methylation % --------------------------------------------------

ggplot(data = as.data.frame(mcols(CGI.pass))) +
  stat_bin(mapping = aes(x = MethTotal), bins = 20) +
  geom_vline(xintercept = 1/3, colour = "red") +
  geom_vline(xintercept = 2/3, colour = "red")

