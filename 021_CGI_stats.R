
library(bsseq)
library(ggplot2)

# Set parameters ----------------------------------------------------------

genomeDir <- 'bostaurus'
outdir <- 'bsseq'

# Load unstranded CG methylation calls
# BS <- readRDS(file.path(outdir, "BS.unstranded.rds"))
BS <- readRDS(file.path(outdir, "BS.rmZero.rds"))

# Load CGI genomic coordinates
CGI <- readRDS(file.path(outdir, "CpG.gr.rds"))

# Load genomic coordinates of genes
genes <- readRDS(file.path(genomeDir, "biomart_ensembl_genes.rds"))

seqlevels(CGI) <- gsub("^chr(Un_)?", "", seqlevels(CGI))
seqlevels(BS) <- gsub("(\\.1)?", "", seqlevels(BS))

# Define cutoffs to keep CGIs ---------------------------------------------

# TODO: use collapseBSseq to avoid rowSums every time
mcols(CGI)[,"sumCovTotal"] <- rowSums(getCoverage(BSseq = BS, regions = CGI, type = "Cov", what = "perRegionTotal"))
mcols(CGI)[,"sumMTotal"] <- rowSums(getCoverage(BSseq = BS, regions = CGI, type = "M", what = "perRegionTotal"))
mcols(CGI)[,"MethTotal"] <- mcols(CGI)[,"sumMTotal"] / mcols(CGI)[,"sumCovTotal"]
mcols(CGI)[,"countCG"] <- countOverlaps(query = CGI, subject = BS)

ggplot(
  data = as.data.frame(mcols(CGI)),
  mapping = aes(x = seq_along(CGI), y = log2(countCG))) +
  geom_point() +
  geom_hline(yintercept = log2(50), colour = "red") +
  theme(
    axis.text = element_text(size = rel(1.5)),
    axis.title = element_text(size = rel(1.5)),
    legend.text = element_text(size = rel(1.5)),
    legend.title = element_text(size = rel(1.5))
  ) +
  ylab(expression(log[2]*" (CG loci)")) +
  xlab(expression("CGI index"))

ggsave(
  filename = file.path(outdir, "CGIs_counts.pdf"),
  width = 9, height = 7)

# Define cutoffs to keep CGIs
# e.g. at least 500 calls over at least 50 CG
ggplot(
  data = as.data.frame(mcols(CGI)),
  mapping = aes(x = log2(countCG), y = log2(sumCovTotal))) +
  geom_point() +
  geom_vline(xintercept = log2(50), colour = "red") +
  geom_hline(yintercept = log2(500), colour = "red") +
  theme(
    axis.text = element_text(size = rel(1.5)),
    axis.title = element_text(size = rel(1.5)),
    legend.text = element_text(size = rel(1.5)),
    legend.title = element_text(size = rel(1.5))
  ) +
  xlab(expression(log[2]*" (CG loci)")) +
  ylab(expression(log[2]*" (total CG)"))

ggsave(
  filename = file.path(outdir, "CGIs_loci-calls.pdf"),
  width = 9, height = 7)

CGI.pass <- CGI[which(
  mcols(CGI)[,"countCG"] >= 50 &
    mcols(CGI)[,"sumCovTotal"] >= 500)]

# Plot CGI methylation % --------------------------------------------------

ggplot(data = as.data.frame(mcols(CGI.pass))) +
  stat_bin(mapping = aes(x = MethTotal), bins = 20) +
  geom_vline(xintercept = 1/3, colour = "red") +
  geom_vline(xintercept = 2/3, colour = "red") +
  theme(
    axis.text = element_text(size = rel(1.5)),
    axis.title = element_text(size = rel(1.5)),
    legend.text = element_text(size = rel(1.5)),
    legend.title = element_text(size = rel(1.5))
  ) +
  xlab("Methylation %")

ggsave(
  filename = file.path(outdir, "CGIs_meth.pdf"),
  width = 9, height = 7)