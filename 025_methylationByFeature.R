# Version 1.9.2
# devtools::install_github("Bioconductor-mirror/bsseq", ref = "bd60d13")

library(bsseq)
library(GenomicRanges)
library(BSgenome.Btaurus.UCSC.bosTau6)
library(DelayedArray)

options(DelayedArray.block.size=getOption("DelayedArray.block.size") * 10L)


# Parameters of analysis ----

outdir <- 'bsseq'
genomeDir <- 'bostaurus'

# Load source data ----

# Load CGI coordinates
CGI.gr <- readRDS(file = file.path(outdir, 'CpG.gr.rds'))
seqlevels(CGI.gr) <- gsub("^chr(Un_)?", "", seqlevels(CGI.gr))

# Load transcript coordinates
transcripts.gr <- readRDS(
  file.path(genomeDir, "biomart_ensembl_transcripts.rds"))
seqlevels(transcripts.gr) <- gsub("\\.[[:digit:]]$", "", seqlevels(transcripts.gr))
# head(genes.gr)

# Load gene coordinates
genes.gr <- readRDS(file.path(genomeDir, "biomart_ensembl_genes.rds"))
seqlevels(genes.gr) <- gsub("\\.[[:digit:]]$", "", seqlevels(genes.gr))
# head(genes.gr)

# Load BSgenome for chromosome lengths
genome <- BSgenome.Btaurus.UCSC.bosTau6
seqlevels(genome) <- gsub("chr", "", seqlevels(genome))
seqlevels(genome) <- gsub("^M$", "MT", seqlevels(genome))
#seqlengths(genome)

# Generate intergenic probes ----

# Generate 5 kb tiles of the genome
# Truncate the last tile rather than letting tiles span across chromosomes
tiles <- unlist(tileGenome(
  seqlengths = seqinfo(genome),
  tilewidth = 5E3,
  cut.last.tile.in.chrom = TRUE))
# summary(width(tiles))
# sum(width(tiles) < 5E3)

# Keep tiles over 2 kb from nearest gene (i.e. intergenic)
tilesDistanceToGene <- distanceToNearest(x = tiles, subject = genes.gr)
# head(tilesDistanceToGene)
intergenic.gr <- tiles[which(mcols(tilesDistanceToGene)[,"distance"] > 2E3)]

# rm(tilesDistanceToGene, tiles)

# Generate gene body probes ----

tilesGenes <- unlist(tile(x = genes.gr, width = 5E3))

# Generate promoters ranges (wrt. CGI) ----

# Generate promoter coordinates
proms.gr <- promoters(x = transcripts.gr, upstream = 1500, downstream = 500)
# Remove duplicated promoters (shared by multiple transcripts)
proms.gr <- proms.gr[!duplicated(proms.gr)]
# which(duplicated(proms.gr))[1]
# head(proms.gr)

# Promoters overlapping CGI, or not
promOverCGI.bool <- overlapsAny(query = proms.gr, subject = CGI.gr)
promsOverCGI.gr <- proms.gr[promOverCGI.bool]
promsNotOverCGI.gr <- proms.gr[!promOverCGI.bool]
rm(promOverCGI.bool)

# Generate CGI ranges (wrt. promoters) ----

# CGI overlapping promoters, or not
CGIOverProm.bool <- overlapsAny(query = CGI.gr, subject = proms.gr)
CGIOverProm.gr <- CGI.gr[CGIOverProm.bool]
CGINotOverProm.gr <- CGI.gr[!CGIOverProm.bool]
rm(CGIOverProm.bool)

# Load methylation calls ----

if (packageVersion("bsseq") == "1.9.2"){
  BS.unstranded <- readRDS(file = file.path(outdir, "BS.unstranded.rds"))
  # BS.2 <- updateObject(BS.unstranded)
  # saveRDS(BS.2, "bsseq/BS.unstranded_DelayedMatrix.rds")
} else {
  BS.unstranded <- readRDS(file = file.path(outdir, "BS.unstranded_DelayedMatrix.rds"))
}

# Collapse by infection group
BS.infection <- collapseBSseq(
  BSseq = BS.unstranded,
  columns = ifelse(grepl("C", colnames(BS.unstranded)), "Control", "M. bovis"))

rm(BS.unstranded)

# saveRDS(BS.infection, file.path(outdir, "BS.infection_DelayedMatrix.rds"))
# BS.infection <- readRDS(file.path(outdir, "BS.infection_DelayedMatrix.rds"))

# Coverage of each base
baseCoverage <- getCoverage(
  BSseq = BS.infection, type = "Cov", what = "perBase")
# Identify loci >=2 calls in Control
# BScontrol <- BS.infection[as.vector(baseCoverage[,"Control"])]
# Identify probes >= 10 loci >= 2 calls in M. bovis
# BSbovis <- BS.infection[as.vector(baseCoverage[,"M. bovis"])]

# rm(baseCoverage)

# Initalise data.frame for ggplot ----

ggData <- data.frame(
  Region = character(), # type as factor when data.frame complete
  Infection = character(), # type as factor when data.frame complete
  Methylation = numeric()
)

# Process intergenic probes ----

# Identify loci that meet inclusion criteria

BS.2 <- BS.infection[rowSums(baseCoverage >= 2) == 2,]
dim(BS.infection); dim(BS.2)

BS.10 <- BS.infection[rowSums(baseCoverage >= 10) == 2,]
dim(BS.infection); dim(BS.10)

# Count loci >= 2 calls in Control
# mcols(intergenic.gr)[,"lociControl.2"] <- countOverlaps(
#   intergenic.gr,
#   BS.2)
# # Count loci >= 2 calls in M. bovis
# mcols(intergenic.gr)[,"lociMbovis.2"] <- countOverlaps(
#   intergenic.gr,
#   BS.2)
mcols(intergenic.gr)[,"loci.2"] <- countOverlaps(
  intergenic.gr,
  BS.2)
mcols(intergenic.gr)[,"loci.10"] <- countOverlaps(
  intergenic.gr,
  BS.10)

# intergenic.gr <- intergenic.gr[
#   mcols(intergenic.gr)[,"lociControl.2"] >= 10 &
#     mcols(intergenic.gr)[,"lociMbovis.2"] >= 10,]

# intergenic.gr <- intergenic.gr[mcols(intergenic.gr)[,"loci.2"] >= 10,]
# intergenic.gr <- intergenic.gr[mcols(intergenic.gr)[,"loci.10"] >= 10,]

# Calculate % methylation in each region for Control & M. bovis
# NOTE: only loci covered >= 2 reads
methIntergenic <- getMeth(
  BSseq = BS.2,
  regions = intergenic.gr,
  type = "raw",
  what = "perRegion")
colnames(methIntergenic) <- c("meth_Control.2", "meth_Mbovis.2")
values(intergenic.gr) <- cbind(
  values(intergenic.gr),
  DataFrame(as.matrix(methIntergenic)))

methIntergenic <- getMeth(
  BSseq = BS.10,
  regions = intergenic.gr,
  type = "raw",
  what = "perRegion")
colnames(methIntergenic) <- c("meth_Control.10", "meth_Mbovis.10")
values(intergenic.gr) <- cbind(
  values(intergenic.gr),
  DataFrame(as.matrix(methIntergenic)))

rm(methIntergenic)

# Don't forget to subset to loci sufficiently covered
summary(intergenic.gr$loci.2)
plot(density(intergenic.gr$loci.2))
pdf("2017-10-04_test_regions/intergenic_hist.pdf", width = 6, height = 4)
hist(intergenic.gr$loci.2, breaks = seq(0, max(intergenic.gr$loci.2)+10, 10))
dev.off()

sink("2017-10-04_test_regions/intergenic_paired.txt")
with(
  subset(
    as.data.frame(intergenic.gr),
    loci.2 >= 10
  ),
  t.test(
    meth_Control.2,
    meth_Mbovis.2,
    paired = TRUE
  )
)
sink()

sink("2017-10-04_test_regions/intergenic_unpaired.txt")
with(
  subset(
    as.data.frame(intergenic.gr),
    loci.2 >= 10
  ),
  t.test(
    meth_Control.2,
    meth_Mbovis.2,
    paired = FALSE
  )
)
sink()

with(
  subset(
    as.data.frame(intergenic.gr),
    loci.2 >= 10 & seqnames != "X"
  ),
  t.test(
    meth_Control.2,
    meth_Mbovis.2,
    paired = FALSE
  )
)

with(
  subset(
    as.data.frame(intergenic.gr),
    loci.10 >= 10 & seqnames != "X"
  ),
  t.test(
    meth_Control.2,
    meth_Mbovis.2,
    paired = FALSE
  )
)

# Remove X chromosome
# intergenic.gr <- intergenic.gr[seqnames(intergenic.gr) != "X"]

# Count of loci at 2x coverage in each region
pdf("2017-10-04_test_regions/intergenic_coverage_loci.pdf")
plot(sort((intergenic.gr$loci.2), decreasing = TRUE), col="blue", type="l")
# same for 10x coverage
lines(sort((intergenic.gr$loci.10), decreasing = TRUE), col="red", type="l")
abline(h=10, lty="dashed")
legend(
  x = "topright",
  legend = c("2x","10x","cutoff"), col=c("blue","red","black"),
  lty= c(1,1,2))
dev.off()

pdf("2017-10-04_test_regions/intergenic_scatter.pdf", width = 6, height = 6)
with(
  subset(
    as.data.frame(intergenic.gr),
    loci.2 >= 10
  ),
  plot(
    meth_Control.2,
    meth_Mbovis.2,
    cex = 0.1
  )
)
abline(a=0,b=1,col="red")
dev.off()

pdf(
  "2017-10-04_test_regions/intergenic_density_difference.pdf",
  width = 6, height = 6)
with(
  subset(
    as.data.frame(intergenic.gr),
    loci.2 >= 10
  ),
  plot(density(
    meth_Control.2 - meth_Mbovis.2
  ))
)
abline(v = 0, col="red")
dev.off()

pdf(
  "2017-10-04_test_regions/intergenic_density_overlay.pdf",
  width = 6, height = 6)
with(
  subset(
    as.data.frame(intergenic.gr),
    loci.2 >= 10
  ),
  plot(density(
    meth_Control.2
  ), col="blue", lty="dotted")
)
with(
  subset(
    as.data.frame(intergenic.gr),
    loci.2 >= 10
  ),
  abline(density(
    meth_Mbovis.2
  ), col="red", lty="dotted")
)
dev.off()

# Save values for plot (Control)
ggData <- rbind(
  ggData,
  data.frame(
    Region = "Intergenic",
    Infection = "Control",
    Value = subset(
      x = values(intergenic.gr),
      subset = loci.2 >= 10,
      select = "meth_Control.2",
      drop = TRUE)
  )
)

# Backup
saveRDS(object = ggData, file = file.path(outdir, "Peat_et_al-ggData.rds"))

# Save values for plot (M. bovis)
ggData <- rbind(
  ggData,
  data.frame(
    Region = "Intergenic",
    Infection = "M. bovis",
    Value = subset(
      x = values(intergenic.gr),
      subset = loci.2 >= 10,
      select = "meth_Mbovis.2",
      drop = TRUE)
  )
)

# Backup
saveRDS(object = ggData, file = file.path(outdir, "Peat_et_al-ggData.rds"))

# rm(intergenic.gr)

# Process gene body probes ----

# Count loci >= 2 calls in Control
# mcols(genes.gr)[,"lociControl"] <- countOverlaps(
#   query = genes.gr, subject = BScontrol)
# Count loci >= 2 calls in M. bovis
# mcols(genes.gr)[,"lociMbovis"] <- countOverlaps(
#   query = genes.gr, subject = BSbovis)

mcols(tilesGenes)[,"loci.2"] <- countOverlaps(
  tilesGenes,
  BS.2)
mcols(tilesGenes)[,"loci.10"] <- countOverlaps(
  tilesGenes,
  BS.10)

# Calculate % methylation in each region for Control & M. bovis
# NOTE: only loci covered >= 2 reads
methGene <- getMeth(
  BSseq = BS.2,
  regions = tilesGenes,
  type = "raw",
  what = "perRegion")
colnames(methGene) <- c("meth_Control.2", "meth_Mbovis.2")
values(tilesGenes) <- cbind(values(tilesGenes), DataFrame(methGene))

rm(methGene)

# Don't forget to subset to loci sufficiently covered
summary(tilesGenes$loci.2)
plot(density(tilesGenes$loci.2))
pdf("2017-10-04_test_regions/genebody_hist.pdf", width = 6, height = 4)
hist(tilesGenes$loci.2, breaks = seq(0, max(tilesGenes$loci.2)+10, 10))
dev.off()

sink("2017-10-04_test_regions/genebody_paired.txt")
with(
  subset(
    as.data.frame(tilesGenes),
    loci.2 >= 10
  ),
  t.test(
    meth_Control.2,
    meth_Mbovis.2,
    paired = TRUE
  )
)
sink()

sink("2017-10-04_test_regions/genebody_unpaired.txt")
with(
  subset(
    as.data.frame(tilesGenes),
    loci.2 >= 10
  ),
  t.test(
    meth_Control.2,
    meth_Mbovis.2,
    paired = FALSE
  )
)
sink()

# Count of loci at 2x coverage in each region
pdf("2017-10-04_test_regions/genes_coverage_loci.pdf")
plot(sort((tilesGenes$loci.2), decreasing = TRUE), col="blue", type="l")
# same for 10x coverage
lines(sort((tilesGenes$loci.10), decreasing = TRUE), col="red", type="l")
abline(h=10, lty="dashed")
legend(
  x = "topright",
  legend = c("2x","10x","cutoff"), col=c("blue","red","black"),
  lty= c(1,1,2))
dev.off()

pdf("2017-10-04_test_regions/genes_scatter.pdf", width = 6, height = 6)
with(
  subset(
    as.data.frame(tilesGenes),
    loci.2 >= 10
  ),
  plot(
    meth_Control.2,
    meth_Mbovis.2,
    cex = 0.1
  )
)
abline(a=0,b=1,col="red")
dev.off()

# Save values for plot (Control)
ggData <- rbind(
  ggData,
  data.frame(
    Region = "Gene body",
    Infection = "Control",
    Value = subset(
      x = values(tilesGenes),
      subset = loci.2 >= 10,
      select = "meth_Control",
      drop = TRUE)
  )
)

# Backup
saveRDS(object = ggData, file = file.path(outdir, "Peat_et_al-ggData.rds"))

# Save values for plot (M. bovis)
ggData <- rbind(
  ggData,
  data.frame(
    Region = "Gene body",
    Infection = "M. bovis",
    Value = subset(
      x = values(tilesGenes),
      subset = loci.2 >= 10,
      select = "meth_Mbovis",
      drop = TRUE)
  )
)

# Backup
saveRDS(object = ggData, file = file.path(outdir, "Peat_et_al-ggData.rds"))

# rm(genes.gr)

# Process promoters overlapping CGI ----

# Count loci >= 2 calls in Control
# mcols(promsOverCGI.gr)[,"lociControl"] <- countOverlaps(
#   query = promsNotOverCGI.gr, subject = BScontrol)
# Count loci >= 2 calls in M. bovis
# mcols(promsOverCGI.gr)[,"lociMbovis"] <- countOverlaps(
#   query = promsOverCGI.gr, subject = BSbovis)
mcols(promsOverCGI.gr)[,"loci.2"] <- countOverlaps(
  query = promsOverCGI.gr, subject = BS.2)
mcols(promsOverCGI.gr)[,"loci.10"] <- countOverlaps(
  query = promsOverCGI.gr, subject = BS.10)

# Calculate % methylation in each region for Control & M. bovis
# NOTE: only loci covered >= 2 reads
methPromOverCGI <- getMeth(
  BSseq = BS.2,
  regions = promsOverCGI.gr,
  type = "raw",
  what = "perRegion")
colnames(methPromOverCGI) <- c("meth_Control.2", "meth_Mbovis.2")
values(promsOverCGI.gr) <- cbind(
  values(promsOverCGI.gr), DataFrame(methPromOverCGI))

rm(methPromOverCGI)

# Don't forget to subset to loci sufficiently covered
summary(promsOverCGI.gr$loci.2)
plot(density(promsOverCGI.gr$loci.2))
pdf("2017-10-04_test_regions/promsOverCGI_hist.pdf", width = 6, height = 4)
hist(promsOverCGI.gr$loci.2, breaks = seq(0, max(promsOverCGI.gr$loci.2)+10, 10))
dev.off()

sink("2017-10-04_test_regions/promsOverCGI_paired.txt")
with(
  subset(
    as.data.frame(promsOverCGI.gr),
    loci.2 >= 5
  ),
  t.test(
    meth_Control.2,
    meth_Mbovis.2,
    paired = TRUE
  )
)
sink()

sink("2017-10-04_test_regions/promsOverCGI_unpaired.txt")
with(
  subset(
    as.data.frame(promsOverCGI.gr),
    loci.2 >= 5
  ),
  t.test(
    meth_Control.2,
    meth_Mbovis.2,
    paired = FALSE
  )
)
sink()

# Count of loci at 2x coverage in each region
pdf("2017-10-04_test_regions/promsOverCGI_coverage_loci.pdf")
plot(sort((promsOverCGI.gr$loci.2), decreasing = TRUE), col="blue", type="l")
# same for 10x coverage
lines(sort((promsOverCGI.gr$loci.10), decreasing = TRUE), col="red", type="l")
abline(h=10, lty="dashed")
legend(
  x = "topright",
  legend = c("2x","10x","cutoff"), col=c("blue","red","black"),
  lty= c(1,1,2))
dev.off()

pdf("2017-10-04_test_regions/promsOverCGI_scatter.pdf", width = 6, height = 6)
with(
  subset(
    as.data.frame(promsOverCGI.gr),
    loci.2 >= 10
  ),
  plot(
    meth_Control.2,
    meth_Mbovis.2,
    cex = 0.1
  )
)
abline(a=0,b=1,col="red")
dev.off()


# Save values for plot (Control)
ggData <- rbind(
  ggData,
  data.frame(
    Region = "CGI Promoters",
    Infection = "Control",
    Value = subset(
      x = values(promsOverCGI.gr),
      subset = loci.2 >= 5,
      select = "meth_Control.2",
      drop = TRUE)
  )
)

# Backup
saveRDS(object = ggData, file = file.path(outdir, "Peat_et_al-ggData.rds"))

# Save values for plot (M. bovis)
ggData <- rbind(
  ggData,
  data.frame(
    Region = "CGI Promoters",
    Infection = "M. bovis",
    Value = subset(
      x = values(promsOverCGI.gr),
      subset = loci.2 >= 5,
      select = "meth_Mbovis.2",
      drop = TRUE)
  )
)

# Backup
saveRDS(object = ggData, file = file.path(outdir, "Peat_et_al-ggData.rds"))

# rm(promsOverCGI.gr)

# Process promoters *not* overlapping CGI ----

# Count loci >= 2 calls in Control
# mcols(promsNotOverCGI.gr)[,"lociControl"] <- countOverlaps(
#   query = promsNotOverCGI.gr, subject = BScontrol)
# Count loci >= 2 calls in M. bovis
# mcols(promsNotOverCGI.gr)[,"lociMbovis"] <- countOverlaps(
#   query = promsNotOverCGI.gr, subject = BSbovis)
mcols(promsNotOverCGI.gr)[,"loci.2"] <- countOverlaps(
  query = promsNotOverCGI.gr, subject = BS.2)
mcols(promsNotOverCGI.gr)[,"loci.10"] <- countOverlaps(
  query = promsNotOverCGI.gr, subject = BS.10)

# Calculate % methylation in each region for Control & M. bovis
# NOTE: including loci covered < 2 reads
methPromNotOverCGI <- getMeth(
  BSseq = BS.2,
  regions = promsNotOverCGI.gr,
  type = "raw",
  what = "perRegion")
colnames(methPromNotOverCGI) <- c("meth_Control.2", "meth_Mbovis.2")
values(promsNotOverCGI.gr) <- cbind(
  values(promsNotOverCGI.gr), DataFrame(methPromNotOverCGI))

rm(methPromNotOverCGI)

# Don't forget to subset to loci sufficiently covered
summary(promsNotOverCGI.gr$loci.2)
plot(density(promsNotOverCGI.gr$loci.2))
pdf("2017-10-04_test_regions/promsNotOverCGI_hist.pdf", width = 6, height = 4)
hist(promsNotOverCGI.gr$loci.2, breaks = seq(0, max(promsNotOverCGI.gr$loci.2)+10, 10))
dev.off()

sink("2017-10-04_test_regions/promsNotOverCGI_paired.txt")
with(
  subset(
    as.data.frame(promsNotOverCGI.gr),
    loci.2 >= 5
  ),
  t.test(
    meth_Control.2,
    meth_Mbovis.2,
    paired = TRUE
  )
)
sink()

sink("2017-10-04_test_regions/promsNotOverCGI_unpaired.txt")
with(
  subset(
    as.data.frame(promsNotOverCGI.gr),
    loci.2 >= 5
  ),
  t.test(
    meth_Control.2,
    meth_Mbovis.2,
    paired = FALSE
  )
)
sink()

# Count of loci at 2x coverage in each region
pdf("2017-10-04_test_regions/promsNotOverCGI_coverage_loci.pdf")
plot(sort((promsNotOverCGI.gr$loci.2), decreasing = TRUE), col="blue", type="l")
# same for 10x coverage
lines(sort((promsNotOverCGI.gr$loci.10), decreasing = TRUE), col="red", type="l")
abline(h=10, lty="dashed")
legend(
  x = "topright",
  legend = c("2x","10x","cutoff"), col=c("blue","red","black"),
  lty= c(1,1,2))
dev.off()

pdf("2017-10-04_test_regions/promsNotOverCGI_scatter.pdf", width = 6, height = 6)
with(
  subset(
    as.data.frame(promsNotOverCGI.gr),
    loci.2 >= 10
  ),
  plot(
    meth_Control.2,
    meth_Mbovis.2,
    cex = 0.1
  )
)
abline(a=0,b=1,col="red")
dev.off()


# Save values for plot
ggData <- rbind(
  ggData,
  data.frame(
    Region = "Non-CGI Promoters",
    Infection = "Control",
    Value = subset(
      x = values(promsNotOverCGI.gr),
      subset = loci.2 >= 5,
      select = "meth_Control.2",
      drop = TRUE)
  )
)

# Backup
saveRDS(object = ggData, file = file.path(outdir, "Peat_et_al-ggData.rds"))

# Save values for plot
ggData <- rbind(
  ggData,
  data.frame(
    Region = "Non-CGI Promoters",
    Infection = "M. bovis",
    Value = subset(
      x = values(promsNotOverCGI.gr),
      subset = loci.2 >= 5,
      select = "meth_Mbovis.2",
      drop = TRUE)
  )
)

# Backup
saveRDS(object = ggData, file = file.path(outdir, "Peat_et_al-ggData.rds"))

# rm(promsOverCGI.gr)

# Process CGI overlapping promoters ----

# Count loci >= 2 calls in Control
# mcols(CGIOverProm.gr)[,"lociControl"] <- countOverlaps(
#   query = CGIOverProm.gr, subject = BScontrol)
# Count loci >= 2 calls in M. bovis
# mcols(CGIOverProm.gr)[,"lociMbovis"] <- countOverlaps(
#   query = CGIOverProm.gr, subject = BSbovis)
mcols(CGIOverProm.gr)[,"loci.2"] <- countOverlaps(
  query = CGIOverProm.gr, subject = BS.2)
mcols(CGIOverProm.gr)[,"loci.10"] <- countOverlaps(
  query = CGIOverProm.gr, subject = BS.10)

# Calculate % methylation in each region for Control & M. bovis
# NOTE: only loci covered >= 2 reads
methCGIOverProm <- getMeth(
  BSseq = BS.2,
  regions = CGIOverProm.gr,
  type = "raw",
  what = "perRegion")
colnames(methCGIOverProm) <- c("meth_Control.2", "meth_Mbovis.2")
values(CGIOverProm.gr) <- cbind(
  values(CGIOverProm.gr), DataFrame(methCGIOverProm))

rm(methCGIOverProm)

# Don't forget to subset to loci sufficiently covered
summary(CGIOverProm.gr$loci.2)
plot(density(CGIOverProm.gr$loci.2))
pdf("2017-10-04_test_regions/CGIsOverProm_hist.pdf", width = 6, height = 4)
hist(CGIOverProm.gr$loci.2, breaks = seq(0, max(CGIOverProm.gr$loci.2)+10, 10))
dev.off()

sink("2017-10-04_test_regions/CGIsOverProm_paired.txt")
with(
  subset(
    as.data.frame(CGIOverProm.gr, row.names = NULL),
    loci.2 >= 4
  ),
  t.test(
    meth_Control.2,
    meth_Mbovis.2,
    paired = TRUE
  )
)
sink()

sink("2017-10-04_test_regions/CGIsOverProm_unpaired.txt")
with(
  subset(
    as.data.frame(CGIOverProm.gr, row.names = NULL),
    loci.2 >= 4
  ),
  t.test(
    meth_Control.2,
    meth_Mbovis.2,
    paired = FALSE
  )
)
sink()

# Count of loci at 2x coverage in each region
pdf("2017-10-04_test_regions/CGIOverProm_coverage_loci.pdf")
plot(sort((CGIOverProm.gr$loci.2), decreasing = TRUE), col="blue", type="l")
# same for 10x coverage
lines(sort((CGIOverProm.gr$loci.10), decreasing = TRUE), col="red", type="l")
abline(h=10, lty="dashed")
legend(
  x = "topright",
  legend = c("2x","10x","cutoff"), col=c("blue","red","black"),
  lty= c(1,1,2))
dev.off()

pdf("2017-10-04_test_regions/CGIOverProm_scatter.pdf", width = 6, height = 6)
with(
  subset(
    as.data.frame(CGIOverProm.gr, row.names = NULL),
    loci.2 >= 4
  ),
  plot(
    meth_Control.2,
    meth_Mbovis.2,
    cex = 0.1
  )
)
abline(a=0,b=1,col="red")
dev.off()


# Save values for plot (Control)
ggData <- rbind(
  ggData,
  data.frame(
    Region = "Promoter CGIs",
    Infection = "Control",
    Value = subset(
      x = values(CGIOverProm.gr),
      subset = loci.2 >= 4,
      select = "meth_Control.2",
      drop = TRUE)
  )
)

# Backup
saveRDS(object = ggData, file = file.path(outdir, "Peat_et_al-ggData.rds"))

# Save values for plot (M. bovis)
ggData <- rbind(
  ggData,
  data.frame(
    Region = "Promoter CGIs",
    Infection = "M. bovis",
    Value = subset(
      x = values(CGIOverProm.gr),
      subset = loci.2 >= 4,
      select = "meth_Mbovis.2",
      drop = TRUE)
  )
)

# Backup
saveRDS(object = ggData, file = file.path(outdir, "Peat_et_al-ggData.rds"))

# Process CGI *not* overlapping promoters ----

# Count loci >= 2 calls in Control
# mcols(CGINotOverProm.gr)[,"lociControl"] <- countOverlaps(
#   query = CGINotOverProm.gr, subject = BScontrol)
# Count loci >= 2 calls in M. bovis
# mcols(CGINotOverProm.gr)[,"lociMbovis"] <- countOverlaps(
#   query = CGINotOverProm.gr, subject = BSbovis)
mcols(CGINotOverProm.gr)[,"loci.2"] <- countOverlaps(
  query = CGINotOverProm.gr, subject = BS.2)
mcols(CGINotOverProm.gr)[,"loci.10"] <- countOverlaps(
  query = CGINotOverProm.gr, subject = BS.10)

# Calculate % methylation in each region for Control & M. bovis
# NOTE: including loci covered < 2 reads
methCGINotOverProm <- getMeth(
  BSseq = BS.2,
  regions = CGINotOverProm.gr,
  type = "raw",
  what = "perRegion")
colnames(methCGINotOverProm) <- c("meth_Control.2", "meth_Mbovis.2")
values(CGINotOverProm.gr) <- cbind(
  values(CGINotOverProm.gr), DataFrame(as.matrix(methCGINotOverProm)))

rm(methCGINotOverProm)

# Don't forget to subset to loci sufficiently covered
summary(CGINotOverProm.gr$loci.2)
plot(density(CGINotOverProm.gr$loci.2))
pdf("2017-10-04_test_regions/CGINotOverProm_hist.pdf", width = 6, height = 4)
hist(CGINotOverProm.gr$loci.2, breaks = seq(0, max(CGINotOverProm.gr$loci.2)+10, 10))
dev.off()

sink("2017-10-04_test_regions/CGINotOverProm_paired.txt")
with(
  subset(
    as.data.frame(CGINotOverProm.gr, row.names = NULL),
    loci.2 >= 4
  ),
  t.test(
    meth_Control.2,
    meth_Mbovis.2,
    paired = TRUE
  )
)
sink()

sink("2017-10-04_test_regions/CGINotOverProm_unpaired.txt")
with(
  subset(
    as.data.frame(CGINotOverProm.gr, row.names = NULL),
    loci.2 >= 4
  ),
  t.test(
    meth_Control.2,
    meth_Mbovis.2,
    paired = FALSE
  )
)
sink()

# Count of loci at 2x coverage in each region
pdf("2017-10-04_test_regions/CGINotOverProm_coverage_loci.pdf")
plot(sort((CGINotOverProm.gr$loci.2), decreasing = TRUE), col="blue", type="l")
# same for 10x coverage
lines(sort((CGINotOverProm.gr$loci.10), decreasing = TRUE), col="red", type="l")
abline(h=10, lty="dashed")
legend(
  x = "topright",
  legend = c("2x","10x","cutoff"), col=c("blue","red","black"),
  lty= c(1,1,2))
dev.off()

pdf("2017-10-04_test_regions/CGINotOverProm_scatter.pdf", width = 6, height = 6)
with(
  subset(
    as.data.frame(CGINotOverProm.gr, row.names = NULL),
    loci.2 >= 4
  ),
  plot(
    meth_Control.2,
    meth_Mbovis.2,
    cex = 0.1
  )
)
abline(a=0,b=1,col="red")
dev.off()



# Save values for plot (Control)
ggData <- rbind(
  ggData,
  data.frame(
    Region = "Non-Promoter CGIs",
    Infection = "Control",
    Value = subset(
      x = values(CGINotOverProm.gr),
      subset = loci.2 >= 4,
      select = "meth_Control.2",
      drop = TRUE)
  )
)

# Backup
saveRDS(object = ggData, file = file.path(outdir, "Peat_et_al-ggData.rds"))

# Save values for plot (M. bovis)
ggData <- rbind(
  ggData,
  data.frame(
    Region = "Non-Promoter CGIs",
    Infection = "M. bovis",
    Value = subset(
      x = values(CGINotOverProm.gr),
      subset = loci.2 >= 4,
      select = "meth_Mbovis.2",
      drop = TRUE)
  )
)

# Backup
saveRDS(object = ggData, file = file.path(outdir, "Peat_et_al-ggData.rds"))

ggData <- readRDS(file = file.path(outdir, "Peat_et_al-ggData.rds"))

# tables count of hemi-methylated CGIs not over promoter

table(subset(as.data.frame(mcols(CGINotOverProm.gr)), lociControl > 10)$meth_Control > 2/3)
table(subset(as.data.frame(mcols(CGINotOverProm.gr)), lociControl > 10)$meth_Control < 1/3)
with(subset(as.data.frame(mcols(CGINotOverProm.gr)), lociControl > 10), table(meth_Control > 1/3 & meth_Control < 2/3))
table(subset(as.data.frame(mcols(CGINotOverProm.gr)), lociMbovis > 10)$meth_Mbovis > 2/3)
table(subset(as.data.frame(mcols(CGINotOverProm.gr)), lociMbovis > 10)$meth_Mbovis < 1/3)
with(subset(as.data.frame(mcols(CGINotOverProm.gr)), lociMbovis > 10), table(meth_Mbovis > 1/3 & meth_Mbovis < 2/3))

table(subset(as.data.frame(mcols(CGINotOverProm.gr)), lociControl > 10)$meth_Control < 1/2)
table(subset(as.data.frame(mcols(CGINotOverProm.gr)), lociControl > 10)$meth_Control > 1/2)
table(subset(as.data.frame(mcols(CGINotOverProm.gr)), lociMbovis > 10)$meth_Mbovis < 1/2)
table(subset(as.data.frame(mcols(CGINotOverProm.gr)), lociMbovis > 10)$meth_Mbovis > 1/2)
# ggplot ----

ggData$Region <- factor(
  x = ggData$Region,
  levels = c(
    "Intergenic",
    "Gene body",
    "CGI Promoters",
    "Non-CGI Promoters",
    "Promoter CGIs",
    "Non-Promoter CGIs"))

ggData$Infection <- factor(
  x = ggData$Infection,
  levels = c("Control", "M. bovis"))

# Backup
saveRDS(object = ggData, file = file.path(outdir, "Peat_et_al-ggData.rds"))
# ggData <- readRDS(file = file.path(outdir, "Peat_et_al-ggData.rds"))

library(ggplot2)
library(RColorBrewer)

colour.infection <- brewer.pal(6, "Paired")[c(2,6)]
names(colour.infection) <- c("Control", "M. bovis")

fill.infection <- brewer.pal(6, "Paired")[c(1,5)]
names(colour.infection) <- c("Control", "M. bovis")

gg <- ggplot(
  data = ggData,
  mapping = aes(
    x = Infection,
    y = Value,
    colour = Infection,
    fill = Infection)
) +
  facet_grid(. ~ Region) +
  scale_fill_manual(values = fill.infection) +
  scale_colour_manual(values = colour.infection)
# + scale_colour_discrete(l = 50)

gg + geom_boxplot() + theme_bw()

ggsave(
  filename = file.path(outdir, "Peat_al_al-boxPlot_paired.pdf"),
  plot = gg + geom_boxplot() + theme_bw(),
  width = 10,
  height = 5)

gg + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5) +
  theme_bw()

ggsave(
  filename = file.path(outdir, "Peat_al_al-violinPlot_paired.pdf"),
  plot = gg + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5) + theme_bw(),
  width = 10,
  height = 5)

# Count of regions in each box ----

regionsCount <- table(ggData[,1:2])

library(reshape2)

regionsCount <- dcast(
  data = as.data.frame(regionsCount),
  formula = Region ~ Infection,
  value.var = "Freq")

write.csv(
  x = as.data.frame(x = regionsCount),
  file = file.path(outdir, "Peat_et_al-regionCount.csv"),
  row.names = FALSE)

# Test if distribution is different in each feature ----

IR.aov <- aov(Value ~ Region*Infection, subset(ggData, Region != "Intergenic"))
summary(IR.aov)
coefficients(IR.aov)
anova(IR.aov)

summary(ggData)
wilcox.test(Value ~ Infection, ggData, Region == "Intergenic", conf.int=TRUE)
wilcox.test(Value ~ Infection, ggData, Region == "Gene body", conf.int=TRUE)
wilcox.test(Value ~ Infection, ggData, Region == "CGI Promoters", conf.int=TRUE)
wilcox.test(Value ~ Infection, ggData, Region == "Non-CGI Promoters", conf.int=TRUE)
wilcox.test(Value ~ Infection, ggData, Region == "Promoter CGIs", conf.int=TRUE, conf.level=0.99)
wilcox.test(Value ~ Infection, ggData, Region == "Non-Promoter CGIs", conf.int=TRUE)

t.test(Value ~ Infection, ggData, Region == "Intergenic")
t.test(Value ~ Infection, ggData, Region == "Gene body")
t.test(Value ~ Infection, ggData, Region == "CGI Promoters")
t.test(Value ~ Infection, ggData, Region == "Non-CGI Promoters")
t.test(Value ~ Infection, ggData, Region == "Promoter CGIs")
t.test(Value ~ Infection, ggData, Region == "Non-Promoter CGIs")

table(subset(ggData, Region=="Intergenic", "Infection"))
