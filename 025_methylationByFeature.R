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

rm(tilesDistanceToGene, tiles)

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
  BSseq = BS.infection, type = "Cov", what = "perBase") > 2
# Identify loci >=2 calls in Control
BScontrol <- BS.infection[as.vector(baseCoverage[,"Control"])]
# Identify probes >= 10 loci >= 2 calls in M. bovis
BSbovis <- BS.infection[as.vector(baseCoverage[,"M. bovis"])]

rm(baseCoverage)

# Initalise data.frame for ggplot ----

ggData <- data.frame(
  Region = character(), # type as factor when data.frame complete
  Infection = character(), # type as factor when data.frame complete
  Methylation = numeric()
)

# Process intergenic probes ----

# Count loci >= 2 calls in Control
mcols(intergenic.gr)[,"lociControl"] <- countOverlaps(
  query = intergenic.gr, subject = BScontrol)
# Count loci >= 2 calls in M. bovis
mcols(intergenic.gr)[,"lociMbovis"] <- countOverlaps(
  query = intergenic.gr, subject = BSbovis)

# Calculate % methylation in each region for Control & M. bovis
# NOTE: including loci covered < 2 reads
methIntergenic <- getMeth(
  BSseq = BS.infection,
  regions = intergenic.gr,
  type = "raw",
  what = "perRegion")
colnames(methIntergenic) <- c("meth_Control", "meth_Mbovis")
values(intergenic.gr) <- cbind(values(intergenic.gr), DataFrame(as.matrix(methIntergenic)))

rm(methIntergenic)

# Save values for plot (Control)
ggData <- rbind(
  ggData,
  data.frame(
    Region = "Intergenic",
    Infection = "Control",
    Value = subset(
      x = values(intergenic.gr),
      subset = lociControl >= 10,
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
    Region = "Intergenic",
    Infection = "M. bovis",
    Value = subset(
      x = values(intergenic.gr),
      subset = lociMbovis >= 10,
      select = "meth_Mbovis",
      drop = TRUE)
  )
)

# Backup
saveRDS(object = ggData, file = file.path(outdir, "Peat_et_al-ggData.rds"))

# rm(intergenic.gr)

# Process gene body probes ----

# Count loci >= 2 calls in Control
mcols(genes.gr)[,"lociControl"] <- countOverlaps(
  query = genes.gr, subject = BScontrol)
# Count loci >= 2 calls in M. bovis
mcols(genes.gr)[,"lociMbovis"] <- countOverlaps(
  query = genes.gr, subject = BSbovis)

# Calculate % methylation in each region for Control & M. bovis
# NOTE: including loci covered < 2 reads
methGene <- getMeth(
  BSseq = BS.infection,
  regions = genes.gr,
  type = "raw",
  what = "perRegion")
colnames(methGene) <- c("meth_Control", "meth_Mbovis")
values(genes.gr) <- cbind(values(genes.gr), DataFrame(methGene))

rm(methGene)

# Save values for plot (Control)
ggData <- rbind(
  ggData,
  data.frame(
    Region = "Gene body",
    Infection = "Control",
    Value = subset(
      x = values(genes.gr),
      subset = lociControl >= 10,
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
      x = values(genes.gr),
      subset = lociMbovis >= 10,
      select = "meth_Mbovis",
      drop = TRUE)
  )
)

# Backup
saveRDS(object = ggData, file = file.path(outdir, "Peat_et_al-ggData.rds"))

# rm(genes.gr)

# Process promoters overlapping CGI ----

# Count loci >= 2 calls in Control
mcols(promsOverCGI.gr)[,"lociControl"] <- countOverlaps(
  query = promsOverCGI.gr, subject = BScontrol)
# Count loci >= 2 calls in M. bovis
mcols(promsOverCGI.gr)[,"lociMbovis"] <- countOverlaps(
  query = promsOverCGI.gr, subject = BSbovis)

# Calculate % methylation in each region for Control & M. bovis
# NOTE: including loci covered < 2 reads
methPromOverCGI <- getMeth(
  BSseq = BS.infection,
  regions = promsOverCGI.gr,
  type = "raw",
  what = "perRegion")
colnames(methPromOverCGI) <- c("meth_Control", "meth_Mbovis")
values(promsOverCGI.gr) <- cbind(
  values(promsOverCGI.gr), DataFrame(methPromOverCGI))

rm(methPromOverCGI)

# Save values for plot (Control)
ggData <- rbind(
  ggData,
  data.frame(
    Region = "CGI Promoters",
    Infection = "Control",
    Value = subset(
      x = values(promsOverCGI.gr),
      subset = lociControl >= 5,
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
    Region = "CGI Promoters",
    Infection = "M. bovis",
    Value = subset(
      x = values(promsOverCGI.gr),
      subset = lociMbovis >= 5,
      select = "meth_Mbovis",
      drop = TRUE)
  )
)

# Backup
saveRDS(object = ggData, file = file.path(outdir, "Peat_et_al-ggData.rds"))

# rm(promsOverCGI.gr)

# Process promoters *not* overlapping CGI ----

# Count loci >= 2 calls in Control
mcols(promsNotOverCGI.gr)[,"lociControl"] <- countOverlaps(
  query = promsNotOverCGI.gr, subject = BScontrol)
# Count loci >= 2 calls in M. bovis
mcols(promsNotOverCGI.gr)[,"lociMbovis"] <- countOverlaps(
  query = promsNotOverCGI.gr, subject = BSbovis)

# Calculate % methylation in each region for Control & M. bovis
# NOTE: including loci covered < 2 reads
methPromNotOverCGI <- getMeth(
  BSseq = BS.infection,
  regions = promsNotOverCGI.gr,
  type = "raw",
  what = "perRegion")
colnames(methPromNotOverCGI) <- c("meth_Control", "meth_Mbovis")
values(promsNotOverCGI.gr) <- cbind(
  values(promsNotOverCGI.gr), DataFrame(methPromNotOverCGI))

rm(methPromNotOverCGI)

# Save values for plot
ggData <- rbind(
  ggData,
  data.frame(
    Region = "Non-CGI Promoters",
    Infection = "Control",
    Value = subset(
      x = values(promsNotOverCGI.gr),
      subset = lociControl >= 5,
      select = "meth_Control",
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
      subset = lociMbovis >= 5,
      select = "meth_Mbovis",
      drop = TRUE)
  )
)

# Backup
saveRDS(object = ggData, file = file.path(outdir, "Peat_et_al-ggData.rds"))

# rm(promsOverCGI.gr)

# Process CGI overlapping promoters ----

# Count loci >= 2 calls in Control
mcols(CGIOverProm.gr)[,"lociControl"] <- countOverlaps(
  query = CGIOverProm.gr, subject = BScontrol)
# Count loci >= 2 calls in M. bovis
mcols(CGIOverProm.gr)[,"lociMbovis"] <- countOverlaps(
  query = CGIOverProm.gr, subject = BSbovis)

# Calculate % methylation in each region for Control & M. bovis
# NOTE: including loci covered < 2 reads
methCGIOverProm <- getMeth(
  BSseq = BS.infection,
  regions = CGIOverProm.gr,
  type = "raw",
  what = "perRegion")
colnames(methCGIOverProm) <- c("meth_Control", "meth_Mbovis")
values(CGIOverProm.gr) <- cbind(
  values(CGIOverProm.gr), DataFrame(methCGIOverProm))

rm(methCGIOverProm)

# Save values for plot (Control)
ggData <- rbind(
  ggData,
  data.frame(
    Region = "Promoter CGIs",
    Infection = "Control",
    Value = subset(
      x = values(CGIOverProm.gr),
      subset = lociControl >= 4,
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
    Region = "Promoter CGIs",
    Infection = "M. bovis",
    Value = subset(
      x = values(CGIOverProm.gr),
      subset = lociMbovis >= 4,
      select = "meth_Mbovis",
      drop = TRUE)
  )
)

# Backup
saveRDS(object = ggData, file = file.path(outdir, "Peat_et_al-ggData.rds"))

# Process CGI *not* overlapping promoters ----

# Count loci >= 2 calls in Control
mcols(CGINotOverProm.gr)[,"lociControl"] <- countOverlaps(
  query = CGINotOverProm.gr, subject = BScontrol)
# Count loci >= 2 calls in M. bovis
mcols(CGINotOverProm.gr)[,"lociMbovis"] <- countOverlaps(
  query = CGINotOverProm.gr, subject = BSbovis)

# Calculate % methylation in each region for Control & M. bovis
# NOTE: including loci covered < 2 reads
methCGINotOverProm <- getMeth(
  BSseq = BS.infection,
  regions = CGINotOverProm.gr,
  type = "raw",
  what = "perRegion")
colnames(methCGINotOverProm) <- c("meth_Control", "meth_Mbovis")
values(CGINotOverProm.gr) <- cbind(
  values(CGINotOverProm.gr), DataFrame(as.matrix(methCGINotOverProm)))

rm(methCGINotOverProm)

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

# Save values for plot (Control)
ggData <- rbind(
  ggData,
  data.frame(
    Region = "Non-Promoter CGIs",
    Infection = "Control",
    Value = subset(
      x = values(CGINotOverProm.gr),
      subset = lociControl >= 4,
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
    Region = "Non-Promoter CGIs",
    Infection = "M. bovis",
    Value = subset(
      x = values(CGINotOverProm.gr),
      subset = lociMbovis >= 4,
      select = "meth_Mbovis",
      drop = TRUE)
  )
)

# Backup
saveRDS(object = ggData, file = file.path(outdir, "Peat_et_al-ggData.rds"))

ggData <- readRDS(file = file.path(outdir, "Peat_et_al-ggData.rds"))

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

colour.infection <- brewer.pal(3, "Set1")[1:2]
names(colour.infection) <- c("M. bovis", "Control")

gg <- ggplot(
  data = ggData,
  mapping = aes(
    x = Infection,
    y = Value,
    colour = Infection,
    fill = Infection)
) +
  facet_grid(. ~ Region) +
  scale_fill_manual(values = colour.infection) +
  scale_colour_manual(values = colour.infection)
# + scale_colour_discrete(l = 50)

gg + geom_boxplot()

ggsave(
  filename = file.path(outdir, "Peat_al_al-boxPlot.pdf"),
  plot = gg + geom_boxplot(),
  width = 10,
  height = 5)

gg + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5) +
  theme_bw()

ggsave(
  filename = file.path(outdir, "Peat_al_al-violinPlot_v2.pdf"),
  plot = gg + geom_violin(draw_quantiles = c(0.25, 0.5, 0.75), alpha = 0.5),
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
